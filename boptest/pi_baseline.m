%% pi_baseline.m: The BOPTEST built-in PI controller (Section 4.4)
%
%  Runs the seven-day scenario with no overrides at all, so the test case's own
%  PI bank drives the heating. It logs the same signals, prints the same KPI
%  report and draws the same comfort figure as the MPC runners, so the runs
%  compare directly.
%
%  BEFORE RUNNING:
%    1. startup
%    2. identification/ekf_sysid_2C     → models/rc_params_2C6z.mat
%    3. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat
%    4. BOPTEST server reachable at http://127.0.0.1:80
%
%  Saves  models/pi[<tag>]_simdata.mat

clearvars -except BOPTEST_TUNC BOPTEST_SUNC BOPTEST_OUT_TAG BOPTEST_START_TIME; clc;

% Optional hooks, set by the sweep drivers:
%   BOPTEST_TUNC        temperature forecast noise: 'medium' | 'high' | ''
%   BOPTEST_SUNC        the same options for solar irradiation
%   BOPTEST_OUT_TAG     suffix appended to the output filename, e.g. '_medium'
%   BOPTEST_START_TIME  scenario start [s]; 0 is the week the thesis reports on,
%                       so any other value runs the controller out of sample
if ~exist('BOPTEST_TUNC','var')       || isempty(BOPTEST_TUNC),       BOPTEST_TUNC       = ''; end
if ~exist('BOPTEST_SUNC','var')       || isempty(BOPTEST_SUNC),       BOPTEST_SUNC       = ''; end
if ~exist('BOPTEST_OUT_TAG','var')    || isempty(BOPTEST_OUT_TAG),    BOPTEST_OUT_TAG    = ''; end
if ~exist('BOPTEST_START_TIME','var') || isempty(BOPTEST_START_TIME), BOPTEST_START_TIME = 0;  end

%% 1. Zone metadata and signal names
load('models/rc_params_2C6z.mat');             % ZONE_NAMES, nZ, Ts
p = load('models/tracking_mpc_params.mat');    % occ_sp_threshold, only for the KPI report

nU = nZ - 1;   % 5 controlled zones, Hal is uncontrolled

SIG_TZ = {'conHeaLiv_reaTZon_y','conHeaRo1_reaTZon_y','conHeaRo2_reaTZon_y', ...
          'conHeaRo3_reaTZon_y','conHeaBth_reaTZon_y','reaTHal_y'};
SIG_QH = {'reaHeaLiv_y','reaHeaRo1_y','reaHeaRo2_y','reaHeaRo3_y','reaHeaBth_y'};

%% 2. Scenario (identical across every controller)
start_time = BOPTEST_START_TIME;
warm_up    = 7*24*3600;
n_days     = 7;
n_steps    = round(n_days * 86400 / Ts);   % 672 steps

%% 3. Deploy BOPTEST
handler = RequestHandler('http://127.0.0.1', 80);
tid     = handler.deploy_test('multizone_residential_hydronic');
handler.set_scenario(start_time, warm_up, BOPTEST_TUNC, BOPTEST_SUNC);
handler.set_step(Ts);
handler.advance(struct());   % first step, the PI is already running

% Under forecast uncertainty BOPTEST caps the /forecast horizon at 48 h, so the
% whole-episode prefetch is replaced by a per-step fetch inside the loop.
USE_NOISY_FCAST = ~isempty(BOPTEST_TUNC) || ~isempty(BOPTEST_SUNC);

%% 4. Setpoint schedule
%  Logged raw, exactly as BOPTEST publishes it, so the figures and the local
%  diagnostics sit on the same band the official KPI is computed against.
n_total = n_steps + 10;

setp_sigs = cell(1, 2*nU);
for j = 1:nU
    setp_sigs{2*j-1} = sprintf('LowerSetp[%s]', ZONE_NAMES{j});
    setp_sigs{2*j  } = sprintf('UpperSetp[%s]', ZONE_NAMES{j});
end

if USE_NOISY_FCAST
    lb_full = zeros(n_steps, nU);   % filled step by step inside the loop
    ub_full = zeros(n_steps, nU);
    Tamb_fc = zeros(n_steps, 1);
    fprintf('[pi_baseline] noisy fcast: per-step refetch (T_unc=%s, S_unc=%s)\n', ...
            BOPTEST_TUNC, BOPTEST_SUNC);
else
    lb_full = zeros(n_total, nU);
    ub_full = zeros(n_total, nU);
    fc_setp = handler.get_forecast(setp_sigs, n_total*Ts, Ts);
    for j = 1:nU
        lb_full(:,j) = fc_setp.(sprintf('LowerSetp_%s_', ZONE_NAMES{j}))(1:n_total);
        ub_full(:,j) = fc_setp.(sprintf('UpperSetp_%s_', ZONE_NAMES{j}))(1:n_total);
    end
    fprintf('LowerSetp sample: '); fprintf('%.1fC ', lb_full(1:4,1)-273.15); fprintf('\n');
    fprintf('UpperSetp sample: '); fprintf('%.1fC ', ub_full(1:4,1)-273.15); fprintf('\n\n');

    fc_tamb = handler.get_forecast({'TDryBul'}, n_steps*Ts, Ts);
    Tamb_fc = fc_tamb.TDryBul(1:n_steps);
end

%% 5. Logging buffers
log_Tz   = zeros(n_steps, nZ);
log_Qhea = zeros(n_steps, nU);
log_ref  = zeros(n_steps, nU);   % LowerSetp [K]
log_ub   = zeros(n_steps, nU);   % UpperSetp [K]

fprintf('PI baseline running: %d zones | %g day(s) | %d steps\n\n', nU, n_days, n_steps);

%% 6. Simulation loop (no overrides, the PI runs inside BOPTEST)
for i = 1:n_steps

    if USE_NOISY_FCAST
        fc_setp_step = handler.get_forecast(setp_sigs, 0, Ts);   % current sample only
        for j = 1:nU
            lb_full(i, j) = fc_setp_step.(sprintf('LowerSetp_%s_', ZONE_NAMES{j}))(1);
            ub_full(i, j) = fc_setp_step.(sprintf('UpperSetp_%s_', ZONE_NAMES{j}))(1);
        end
        fc_tamb_step = handler.get_forecast({'TDryBul'}, 0, Ts);
        Tamb_fc(i)   = fc_tamb_step.TDryBul(1);
    end

    res = handler.advance(struct());   % zero overrides

    for zi = 1:nZ, log_Tz(i,zi)  = res.(SIG_TZ{zi}); end
    for j  = 1:nU, log_Qhea(i,j) = res.(SIG_QH{j});  end
    log_ref(i,:) = lb_full(min(i, end), :);
    log_ub(i,:)  = ub_full(min(i, end), :);

    if mod(i,16)==0
        fprintf('  step %4d/%d  Tamb=%+.1fC  T_air=[', ...
                i, n_steps, Tamb_fc(min(i,end))-273.15);
        fprintf('%.1f ', log_Tz(i,1:nU)-273.15);
        fprintf(']\n');
    end
end

%% 7. KPIs and figure
kpi = handler.get_kpi();
handler.stop_test(tid);

log_Tamb = Tamb_fc(:);
metrics  = kpi_report('pi_baseline', kpi, log_Tz, log_Qhea, log_ref, log_ub, ...
                      ZONE_NAMES, Ts, p.occ_sp_threshold);
comfort_plot('pi_baseline', log_Tz, log_ref, log_ub, log_Tamb, kpi, ZONE_NAMES, Ts);

%% 8. Save
sim_pi.t       = (0:n_steps-1)' * Ts;
sim_pi.Tz      = log_Tz;
sim_pi.Qhea    = log_Qhea;
sim_pi.ref     = log_ref;      % raw lower setpoint schedule [K]
sim_pi.ref_ub  = log_ub;       % raw upper setpoint schedule [K]
sim_pi.Tamb    = log_Tamb;
sim_pi.kpi     = kpi;
sim_pi.metrics = metrics;
sim_pi.BOPTEST_TUNC = BOPTEST_TUNC;
sim_pi.BOPTEST_SUNC = BOPTEST_SUNC;

out_fname = sprintf('pi%s_simdata.mat', BOPTEST_OUT_TAG);
save(fullfile('models', out_fname), 'sim_pi', 'kpi', 'metrics', 'Ts', ...
     'BOPTEST_TUNC', 'BOPTEST_SUNC', 'start_time');
fprintf('[pi_baseline] Saved  models/%s\n', out_fname);
