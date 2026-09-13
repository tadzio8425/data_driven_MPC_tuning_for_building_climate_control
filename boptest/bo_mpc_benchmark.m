%% bo_mpc_benchmark.m: Evaluate the BO-tuned MPC on a fresh episode (Chapter 6)
%
%  Loads θ_best from bo_mpc_tune and runs one seven-day episode with the cost
%  parameters frozen. The scenario is the one every other runner uses, so the
%  KPI tables compare directly.
%
%  BEFORE RUNNING:
%    1. startup
%    2. bo_mpc_tune  → models/bo_mpc[_bias]_tune_data.mat
%    3. BOPTEST server reachable at http://127.0.0.1:80
%
%  Saves  models/bo_mpc[_bias][<tag>]_simdata.mat

clearvars -except use_bias BOPTEST_TUNC BOPTEST_SUNC BOPTEST_OUT_TAG ...
                  BO_THETA_TAG BOPTEST_START_TIME; clc;

% Optional hooks, set by the sweep drivers:
%   BOPTEST_TUNC        temperature forecast noise: 'medium' | 'high' | ''
%   BOPTEST_SUNC        the same options for solar irradiation
%   BOPTEST_OUT_TAG     suffix appended to the output filename
%   BO_THETA_TAG        suffix on the input tune-data, i.e. which θ to load
%   BOPTEST_START_TIME  scenario start [s]; θ was tuned at 0, so any other
%                       value evaluates it out of sample
if ~exist('BOPTEST_TUNC','var')       || isempty(BOPTEST_TUNC),       BOPTEST_TUNC       = ''; end
if ~exist('BOPTEST_SUNC','var')       || isempty(BOPTEST_SUNC),       BOPTEST_SUNC       = ''; end
if ~exist('BOPTEST_OUT_TAG','var')    || isempty(BOPTEST_OUT_TAG),    BOPTEST_OUT_TAG    = ''; end
if ~exist('BO_THETA_TAG','var')       || isempty(BO_THETA_TAG),       BO_THETA_TAG       = ''; end
if ~exist('BOPTEST_START_TIME','var') || isempty(BOPTEST_START_TIME), BOPTEST_START_TIME = 0;  end
if ~exist('use_bias','var'),           use_bias = false; end

%% 1. Load the model and the precomputed QP matrices
load('models/rc_params_2C6z.mat');            % ZONE_NAMES, nZ, Ts
nZ_boptest = nZ;
p  = load('models/tracking_mpc_params.mat');
nU = p.nU;  nYc = p.nYc;
clear rl_mpc_step

%% 2. Load θ_best
if use_bias
    tune_fname = sprintf('bo_mpc_bias%s_tune_data.mat', BO_THETA_TAG);
else
    tune_fname = sprintf('bo_mpc%s_tune_data.mat',      BO_THETA_TAG);
end
tune_path = fullfile('models', tune_fname);
if ~exist(tune_path, 'file')
    error('[bo_bench] %s not found. Run bo_mpc_tune first.', tune_path);
end

S        = load(tune_path, 'theta_best', 'x_best', 'L_best', 'use_bias');
theta    = S.theta_best;
x_best   = S.x_best;
L_best   = S.L_best;
use_bias = S.use_bias;   % must match the configuration the search was run under

fprintf('[bo_bench] θ_best from %s  (use_bias=%d)\n', tune_path, use_bias);
fprintf('[bo_bench] J_real at that point during the search = %.4f\n', L_best);
fprintf('[bo_bench] q_occ    = %s\n', mat2str(theta.q_diag(:,1).', 3));
fprintf('[bo_bench] r_u      = %.4g\n', theta.r_u(1));
fprintf('[bo_bench] w_offset = %s\n', mat2str(theta.w_offset.', 3));

%% 3. Scenario and BOPTEST deployment
start_time = BOPTEST_START_TIME;
warm_up    = 7*24*3600;
n_days     = 7;
n_steps    = round(n_days * 86400 / Ts);

SIG_TZ  = {'conHeaLiv_reaTZon_y','conHeaRo1_reaTZon_y','conHeaRo2_reaTZon_y', ...
           'conHeaRo3_reaTZon_y','conHeaBth_reaTZon_y','reaTHal_y'};
SIG_QH  = {'reaHeaLiv_y','reaHeaRo1_y','reaHeaRo2_y','reaHeaRo3_y','reaHeaBth_y'};
SIG_CON = {'InternalGainsCon[Liv]','InternalGainsCon[Ro1]','InternalGainsCon[Ro2]', ...
           'InternalGainsCon[Ro3]','InternalGainsCon[Bth]','InternalGainsCon[Hal]'};
SIG_RAD = {'InternalGainsRad[Liv]','InternalGainsRad[Ro1]','InternalGainsRad[Ro2]', ...
           'InternalGainsRad[Ro3]','InternalGainsRad[Bth]','InternalGainsRad[Hal]'};

handler = RequestHandler('http://127.0.0.1', 80);
tid     = handler.deploy_test('multizone_residential_hydronic');
handler.set_scenario(start_time, warm_up, BOPTEST_TUNC, BOPTEST_SUNC);
handler.set_step(Ts);
res = handler.advance(struct());

% See fetch_forecast_window.m for why a noisy forecast needs a per-step fetch.
USE_NOISY_FCAST = ~isempty(BOPTEST_TUNC) || ~isempty(BOPTEST_SUNC);

%% 4. Forecast prefetch
n_total = n_steps + p.N_p;

setp_sigs = cell(1, 2*nU);
for j = 1:nU
    setp_sigs{2*j-1} = sprintf('LowerSetp[%s]', ZONE_NAMES{j});
    setp_sigs{2*j  } = sprintf('UpperSetp[%s]', ZONE_NAMES{j});
end
dist_sigs = [{'TDryBul','HGloHor'}, SIG_CON(1:nU), SIG_RAD(1:nU)];

if USE_NOISY_FCAST
    Tamb_full = zeros(p.N_p, 1);      % refilled each step by fetch_forecast_window
    HGlo_full = zeros(p.N_p, 1);
    Qint_full = zeros(nU, p.N_p);
    lb_mpc    = zeros(p.N_p, nU);
    ub_mpc    = zeros(p.N_p, nU);
    fprintf('[bo_bench] noisy fcast: per-step refetch (T_unc=%s, S_unc=%s)\n', ...
            BOPTEST_TUNC, BOPTEST_SUNC);
else
    lb_mpc  = zeros(n_total, nU);
    ub_mpc  = zeros(n_total, nU);
    fc_setp = handler.get_forecast(setp_sigs, n_total*Ts, Ts);
    for j = 1:nU
        lb_mpc(:,j) = fc_setp.(sprintf('LowerSetp_%s_', ZONE_NAMES{j}))(1:n_total);
        ub_mpc(:,j) = fc_setp.(sprintf('UpperSetp_%s_', ZONE_NAMES{j}))(1:n_total);
    end

    fc_dist   = handler.get_forecast(dist_sigs, n_total*Ts, Ts);
    Tamb_full = fc_dist.TDryBul(1:n_total);
    HGlo_full = fc_dist.HGloHor(1:n_total);
    Qint_full = zeros(nU, n_total);
    for j = 1:nU
        qc = fc_dist.(regexprep(SIG_CON{j},'[\[\]]','_'))(1:n_total);
        qr = fc_dist.(regexprep(SIG_RAD{j},'[\[\]]','_'))(1:n_total);
        Qint_full(j,:) = qc(:).' + qr(:).';
    end
end

%% 5. Initial state and logging buffers
Tz_all = zeros(nZ_boptest,1);
for i = 1:nZ_boptest, Tz_all(i) = res.(SIG_TZ{i}); end
Tair_init = Tz_all(1:nU);
z_aug     = [Tair_init; Tair_init + 1.5; zeros(p.nd_bias, 1)];
u_prev    = zeros(nU,1);

log_Tz   = zeros(n_steps, nZ_boptest);
log_v    = zeros(n_steps, nU);
log_Qhea = zeros(n_steps, nU);
log_ref  = zeros(n_steps, nU);
log_ub   = zeros(n_steps, nU);
log_Tamb = zeros(n_steps, 1);
log_V    = zeros(n_steps, 1);

fprintf('\n[bo_bench] %d steps  (%g day(s), Ts=%.0fs)  greedy, θ_best frozen\n\n', ...
        n_steps, n_days, Ts);

%% 6. Evaluation loop
for i = 1:n_steps

    if USE_NOISY_FCAST
        [Tamb_full, HGlo_full, Qint_w, lb_mpc, ub_mpc] = fetch_forecast_window( ...
            handler, dist_sigs, setp_sigs, ZONE_NAMES, nU, p.N_p, Ts);
        Qint_full(1:nU, :) = Qint_w;
        i_loc = 1;
        n_loc = p.N_p;
    else
        i_loc = i;
        n_loc = n_total;
    end

    [Tz_meas, dvec, lb_h] = build_step_inputs( ...
        res, lb_mpc, i_loc, p, Tamb_full, HGlo_full, Qint_full, n_loc, SIG_TZ, nU, nYc);

    [v_cmd, ~, z_aug_k, ~, ~, V_k, ~] = ...
        rl_mpc_step(Tz_meas, dvec, lb_h, z_aug, u_prev, theta, [], use_bias);

    ov = struct('oveTSetPumBoi_activate', 1, 'oveTSetPumBoi_u', p.w_occ);
    for j = 1:nU
        zn = ZONE_NAMES{j};
        ov.(sprintf('conHea%s_oveActHea_activate',  zn)) = 1;
        ov.(sprintf('conHea%s_oveActHea_u',         zn)) = v_cmd(j);
        ov.(sprintf('conHea%s_oveTSetHea_activate', zn)) = 0;
    end
    res_next = handler.advance(ov);

    Q_heat  = zeros(nU,1);
    Tz_next = zeros(nU,1);
    for j = 1:nU
        Q_heat(j)  = res_next.(SIG_QH{j});
        Tz_next(j) = res_next.(SIG_TZ{j});
    end

    log_Tz(i,1:nU) = Tz_next.';
    log_Tz(i,end)  = res_next.(SIG_TZ{end});
    log_v(i,:)     = v_cmd.';
    log_Qhea(i,:)  = Q_heat.';
    log_ref(i,:)   = lb_mpc(min(i_loc,end),:);
    log_ub(i,:)    = ub_mpc(min(i_loc,end),:);   % raw schedule, for plots and KPIs
    log_Tamb(i)    = dvec(1);
    log_V(i)       = V_k;

    z_aug  = z_aug_k;
    u_prev = Q_heat;
    res    = res_next;

    if mod(i,16) == 0
        fprintf('  step %4d/%d  V=%8.1f  Tz=[', i, n_steps, V_k);
        fprintf('%5.2f ', Tz_next - 273.15);
        fprintf(']  v=['); fprintf('%4.2f ', v_cmd); fprintf(']\n');
    end
end

%% 7. KPIs, figure and save
kpi = handler.get_kpi();
handler.stop_test(tid);

if use_bias, bias_tag = '_bias'; else, bias_tag = ''; end
ctrl_name = sprintf('bo_mpc%s', bias_tag);

metrics = kpi_report(ctrl_name, kpi, log_Tz, log_Qhea, log_ref, log_ub, ...
                     ZONE_NAMES, Ts, p.occ_sp_threshold);
comfort_plot(ctrl_name, log_Tz, log_ref, log_ub, log_Tamb, kpi, ZONE_NAMES, Ts);

out_fname = sprintf('bo_mpc%s%s_simdata.mat', bias_tag, BOPTEST_OUT_TAG);
save(fullfile('models', out_fname), ...
     'log_Tz','log_v','log_Qhea','log_ref','log_ub','log_Tamb','log_V', ...
     'theta','x_best','L_best','kpi','metrics','Ts','use_bias', ...
     'BOPTEST_TUNC','BOPTEST_SUNC','BO_THETA_TAG','start_time');
fprintf('\n  Saved: models/%s\n', out_fname);
