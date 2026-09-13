%% ekf_mpc.m: The offset-free MPC (Section 4.3)
%
%  Same controller as mpc_baseline.m, with the disturbance estimation switched
%  on. Five bias states are appended to the model, the Kalman filter estimates
%  them from the measured outputs, and the estimate d̂ shifts the steady-state
%  target. That is the integral action the plain MPC lacks.
%
%    Model:       2C RC per zone (air node T_z + internal mass T_i)
%    Estimator:   15-state KF on [x(10); d_bias(5)], fixed DARE gain
%    Controller:  setpoint-tracking QP (Magni & Scattolini §12.2.3), direct
%                 valve commands, BOPTEST's inner PI disabled
%
%  BEFORE RUNNING:
%    1. startup
%    2. identification/ekf_sysid_2C     → models/rc_params_2C6z.mat
%    3. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat
%    4. BOPTEST server reachable at http://127.0.0.1:80
%
%  Saves  models/mpc_simdata.mat, or models/ekf_mpc<tag>_simdata.mat when tagged

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

%% 1. Load the model and the precomputed QP matrices
load('models/rc_params_2C6z.mat');            % ZONE_NAMES, nZ, Ts, nU_ss
nZ_boptest = nZ;      % 6 zones in BOPTEST
nZ         = nU_ss;   % 5 zones in the MPC model

if ~exist('models/tracking_mpc_params.mat', 'file')
    error('[ekf_mpc] Run mpc/tracking_mpc_params_build first.');
end
p = load('models/tracking_mpc_params.mat');
clear tracking_mpc_step

nU = p.nU;  nYc = p.nYc;  nx = p.nx;  nYm = p.nYm;

fprintf('[ekf_mpc] ||Lkf||_F=%.4f   nx=%d   nz=%d   Np=%d   Ts=%.0fs\n', ...
        norm(p.Lkf,'fro'), p.nx, p.nz, p.N_p, Ts);
fprintf('[ekf_mpc] Disturbance estimation active: nd_bias=%d (offset-free MPC)\n', p.nd_bias);

%% 2. Signal names
SIG_TZ  = {'conHeaLiv_reaTZon_y','conHeaRo1_reaTZon_y','conHeaRo2_reaTZon_y', ...
           'conHeaRo3_reaTZon_y','conHeaBth_reaTZon_y','reaTHal_y'};
SIG_QH  = {'reaHeaLiv_y','reaHeaRo1_y','reaHeaRo2_y','reaHeaRo3_y','reaHeaBth_y'};
SIG_CON = {'InternalGainsCon[Liv]','InternalGainsCon[Ro1]','InternalGainsCon[Ro2]', ...
           'InternalGainsCon[Ro3]','InternalGainsCon[Bth]','InternalGainsCon[Hal]'};
SIG_RAD = {'InternalGainsRad[Liv]','InternalGainsRad[Ro1]','InternalGainsRad[Ro2]', ...
           'InternalGainsRad[Ro3]','InternalGainsRad[Bth]','InternalGainsRad[Hal]'};

%% 3. Scenario (identical across every controller)
start_time = BOPTEST_START_TIME;
warm_up    = 7*24*3600;
n_days     = 7;
n_steps    = round(n_days * 86400 / Ts);

%% 4. Deploy BOPTEST
handler = RequestHandler('http://127.0.0.1', 80);
tid     = handler.deploy_test('multizone_residential_hydronic');
handler.set_scenario(start_time, warm_up, BOPTEST_TUNC, BOPTEST_SUNC);
handler.set_step(Ts);
res = handler.advance(struct());

% See fetch_forecast_window.m for why a noisy forecast needs a per-step fetch.
USE_NOISY_FCAST = ~isempty(BOPTEST_TUNC) || ~isempty(BOPTEST_SUNC);

%% 5. Forecast prefetch
n_total = n_steps + p.N_p;

setp_sigs = cell(1, 2*nU);
for j = 1:nU
    setp_sigs{2*j-1} = sprintf('LowerSetp[%s]', ZONE_NAMES{j});
    setp_sigs{2*j  } = sprintf('UpperSetp[%s]', ZONE_NAMES{j});
end
dist_sigs = [{'TDryBul','HGloHor'}, SIG_CON(1:nZ), SIG_RAD(1:nZ)];

if USE_NOISY_FCAST
    Tamb_full = zeros(p.N_p, 1);      % refilled each step by fetch_forecast_window
    HGlo_full = zeros(p.N_p, 1);
    Qint_full = zeros(nZ, p.N_p);
    lb_mpc    = zeros(p.N_p, nU);
    ub_mpc    = zeros(p.N_p, nU);
    fprintf('[ekf_mpc] noisy fcast: per-step refetch (T_unc=%s, S_unc=%s)\n', ...
            BOPTEST_TUNC, BOPTEST_SUNC);
else
    % Schedule and disturbances, one bulk call each instead of one per step.
    % The QP sees the raw schedule; y_max in tracking_mpc_step is the only ceiling.
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
    Qint_full = zeros(nZ, n_total);
    for j = 1:nZ
        qc = fc_dist.(regexprep(SIG_CON{j},'[\[\]]','_'))(1:n_total);
        qr = fc_dist.(regexprep(SIG_RAD{j},'[\[\]]','_'))(1:n_total);
        Qint_full(j,:) = qc(:).' + qr(:).';
    end

    fprintf('LowerSetp sample: '); fprintf('%.1fC ', lb_mpc(1:4,1)-273.15); fprintf('\n');
    fprintf('UpperSetp sample: '); fprintf('%.1fC ', ub_mpc(1:4,1)-273.15); fprintf('\n\n');
end

%% 6. Initialise the augmented KF state from the first measurement
Tz_all = zeros(nZ_boptest,1);
for i = 1:nZ_boptest, Tz_all(i) = res.(SIG_TZ{i}); end
Tair_init = Tz_all(1:nZ);
z_aug     = [Tair_init; Tair_init + 1.5; zeros(p.nd_bias, 1)];   % [x; d_bias=0]
u_prev    = zeros(nU,1);

fprintf('KF init from measurement.  T_air = '); fprintf('%.1fC ', Tz_all-273.15); fprintf('\n\n');

%% 7. Logging buffers
log_Tz    = zeros(n_steps, nZ_boptest);
log_u     = zeros(n_steps, nU);
log_Qhea  = zeros(n_steps, nU);
log_ref   = zeros(n_steps, nU);
log_ub    = zeros(n_steps, nU);
log_Tamb  = zeros(n_steps, 1);
log_d     = zeros(n_steps, p.nD);
log_xhat  = zeros(n_steps, p.nz);
log_setpt = zeros(n_steps, nU);
log_innov = zeros(n_steps, nU);
log_d_hat = zeros(n_steps, nU);

fprintf('MPC running [tracking, offset-free]: %d zones | %g day(s) | Np=%d (%.0fh)\n', ...
        nU, n_days, p.N_p, p.N_p*Ts/3600);
fprintf('  KF: z=[T_z(%d); T_i(%d); d_bias(%d)]  nz=%d\n\n', nZ, nZ, p.nd_bias, p.nz);

%% 8. MPC loop
for i = 1:n_steps

    % 1. Measure
    Tz_meas = zeros(nYm,1);
    for zi = 1:nZ, Tz_meas(zi) = res.(SIG_TZ{zi}); end
    T_Hal_now = res.(SIG_TZ{end});

    % 2. Forecast window: a slice of the bulk arrays, or a fresh fetch under noise
    if USE_NOISY_FCAST
        [Tamb_full, HGlo_full, Qint_w, lb_mpc, ub_mpc] = fetch_forecast_window( ...
            handler, dist_sigs, setp_sigs, ZONE_NAMES, nU, p.N_p, Ts);
        Qint_full(1:nU, :) = Qint_w;
        i_loc = 1;              % the fetched window starts at the current step
        n_loc = p.N_p;
    else
        i_loc = i;              % the bulk arrays are indexed by absolute step
        n_loc = n_total;
    end

    % 3. Disturbance vector over the horizon
    dvec = zeros(p.nD * p.N_p, 1);
    for k = 1:p.N_p
        ii = min(i_loc + k - 1, n_loc);
        dvec((k-1)*p.nD + (1:p.nD)) = [Tamb_full(ii); HGlo_full(ii);
                                        Qint_full(:,ii); T_Hal_now];
    end
    dvec(~isfinite(dvec)) = 0;

    % 4. Comfort band, current step and N_p-step look-ahead
    lb_now_K = lb_mpc(min(i_loc, size(lb_mpc,1)), :)';
    ub_now_K = ub_mpc(min(i_loc, size(ub_mpc,1)), :)';
    w_target = p.w_occ * ones(nU, 1);

    lb_horizon = zeros(nYc * p.N_p, 1);
    for k = 1:p.N_p
        ii    = min(i_loc + k - 1, size(lb_mpc, 1));
        lb_k  = lb_mpc(ii,:)';
        occ_k = double(lb_k > p.occ_sp_threshold);
        lb_horizon((k-1)*nYc+1 : k*nYc) = occ_k .* w_target + (1-occ_k) .* lb_k;
    end

    % 5. MPC step: KF update and QP.
    %    The returned u_prev is discarded: step 7 replaces it with the heat
    %    BOPTEST actually delivered, which is what the KF must propagate.
    [v_cmd, u_apply, z_aug, ~, innov_kf] = tracking_mpc_step( ...
        Tz_meas, dvec, lb_horizon, z_aug, u_prev);

    % 6. One-step model prediction, bias-corrected, logged as the implied setpoint
    x_hat        = z_aug(1:nx);
    d_hat_kf     = z_aug(nx+1:nx+p.nd_bias);
    Tamb_now     = dvec(1);
    d_now        = dvec(1:p.nD);
    T_model_next = p.C_ctrl_m * (p.A * x_hat + p.Bu * u_apply + p.Bd * d_now) + d_hat_kf;
    T_sp_mpc     = max(lb_now_K, min(ub_now_K, T_model_next));

    % 7. Apply the valve commands directly; BOPTEST's inner PI stays disabled
    ov = struct('oveTSetPumBoi_activate', 1, 'oveTSetPumBoi_u', p.w_occ);
    for j = 1:nU
        zn = ZONE_NAMES{j};
        ov.(sprintf('conHea%s_oveActHea_activate',  zn)) = 1;
        ov.(sprintf('conHea%s_oveActHea_u',         zn)) = v_cmd(j);
        ov.(sprintf('conHea%s_oveTSetHea_activate', zn)) = 0;
    end
    res = handler.advance(ov);

    % 8. Log
    for zi = 1:nZ_boptest, log_Tz(i,zi)  = res.(SIG_TZ{zi}); end
    for j  = 1:nU,         log_Qhea(i,j) = res.(SIG_QH{j});  end

    u_prev = log_Qhea(i,:)';
    log_u(i,:)     = v_cmd';
    log_setpt(i,:) = T_sp_mpc';
    log_ref(i,:)   = lb_mpc(min(i_loc,end),:);
    log_ub(i,:)    = ub_mpc(min(i_loc,end),:);   % raw schedule, for plots and KPIs
    log_Tamb(i)    = Tamb_now;
    log_d(i,:)     = d_now';
    log_xhat(i,:)  = z_aug';
    log_d_hat(i,:) = d_hat_kf';
    log_innov(i,:) = innov_kf';

    if mod(i,16)==0
        fprintf('  step %4d/%d  Tamb=%+.1fC  T_air=[', i, n_steps, Tamb_now-273.15);
        fprintf('%.1f ', log_Tz(i,1:nU)-273.15);
        fprintf(']  v=['); fprintf('%.2f ', v_cmd');
        fprintf(']  d_hat=['); fprintf('%+.2f ', d_hat_kf'); fprintf(']K\n');
    end
end

%% 9. KPIs and figure
kpi = handler.get_kpi();
handler.stop_test(tid);

metrics = kpi_report('ekf_mpc', kpi, log_Tz, log_Qhea, log_ref, log_ub, ...
                     ZONE_NAMES, Ts, p.occ_sp_threshold);
comfort_plot('ekf_mpc', log_Tz, log_ref, log_ub, log_Tamb, kpi, ZONE_NAMES, Ts);

%% 10. Estimated bias
%  The peaks quoted in Section 4.3: how much of the mismatch the estimator had
%  to absorb, per zone, over the cold morning warm-ups.
fprintf('  Peak estimated bias d̂ [K]:  ');
for j = 1:nU
    fprintf('%s=%+.2f  ', ZONE_NAMES{j}, max(log_d_hat(:,j)));
end
fprintf('\n\n');

%% 11. Save
ss_model.A          = p.A;
ss_model.Bu         = p.Bu;
ss_model.Bd         = p.Bd;
ss_model.C_ctrl     = p.C_ctrl_m;
ss_model.C_full     = p.C_full_m;
ss_model.Ts         = Ts;
ss_model.Qdot_max   = p.Qdot_max_vec(:);
ss_model.ZONE_NAMES = ZONE_NAMES;
ss_model.nx  = nx;   ss_model.nU  = nU;   ss_model.nD  = p.nD;
ss_model.nYc = nYc;  ss_model.nYm = nYm;

sim.t       = (0:n_steps-1)' * Ts;
sim.Tz      = log_Tz;
sim.x_hat   = log_xhat;
sim.u_mpc   = log_u;
sim.setpt   = log_setpt;
sim.Qhea    = log_Qhea;
sim.ref     = log_ref;         % raw lower setpoint schedule [K]
sim.ref_ub  = log_ub;          % raw upper setpoint schedule [K]
sim.d       = log_d;
sim.Tamb    = log_Tamb;
sim.innov   = log_innov;
sim.d_hat   = log_d_hat;
sim.kpi     = kpi;
sim.metrics = metrics;
sim.BOPTEST_TUNC = BOPTEST_TUNC;
sim.BOPTEST_SUNC = BOPTEST_SUNC;

if isempty(BOPTEST_OUT_TAG)
    out_fname = 'mpc_simdata.mat';   % the deterministic offset-free run
else
    out_fname = sprintf('ekf_mpc%s_simdata.mat', BOPTEST_OUT_TAG);
end
save(fullfile('models', out_fname), 'ss_model', 'sim', 'kpi', 'metrics', 'Ts', ...
     'BOPTEST_TUNC', 'BOPTEST_SUNC', 'start_time');
fprintf('[ekf_mpc] Saved  models/%s\n', out_fname);
