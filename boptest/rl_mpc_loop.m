%% rl_mpc_loop.m: RL tuning of the MPC cost parameters (Chapter 5)
%
%  On-policy Q-learning over one closed-loop episode, following the projected
%  update of Gros & Zanon (arXiv:1904.04152, eq. 32). It learns exactly the
%  eleven parameters the BO of Chapter 6 tunes, so the two methods are compared
%  on one and the same search space:
%
%      q_occ[5]     per-zone occupied tracking weights   ∈ [1e-2, 1e4]
%      r_u          one shared input penalty             ∈ [1e-4, 1e4]
%      w_offset[5]  reference offsets [K]                ∈ [-2, 2]
%
%  Every other entry of θ (q_unocc and the three slack penalties) is held at the
%  default the BO also fixes it at.
%
%  One iteration:
%      s_k → rl_mpc_step (V mode) → V_θ(s_k), ∇_θ, policy
%            ε-greedy: a_k = policy + noise on 10% of the steps
%            apply a_k to BOPTEST, observe Q_heat, T_z(s_{k+1}) and L_real
%      τ_k = L_real + γ·V_θ(s_{k+1}) − Q_θ(s_k, a_k)
%      θ   ← P_Φ( θ + α · clip(τ_k − baseline) · ∇_θ Q )
%
%  Three departures from the literal eq. 32, all needed for stability here and
%  all listed in Section 5.3:
%    * an EMA baseline on τ. V_θ is in tracking units while L_real is economic,
%      so the raw TD error carries a large constant offset. Without the baseline
%      the sign-definite q-gradient drives every weight to its lower bound.
%    * a small step size, the same for all three parameter groups.
%    * clipping of the centred τ, against the spikes produced when a constraint
%      becomes active or inactive between two steps.
%
%  BEFORE RUNNING:
%    1. startup
%    2. identification/ekf_sysid_2C     → models/rc_params_2C6z.mat
%    3. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat
%    4. BOPTEST server reachable at http://127.0.0.1:80
%
%  Saves  models/rl_mpc[_bias]_simdata.mat

clearvars -except use_bias BOPTEST_TUNC BOPTEST_SUNC BOPTEST_OUT_TAG ...
                  RL_SEED RL_START_TIME; clc;

% Optional hooks, set by the sweep drivers:
%   BOPTEST_TUNC     temperature forecast noise: 'medium' | 'high' | ''
%   BOPTEST_SUNC     the same options for solar irradiation
%   BOPTEST_OUT_TAG  suffix appended to the output filename, e.g. '_medium'
%   RL_SEED          seed of the ε-greedy exploration
%   RL_START_TIME    scenario start [s]; 0 is the week the thesis trains on
if ~exist('BOPTEST_TUNC','var')    || isempty(BOPTEST_TUNC),    BOPTEST_TUNC    = ''; end
if ~exist('BOPTEST_SUNC','var')    || isempty(BOPTEST_SUNC),    BOPTEST_SUNC    = ''; end
if ~exist('BOPTEST_OUT_TAG','var') || isempty(BOPTEST_OUT_TAG), BOPTEST_OUT_TAG = ''; end
if ~exist('RL_SEED','var')         || isempty(RL_SEED),         RL_SEED         = 47; end
if ~exist('RL_START_TIME','var')   || isempty(RL_START_TIME),   RL_START_TIME   = 0;  end

%% 1. Load the model and the precomputed QP matrices
load('models/rc_params_2C6z.mat');            % ZONE_NAMES, nZ, Ts
nZ_boptest = nZ;
p  = load('models/tracking_mpc_params.mat');
nU = p.nU;  nYc = p.nYc;
clear rl_mpc_step

%% 2. Initial θ
%  rl_mpc_step needs the whole struct, but only q_diag(:,1), r_u and w_offset
%  are learned. The rest holds the values the BO also keeps fixed.
q_init = diag(p.Q_tilde(1:nYc, 1:nYc));
theta.q_diag   = repmat(q_init, 1, 2);   % col 1 = q_occ (learned), col 2 = q_unocc
theta.r_u      = repmat(p.R_u,  2, 1);   % both entries kept equal, one shared value
theta.w_offset = zeros(nYc, 1);
theta.g1_over  = p.g1_over;
theta.g1_under = p.g1_under;
theta.g2       = p.g2;

fprintf('[rl_mpc] learning 11 parameters: q_occ[5], shared r_u, w_offset[5]\n');
fprintf('  q_occ init = %s\n', mat2str(theta.q_diag(:,1).', 4));
fprintf('  r_u init = %.3g   w_offset init = 0\n', theta.r_u(1));

%% 3. Learning hyperparameters (Table B.2)
%  use_bias = false puts the learning on the MPC of Section 4.2, where θ is the
%  only compensator of the mismatch; use_bias = true puts it on the offset-free
%  MPC of Section 4.3.
if ~exist('use_bias', 'var'), use_bias = false; end
gamma_rl  = 1.0;        % undiscounted, matching J_real
beta_base = 0.02;       % EMA rate of the TD baseline, about a 50-step memory
eps_expl  = 0.10;       % exploration probability
expl_std  = 0.05;       % uniform valve perturbation when exploring
tau_clip  = 1e3;        % TD-error clip

alpha.q_diag   = 1e-7;
alpha.r_u      = 1e-7;
alpha.w_offset = 1e-7;

w_E = 1.0; w_d = 1.0; w_o = 1.0;   % J_real = w_E·E + w_d·under + w_o·over

% The admissible box Φ, identical to the BO search space (Table B.1).
theta_lb.q_diag   =  1e-2 * ones(nYc,1);  theta_ub.q_diag   = 1e4 * ones(nYc,1);
theta_lb.r_u      =  1e-4;                theta_ub.r_u      = 1e4;
theta_lb.w_offset = -2.0  * ones(nYc,1);  theta_ub.w_offset = 2.0 * ones(nYc,1);

%% 4. Scenario and BOPTEST deployment
start_time = RL_START_TIME;
warm_up    = 7*24*3600;
n_days     = 7;
n_steps    = round(n_days * 86400 / Ts);   % 672 steps, one training episode

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

%% 5. Forecast prefetch
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
    fprintf('[rl_mpc] noisy fcast: per-step refetch (T_unc=%s, S_unc=%s)\n', ...
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

%% 6. Initial state and logging buffers
Tz_all = zeros(nZ_boptest,1);
for i = 1:nZ_boptest, Tz_all(i) = res.(SIG_TZ{i}); end
Tair_init = Tz_all(1:nU);
z_aug     = [Tair_init; Tair_init + 1.5; zeros(p.nd_bias, 1)];
u_prev    = zeros(nU,1);

log_Tz      = zeros(n_steps, nZ_boptest);
log_v       = zeros(n_steps, nU);
log_Qhea    = zeros(n_steps, nU);
log_ref     = zeros(n_steps, nU);
log_ub      = zeros(n_steps, nU);
log_Tamb    = zeros(n_steps, 1);
log_L       = zeros(n_steps, 1);
log_V       = zeros(n_steps, 1);
log_V_next  = zeros(n_steps, 1);
log_tau     = zeros(n_steps, 1);
log_tau_bar = zeros(n_steps, 1);
log_theta   = zeros(n_steps, 11);     % [q_occ(5); r_u; w_offset(5)]
log_explore = false(n_steps, 1);

fprintf('[rl_mpc] %d steps  (%g day(s), Ts=%.0fs)  use_bias=%d  seed=%d\n', ...
        n_steps, n_days, Ts, use_bias, RL_SEED);
fprintf('  centred τ, clip=%g, alpha=%.0e, eps=%.2f, gamma=%.2f\n\n', ...
        tau_clip, alpha.q_diag, eps_expl, gamma_rl);

%% 7. Training loop
%  A BOPTEST failure breaks the loop, and the script then reports on the partial
%  episode rather than losing the run.
n_done   = n_steps;
loop_err = '';
tau_bar  = 0;            % EMA baseline of the TD error
rng(RL_SEED);
t_run_start = tic;

for i = 1:n_steps

    % 1. Refresh the forecast window under BOPTEST noise
    if USE_NOISY_FCAST
        [Tamb_full, HGlo_full, Qint_w, lb_mpc, ub_mpc] = fetch_forecast_window( ...
            handler, dist_sigs, setp_sigs, ZONE_NAMES, nU, p.N_p, Ts);
        Qint_full(1:nU, :) = Qint_w;
        i_loc = 1;          % the fetched window starts at the current step
        n_loc = p.N_p;
    else
        i_loc = i;
        n_loc = n_total;
    end

    % 2. Assemble s_k from the current advance response
    [Tz_meas_k, dvec_k, lb_h_k] = build_step_inputs( ...
        res, lb_mpc, i_loc, p, Tamb_full, HGlo_full, Qint_full, n_loc, SIG_TZ, nU, nYc);

    % 3. V mode: policy, value and gradient at s_k
    [v_cmd, ~, z_aug_k, ~, ~, V_k, grad_V_k] = ...
        rl_mpc_step(Tz_meas_k, dvec_k, lb_h_k, z_aug, u_prev, theta, [], use_bias);

    % 4. ε-greedy exploration
    explore = rand() < eps_expl;
    if explore
        v_apply = min(max(v_cmd + expl_std * (2*rand(nU,1)-1), 0), 1);
    else
        v_apply = v_cmd;
    end

    % 5. Apply to BOPTEST
    ov = struct('oveTSetPumBoi_activate', 1, 'oveTSetPumBoi_u', p.w_occ);
    for j = 1:nU
        zn = ZONE_NAMES{j};
        ov.(sprintf('conHea%s_oveActHea_activate',  zn)) = 1;
        ov.(sprintf('conHea%s_oveActHea_u',         zn)) = v_apply(j);
        ov.(sprintf('conHea%s_oveTSetHea_activate', zn)) = 0;
    end
    try
        res_next = handler.advance(ov);
    catch ME
        n_done   = i - 1;
        loop_err = ME.message;
        fprintf('\n[rl_mpc] BOPTEST advance failed at step %d/%d: %s\n', i, n_steps, loop_err);
        fprintf('         Truncating the logs to %d completed steps.\n\n', n_done);
        break
    end

    % 6. Realised stage cost l_real(s_k, v_apply)
    Q_heat  = zeros(nU,1);
    Tz_next = zeros(nU,1);
    for j = 1:nU
        Q_heat(j)  = res_next.(SIG_QH{j});
        Tz_next(j) = res_next.(SIG_TZ{j});
    end
    lb_now = lb_mpc(min(i_loc,end),:).';
    ub_now = ub_mpc(min(i_loc,end),:).';
    occ    = double(lb_now > p.occ_sp_threshold);
    L_k    = w_E * sum(Q_heat) * Ts/3.6e6 ...
           + w_d * sum(occ .* max(0, lb_now - Tz_next)) * Ts/3600 ...
           + w_o * sum(occ .* max(0, Tz_next - ub_now)) * Ts/3600;

    % 7. Q mode at (s_k, v_apply). On a greedy step it coincides with V mode,
    %    so the extra QP solve is only paid when exploring.
    if explore
        [~, ~, ~, ~, ~, Q_k, grad_Q_k] = ...
            rl_mpc_step(Tz_meas_k, dvec_k, lb_h_k, z_aug, u_prev, theta, v_apply, use_bias);
    else
        Q_k      = V_k;
        grad_Q_k = grad_V_k;
    end

    % 8. V at s_{k+1}, the second half of the TD target
    [Tz_meas_n, dvec_n, lb_h_n] = build_step_inputs( ...
        res_next, lb_mpc, i_loc+1, p, Tamb_full, HGlo_full, Qint_full, n_loc, SIG_TZ, nU, nYc);
    [~, ~, ~, ~, ~, V_next, ~] = ...
        rl_mpc_step(Tz_meas_n, dvec_n, lb_h_n, z_aug_k, Q_heat, theta, [], use_bias);

    % 9. TD update on the eleven learned parameters, centred then clipped
    tau = L_k + gamma_rl * V_next - Q_k;
    if isfinite(tau)
        tau_bar = (1 - beta_base) * tau_bar + beta_base * tau;
        tau_upd = max(min(tau - tau_bar, tau_clip), -tau_clip);
        theta   = theta_step_sgd(theta, alpha, tau_upd, grad_Q_k);
        theta   = theta_project(theta, theta_lb, theta_ub);
    end

    % 10. Log and roll the state forward
    log_Tz(i,1:nU) = Tz_next.';
    log_Tz(i,end)  = res_next.(SIG_TZ{end});
    log_v(i,:)     = v_apply.';
    log_Qhea(i,:)  = Q_heat.';
    log_ref(i,:)   = lb_mpc(min(i_loc,end),:);
    log_ub(i,:)    = ub_mpc(min(i_loc,end),:);   % raw schedule, for plots and KPIs
    log_Tamb(i)    = dvec_k(1);
    log_L(i)       = L_k;
    log_V(i)       = V_k;
    log_V_next(i)  = V_next;
    log_tau(i)     = tau;
    log_tau_bar(i) = tau_bar;
    log_theta(i,:) = [theta.q_diag(:,1); theta.r_u(1); theta.w_offset].';
    log_explore(i) = explore;

    z_aug  = z_aug_k;     % the posterior at s_k becomes the next prior
    u_prev = Q_heat;      % the heat actually delivered at step k
    res    = res_next;

    if mod(i, 16) == 0
        fprintf('  step %4d/%d  L=%.3f  V=%.2f  V_+=%.2f  tau=%+.3f  q1=%.2g  e=%d\n', ...
                i, n_steps, L_k, V_k, V_next, tau, theta.q_diag(1), explore);
    end
end

%% 8. Trim the logs and collect the KPIs
log_Tz      = log_Tz(1:n_done, :);
log_v       = log_v(1:n_done, :);
log_Qhea    = log_Qhea(1:n_done, :);
log_ref     = log_ref(1:n_done, :);
log_ub      = log_ub(1:n_done, :);
log_Tamb    = log_Tamb(1:n_done);
log_L       = log_L(1:n_done);
log_V       = log_V(1:n_done);
log_V_next  = log_V_next(1:n_done);
log_tau     = log_tau(1:n_done);
log_tau_bar = log_tau_bar(1:n_done);
log_theta   = log_theta(1:n_done, :);
log_explore = log_explore(1:n_done);

% The KPIs are only meaningful over a complete episode.
if isempty(loop_err)
    kpi = handler.get_kpi();
else
    kpi = struct('tdis_tot', nan, 'ener_tot', nan, 'cost_tot', nan, 'emis_tot', nan);
end
try
    handler.stop_test(tid);
catch
end

t_elapsed = toc(t_run_start);
fprintf('\n  Completed %d / %d steps.  Sum L_real = %.3f  (mean %.4f per step)\n', ...
        n_done, n_steps, sum(log_L), mean(log_L));
fprintf('  Wall-clock = %.1f s  (%.2f s/step)   seed=%d\n', ...
        t_elapsed, t_elapsed/max(n_done,1), RL_SEED);
fprintf('  Learned θ:  q_occ=%s  r_u=%.4g  w_offset=%s\n', ...
        mat2str(theta.q_diag(:,1).',4), theta.r_u(1), mat2str(theta.w_offset.',3));

if use_bias, bias_tag = '_bias'; else, bias_tag = ''; end
ctrl_name = sprintf('rl_mpc%s', bias_tag);
metrics   = kpi_report(ctrl_name, kpi, log_Tz, log_Qhea, log_ref, log_ub, ...
                       ZONE_NAMES, Ts, p.occ_sp_threshold);

%% 9. Save before plotting, so a figure error cannot lose the episode
out_fname = sprintf('rl_mpc%s%s_simdata.mat', bias_tag, BOPTEST_OUT_TAG);
save(fullfile('models', out_fname), ...
     'log_Tz','log_v','log_Qhea','log_ref','log_ub','log_Tamb', ...
     'log_L','log_V','log_V_next','log_tau','log_tau_bar','log_theta','log_explore', ...
     'theta','kpi','metrics','Ts','use_bias','t_elapsed', ...
     'RL_SEED','RL_START_TIME','BOPTEST_TUNC','BOPTEST_SUNC');
fprintf('  Saved: models/%s\n', out_fname);

%% 10. Figure
%  The parameter-evolution figure of Chapter 5 is drawn separately, by
%  generate_controller_figures, straight from log_theta.
comfort_plot(ctrl_name, log_Tz, log_ref, log_ub, log_Tamb, kpi, ZONE_NAMES, Ts);


%% Helpers
function theta = theta_step_sgd(theta, alpha, tau, grad)
% Gradient step (eq. 32) on q_occ, the shared r_u and w_offset only.
    theta.q_diag(:,1) = theta.q_diag(:,1) + alpha.q_diag * tau * grad.q_diag(:,1);
    theta.r_u         = theta.r_u + alpha.r_u * tau * (grad.r_u(1) + grad.r_u(2)); % shared
    theta.w_offset    = theta.w_offset + alpha.w_offset * tau * grad.w_offset;
end

function theta = theta_project(theta, lb, ub)
% Projection onto the admissible box Φ (eq. 33).
    theta.q_diag(:,1) = min(max(theta.q_diag(:,1), lb.q_diag),   ub.q_diag);
    theta.r_u         = min(max(theta.r_u,         lb.r_u),      ub.r_u);
    theta.w_offset    = min(max(theta.w_offset,    lb.w_offset), ub.w_offset);
end
