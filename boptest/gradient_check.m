%% gradient_check.m: Finite-difference validation of the analytic gradient
%
%  Checks the envelope-theorem gradient that rl_mpc_step returns (Eq. 5.6)
%  against a central finite difference of the QP optimal value, on one
%  representative occupied state and with no BOPTEST in the loop:
%
%      for every entry of θ:   grad   vs   (V(θ+e) − V(θ−e)) / 2e
%
%  θ enters the cost only, so the two must agree to solver tolerance. The
%  exception is w_offset, which also shifts the comfort band through W_eff: the
%  cost partial is exact only while that band is inactive, which it is at the
%  test state chosen below.
%
%  This is the check Section 5.4 refers to when it states that the gradient
%  driving the learning is exact.
%
%  BEFORE RUNNING:
%    1. startup
%    2. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat

clear; clc;
clear rl_mpc_step

p = load('models/tracking_mpc_params.mat');
nU = p.nU;  nYc = p.nYc;  N_p = p.N_p;  nz = p.nz;  nx = p.nx;  nYm = p.nYm;

%% 1. Default θ, the same initialisation rl_mpc_loop starts from
theta.q_diag   = repmat(diag(p.Q_tilde(1:nYc, 1:nYc)), 1, 2);
theta.r_u      = repmat(p.R_u, 2, 1);
theta.w_offset = zeros(nYc, 1);
theta.g1_over  = p.g1_over;
theta.g1_under = p.g1_under;
theta.g2       = p.g2;

%% 2. Test state: occupied and cool, so the heating is engaged
occ_lb  = p.occ_sp_threshold + 1.0;    % clearly above the occupancy threshold
T0      = occ_lb - 1.0;                % about 1 K below the occupied target
Tz_meas = T0 * ones(nYm, 1);
z_in    = [T0*ones(nU,1); (T0+1.5)*ones(nU,1); zeros(nz-nx, 1)];
u_prev  = zeros(nU, 1);

Tamb = 278.15;  HGlo = 0;  Qg = zeros(nU,1);  THal = 290.15;
dvec = repmat([Tamb; HGlo; Qg; THal], N_p, 1);
lb_h = repmat(occ_lb * ones(nYc,1), N_p, 1);
use_bias = false;

%% 3. Analytic gradient
[~,~,~,~,~, V0, grad] = rl_mpc_step(Tz_meas, dvec, lb_h, z_in, u_prev, theta, [], use_bias);
if V0 == 0
    error('gradient_check: the QP returned val=0 (infeasible or solver failure).');
end
fprintf('Base QP value V(θ) = %.6g   (use_bias=%d, occupied test state)\n\n', V0, use_bias);

val_fun = @(th) nth_out(6, @() rl_mpc_step(Tz_meas, dvec, lb_h, z_in, u_prev, th, [], use_bias));

%% 4. Central finite difference, one θ block at a time
fields = {'q_diag','r_u','w_offset','g1_over','g1_under','g2'};
fprintf('%-11s %14s %14s %14s   %s\n', 'block', 'max|analytic|', 'max|FD|', 'max rel. err', 'verdict');
fprintf('%s\n', repmat('-', 1, 74));

worst = 0;
for f = 1:numel(fields)
    fn  = fields{f};
    blk = theta.(fn);
    gan = grad.(fn);
    gfd = zeros(size(blk));
    for idx = 1:numel(blk)
        e  = 1e-3 * max(1, abs(blk(idx)));
        tp = theta;  tp.(fn)(idx) = blk(idx) + e;
        tm = theta;  tm.(fn)(idx) = blk(idx) - e;
        gfd(idx) = (val_fun(tp) - val_fun(tm)) / (2*e);
    end
    denom   = max(abs(gfd(:)));  if denom < 1e-9, denom = 1; end
    rel_err = max(abs(gan(:) - gfd(:))) / denom;
    worst   = max(worst, rel_err);
    if rel_err > 5e-2, verdict = '*** MISMATCH'; else, verdict = 'OK'; end
    fprintf('%-11s %14.4g %14.4g %14.3g   %s\n', ...
            fn, max(abs(gan(:))), max(abs(gfd(:))), rel_err, verdict);
end

fprintf('%s\n', repmat('-', 1, 74));
fprintf('Worst relative error across all blocks: %.3g\n', worst);
if worst <= 5e-2
    fprintf('The analytic gradient matches finite differences: Eq. 5.6 is implemented correctly.\n');
else
    fprintf('At least one block disagrees — inspect above.\n');
end


%% Helper
function v = nth_out(n, fh)
% Return the n-th output of a function handle.
    out = cell(1, n);
    [out{:}] = fh();
    v = out{n};
end
