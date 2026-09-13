function [v_cmd, u_apply, z_out, u_prev_out, innov_out, val, grad] = ...
        rl_mpc_step(Tz_meas, dvec, lb_horizon, z_in, u_prev_in, theta, a_fixed, use_bias)
%RL_MPC_STEP  One step of the parametrised MPC, with value and gradient.
%
%  Same KF, prediction model, heating curve and constraints as
%  tracking_mpc_step. What this function adds is the cost parametrisation θ and
%  the sensitivities the RL update of Chapter 5 needs. The BO of Chapter 6 uses
%  the same function and simply ignores the gradient.
%
%  Parametrised stage cost (Gros & Zanon, eq. 21a):
%    l_θ(x,u) = Σ_j q_j · M_occ · (y_j − w_j − Δȳ_j)²  +  r_u · ‖V − V̄‖²
%               + slack penalties g1_over, g1_under, g2
%
%  θ fields:
%    q_diag    nYc×2   tracking weights, col 1 occupied / col 2 unoccupied
%    r_u       2×1     input penalty,    entry 1 occupied / entry 2 unoccupied
%    w_offset  nYc×1   reference offsets Δȳ [K], occupied steps only
%    g1_over, g1_under, g2   slack penalties, held at their defaults
%  Only q_diag(:,1), r_u and w_offset are tuned: the eleven parameters of
%  Appendix B, shared by RL and BO.
%
%  Modes:
%    a_fixed = []     → V mode, val = V_θ(s)
%    a_fixed = a_k    → Q mode, val = Q_θ(s,a) with V_b(1:nU) fixed to a_k
%
%  Sensitivities: θ enters the cost only, so by the envelope theorem grad is
%  the partial derivative of the QP cost at z*. gradient_check.m verifies it
%  against finite differences.

%  The parameter file is constant over a run, so it is read once and kept.
%  The runners call `clear rl_mpc_step` before an episode to drop it.
persistent p
if isempty(p)
    p = load('models/tracking_mpc_params.mat');
end

nx  = p.nx;    nU     = p.nU;      nYc   = p.nYc;
nz  = p.nz;    nD     = p.nD;
N_p = p.N_p;   N_b    = p.N_b;     b_block = p.b_block;
nU_b = p.nU_b; nY_tot = p.nY_tot;  nZ_t  = p.nZ_t;

%% 1. KF update (15-state augmented)
%  use_bias toggles the offset-free augmentation:
%    true  → d̂ enters Y_free and the steady-state target (offset-free MPC)
%    false → bias entries held at 0, so θ is the only compensator of the
%            model-plant mismatch
if nargin < 8 || isempty(use_bias), use_bias = true; end

z_aug = z_in(1:nz);
if ~use_bias
    z_aug(nx+1:nz) = 0;          % drop bias before prediction
end
u_prev = u_prev_in(:);
d_dist = dvec(1:nD);
z_pred = p.A_aug * z_aug + p.Bu_aug * u_prev + p.Bd_aug * d_dist;
innov  = Tz_meas(:) - p.C_obs * z_pred;
z_hat  = z_pred + p.Lkf * innov;
if ~all(isfinite(z_hat)), z_hat = z_pred; end
if ~use_bias
    z_hat(nx+1:nz) = 0;          % do not carry bias to the next step
end
x_hat = z_hat(1:nx);
d_hat = z_hat(nx+1:nz);

%% 2. Radiator capacity from the heating curve (Wetter 2009)
Tamb_now = dvec(1);
qRel     = max(0, (p.HC_TRoo - Tamb_now) / (p.HC_TRoo - p.HC_TOut_nom));
T_supply = p.HC_TRoo + ((p.HC_TSup_nom + p.HC_TRet_nom)/2 - p.HC_TRoo) * qRel^(1/p.HC_m) ...
           + (p.HC_TSup_nom - p.HC_TRet_nom)/2 * qRel;
T_supply = max(p.HC_min, min(p.HC_max, T_supply));
dT_ref   = p.T_supply_ref - p.HC_TRoo;
Qdot_now = p.Qdot_max_vec(:) .* max(0.1, (T_supply - x_hat(1:nU)) / dT_ref);
ub_cap_b = Qdot_now ./ p.Qdot_max_vec(:);
ub_Z_b   = repmat(ub_cap_b, N_b, 1);

%% 3. Free response over the horizon
dvec_c = dvec(:); dvec_c(~isfinite(dvec_c)) = 0;
Y_free = p.Psi * x_hat + p.Theta_d * dvec_c + p.D_ext * d_hat;

%% 4. Reference and occupancy mask (same logic as tracking_mpc_step)
W     = zeros(nY_tot, 1);
M_occ = zeros(nY_tot, 1);
lb_h  = lb_horizon(:);
if numel(lb_h) < nY_tot, lb_h = repmat(lb_h(1:nYc), N_p, 1); end
for k = 1:N_p
    idx   = (k-1)*nYc+1 : k*nYc;
    occ_k = double(lb_h(idx) > p.occ_sp_threshold);
    W(idx)     = lb_h(idx) .* occ_k;
    M_occ(idx) = occ_k;
end
% Pre-heat window: track only if occupancy starts inside the horizon.
if M_occ(1) < 0.5
    first_occ_idx = 0;
    for k = 1:N_p
        if M_occ((k-1)*nYc + 1) > 0.5, first_occ_idx = k; break; end
    end
    fb = ceil(first_occ_idx / b_block);
    if fb == 0 || fb > N_b, M_occ(:) = 0; W(:) = 0; end
end
% Release the lower bound to y_min wherever the mask was cleared, so the zones
% may drift freely through a deep setback.
for k = 1:N_p
    idx = (k-1)*nYc+1 : k*nYc;
    if M_occ(idx(1)) < 0.5 && lb_h(idx(1)) > p.occ_sp_threshold
        lb_h(idx) = p.y_min * ones(nYc, 1);
    end
end

%% 5. Apply θ to the reference
%  Per-zone shift of the occupied reference: w_eff = w + Δȳ.
W_eff = W;
for k = 1:N_p
    idx = (k-1)*nYc+1 : k*nYc;
    W_eff(idx) = W(idx) + theta.w_offset(:) .* M_occ(idx);
end

%% 6. Steady-state input target
%  Built on W, not W_eff, so V̄ stays independent of θ and the gradient below
%  keeps its closed form.
if any(M_occ > 0.5)
    w_target = W(1:nYc);
    for k = 1:nYc
        if M_occ(k) < 0.5, w_target(k) = p.y_min; end
    end
    u_bar = p.K_u_d * d_dist + p.K_u_w * (w_target - d_hat);
else
    u_bar = zeros(nU, 1);
end
v_bar   = min(max(u_bar ./ p.Qdot_max_vec(:), 0), ub_cap_b);
V_bar_b = p.V_rep * v_bar;

%% 7. θ-dependent QP cost
%  The weights are binned by occupancy: per (zone, step) for the tracking
%  weight, and per block for the input penalty, where a block counts as
%  occupied if any of its steps is. That matches the pre-heat convention above.
M_mat_zk = reshape(M_occ, nYc, N_p);

block_occ = zeros(N_b, 1);
for b = 1:N_b
    step_lo = (b-1)*b_block + 1;
    step_hi = min(b*b_block, N_p);
    if step_hi >= step_lo
        block_occ(b) = double(any(any(M_mat_zk(:, step_lo:step_hi) > 0.5)));
    end
end

q_occ   = theta.q_diag(:, 1);
q_unocc = theta.q_diag(:, 2);
q_full  = zeros(nY_tot, 1);
for k = 1:N_p
    idx         = (k-1)*nYc + (1:nYc);
    occ_k       = M_mat_zk(:, k);
    q_full(idx) = occ_k .* q_occ + (1 - occ_k) .* q_unocc;
end
Q_tilde  = diag(q_full);
Q_masked = diag(M_occ) * Q_tilde;

r_u_b  = block_occ * theta.r_u(1) + (1 - block_occ) * theta.r_u(2);
R_Vb   = diag(kron(r_u_b, ones(nU, 1)));

% Factor 2 for the OSQP standard form min ½z'Pz + q'z.
H_vb = 2 * (p.Gamma_u_norm_b' * Q_tilde * p.Gamma_u_norm_b + R_Vb);
H_t  = blkdiag(H_vb, ...
               2*theta.g1_over  * eye(nY_tot), ...
               2*theta.g1_under * eye(nY_tot), ...
               2*theta.g2       * eye(nY_tot));
H_t  = (H_t + H_t') / 2;

e_free = Y_free - W_eff;
f_vb   = 2 * (p.Gamma_u_norm_b' * Q_masked * e_free - R_Vb * V_bar_b);
f_t    = [f_vb; zeros(nZ_t - nU_b, 1)];

%% 8. Constraints
b_com_upper = W_eff   - Y_free + p.Delta;
b_com_lower = Y_free  - lb_h   + p.Delta;
b_com_upper(M_occ < 0.5) = 1e6;
b_prt_upper = p.y_max - Y_free  + p.Delta;
b_prt_lower = Y_free  - p.y_min + p.Delta;
b_ineq_t    = [b_com_upper; b_com_lower; b_prt_upper; b_prt_lower];

n_ineq = 4 * nY_tot;
A_osqp = p.A_osqp_t;
l_osqp = [-inf(n_ineq,1); p.lb_z_t];
u_osqp = [b_ineq_t; ub_Z_b; inf(3*nY_tot,1)];

% Q mode: append the equality V_b(1:nU) = a_fixed. The action is clipped to the
% capacity bound first, so exploration noise cannot make the equality infeasible.
if ~isempty(a_fixed)
    a_clip = min(max(a_fixed(:), 0), ub_cap_b);
    A_osqp = [A_osqp; eye(nU), zeros(nU, nZ_t - nU)];
    l_osqp = [l_osqp; a_clip];
    u_osqp = [u_osqp; a_clip];
end

%% 9. Solve the QP
v_cmd      = zeros(nU,1);
u_apply    = zeros(nU,1);
z_out      = z_hat;
u_prev_out = zeros(nU,1);
innov_out  = innov;
val        = 0;
grad       = zeros_like_theta(theta);

if ~all(isfinite(f_t)) || ~all(isfinite(b_ineq_t)), return; end

z_sol = osqp_solve_wrapper(H_t, f_t, A_osqp, l_osqp, u_osqp);
if isempty(z_sol), return; end

V_b_sol = z_sol(1:nU_b);
dplus   = z_sol(nU_b + (1:nY_tot));
dminus  = z_sol(nU_b + nY_tot + (1:nY_tot));
sigma   = z_sol(nU_b + 2*nY_tot + (1:nY_tot));

v_opt      = min(max(V_b_sol(1:nU), 0), ub_cap_b);
u_apply    = p.Qdot_max_vec(:) .* v_opt;
v_cmd      = min(max(u_apply ./ Qdot_now, 0), 1);
u_prev_out = u_apply;

%% 10. Value and gradient at the QP optimum (envelope theorem, eq. 29)
%  With v = V_b*, the QP cost reads
%    C(v,θ) = v' G' Q̃ G v + 2 v' G' Q_masked e_free
%           + v' R v − 2 v' R V̄ + slack terms,
%  and V_θ(s) = C(v*,θ) + e_free' Q_masked e_free + V̄' R V̄. Since θ enters the
%  cost only, ∂V/∂θ is read straight off these terms:
%    ∂/∂q_occ_j   = Σ_{k occ}   (Y_v + e_free)²      ∂/∂r_occ   = Σ_{b occ}   ‖V_b − V̄‖²
%    ∂/∂q_unocc_j = Σ_{k unocc} Y_v²                 ∂/∂r_unocc = Σ_{b unocc} ‖V_b − V̄‖²
%    ∂/∂Δȳ_j      = −2 q_occ_j Σ_{k occ} (Y_v + e_free)
%    ∂/∂g         = ‖slack‖²
Y_v = p.Gamma_u_norm_b * V_b_sol;

dplus_sq   = dplus.'  * dplus;
dminus_sq  = dminus.' * dminus;
sigma_sq   = sigma.'  * sigma;
slack_quad = theta.g1_over * dplus_sq + theta.g1_under * dminus_sq + theta.g2 * sigma_sq;

val = Y_v.' * Q_tilde * Y_v ...
    + 2 * Y_v.' * Q_masked * e_free ...
    + e_free.' * Q_masked * e_free ...
    + V_b_sol.' * R_Vb * V_b_sol ...
    - 2 * V_b_sol.' * R_Vb * V_bar_b ...
    + V_bar_b.' * R_Vb * V_bar_b ...
    + slack_quad;

% Step × zone reshapes, so each gradient is a sum over the occupied steps.
Y_v_mat    = reshape(Y_v,    nYc, N_p).';
e_free_mat = reshape(e_free, nYc, N_p).';
M_mat      = reshape(M_occ,  nYc, N_p).';

grad.q_diag(:, 1) = sum(M_mat .* (Y_v_mat + e_free_mat).^2, 1).';
grad.q_diag(:, 2) = sum((1 - M_mat) .* Y_v_mat.^2, 1).';

V_b_mat    = reshape(V_b_sol, nU, N_b);
V_bar_mat  = reshape(V_bar_b, nU, N_b);
block_quad = sum((V_b_mat - V_bar_mat).^2, 1).';
grad.r_u(1) = sum(block_occ .* block_quad);
grad.r_u(2) = sum((1 - block_occ) .* block_quad);

grad.w_offset = -2 * q_occ .* sum(M_mat .* (Y_v_mat + e_free_mat), 1).';

grad.g1_over  = dplus_sq;
grad.g1_under = dminus_sq;
grad.g2       = sigma_sq;

end


%% Helpers
function g = zeros_like_theta(theta)
    g = struct('q_diag',   zeros(size(theta.q_diag)), ...
               'r_u',      zeros(size(theta.r_u)), ...
               'w_offset', zeros(size(theta.w_offset)), ...
               'g1_over',  0, ...
               'g1_under', 0, ...
               'g2',       0);
end
