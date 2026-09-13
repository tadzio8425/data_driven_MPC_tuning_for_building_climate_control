function [v_cmd, u_apply, z_out, u_prev_out, innov_out] = ...
        tracking_mpc_step(Tz_meas, dvec, lb_horizon, z_in, u_prev_in, use_bias)
%TRACKING_MPC_STEP  One step of the KF and the setpoint-tracking MPC.
%
%  Offset-free MPC with disturbance estimation (Magni & Scattolini §12.2.3).
%  Augmented state z_aug = [x(nx); d_bias(nd_bias)] with M = 0, F = I:
%    x(k+1) = A x(k) + Bu u(k) + Bd d_ext(k)
%    d(k+1) = d(k)
%    y(k)   = C x(k) + d(k)
%
%  The steady-state target (x̄, ū) is solved at every step from the current d̂
%  and the setpoint,
%    [I-A  -Bu] [x̄]   [Bd·d_now]
%    [ C    0 ] [ū] = [ w − d̂  ]
%  and the QP minimises
%    J = Σ_k ‖y_k − w_k‖²_Q + ‖u_k − ū‖²_R
%        + δ⁺' G1_over δ⁺ + δ⁻' G1_under δ⁻ + σ' G2 σ
%  over z = [V_b(nU_b); δ⁺(nY_tot); δ⁻(nY_tot); σ(nY_tot)].
%
%  use_bias (optional, default true):
%    true  → d̂ enters the target and Y_free, giving the integral action of the
%            offset-free MPC of Section 4.3
%    false → bias entries held at 0 throughout. This is the MPC of Section 4.2,
%            which carries the model-plant mismatch uncompensated.

%  The parameter file is constant over a run, so it is read once and kept.
%  The runners call `clear tracking_mpc_step` before an episode to drop it.
persistent p
if isempty(p)
    p = load('models/tracking_mpc_params.mat');
end

if nargin < 6 || isempty(use_bias)
    use_bias = true;
end

nx  = p.nx;    nU     = p.nU;      nYc  = p.nYc;
nz  = p.nz;    nD     = p.nD;
N_p = p.N_p;   N_b    = p.N_b;     b_block = p.b_block;
nU_b = p.nU_b; nY_tot = p.nY_tot;  nZ_t = p.nZ_t;

z_aug = z_in(1:nz);
if ~use_bias
    z_aug(nx+1:nz) = 0;         % drop bias before prediction
end
u_prev = u_prev_in(:);

%% 1. KF update (fixed DARE gain, 15-state augmented)
d_dist    = dvec(1:nD);
z_pred    = p.A_aug * z_aug + p.Bu_aug * u_prev + p.Bd_aug * d_dist;
innov     = Tz_meas(:) - p.C_obs * z_pred;
innov_out = innov;
z_hat     = z_pred + p.Lkf * innov;
if ~all(isfinite(z_hat)), z_hat = z_pred; end
if ~use_bias
    z_hat(nx+1:nz) = 0;         % do not carry bias to the next step
end

x_hat = z_hat(1:nx);
d_hat = z_hat(nx+1:nz);

%% 2. Radiator capacity from the heating curve (Wetter 2009)
Tamb_now = dvec(1);
qRel     = max(0, (p.HC_TRoo - Tamb_now) / (p.HC_TRoo - p.HC_TOut_nom));
T_supply = p.HC_TRoo ...
         + ((p.HC_TSup_nom + p.HC_TRet_nom)/2 - p.HC_TRoo) * qRel^(1/p.HC_m) ...
         + (p.HC_TSup_nom - p.HC_TRet_nom)/2 * qRel;
T_supply = max(p.HC_min, min(p.HC_max, T_supply));
dT_ref   = p.T_supply_ref - p.HC_TRoo;
Qdot_now = p.Qdot_max_vec(:) .* max(0.1, (T_supply - x_hat(1:nU)) / dT_ref);
ub_cap_b = Qdot_now ./ p.Qdot_max_vec(:);    % [nU x 1], <= 1
ub_Z_b   = repmat(ub_cap_b, N_b, 1);         % [nU_b x 1]

%% 3. Free response (the output includes the bias estimate)
dvec_c = dvec(:);  dvec_c(~isfinite(dvec_c)) = 0;
Y_free = p.Psi * x_hat + p.Theta_d * dvec_c + p.D_ext * d_hat;

%% 4. Reference and occupancy mask
W     = zeros(nY_tot, 1);
M_occ = zeros(nY_tot, 1);
lb_h  = lb_horizon(:);
if numel(lb_h) < nY_tot
    lb_h = repmat(lb_h(1:nYc), N_p, 1);
end
for k = 1:N_p
    idx   = (k-1)*nYc+1 : k*nYc;
    occ_k = double(lb_h(idx) > p.occ_sp_threshold);
    W(idx)     = lb_h(idx) .* occ_k;
    M_occ(idx) = occ_k;
end

% Pre-heat window: track only if occupancy starts within the blocked horizon,
% otherwise release the reference and let the zones drift through the setback.
if M_occ(1) < 0.5
    first_occ_idx = 0;
    for k = 1:N_p
        if M_occ((k-1)*nYc + 1) > 0.5
            first_occ_idx = k;
            break;
        end
    end
    fb = ceil(first_occ_idx / b_block);
    if fb == 0 || fb > N_b
        M_occ(:) = 0;
        W(:)     = 0;
    end
end

% Release the lower bound to y_min wherever the mask was cleared.
for k = 1:N_p
    idx = (k-1)*nYc+1 : k*nYc;
    if M_occ(idx(1)) < 0.5 && lb_h(idx(1)) > p.occ_sp_threshold
        lb_h(idx) = p.y_min * ones(nYc, 1);
    end
end

%% 5. Steady-state input target (absorbs d̂ and the current disturbance)
%  Uses the current-step setpoint, i.e. the constant-reference assumption.
if any(M_occ > 0.5)
    w_target = W(1:nYc);
    for k = 1:nYc
        if M_occ(k) < 0.5, w_target(k) = p.y_min; end   % unoccupied: protective bound
    end
    u_bar = p.K_u_d * d_dist + p.K_u_w * (w_target - d_hat);
else
    u_bar = zeros(nU, 1);   % nothing occupied in the horizon
end

v_bar   = min(max(u_bar ./ p.Qdot_max_vec(:), 0), ub_cap_b);
V_bar_b = p.V_rep * v_bar;                % [nU_b x 1]

%% 6. Time-varying QP terms
Q_masked = diag(M_occ) * p.Q_tilde;
e_free   = Y_free - W;

% Linear term: tracking 2·Γ'·Q·e_free, input cost −2·R_u·V̄_b.
f_vb = 2 * (p.Gamma_u_norm_b' * Q_masked * e_free - p.R_u * V_bar_b);
f_t  = [f_vb; zeros(nZ_t - nU_b, 1)];

% Constraint right-hand side, in absolute output form.
b_com_upper = W      - Y_free  + p.Delta;
b_com_lower = Y_free - lb_h    + p.Delta;
b_com_upper(M_occ < 0.5) = 1e6;
b_prt_upper = p.y_max - Y_free  + p.Delta;
b_prt_lower = Y_free  - p.y_min + p.Delta;
b_ineq_t    = [b_com_upper; b_com_lower; b_prt_upper; b_prt_lower];

%% 7. Solve the QP
%  A_osqp_t = [A_ineq_t; eye(nZ_t)]: V_b in [0, ub_Z_b], slacks >= 0.
u_apply = zeros(nU, 1);
if all(isfinite(f_t)) && all(isfinite(b_ineq_t))
    n_ineq = 4 * nY_tot;
    l_osqp = [-inf(n_ineq, 1);        % comfort and protection are one-sided
              p.lb_z_t];              % V_b >= 0, slacks >= 0
    u_osqp = [b_ineq_t;               % comfort and protection upper bounds
              ub_Z_b;                 % V_b <= ub_cap
              inf(3*nY_tot, 1)];      % slacks unbounded above

    z_sol = osqp_solve_wrapper(p.H_t, f_t, p.A_osqp_t, l_osqp, u_osqp);

    if ~isempty(z_sol)
        v_opt   = min(max(z_sol(1:nU), 0), ub_cap_b);
        u_apply = p.Qdot_max_vec(:) .* v_opt;
    end
end

%% 8. Valve command and state hand-over
v_cmd      = min(max(u_apply ./ Qdot_now, 0), 1);
z_out      = z_hat;
u_prev_out = u_apply;
