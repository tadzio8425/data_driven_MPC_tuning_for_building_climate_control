%% tracking_mpc_params_build.m: Precompute the constant MPC matrices
%
%  Assembles, once and offline, everything in the QP of Appendix A that does not
%  change during a run: the condensed prediction maps, the Hessian H, the
%  constraint matrix G and the steady-state target solver. The step functions
%  then only rebuild the linear term and the bounds.
%
%  Offset-free setpoint-tracking MPC with disturbance estimation
%  (Magni & Scattolini §12.2.3), built on the identified 2C model.
%
%  Augmented system (M = 0, F = I, additive output bias per zone):
%      x(k+1) = A x(k) + Bu u(k) + Bd d_ext(k)
%      d(k+1) = d(k)
%      y(k)   = C x(k) + d(k)
%
%  Decision variable, with move-blocking and absolute valve positions:
%      z = [V_b(nU_b);  δ⁺(nY_tot);  δ⁻(nY_tot);  σ(nY_tot)]
%  where V_b ∈ [0,1] is the normalised valve opening, u = P_nom · V.
%
%  BEFORE RUNNING:
%    1. startup
%    2. identification/ekf_sysid_2C   → models/rc_params_2C6z.mat
%
%  Saves  models/tracking_mpc_params.mat

clear; clc;

%% 1. Load the identified model
load('models/rc_params_2C6z.mat');

nZ  = 5;             % Hal is a disturbance, not a state of the MPC model
nU  = nU_ss;         %  5  controlled zones
nD  = nD_ss;         %  8  [Tamb; HGlo; Qint(5); T_Hal]
nx  = nx_ss;         % 10  [T_air(5); T_wall(5)]
nYm = nZ;            %  5  KF measurements
nYc = nU;            %  5  MPC outputs
nd_bias = nYc;       %  5  additive output-bias states (F = I)
nz  = nx + nd_bias;  % 15  augmented KF state

A        = Ad;
Bu       = Bud;
Bd       = Bdd;
C_full_m = C_full;
C_ctrl_m = C_ctrl;

% Nominal radiator power per zone, taken as the peak heating power the baseline
% PI delivers there (Table 4.3).
Qdot_max_vec = [1800, 1000, 900, 1100, 800]';   % [W]

%% 2. Cost parameters (Table 4.2)
% N_p = 24 at Ts = 900 s gives a 6 h prediction horizon.
N_p     = 24;
b_block = 1;
N_b     = N_p / b_block;

gamma    = 100;
R_u      = 1;
% The comfort slacks are asymmetric on purpose. The system is heating-only, so
% undershoot during occupancy is the violation that actually happens: penalising
% it three times harder makes the QP prefer a slight overshoot to a cold zone.
g1_over  = 10;
g1_under = 30;
g2       = 100;
Delta    = 0.5;               % band tolerance [K]
y_min    = 15 + 273.15;       % protection band, fixed
y_max    = 25 + 273.15;

w_occ            = 19 + 273.15;   % occupied comfort reference
occ_sp_threshold = 17 + 273.15;   % a scheduled lb above this means "occupied"

%% 3. Heating curve (BOPTEST, Wetter 2009)
HC_TRoo      = w_occ;
HC_TSup_nom  = 60 + 273.15;
HC_TRet_nom  = 50 + 273.15;
HC_TOut_nom  = -5 + 273.15;
HC_m         = 1.3;            % heat-emission exponent
HC_min       = 35 + 273.15;    % supply-temperature clip
HC_max       = 80 + 273.15;
T_supply_ref = 60 + 273.15;

%% 4. Kalman filter
q_Tair  = (0.05)^2;
q_Twall = (0.10)^2;
q_bias  = (0.4)^2;      % random-walk variance of the bias states: the single
                        % knob that sets how fast the integral action responds

Q_kf     = blkdiag(q_Tair * eye(nZ), q_Twall * eye(nZ));
R_kf     = (0.01)^2 * eye(nYm);

A_aug    = blkdiag(A, eye(nd_bias));       % [A 0; 0 I]   15x15
Bu_aug   = [Bu; zeros(nd_bias, nU)];       % [Bu; 0]      15x5
Bd_aug   = [Bd; zeros(nd_bias, nD)];       % [Bd; 0]      15xnD
C_obs    = [C_ctrl_m, eye(nd_bias)];       % [C  I]        5x15
Q_kf_aug = blkdiag(Q_kf, q_bias * eye(nd_bias));

% Observability of the augmented pair, PBH test at z = 1 (Theorem 12.2.1).
obs_rank = rank([A - eye(nx), zeros(nx, nd_bias); C_ctrl_m, eye(nd_bias)]);
assert(obs_rank == nz, 'Augmented system not observable (Theorem 12.2.1)');
fprintf('Observability check passed (rank = %d = nz)\n', obs_rank);

[P_ss, ~, ~] = dare(A_aug', C_obs', Q_kf_aug, R_kf);
Lkf = P_ss * C_obs' / (C_obs * P_ss * C_obs' + R_kf);

cl_eigs = eig((eye(nz) - Lkf*C_obs) * A_aug);
assert(all(abs(cl_eigs) < 1), 'KF is unstable — check the identified model.');
fprintf('DARE solved.  nz=%d  ||Lkf||_F = %.4f   max|cl_eig| = %.6f\n', ...
        nz, norm(Lkf,'fro'), max(abs(cl_eigs)));

%% 5. Steady-state target solver
%  T_ss [x̄; ū] = [Bd·d_now; w − d̂].  Only the rows that yield ū are needed at
%  run time, so they are precomputed here and the step function evaluates
%  ū = K_u_d·d_now + K_u_w·(w − d̂).
T_ss = [eye(nx) - A, -Bu; C_ctrl_m, zeros(nYc, nU)];
assert(rank(T_ss) == nx + nU, 'T_ss singular — transmission zero at z = 1');

N_d   = T_ss \ [Bd;             zeros(nYc, nD)];
N_w   = T_ss \ [zeros(nx, nYc); eye(nYc)      ];
K_u_d = N_d(nx+1:end, :);   % nU x nD
K_u_w = N_w(nx+1:end, :);   % nU x nYc

%% 6. Condensed prediction (Appendix A.1)
%  Y_free = Psi·x̂ + Theta_d·D + D_ext·d̂,  and the forced response is Gamma_u·U.
Psi     = zeros(nYc*N_p, nx);
Gamma_u = zeros(nYc*N_p, nU*N_p);
Theta_d = zeros(nYc*N_p, nD*N_p);
Ak = eye(nx);
for k = 1:N_p
    Ak   = A * Ak;
    rows = (k-1)*nYc+1 : k*nYc;
    Psi(rows, :) = C_ctrl_m * Ak;
    Akj = eye(nx);
    for j = 1:k
        Gamma_u(rows, (j-1)*nU+1:j*nU) = C_ctrl_m * Akj * Bu;
        Theta_d(rows, (j-1)*nD+1:j*nD) = C_ctrl_m * Akj * Bd;
        Akj = A * Akj;
    end
end
D_ext = kron(ones(N_p, 1), eye(nYc));   % broadcasts the constant bias estimate

%% 7. Move-blocking and the valve normalisation
nU_tot = nU * N_p;
nY_tot = nYc * N_p;
nU_b   = nU * N_b;
nZ_t   = nU_b + 3*nY_tot;

W_u            = kron(eye(N_p), diag(Qdot_max_vec));
Gamma_u_norm   = Gamma_u * W_u;                    % input expressed as valve opening
M_block        = kron(eye(N_b), kron(ones(b_block,1), eye(nU)));
Gamma_u_norm_b = Gamma_u_norm * M_block;           % V_b → ΔY
V_rep          = kron(ones(N_b, 1), eye(nU));      % ū (nU) → V̄_b (nU_b)

%% 8. Tracking weight from the input sensitivities (Eq. 4.9)
%  Q is not chosen by hand. Each zone is weighted by the inverse square of the
%  temperature rise one step of its own radiator produces, so that a kelvin of
%  error costs the same everywhere regardless of how strong that radiator is.
CBu_diag    = diag(C_ctrl_m * Bu);
assert(all(CBu_diag > 0), 'Non-positive C·Bu diagonal — check the model matrices.');
sensitivity = CBu_diag .* Qdot_max_vec;            % eta_i [K]

F       = diag(1 ./ sensitivity.^2);
Q_zone  = gamma * R_u * F;
Q_tilde = kron(eye(N_p), Q_zone);

ZN_MPC  = ZONE_NAMES(1:nU);
F_ratio = diag(F) / min(diag(F));
fprintf('\n  ──────── Per-zone tracking weight ────────\n');
fprintf('  %-18s ', 'Zone:');            fprintf('%10s ',   ZN_MPC{:});      fprintf('\n');
fprintf('  %-18s ', 'P_nom [W]:');       fprintf('%10.0f ', Qdot_max_vec');  fprintf('\n');
fprintf('  %-18s ', 'sensitivity [K]:'); fprintf('%10.3f ', sensitivity');   fprintf('\n');
fprintf('  %-18s ', 'Q_i / min(Q_i):');  fprintf('%10.2f ', F_ratio');       fprintf('\n');
[~, jmax] = max(F_ratio);
[~, jmin] = min(F_ratio);
fprintf('  → the QP weights %s %.1fx above %s, so that zone leads the tracking.\n\n', ...
        ZN_MPC{jmax}, F_ratio(jmax), ZN_MPC{jmin});

%% 9. Hessian (Appendix A.2)
%  Constant, because the target V̄_b enters only through the linear term.
H_vb = 2 * (Gamma_u_norm_b' * Q_tilde * Gamma_u_norm_b + R_u * eye(nU_b));
H_t  = blkdiag(H_vb, 2*g1_over*eye(nY_tot), 2*g1_under*eye(nY_tot), 2*g2*eye(nY_tot));
H_t  = (H_t + H_t') / 2;

%% 10. Constraint matrix (Appendix A.3)
%  Four band blocks on the absolute output Y = Y_free + Γ·V_b, then the box on z.
O = zeros(nY_tot, nY_tot);
A_ineq_t = [ Gamma_u_norm_b, -eye(nY_tot),  O,            O            ;   % comfort, from above
            -Gamma_u_norm_b,  O,           -eye(nY_tot),  O            ;   % comfort, from below
             Gamma_u_norm_b,  O,            O,           -eye(nY_tot)  ;   % protection, from above
            -Gamma_u_norm_b,  O,            O,           -eye(nY_tot) ];   % protection, from below

lb_z_t   = zeros(nZ_t, 1);              % V_b >= 0 and all slacks >= 0
A_osqp_t = [A_ineq_t; eye(nZ_t)];       % the full OSQP constraint matrix

%% 11. Save
z_aug_init = [(19 + 273.15) * ones(nx, 1); zeros(nd_bias, 1)];

save('models/tracking_mpc_params.mat', ...
    'A_aug','Bu_aug','Bd_aug','C_obs','Lkf', ...
    'A','Bu','Bd','C_ctrl_m','C_full_m', ...
    'Gamma_u_norm_b','Psi','Theta_d','D_ext','V_rep','M_block', ...
    'H_t','A_ineq_t','Q_tilde','lb_z_t','A_osqp_t', ...
    'K_u_d','K_u_w','T_ss', ...
    'Qdot_max_vec', ...
    'HC_TRoo','HC_TSup_nom','HC_TRet_nom','HC_TOut_nom','HC_m','HC_min','HC_max','T_supply_ref', ...
    'N_p','nU','nD','nYc','nYm','nx','nZ','nz','nd_bias', ...
    'nU_tot','nY_tot','nU_b','nZ_t','b_block','N_b', ...
    'gamma','R_u','g1_over','g1_under','g2','Delta','y_min','y_max', ...
    'w_occ','occ_sp_threshold', ...
    'z_aug_init','Ts');

fprintf('Saved → models/tracking_mpc_params.mat\n');
fprintf('  Horizon: Np=%d  Ts=%.0fs  |  gamma=%.1f  R_u=%g\n', N_p, Ts, gamma, R_u);
fprintf('  g1_over=%g  g1_under=%g  g2=%g  Delta=%.2fK  protection=[%.0f,%.0f]C\n', ...
        g1_over, g1_under, g2, Delta, y_min-273.15, y_max-273.15);
fprintf('  Decision var: z=[V_b(%d); d+(%d); d-(%d); sigma(%d)]  n_z=%d\n', ...
        nU_b, nY_tot, nY_tot, nY_tot, nZ_t);

clear tracking_mpc_step rl_mpc_step
