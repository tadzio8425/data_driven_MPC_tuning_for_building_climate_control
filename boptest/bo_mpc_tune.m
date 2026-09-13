%% bo_mpc_tune.m: Bayesian optimisation of the MPC cost parameters (Chapter 6)
%
%  Minimises the realised closed-loop cost J_real directly, treating one
%  seven-day BOPTEST episode as a single black-box evaluation. The search never
%  reads the model or its gradient: the only thing it sees is what the plant
%  actually cost. MATLAB's bayesopt (Statistics and Machine Learning Toolbox)
%  drives it, and every evaluation is a call to run_mpc_episode.
%
%  Optimised parameters, the same eleven the RL of Chapter 5 learns:
%    q_Liv q_Ro1 q_Ro2 q_Ro3 q_Bth   occupied tracking weights   [log scale]
%    r_u                             shared input penalty        [log scale]
%    w_Liv w_Ro1 w_Ro2 w_Ro3 w_Bth   reference offsets [K]       [linear]
%  Bounds in Table B.1. Everything else in θ stays at the default of
%  tracking_mpc_params.mat.
%
%  Set before running:
%    use_bias = false   the MPC of Section 4.2, which is the search reported in
%                       Chapter 6 (the default)
%    use_bias = true    the same search on the offset-free MPC, the control
%                       experiment of Section 6.4
%
%  BEFORE RUNNING:
%    1. startup
%    2. identification/ekf_sysid_2C     → models/rc_params_2C6z.mat
%    3. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat
%    4. BOPTEST server reachable at http://127.0.0.1:80
%
%  Saves  models/bo_mpc[_bias]_tune_data.mat
%    bo_results   the full bayesopt results object, including the trace that
%                 generate_controller_figures draws as Figure 6.1
%    theta_best   θ struct at the best observed point
%    x_best       the corresponding 1x1 table of BO variables
%    L_best       the lowest J_real found
%
%  Twenty evaluations, roughly five hours, because each one is a complete
%  seven-day closed-loop run. Afterwards run bo_mpc_benchmark to evaluate
%  θ_best on a fresh episode.
%
%  Reference: Lu, Kumar & Zavala, "MPC controller tuning using Bayesian
%  optimization techniques", arXiv:2009.14175, 2020.

if ~exist('use_bias', 'var'), use_bias = false; end
clearvars -except use_bias BOPTEST_TUNC BOPTEST_SUNC BOPTEST_OUT_TAG; clc;

% Optional forecast-uncertainty hooks, set by the sweep drivers:
%   BOPTEST_TUNC     temperature forecast noise: 'medium' | 'high' | ''
%   BOPTEST_SUNC     the same options for solar irradiation
%   BOPTEST_OUT_TAG  suffix appended to the tune-data filename
if ~exist('BOPTEST_TUNC','var')    || isempty(BOPTEST_TUNC),    BOPTEST_TUNC    = ''; end
if ~exist('BOPTEST_SUNC','var')    || isempty(BOPTEST_SUNC),    BOPTEST_SUNC    = ''; end
if ~exist('BOPTEST_OUT_TAG','var') || isempty(BOPTEST_OUT_TAG), BOPTEST_OUT_TAG = ''; end

%% 1. Defaults, which are also the warm start
p   = load('models/tracking_mpc_params.mat');
nU  = p.nU;
nYc = p.nYc;
q_default = diag(p.Q_tilde(1:nYc, 1:nYc));

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  Bayesian optimisation of the MPC cost parameters\n');
fprintf('  use_bias = %d   |   nYc = %d   |   Ts = %d s\n', use_bias, nYc, p.Ts);
fprintf('══════════════════════════════════════════════════════\n\n');

%% 2. Search space (Table B.1)
%  The weights and the input penalty span several orders of magnitude, so they
%  are searched on a logarithmic scale; the offsets are searched linearly.
q_lo = 1e-2;  q_hi = 1e4;
r_lo = 1e-4;  r_hi = 1e4;
w_lo = -2;    w_hi = 2;

vars = [
    optimizableVariable('q_Liv', [q_lo, q_hi], 'Transform', 'log')
    optimizableVariable('q_Ro1', [q_lo, q_hi], 'Transform', 'log')
    optimizableVariable('q_Ro2', [q_lo, q_hi], 'Transform', 'log')
    optimizableVariable('q_Ro3', [q_lo, q_hi], 'Transform', 'log')
    optimizableVariable('q_Bth', [q_lo, q_hi], 'Transform', 'log')
    optimizableVariable('r_u',   [r_lo, r_hi], 'Transform', 'log')
    optimizableVariable('w_Liv', [w_lo, w_hi])
    optimizableVariable('w_Ro1', [w_lo, w_hi])
    optimizableVariable('w_Ro2', [w_lo, w_hi])
    optimizableVariable('w_Ro3', [w_lo, w_hi])
    optimizableVariable('w_Bth', [w_lo, w_hi])
];

fprintf('Variable ranges:\n');
fprintf('  q_*  : [%.2g, %.2g]  (log)\n',    q_lo, q_hi);
fprintf('  r_u  : [%.2g, %.2g]  (log)\n',    r_lo, r_hi);
fprintf('  w_*  : [%.2g, %.2g]  (linear, K)\n\n', w_lo, w_hi);

%% 3. Warm start from the hand-chosen defaults of Table 4.2
InitialX = table(q_default(1), q_default(2), q_default(3), q_default(4), q_default(5), ...
                 p.R_u, 0, 0, 0, 0, 0, ...
    'VariableNames', {'q_Liv','q_Ro1','q_Ro2','q_Ro3','q_Bth', ...
                      'r_u','w_Liv','w_Ro1','w_Ro2','w_Ro3','w_Bth'});

%% 4. Run the search
n_evals  = 20;
rng_seed = 42;

% Under forecast uncertainty the objective is stochastic: the same θ returns a
% different J_real on every run. Declaring it deterministic would make the GP
% fit that noise, so the flag follows the noise setting.
is_deterministic = isempty(BOPTEST_TUNC) && isempty(BOPTEST_SUNC);

rng(rng_seed);
fprintf('Starting bayesopt:  MaxObjectiveEvaluations=%d   seed=%d\n\n', n_evals, rng_seed);

bo_results = bayesopt(@(x) run_mpc_episode(x, use_bias, BOPTEST_TUNC, BOPTEST_SUNC), vars, ...
    'MaxObjectiveEvaluations', n_evals, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'IsObjectiveDeterministic', is_deterministic, ...
    'InitialX',    InitialX, ...
    'Verbose',     1, ...
    'UseParallel', false);

%% 5. Best point
x_best = bo_results.XAtMinObjective;
L_best = bo_results.MinObjective;

theta_best.q_diag   = [[x_best.q_Liv; x_best.q_Ro1; x_best.q_Ro2; ...
                        x_best.q_Ro3; x_best.q_Bth], q_default];
theta_best.r_u      = [x_best.r_u; x_best.r_u];
theta_best.w_offset = [x_best.w_Liv; x_best.w_Ro1; x_best.w_Ro2; ...
                       x_best.w_Ro3; x_best.w_Bth];
theta_best.g1_over  = p.g1_over;
theta_best.g1_under = p.g1_under;
theta_best.g2       = p.g2;

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  BO complete  |  %d evaluations  |  J_real,best = %.4f\n', ...
        bo_results.NumObjectiveEvaluations, L_best);
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  q   = [%8.3g %8.3g %8.3g %8.3g %8.3g]\n', ...
        x_best.q_Liv, x_best.q_Ro1, x_best.q_Ro2, x_best.q_Ro3, x_best.q_Bth);
fprintf('  r_u = %.4g\n', x_best.r_u);
fprintf('  w   = [%+.3f %+.3f %+.3f %+.3f %+.3f]  K\n', ...
        x_best.w_Liv, x_best.w_Ro1, x_best.w_Ro2, x_best.w_Ro3, x_best.w_Bth);

%% 6. Save
if use_bias, bias_tag = '_bias'; else, bias_tag = ''; end
out_fname = sprintf('bo_mpc%s%s_tune_data.mat', bias_tag, BOPTEST_OUT_TAG);
save(fullfile('models', out_fname), ...
     'bo_results', 'theta_best', 'x_best', 'L_best', 'use_bias', ...
     'BOPTEST_TUNC', 'BOPTEST_SUNC');

fprintf('\n  Saved: models/%s\n', out_fname);
fprintf('  Next:  bo_mpc_benchmark, to evaluate θ_best on a fresh episode.\n\n');
