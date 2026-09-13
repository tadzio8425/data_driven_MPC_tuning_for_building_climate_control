%% run_baselines.m: Menu launcher for every controller, sweep and figure
%
%  Each controller runs the same seven-day BOPTEST scenario and saves its own
%  .mat under models/, with identical log fields, so any two runs can be
%  overlaid directly. Every entry can also be called on its own from the command
%  line, which is what the sweeps do.
%
%  BEFORE RUNNING:
%    1. startup
%    2. identification/ekf_sysid_2C     → models/rc_params_2C6z.mat
%    3. mpc/tracking_mpc_params_build   → models/tracking_mpc_params.mat
%
%  Options 1 to 17 need a BOPTEST server at http://127.0.0.1:80 serving the
%  multizone_residential_hydronic test case. Options 8, 10 and 12 also need the
%  Statistics and Machine Learning Toolbox, for bayesopt.
%
%  Notes
%    * Train before evaluating: 4 before 5, 6 before 7, 8 before 9, 10 before 11.
%      Each benchmark reads what its tuner saved.
%    * Options 1, 2, 3, 7 and 9 are the five controllers the thesis reports, and
%      they are the clean column of every Chapter 7 table. Run them first.
%    * Options 10 and 11 are the control experiment of Section 6.4: the same
%      twenty-episode search on top of the offset-free MPC.
%    * Nothing is retuned in Chapter 7. The sweeps deploy the nominal θ on a
%      degraded forecast (12) or on a week nobody tuned on (14-17), and every
%      out-of-sample output carries a _wk<N> tag, so no reported run is touched.
%    * Options 13, 18, 19 and 20 read .mat files only and need no server.

fprintf('\n──────────────────────────────────────────────────────────\n');
fprintf('  Pick what to run:\n');
fprintf('   ── The five reported controllers (Chapters 4-6) ──\n');
fprintf('    1.   pi_baseline       (BOPTEST built-in PI)\n');
fprintf('    2.   mpc_baseline      (MPC with mismatch, no bias states)\n');
fprintf('    3.   ekf_mpc           (offset-free MPC, augmented KF)\n');
fprintf('    4.   rl_mpc_loop       (RL training on the MPC)\n');
fprintf('    5.   rl_mpc_benchmark  (greedy eval of that theta)\n');
fprintf('    6.   rl_mpc_loop       (RL training on the offset-free MPC)\n');
fprintf('    7.   rl_mpc_benchmark  (greedy eval, offset-free)\n');
fprintf('    8.   bo_mpc_tune       (Bayesian optimisation)          [~5 h]\n');
fprintf('    9.   bo_mpc_benchmark  (eval of the BO-best theta)\n');
fprintf('   ── BO on the offset-free MPC (Section 6.4) ──\n');
fprintf('   10.   bo_mpc_tune       (offset-free base)               [~5 h]\n');
fprintf('   11.   bo_mpc_benchmark  (eval, offset-free base)\n');
fprintf('   ── Forecast uncertainty (Section 7.3) ──\n');
fprintf('   12.   sweep_uncertainty         (5 controllers x 2 levels) [~2 h]\n');
fprintf('   13.   build_uncertainty_report  (table + figure, no server)\n');
fprintf('   ── Out-of-sample weeks (Section 7.4) ──\n');
fprintf('   14.   sweep_out_of_sample(2)    (identification record)  [~20 min]\n');
fprintf('   15.   sweep_out_of_sample(3)    (validation record)      [~20 min]\n');
fprintf('   16.   sweep_out_of_sample(4)    (validation record)      [~20 min]\n');
fprintf('   17.   sweep_out_of_sample(5)    (untouched week)         [~20 min]\n');
fprintf('   18.   plot_out_of_sample        (week ladder, no server)\n');
fprintf('   ── Figures and checks ──\n');
fprintf('   19.   generate_controller_figures  (Ch. 4-6, no server)\n');
fprintf('   20.   gradient_check               (Section 5.4, no server)\n');
fprintf('    0.   cancel\n');
fprintf('──────────────────────────────────────────────────────────\n');

choice = input('  Enter choice [0-20]: ');

switch choice
    case 1,  pi_baseline
    case 2,  mpc_baseline
    case 3,  ekf_mpc
    case 4,  use_bias = false; rl_mpc_loop
    case 5,  use_bias = false; rl_mpc_benchmark
    case 6,  use_bias = true;  rl_mpc_loop
    case 7,  use_bias = true;  rl_mpc_benchmark
    case 8,  use_bias = false; bo_mpc_tune
    case 9,  use_bias = false; bo_mpc_benchmark
    case 10, use_bias = true;  bo_mpc_tune
    case 11, use_bias = true;  bo_mpc_benchmark
    case 12, sweep_uncertainty
    case 13, build_uncertainty_report
    case 14, sweep_out_of_sample(2)
    case 15, sweep_out_of_sample(3)
    case 16, sweep_out_of_sample(4)
    case 17, sweep_out_of_sample(5)
    case 18, plot_out_of_sample
    case 19, generate_controller_figures
    case 20, gradient_check
    case 0,  fprintf('  Cancelled.\n');
    otherwise, fprintf('  Invalid choice.\n');
end
