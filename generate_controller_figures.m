%% generate_controller_figures.m: Export the controller figures of the thesis
%
%  Redraws every figure of Chapters 4, 5 and 6 straight from the simdata under
%  models/. No BOPTEST server is needed: the .mat files have to exist already,
%  and any entry whose file is missing is skipped without an error.
%
%  Run from src/ after startup.
%
%  Output tree, one entry per figure the thesis includes:
%    figures/
%      ch04_baselines/
%        mpc_baseline.png        the MPC that carries the mismatch  (Fig. 4.2)
%        akf_mpc.png             the offset-free MPC                (Fig. 4.3)
%      ch05_rl/
%        comfort_benchmark.png   RL on the MPC without the correction (Fig. 5.2)
%        rl_param_evolution.png  parameter and TD-error traces        (Fig. 5.3)
%        rl_woffset.png          learned offsets, with and without the estimator
%      ch06_bo/
%        bo_convergence.png      the twenty BO episodes             (Fig. 6.1)
%        bo_mpc.png              the BO-tuned MPC                   (Fig. 6.2)
%
%  The Chapter 3 figures come out of ekf_sysid_2C, and the two Chapter 7 figures
%  out of build_uncertainty_report and plot_out_of_sample.

clearvars; clc;

RC     = load('models/rc_params_2C6z.mat');
ALL_ZN = RC.ZONE_NAMES;        % 1x6, including Hal
ZN     = RC.ZONE_NAMES(1:5);   % the controlled zones
Ts     = RC.Ts;

nYc   = numel(ZN);
idx_q = 1:nYc;  idx_r = nYc+1;  idx_w = nYc+2 : 2*nYc+1;   % columns of log_theta

for d = {'figures/ch04_baselines', 'figures/ch05_rl', 'figures/ch06_bo'}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

%% 1. Comfort figures
%  comfort_benchmark.png is drawn from the GREEDY evaluation, not from the
%  training episode, so that its header carries the row Table 5.1 reports.
%
%       simdata file                     output                            controller key
COMFORT = { ...
 'mpc_baseline_simdata.mat',      'ch04_baselines/mpc_baseline.png', 'mpc_baseline'  ; ...
 'mpc_simdata.mat',               'ch04_baselines/akf_mpc.png',      'ekf_mpc'       ; ...
 'rl_mpc_benchmark_simdata.mat',  'ch05_rl/comfort_benchmark.png',   'rl_mpc_bench'  ; ...
 'bo_mpc_simdata.mat',            'ch06_bo/bo_mpc.png',              'bo_mpc'        };

fprintf('=== Comfort figures (%d) ===\n\n', size(COMFORT,1));
for k = 1:size(COMFORT,1)
    fpath = fullfile('models', COMFORT{k,1});
    if ~exist(fpath, 'file')
        fprintf('[skip]  %s\n', COMFORT{k,1});
        continue
    end
    fprintf('[gen]   %s\n', COMFORT{k,1});

    S = load(fpath);
    [Tz, ref, ub, Tamb] = extract_logs(S);
    f = comfort_plot(COMFORT{k,3}, Tz, ref, ub, Tamb, S.kpi, ALL_ZN, Ts, ...
                     'SaveTo', fullfile('figures', COMFORT{k,2}));
    close(f);
end

%% 2. RL parameter evolution (Figure 5.3)
%  The training trace on the MPC without the correction: learned weights,
%  reference offsets and the TD error in one figure. The axis labels use the
%  symbols of the thesis, not the internal field names.
fprintf('\n=== RL parameter evolution ===\n');
if exist('models/rl_mpc_simdata.mat', 'file')
    S   = load('models/rl_mpc_simdata.mat', 'log_theta', 'log_tau');
    t_h = (1:size(S.log_theta,1)).' * Ts / 3600;

    f = figure('Color','w','Position',[60 60 1200 760],'Visible','off');

    subplot(2,2,1);
    plot(t_h, S.log_theta(:,idx_q)); grid on; box on;
    ylabel('$Q_i$', 'Interpreter','latex');
    title('Tracking weights', 'FontWeight','normal');
    legend(ZN, 'Location','best', 'FontSize',7, 'Box','off');

    subplot(2,2,2);
    plot(t_h, S.log_theta(:,idx_r)); grid on; box on;
    ylabel('$R$', 'Interpreter','latex');
    title('Input penalty (shared)', 'FontWeight','normal');

    subplot(2,2,3);
    plot(t_h, S.log_theta(:,idx_w)); grid on; box on;
    yline(0, 'k-', 'LineWidth',0.6, 'HandleVisibility','off');
    ylabel('$\Delta \bar{y}_i$ [K]', 'Interpreter','latex');
    xlabel('Time [h]');
    title('Reference offsets', 'FontWeight','normal');
    legend(ZN, 'Location','best', 'FontSize',7, 'Box','off');

    subplot(2,2,4);
    plot(t_h, S.log_tau, 'Color',[0.2 0.4 0.8]); grid on; box on;
    ylabel('$\tau$', 'Interpreter','latex');
    xlabel('Time [h]');
    title('TD error', 'FontWeight','normal');

    sgtitle('RL-MPC: parameter and TD-error evolution', 'FontWeight','bold');
    save_png(f, 'figures/ch05_rl/rl_param_evolution.png');
else
    fprintf('[skip]  rl_param_evolution — rl_mpc_simdata.mat not found\n');
end

%% 3. Learned reference offsets, with and without the estimator
%  Both bars are tail-averaged over the last five days, which is exactly what
%  rl_mpc_benchmark does to θ before evaluating it, so the figure shows the
%  offsets the benchmarked policies actually use. The directional contrast,
%  negative without the estimator and positive with it, is the chapter's
%  diagnosis.
fprintf('\n=== RL reference-offset comparison ===\n');
if exist('models/rl_mpc_simdata.mat', 'file') && exist('models/rl_mpc_bias_simdata.mat', 'file')
    A = load('models/rl_mpc_simdata.mat',      'log_theta');   % without the estimator
    B = load('models/rl_mpc_bias_simdata.mat', 'log_theta');   % with it
    W = [tail_mean(A.log_theta(:, idx_w), Ts), tail_mean(B.log_theta(:, idx_w), Ts)];

    C_NB = [0.85 0.33 0.10];   % orange: without the estimator
    C_OF = [0.12 0.47 0.71];   % blue:   with the offset-free estimator

    f  = figure('Color','w','Position',[80 80 760 440],'Visible','off');
    hb = bar(W, 'grouped', 'EdgeColor','none'); hold on;
    hb(1).FaceColor = C_NB;  hb(2).FaceColor = C_OF;
    yline(0, 'k-', 'LineWidth',0.8, 'HandleVisibility','off');
    set(gca, 'XTickLabel', ZN, 'FontSize',12); box on; grid on;
    ylabel('$\Delta \bar{y}_i$  [K]', 'Interpreter','latex', 'FontSize',13);
    xlabel('Zone');
    legend({'Without estimator','With offset-free estimator'}, ...
           'Location','northoutside','Orientation','horizontal','Box','off','FontSize',9);
    title('Learned reference offsets per zone', 'FontWeight','bold');
    save_png(f, 'figures/ch05_rl/rl_woffset.png');
else
    fprintf('[skip]  rl_woffset — needs both rl_mpc and rl_mpc_bias simdata\n');
end

%% 4. BO convergence (Figure 6.1)
%  The running best against the realised cost of each individual episode, over
%  the twenty evaluations of the search.
fprintf('\n=== BO convergence ===\n');
if exist('models/bo_mpc_tune_data.mat', 'file')
    S    = load('models/bo_mpc_tune_data.mat', 'bo_results');
    J    = S.bo_results.ObjectiveTrace;
    Jmin = S.bo_results.ObjectiveMinimumTrace;
    ev   = 1:numel(J);
    C_BO = [0 0.4470 0.7410];

    f = figure('Color','w','Position',[80 80 900 520],'Visible','off');
    plot(ev, J, 'o', 'MarkerEdgeColor', [0.35 0.35 0.35], 'MarkerSize', 6, ...
         'LineWidth', 0.8, 'DisplayName', 'Episode J_{real}'); hold on;
    plot(ev, Jmin, '-o', 'Color', C_BO, 'LineWidth', 2, 'MarkerSize', 6, ...
         'MarkerFaceColor', C_BO, 'MarkerEdgeColor', C_BO, 'DisplayName', 'Best so far');
    grid on; box on;
    set(gca, 'FontSize', 12, 'XTick', 2:2:numel(J));
    xlim([0.5, numel(J) + 0.5]);
    xlabel('Evaluation (7-day episode)');
    ylabel('J_{real}');
    legend('Location', 'northeast', 'Box', 'off', 'FontSize', 11);
    save_png(f, 'figures/ch06_bo/bo_convergence.png');
else
    fprintf('[skip]  bo_convergence — bo_mpc_tune_data.mat not found\n');
end

fprintf('\nDone. Figures written to figures/.\n');


%% Local helpers
function [Tz, ref, ub, Tamb] = extract_logs(S)
%EXTRACT_LOGS  Pull the logs out of either save format.
%  The PI and MPC runners nest theirs in a sim / sim_pi struct; the RL and BO
%  runners keep theirs at the top level.
    if isfield(S, 'log_Tz')
        Tz = S.log_Tz;  ref = S.log_ref;  ub = S.log_ub;  Tamb = S.log_Tamb;
    elseif isfield(S, 'sim')
        Tz = S.sim.Tz;  ref = S.sim.ref;  ub = S.sim.ref_ub;  Tamb = S.sim.Tamb;
    else
        Tz = S.sim_pi.Tz;  ref = S.sim_pi.ref;  ub = S.sim_pi.ref_ub;  Tamb = S.sim_pi.Tamb;
    end
end

function m = tail_mean(X, Ts)
%TAIL_MEAN  Average the last five days of a logged trace, one value per column.
    win = min(size(X,1), round(5 * 86400 / Ts));
    m   = mean(X(end-win+1:end, :), 1).';
end

function save_png(f, fpath)
    print(f, fpath, '-dpng', '-r150');
    fprintf('        → %s\n', fpath);
    close(f);
end
