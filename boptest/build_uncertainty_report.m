%% build_uncertainty_report.m: Table 7.3 and Figure 7.1
%
%  Collects the KPIs of the five controllers over the three forecast levels and
%  builds the degradation table Chapter 7 reports, together with the bar chart
%  that accompanies it. Every controller keeps its parameters frozen at the
%  nominal values, so a row is one controller read under three forecasts.
%
%  The clean column is the reported run of each controller; the medium and high
%  columns come from sweep_uncertainty.
%
%  This script only reads .mat files, so it needs no BOPTEST server.
%
%  BEFORE RUNNING:
%    1. startup
%    2. the five reported runs (run_baselines options 1, 2, 3, 7 and 9)
%    3. sweep_uncertainty
%
%  Outputs
%    models/uncertainty_summary.mat     (T_degrade)
%    figures/uncertainty_degradation.png

clear; clc;

if ~exist('figures', 'dir'), mkdir('figures'); end

levels       = {'',      '_medium', '_high'};
level_labels = {'clean', 'medium',  'high' };

%  Each runner inserts its tag before '_simdata.mat'. Two of the clean runs keep
%  a name of their own: the offset-free MPC is saved as mpc_simdata, and the RL
%  arm uses the offset-free configuration throughout, so that its three columns
%  describe one controller rather than two.
%
%  The labels are the ones the published figure carries, so a regenerated
%  Figure 7.1 matches the one in the thesis.
%
%        label                 clean file                           tagged pattern
ctrl = { ...
 'PI',                 'pi_simdata.mat',                    'pi%s_simdata.mat'                    ; ...
 'MPC',                'mpc_baseline_simdata.mat',          'mpc_baseline%s_simdata.mat'          ; ...
 'Offset-free MPC',    'mpc_simdata.mat',                   'ekf_mpc%s_simdata.mat'               ; ...
 'RL-MPC (clean \phi)','rl_mpc_benchmark_bias_simdata.mat', 'rl_mpc_benchmark_bias%s_simdata.mat' ; ...
 'BO-MPC (clean \phi)','bo_mpc_simdata.mat',                'bo_mpc%s_simdata.mat'                };

nC = size(ctrl, 1);
nL = numel(levels);

%% 1. Collect the KPIs
T_degrade = cell2table(repmat({'', '', NaN, NaN, NaN, NaN}, nC*nL, 1), ...
    'VariableNames', {'Controller','Condition','tdis_tot','ener_tot','cost_tot','emis_tot'});

row = 0;
for c = 1:nC
    for L = 1:nL
        if isempty(levels{L})
            fname = ctrl{c,2};
        else
            fname = sprintf(ctrl{c,3}, levels{L});
        end
        kpi = load_kpi(fullfile('models', fname));
        row = row + 1;

        if isempty(kpi)
            fprintf('  [miss] %-40s %s / %s\n', fname, ctrl{c,1}, level_labels{L});
            T_degrade(row, :) = {ctrl{c,1}, level_labels{L}, NaN, NaN, NaN, NaN};
        else
            fprintf('  [ ok ] %-40s %s / %s\n', fname, ctrl{c,1}, level_labels{L});
            T_degrade(row, :) = {ctrl{c,1}, level_labels{L}, ...
                                 kpi.tdis_tot, kpi.ener_tot, kpi.cost_tot, kpi.emis_tot};
        end
    end
end

fprintf('\n══════════════════════════════════════════════════════════════════════════\n');
fprintf('  DEGRADATION TABLE — nominal φ deployed under three forecast levels\n');
fprintf('══════════════════════════════════════════════════════════════════════════\n');
disp(T_degrade);

save('models/uncertainty_summary.mat', 'T_degrade');
fprintf('  Saved: models/uncertainty_summary.mat\n\n');

%% 2. Figure
metrics  = {'tdis_tot',    'ener_tot',  'cost_tot',  'emis_tot'};
m_units  = {'K·h/zone',    'kWh/m²',    '$/m²',      'kgCO_2/m²'};
m_labels = {'Discomfort',  'Energy',    'Cost',      'Emissions'};
m_titles = {'T_{dis,tot}', 'E_{tot}',   'C_{tot}',   'M_{CO_2,tot}'};

fig = figure('Name', 'Forecast-uncertainty degradation', ...
             'Color', 'w', 'Position', [60 60 1300 700]);

for m = 1:numel(metrics)
    subplot(2, 2, m);
    G = reshape(T_degrade.(metrics{m}), nL, nC).';   % rows = controller, cols = level
    bar(G, 'grouped');
    set(gca, 'XTickLabel', ctrl(:,1), 'FontSize', 9);
    xtickangle(20);
    ylabel(sprintf('%s [%s]', m_labels{m}, m_units{m}));
    title(m_titles{m}, 'FontWeight', 'normal');
    % Headroom, so the legend never sits on top of the bars.
    ymax = max(G(:), [], 'omitnan');
    if isfinite(ymax) && ymax > 0, ylim([0, ymax * 1.30]); end
    legend(level_labels, 'Location', 'north', 'Orientation', 'horizontal', ...
           'FontSize', 8, 'Box', 'off');
    grid on; box on;
end
sgtitle('Controller degradation under forecast uncertainty', 'FontWeight', 'bold');

exportgraphics(fig, 'figures/uncertainty_degradation.png', 'Resolution', 150);
fprintf('  Figure saved: figures/uncertainty_degradation.png\n');
