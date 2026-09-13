%% plot_out_of_sample.m: Table 7.4 and Figure 7.2, the week ladder
%
%  Reads the KPI of every controller on every week already on disk and draws one
%  line per controller across the weeks. What the figure is about is the
%  ordering, not the level: a controller whose correction is recomputed online
%  holds its position, while a controller whose correction was fitted on the
%  identification record loses it as soon as the week falls outside that record.
%
%  The second panel exists because the weeks are not equally hard. The PI
%  baseline alone ranges from 10.4 to 13.9 K·h/zone across the five weeks, so on
%  the absolute panel every line rises and falls together and the weather does
%  most of the moving. Dividing by the PI baseline of the same week removes that
%  common mode: PI becomes the flat line at one, anything above it is worse than
%  the controller the building already has, and what is left is the ordering.
%
%  This script only reads .mat files, so it needs no BOPTEST server. It draws
%  whatever weeks are on disk and skips the rest.
%
%  BEFORE RUNNING:
%    1. startup
%    2. the five reported runs (run_baselines options 1, 2, 3, 7 and 9), which
%       supply the week-1 column
%    3. one or more out-of-sample sweeps, for the remaining columns
%
%  Outputs
%    figures/out_of_sample_slope.png
%    models/out_of_sample_slope_data.mat   (T_slope, weeks, Tdis, Ener, Rel)

clear; clc;

WEEK_MAX = 12;   % past here the heating season stops being a heating season
if ~exist('figures', 'dir'), mkdir('figures'); end

%% 1. Controllers, in the order of the reported table
%  Column 2 is the week-1 filename, which carries no tag; column 3 is the
%  pattern the sweep writes for every other week. The colours are the ZCOL
%  palette of mpc/comfort_plot, so this figure sits beside the rest of the
%  thesis: blue for the offset-free MPC, orange-red for the BO-tuned one, green
%  for the RL-tuned one, purple for the uncorrected MPC and a dashed grey for
%  the PI baseline, which is a reference rather than a proposal. Distinct
%  markers keep the five series separable in greyscale.
ctrl = { ...
 'PI baseline',          'pi_simdata.mat',                    'pi_wk%d_simdata.mat',                    [0.50 0.50 0.50], 1.1, '--o' ; ...
 'MPC (uncorrected)',    'mpc_baseline_simdata.mat',          'mpc_baseline_wk%d_simdata.mat',          [0.58 0.40 0.74], 1.3, '-.^' ; ...
 'Offset-free MPC',      'mpc_simdata.mat',                   'ekf_mpc_wk%d_simdata.mat',               [0.12 0.47 0.71], 1.6, '-o'  ; ...
 'RL-tuned offset-free', 'rl_mpc_benchmark_bias_simdata.mat', 'rl_mpc_benchmark_bias_wk%d_simdata.mat', [0.27 0.63 0.27], 1.6, '-s'  ; ...
 'BO-tuned MPC',         'bo_mpc_simdata.mat',                'bo_mpc_wk%d_simdata.mat',                [0.86 0.35 0.19], 1.6, '-d'  };

nC    = size(ctrl, 1);
i_ref = 1;   % the row every other row is normalised against

%% 2. Collect every week that has at least one controller on disk
weeks = [];
for w = 1:WEEK_MAX
    for c = 1:nC
        if exist(fullfile('models', week_filename(ctrl(c,:), w)), 'file')
            weeks(end+1) = w; %#ok<SAGROW>
            break
        end
    end
end
if isempty(weeks)
    error(['No simdata found in models/. Run the five reported controllers ' ...
           'first (run_baselines options 1, 2, 3, 7 and 9).']);
end

Tdis = NaN(nC, numel(weeks));
Ener = NaN(nC, numel(weeks));
for c = 1:nC
    for k = 1:numel(weeks)
        kpi = load_kpi(fullfile('models', week_filename(ctrl(c,:), weeks(k))));
        if ~isempty(kpi)
            Tdis(c,k) = kpi.tdis_tot;
            Ener(c,k) = kpi.ener_tot;
        end
    end
end

% Normalised to the reference controller on the same week. A week whose
% reference run is missing normalises to NaN rather than to something misleading.
Rel = Tdis ./ repmat(Tdis(i_ref,:), nC, 1);

%% 3. Print the tables, so the numbers can be read off for the thesis
vnames  = [{'Controller'}, arrayfun(@(w) sprintf('wk%d', w), weeks, 'Uni', false)];
T_slope = cell2table([ctrl(:,1), num2cell(Tdis)], 'VariableNames', vnames);

fprintf('\n══════════════════════════════════════════════════════════════════\n');
fprintf('  THERMAL DISCOMFORT ACROSS THE WEEK LADDER  [K·h/zone]\n');
fprintf('══════════════════════════════════════════════════════════════════\n');
disp(T_slope);

fprintf('  Week roles: ');
for k = 1:numel(weeks), fprintf('%d=%s   ', weeks(k), week_role(weeks(k))); end

fprintf('\n\n  Relative to the %s on the same week (reference = 1.00):\n', ctrl{i_ref,1});
print_matrix(Rel, ctrl, weeks, '%7.2f');
fprintf('\n  Energy [kWh/m²]:\n');
print_matrix(Ener, ctrl, weeks, '%7.2f');

% The sign of this row is the finding: negative means the tuning stage paid off
% on that week, positive means it cost more than it returned.
i_of = find(strcmp(ctrl(:,1), 'Offset-free MPC'), 1);
fprintf('\n  Gain over the offset-free MPC (negative = the tuning helped):\n');
print_matrix(Tdis - repmat(Tdis(i_of,:), nC, 1), ctrl, weeks, '%+7.2f');

missing = sum(isnan(Tdis(:)));
if missing > 0
    fprintf('\n  %d of %d cells are still empty — run the missing sweeps.\n', ...
            missing, numel(Tdis));
end
fprintf('\n');

%% 4. The figure: absolute on top, normalised below
fig = figure('Name', 'Out-of-sample week ladder', 'Color', 'w', ...
             'Position', [60 60 820 730]);

ax_abs = subplot(2, 1, 1);
draw_ladder(ax_abs, Tdis, weeks, ctrl, 'T_{dis} [K\cdoth/zone]', ...
            struct('annotate_best', true, 'show_xlabels', false, ...
                   'show_legend', false, 'label_bands', true, 'ref_line', NaN));

ax_rel = subplot(2, 1, 2);
draw_ladder(ax_rel, Rel, weeks, ctrl, 'T_{dis} relative to PI [-]', ...
            struct('annotate_best', false, 'show_xlabels', true, ...
                   'show_legend', true, 'label_bands', false, 'ref_line', 1));

exportgraphics(fig, 'figures/out_of_sample_slope.png', 'Resolution', 150);
save('models/out_of_sample_slope_data.mat', 'T_slope', 'weeks', 'Tdis', 'Ener', 'Rel');

fprintf('  Figure saved: figures/out_of_sample_slope.png\n');
fprintf('  Data saved:   models/out_of_sample_slope_data.mat\n\n');


%% Local helpers
function f = week_filename(ctrl_row, w)
% Week 1 is the reported run and carries no tag; every other week is tagged.
    if w == 1
        f = ctrl_row{2};
    else
        f = sprintf(ctrl_row{3}, w);
    end
end

function r = week_role(w)
% What the controllers had already seen of this week when they were built.
    if     w == 1, r = 'θ tuned here';
    elseif w == 2, r = 'model fitted';
    elseif w <= 4, r = 'validation';
    else,          r = 'untouched';
    end
end

function print_matrix(M, ctrl, weeks, fmt)
% One row per controller, one column per week, under a week header.
    fprintf('    %-24s', '');
    for k = 1:numel(weeks), fprintf('%7s', sprintf('wk%d', weeks(k))); end
    fprintf('\n');
    for c = 1:size(M, 1)
        fprintf('    %-24s', ctrl{c,1});
        fprintf(fmt, M(c,:));   % the format cycles over the row
        fprintf('\n');
    end
end

function draw_ladder(ax, Y, weeks, ctrl, ylab, opts)
% One line per controller across the weeks, in the style of the other figures of
% the thesis: grid on, box on, a legend with no frame, reference level as yline.
    nC = size(ctrl, 1);
    n  = numel(weeks);
    x  = 1:n;

    hold(ax, 'on');

    ymax = max(Y(:), [], 'omitnan');
    if ~isfinite(ymax) || ymax <= 0, ymax = 1; end
    ylim(ax, [0, ymax * 1.15]);
    xlim(ax, [0.6, n + 0.4]);

    % Shade the two nested records the fitted corrections depend on. The outer
    % band is the identification record, which fixes the model, and the darker
    % sub-band inside it is the single week θ was tuned on. Everything to the
    % right of the outer edge is a week no stage of the design ever saw, and
    % that edge is where the sign of both fitted corrections changes, so the
    % shading carries the finding rather than decorating the panel.
    yr    = ylim(ax);
    k_tun = find(weeks <= 1, 1, 'last');   % tuned on week 1 alone
    k_id  = find(weeks <= 2, 1, 'last');   % identification record, days 0-14
    if ~isempty(k_tun)
        shade(ax, [0.6, k_tun + 0.5], yr, [0.84 0.87 0.92]);
    end
    if ~isempty(k_id) && k_id > k_tun
        shade(ax, [k_tun + 0.5, k_id + 0.5], yr, [0.92 0.94 0.96]);
        plot(ax, (k_tun + 0.5)*[1 1], yr, ':', 'Color', [0.62 0.66 0.72], ...
             'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
    if opts.label_bands && ~isempty(k_id)
        text(ax, (0.6 + k_id + 0.5)/2, yr(2), 'identification record', ...
             'FontSize', 8, 'FontAngle', 'italic', 'Color', [0.35 0.35 0.35], ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        text(ax, (0.6 + k_tun + 0.5)/2, yr(2) - 0.062*diff(yr), '\theta tuned', ...
             'FontSize', 8, 'FontAngle', 'italic', 'Color', [0.30 0.34 0.42], ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    % On the normalised panel the reference level is the whole point: above it
    % the controller is worse than the baseline the building already has.
    if isfinite(opts.ref_line)
        yline(ax, opts.ref_line, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end

    h = gobjects(nC, 1);
    for c = 1:nC
        h(c) = plot(ax, x, Y(c,:), ctrl{c,6}, 'Color', ctrl{c,4}, ...
                    'LineWidth', ctrl{c,5}, 'MarkerSize', 5, ...
                    'MarkerFaceColor', ctrl{c,4}, 'MarkerEdgeColor', 'none');
    end

    % Call out the best week-1 value: it is the number the thesis headlines. The
    % label sits to the left of the marker, since the best value is low on the
    % axis and anything below it would be clipped by the axis line.
    if opts.annotate_best && any(isfinite(Y(:,1)))
        [vbest, cbest] = min(Y(:,1));
        text(ax, 1 - 0.09, vbest, sprintf('%.1f', vbest), ...
             'Color', ctrl{cbest,4}, 'FontSize', 8, ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    end

    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'XTick', x, 'FontSize', 9, 'Layer', 'top');
    ylabel(ax, ylab);

    if opts.show_xlabels
        set(ax, 'XTickLabel', arrayfun(@(w) sprintf('Week %d', w), weeks, 'Uni', false));
        xlabel(ax, 'Evaluation week');
    else
        set(ax, 'XTickLabel', []);
    end

    if opts.show_legend
        legend(ax, h, ctrl(:,1), 'Location', 'southoutside', ...
               'Orientation', 'horizontal', 'NumColumns', 3, ...
               'FontSize', 8, 'Box', 'off');
    end

    hold(ax, 'off');
end

function shade(ax, xr, yr, col)
% One flat background band between two x limits, drawn behind everything else.
    h = patch(ax, [xr(1) xr(2) xr(2) xr(1)], [yr(1) yr(1) yr(2) yr(2)], col, ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(h, 'bottom');
end
