function sweep_out_of_sample(week_in)
%SWEEP_OUT_OF_SAMPLE  The five controllers on a week nobody tuned on (Section 7.4).
%
%  Re-runs the reported comparison with every parameter frozen and only the
%  start of the scenario moved forward, so that the RL θ and the BO θ meet a
%  week they never saw. Week N starts at (N−1)·7 days, and BOPTEST warms up on
%  the week before it, so the building enters the episode warm.
%
%  Nothing is retuned and nothing reported is overwritten: every output carries
%  a _wk<N> tag. Existing outputs are skipped, so the sweep can be resumed.
%
%  USAGE
%    sweep_out_of_sample(2)    identification record, second half
%    sweep_out_of_sample(3)    validation record  (the default)
%    sweep_out_of_sample(4)    validation record, beside week 3
%    sweep_out_of_sample(5)    no role anywhere in the thesis
%  Week 1 needs no sweep: it is the reported week, run by run_baselines options
%  1, 2, 3, 7 and 9. Run several weeks and plot_out_of_sample draws them all.
%
%  WHICH WEEK
%    BOPTEST counts seconds from 1 January and every dataset starts at t = 0:
%      week 1     days  0-7   identification record, and the week θ was tuned on
%      week 2     days  7-14  identification record, θ never tuned here
%      weeks 3-4  days 14-28  validation record, neither the model nor θ saw it
%      week 5     days 28-35  no role in this work
%    The four of them form a ladder, each one step further from what was fitted.
%
%  BEFORE RUNNING:
%    1. startup
%    2. the five reported runs (run_baselines options 1, 2, 3, 7 and 9), which
%       supply the week-1 column
%    3. BOPTEST server reachable at http://127.0.0.1:80
%
%  Outputs
%    models/*_wk<N>_simdata.mat              one per controller
%    models/out_of_sample_wk<N>_summary.mat  week 1 against week N

%% 1. Which week
if nargin < 1 || isempty(week_in)
    week = 3;
else
    week = week_in;
end
if ~isscalar(week) || ~isnumeric(week) || week ~= fix(week) || week < 2
    error(['week must be an integer >= 2. Week 1 is the reported week, and ' ...
           'run_baselines options 1, 2, 3, 7 and 9 already run it.']);
end

start_time = (week - 1) * 7 * 86400;
tag        = sprintf('_wk%d', week);

%        runner script       use_bias  tagged output                          reported output                      label
ctrl = { ...
 'pi_baseline',      [],     'pi%s_simdata.mat',                    'pi_simdata.mat',                    'PI baseline'     ; ...
 'mpc_baseline',     [],     'mpc_baseline%s_simdata.mat',          'mpc_baseline_simdata.mat',          'MPC'             ; ...
 'ekf_mpc',          [],     'ekf_mpc%s_simdata.mat',               'mpc_simdata.mat',                   'Offset-free MPC' ; ...
 'rl_mpc_benchmark', true,   'rl_mpc_benchmark_bias%s_simdata.mat', 'rl_mpc_benchmark_bias_simdata.mat', 'RL-tuned MPC'    ; ...
 'bo_mpc_benchmark', false,  'bo_mpc%s_simdata.mat',                'bo_mpc_simdata.mat',                'BO-tuned MPC'    };

%% 2. Pre-flight: the two tuned controllers need their θ on disk
if ~exist('models/rl_mpc_bias_simdata.mat', 'file')
    warning('rl_mpc_bias_simdata.mat not found — the RL row will fail. Train it first.');
end
if ~exist('models/bo_mpc_tune_data.mat', 'file')
    warning('bo_mpc_tune_data.mat not found — the BO row will fail. Tune it first.');
end

fprintf('\n══════════════════════════════════════════════════════════\n');
fprintf('  Out-of-sample sweep — week %d\n', week);
fprintf('  start_time = %d s (day %d), warm-up on the week before it\n', ...
        start_time, start_time/86400);
fprintf('  %d controllers, clean forecast, every parameter frozen\n', size(ctrl,1));
fprintf('══════════════════════════════════════════════════════════\n');

%% 3. Run
t_total = tic;

for c = 1:size(ctrl, 1)
    script = ctrl{c,1};
    out    = sprintf(ctrl{c,3}, tag);

    if exist(fullfile('models', out), 'file')
        fprintf('\n[skip] %s — %s exists\n', script, out);
        continue
    end

    fprintf('\n────────────────────────────────────────────────────\n');
    fprintf('  [%d/%d] %s  |  %s on week %d  →  %s\n', ...
            c, size(ctrl,1), script, ctrl{c,5}, week, out);
    fprintf('────────────────────────────────────────────────────\n');

    t_run = tic;
    try
        invoke_runner(script, start_time, tag, ctrl{c,2});
    catch ME
        fprintf(2, '[FAIL] %s on week %d: %s\n', script, week, ME.message);
        close all
        continue
    end
    fprintf('[done] %s in %.1f min\n', script, toc(t_run)/60);
    close all
end

fprintf('\n══════════════════════════════════════════════════════════\n');
fprintf('  SWEEP COMPLETE in %.1f min\n', toc(t_total)/60);
fprintf('══════════════════════════════════════════════════════════\n');

%% 4. Week 1 against week N
vars  = {'Controller', 'tdis_wk1', sprintf('tdis_wk%d', week), 'tdis_delta', ...
         'ener_wk1',   sprintf('ener_wk%d', week)};
T_oos = cell2table(repmat({'', NaN, NaN, NaN, NaN, NaN}, size(ctrl,1), 1), ...
                   'VariableNames', vars);

for c = 1:size(ctrl, 1)
    kpi_in  = load_kpi(fullfile('models', ctrl{c,4}));
    kpi_out = load_kpi(fullfile('models', sprintf(ctrl{c,3}, tag)));
    if isempty(kpi_out)
        fprintf('  [miss] %-40s week %d\n', sprintf(ctrl{c,3}, tag), week);
    end
    t_in  = getkpi(kpi_in,  'tdis_tot');
    t_out = getkpi(kpi_out, 'tdis_tot');
    T_oos(c, :) = {ctrl{c,5}, t_in, t_out, t_out - t_in, ...
                   getkpi(kpi_in, 'ener_tot'), getkpi(kpi_out, 'ener_tot')};
end

fprintf('\n══════════════════════════════════════════════════════════════════\n');
fprintf('  WEEK 1 (reported, in sample) AGAINST WEEK %d (out of sample)\n', week);
fprintf('  Discomfort in K·h/zone, energy in kWh/m²\n');
fprintf('══════════════════════════════════════════════════════════════════\n');
disp(T_oos);

out_mat = fullfile('models', sprintf('out_of_sample_wk%d_summary.mat', week));
save(out_mat, 'T_oos', 'week', 'start_time');
fprintf('  Saved: %s\n', out_mat);
fprintf('  Next:  plot_out_of_sample, once every week has been run.\n\n');
end


%% Helpers
function invoke_runner(script, start_time, out_tag, use_bias_in)
% Run one controller in its own workspace, so its clearvars stays local.
    BOPTEST_TUNC       = '';          %#ok<NASGU>  clean forecast, as in Table 7.2
    BOPTEST_SUNC       = '';          %#ok<NASGU>
    BOPTEST_OUT_TAG    = out_tag;     %#ok<NASGU>
    BOPTEST_START_TIME = start_time;  %#ok<NASGU>
    RL_THETA_TAG       = '';          %#ok<NASGU>  the θ trained on week 1
    BO_THETA_TAG       = '';          %#ok<NASGU>  the θ tuned on week 1
    if ~isempty(use_bias_in)
        use_bias = use_bias_in;       %#ok<NASGU>
    end
    run(script);
end

function v = getkpi(kpi, name)
    if ~isempty(kpi) && isfield(kpi, name), v = kpi.(name); else, v = NaN; end
end
