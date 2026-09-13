function sweep_uncertainty()
%SWEEP_UNCERTAINTY  The five controllers under a degraded forecast (Section 7.3).
%
%  BOPTEST's uncertainty emulator perturbs the forecast handed to the controller
%  while leaving the simulated building untouched. This driver re-runs the five
%  reported controllers at the two degraded levels, medium and high, applied to
%  both the ambient temperature and the solar irradiation.
%
%  Nothing is retuned. Every controller deploys exactly the parameters of the
%  previous chapters, so the three columns of Table 7.3 describe one and the
%  same controller under three forecasts.
%
%  The clean column is not produced here: it is the reported run of each
%  controller, which run_baselines options 1, 2, 3, 7 and 9 already wrote. That
%  keeps a sweep from ever overwriting a reported result.
%
%  Written as a function, not a script, because every runner opens with
%  `clearvars -except BOPTEST_*`, which would wipe the loop state if it ran in
%  the caller's workspace. The helper at the end gives each run its own.
%
%  Existing outputs are skipped, so the sweep is safe to resume.
%
%  BEFORE RUNNING:
%    1. startup
%    2. the five reported runs (run_baselines options 1, 2, 3, 7 and 9)
%    3. BOPTEST server reachable at http://127.0.0.1:80
%
%  Produces the ten degraded-forecast files build_uncertainty_report aggregates.

levels = {'medium', 'high'};

%        runner script       use_bias   tagged output
ctrl = { ...
 'pi_baseline',      [],     'pi%s_simdata.mat'                    ; ...
 'mpc_baseline',     [],     'mpc_baseline%s_simdata.mat'          ; ...
 'ekf_mpc',          [],     'ekf_mpc%s_simdata.mat'               ; ...
 'rl_mpc_benchmark', true,   'rl_mpc_benchmark_bias%s_simdata.mat' ; ...
 'bo_mpc_benchmark', false,  'bo_mpc%s_simdata.mat'                };

fprintf('\n══════════════════════════════════════════════════════════\n');
fprintf('  Forecast-uncertainty sweep\n');
fprintf('  %d controllers × %d degraded levels = %d episodes\n', ...
        size(ctrl,1), numel(levels), size(ctrl,1)*numel(levels));
fprintf('══════════════════════════════════════════════════════════\n');

t_total = tic;

for c = 1:size(ctrl, 1)
    for L = 1:numel(levels)
        lvl = levels{L};
        tag = ['_' lvl];
        out = sprintf(ctrl{c,3}, tag);

        if exist(fullfile('models', out), 'file')
            fprintf('\n[skip] %-18s %-6s — %s exists\n', ctrl{c,1}, lvl, out);
            continue
        end

        fprintf('\n────────────────────────────────────────────────────\n');
        fprintf('  %s  |  %s forecast noise  →  %s\n', ctrl{c,1}, lvl, out);
        fprintf('────────────────────────────────────────────────────\n');

        t_run = tic;
        try
            invoke_runner(ctrl{c,1}, lvl, tag, ctrl{c,2});
        catch ME
            fprintf(2, '[FAIL] %s under %s: %s\n', ctrl{c,1}, lvl, ME.message);
            close all
            continue
        end
        fprintf('[done] %s in %.1f min\n', ctrl{c,1}, toc(t_run)/60);
        close all
    end
end

fprintf('\n══════════════════════════════════════════════════════════\n');
fprintf('  SWEEP COMPLETE in %.1f min\n', toc(t_total)/60);
fprintf('  Next: build_uncertainty_report, to aggregate the KPIs.\n');
fprintf('══════════════════════════════════════════════════════════\n\n');
end


function invoke_runner(script, lvl, out_tag, use_bias_in)
% Run one controller in its own workspace, so its clearvars stays local.
    BOPTEST_TUNC    = lvl;      %#ok<NASGU>  both weather channels are degraded
    BOPTEST_SUNC    = lvl;      %#ok<NASGU>
    BOPTEST_OUT_TAG = out_tag;  %#ok<NASGU>
    RL_THETA_TAG    = '';       %#ok<NASGU>  the θ trained on the clean forecast
    BO_THETA_TAG    = '';       %#ok<NASGU>  the θ tuned on the clean forecast
    if ~isempty(use_bias_in)
        use_bias = use_bias_in; %#ok<NASGU>
    end
    run(script);
end
