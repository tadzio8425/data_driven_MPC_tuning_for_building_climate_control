%% data_collection_baseline.m: Passive data collection with BOPTEST baseline controller
%
%  Collects 4 weeks of data from the BOPTEST multizone_residential_hydronic
%  testcase WITHOUT overriding any actuator.  The built-in PI controllers
%  drive the heating system exactly as they would in a real building.
%
%  This simulates the realistic scenario where an MPC designer only has
%  access to operational (closed-loop) data: no PRBS excitation and no valve
%  overrides, because the building is occupied and comfort must be maintained.
%
%  Output: data/baseline_identification.csv  (weeks 1-2)
%          data/baseline_validation.csv      (weeks 3-4)
%
%  BEFORE RUNNING:
%    Start the BOPTEST Docker container on port 80.

url         = 'http://127.0.0.1';
testcase    = 'multizone_residential_hydronic';
step_period = 900;                           % 15 min per step

%% Timing
%  Aligned with the MLPRS experiment: the same four weeks from the start of the
%  simulated year, days 0 to 28, with a one-week warm-up before them. The two
%  datasets therefore cover the same weather and the same occupancy schedule,
%  which is what makes the excitation comparison of Section 3.3 fair.
start_time  = 0;                            % [s] from 1 January, as in the MLPRS run
warm_up     = 7 * 24 * 3600;                % 1-week warm-up [s]
max_steps   = (4*7*24*3600) / step_period;  % 4 weeks total = 2688 steps

%% Zone definitions
ctrl_zones    = {'Liv', 'Ro1', 'Ro2', 'Ro3', 'Bth'};
all_zones     = {'Liv', 'Ro1', 'Ro2', 'Ro3', 'Bth', 'Hal'};
nZones_ctrl   = numel(ctrl_zones);
nZones_all    = numel(all_zones);

zone_temp_sig = { ...
    'conHeaLiv_reaTZon_y', 'conHeaRo1_reaTZon_y', 'conHeaRo2_reaTZon_y', ...
    'conHeaRo3_reaTZon_y', 'conHeaBth_reaTZon_y', 'reaTHal_y' };

zone_heat_sig = { ...
    'reaHeaLiv_y', 'reaHeaRo1_y', 'reaHeaRo2_y', ...
    'reaHeaRo3_y', 'reaHeaBth_y' };

%% Disturbance forecast signals
forecast_names = { ...
    'InternalGainsCon[Liv]', 'InternalGainsRad[Liv]', ...
    'InternalGainsCon[Ro1]', 'InternalGainsRad[Ro1]', ...
    'InternalGainsCon[Ro2]', 'InternalGainsRad[Ro2]', ...
    'InternalGainsCon[Ro3]', 'InternalGainsRad[Ro3]', ...
    'InternalGainsCon[Bth]', 'InternalGainsRad[Bth]', ...
    'InternalGainsCon[Hal]', 'InternalGainsRad[Hal]'  };
forecast_fields = cellfun(@(s) regexprep(s, '[\[\]]', '_'), ...
                           forecast_names, 'UniformOutput', false);

%% Deploy testcase
handler = RequestHandler(url, 80);
tid     = handler.deploy_test(testcase);
handler.set_scenario(start_time, warm_up);
handler.set_step(step_period);

%% Data collection loop (NO overrides)
all_data = table();

fprintf('Baseline data collection: %d steps (4 weeks)\n', max_steps);
fprintf('  No actuator overrides — BOPTEST PI controllers active.\n\n');

for i = 1:max_steps

    % 1. Fetch one-step-ahead disturbance forecast
    fcast = handler.get_forecast(forecast_names, step_period, step_period);

    % 2. Advance simulation: empty struct = no overrides
    res = handler.advance(struct());

    % 3. Log row
    new_row = struct2table(res, 'AsArray', true);

    % Log delivered heating power per zone
    for k = 1:nZones_ctrl
        if isfield(res, zone_heat_sig{k})
            new_row.(['Qdot_' ctrl_zones{k}]) = res.(zone_heat_sig{k});
        else
            new_row.(['Qdot_' ctrl_zones{k}]) = 0;
        end
    end

    % Disturbance forecast (one-step ahead)
    for k = 1:numel(forecast_names)
        field = forecast_fields{k};
        if isfield(fcast, field)
            val = fcast.(field);
            if isnumeric(val) && ~isempty(val)
                new_row.(['fcast_' field]) = val(1);
            end
        end
    end

    all_data = [all_data; new_row]; %#ok<AGROW>

    if mod(i, 50) == 0
        temps = zeros(1, nZones_all);
        for k = 1:nZones_all
            if isfield(res, zone_temp_sig{k})
                temps(k) = res.(zone_temp_sig{k}) - 273.15;
            end
        end
        fprintf('Step %4d/%d  T=[%s] °C\n', ...
                i, max_steps, sprintf('%.1f ', temps));
    end
end

%% Save
if ~exist('data', 'dir'), mkdir('data'); end
mid = floor(height(all_data) / 2);
writetable(all_data(1:mid,:),     'data/baseline_identification.csv');
writetable(all_data(mid+1:end,:), 'data/baseline_validation.csv');
fprintf('\nDone. %d rows saved (ident: %d  val: %d).\n', ...
        height(all_data), mid, height(all_data)-mid);

handler.stop_test(tid);
