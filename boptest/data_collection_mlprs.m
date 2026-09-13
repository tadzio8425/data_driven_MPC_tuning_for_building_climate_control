%% data_collection_mlprs.m: MLPRS valve excitation for grey-box system ID
%
%  Collects 4 weeks of heating data from the BOPTEST multizone_residential_hydronic
%  testcase using a Multi-Level Pseudo-Random Sequence (MLPRS) signal design.
%
%  MLPRS vs. RANDOM MULTI-LEVEL (data_collection_prbs.m)
% 
%  data_collection_prbs.m  :  Random hold durations, random level selection.
%                             Good coverage but irregular spectrum.
%
%  data_collection_mlprs.m :  Levels determined by a maximal-length Linear
%                             Feedback Shift Register (LFSR) clocked at a
%                             fixed rate.  Key advantages:
%
%    Flat power spectrum    PRBS has near-flat PSD up to the clock
%                           frequency.  Excitation energy is spread evenly
%                           across ALL thermal time constants of interest
%                           (τ_air ≈ 1-3 h, τ_wall ≈ 5-20 h).
%
%    Known autocorrelation : R_uu(τ) ≈ δ(τ) (impulse-like).  This makes
%                            the information matrix well-conditioned and
%                            simplifies bias/variance analysis of the
%                            identified model.
%
%    Reproducible experiment : no rng() dependence.  Given the same LFSR
%                              seed and taps, the exact same valve sequence
%                              is replayed every time.  Peer-reviewable.
%
%    Spatial decorrelation : each zone runs an independent LFSR (same
%                            polynomial, different seed) so cross-zone
%                            valve commands are orthogonal by construction.
%
%  SIGNAL DESIGN
% 
%  LFSR degree   n = 10  →  sequence period = 2^10 - 1 = 1023 clock ticks
%  Clock period  T_c = 4 steps = 1 h  (at Ts = 15 min)
%  Full period   1023 × 1 h ≈ 42.6 days  (> 4-week experiment: no wrap)
%
%  Two LFSR output bits per clock tick → 4 discrete amplitude levels:
%     (0,0) → v = 0.00   (off)
%     (0,1) → v = 0.33   (low)
%     (1,0) → v = 0.67   (medium)
%     (1,1) → v = 1.00   (full)
%
%  Primitive polynomial: x^10 + x^7 + 1  (taps at positions 7 and 10).
%  This polynomial is verified maximal-length over GF(2).
%
%  EXPERIMENT WINDOW (days 0-28 of the simulated year):
%    Four weeks from the start of the year, a period with moderate heating
%    demand in the Bordeaux climate.  Heating is active but unsaturated for
%    most of the day, so the valve signal produces visible temperature
%    responses at all four amplitude levels.  The first two weeks are the
%    identification set and the last two the cross-validation set.
%
%  Output:
%    data/multilevel_identification.csv   (weeks 1-2, days  0-14)
%    data/multilevel_validation.csv       (weeks 3-4, days 14-28)

%% Test-case parameters
url         = 'http://127.0.0.1';
testcase    = 'multizone_residential_hydronic';
step_period = 900;                               % 15 min per step  [s]
n_weeks     = 4;
max_steps   = n_weeks * 7 * 24 * 3600 / step_period;  % 2688 steps

start_time  = 0;                                % [s] from 1 January
warm_up     = 7  * 24 * 3600;                   % 1-week warm-up  [s]

%% MLPRS design parameters
clock_period = 4;                  % steps between LFSR clock ticks (= 1 h)
valve_levels = [0, 1/3, 2/3, 1];  % 4 amplitude levels  ∈ [0, 1]
nLevels      = numel(valve_levels);

%  LFSR: degree 10, primitive polynomial x^10 + x^7 + 1
%  Convention:  reg = [newest, ..., oldest]
%               output bit  = reg(end)
%               new bit     = reg(end) XOR reg(end-3)  (taps at indices end, end-3)
%  Period = 2^10 - 1 = 1023 clock ticks  ≈ 1023 h ≈ 42.6 days  (> 4 weeks)
LFSR_TAPS = [3, 10];   % 1-based positions in reg vector; polynomial x^10 + x^7 + 1

%  Per-zone initial LFSR seeds : all non-zero, mutually different.
%  Distinct seeds produce decorrelated sequences despite sharing the same polynomial.
zone_seeds = { ...
    [1,0,1,0,0,0,1,0,0,0], ...   % Liv
    [0,1,0,1,0,0,0,0,1,0], ...   % Ro1
    [1,1,0,0,0,1,0,0,0,1], ...   % Ro2
    [0,0,1,0,1,0,0,1,0,1], ...   % Ro3
    [1,0,0,1,0,1,0,1,0,0]  };    % Bth

%% Safety and boiler
T_min      = 287.15;   % K = 14°C, safety floor: open valve fully if too cold
T_max      = 296.15;   % K = 23°C, safety ceiling: close valve if too hot
%                                     (intervene early, prevents solver crash)
T_boi      = 303.15;   % K = 30°C, keeps boiler on without PI windup
%                                    (95°C = 368.15 K caused Modelica FMU
%                                     deadlock after ~720 steps via integral
%                                     windup in the boiler PI controller;
%                                     30°C is high enough that zones never
%                                     satisfy it, so the pump runs all the
%                                     time, identical to data_collection_prbs.m)
T_cool_off = 303.15;   % K, cooling disabled throughout

%% Zone and signal definitions
ctrl_zones   = {'Liv', 'Ro1', 'Ro2', 'Ro3', 'Bth'};
all_zones    = {'Liv', 'Ro1', 'Ro2', 'Ro3', 'Bth', 'Hal'};
nZones_ctrl  = numel(ctrl_zones);
nZones_all   = numel(all_zones);

zone_temp_sig = { ...
    'conHeaLiv_reaTZon_y', 'conHeaRo1_reaTZon_y', 'conHeaRo2_reaTZon_y', ...
    'conHeaRo3_reaTZon_y', 'conHeaBth_reaTZon_y', 'reaTHal_y' };

zone_heat_sig = { ...
    'reaHeaLiv_y', 'reaHeaRo1_y', 'reaHeaRo2_y', ...
    'reaHeaRo3_y', 'reaHeaBth_y' };

forecast_names = { ...
    'InternalGainsCon[Liv]', 'InternalGainsRad[Liv]', ...
    'InternalGainsCon[Ro1]', 'InternalGainsRad[Ro1]', ...
    'InternalGainsCon[Ro2]', 'InternalGainsRad[Ro2]', ...
    'InternalGainsCon[Ro3]', 'InternalGainsRad[Ro3]', ...
    'InternalGainsCon[Bth]', 'InternalGainsRad[Bth]', ...
    'InternalGainsCon[Hal]', 'InternalGainsRad[Hal]'  };
forecast_fields = cellfun(@(s) regexprep(s,'[\[\]]','_'), ...
                           forecast_names, 'UniformOutput', false);

%% Initialise LFSR state
lfsr_regs   = zone_seeds;             % one 10-bit register per zone
%  Stagger zone clocks so they never all tick on the same step.
%  Simultaneous valve changes in all 5 zones create a burst of Modelica
%  events that can time-out the solver.  Offset by 1 step per zone:
%    Liv ticks at steps 4,8,12,...   Ro1 at 5,9,13,...  etc.
clock_count = [0, 1, 2, 3, 0];       % staggered starting offsets
valve_idx   = ones(1, nZones_ctrl);   % current level index (1 = valve off)

%  Advance each zone's LFSR by its index × 1 extra tick so they start
%  at different points in their sequence even before the phase shift.
for k = 1:nZones_ctrl
    for t = 1:k
        [~, lfsr_regs{k}] = lfsr_next(lfsr_regs{k}, LFSR_TAPS);
        [~, lfsr_regs{k}] = lfsr_next(lfsr_regs{k}, LFSR_TAPS);
    end
end

%% Deploy testcase
handler = RequestHandler(url, 80);
tid     = handler.deploy_test(testcase);
handler.set_scenario(start_time, warm_up);
handler.set_step(step_period);

%% Logging
all_data     = table();
override_log = zeros(max_steps, nZones_ctrl);
prev_temps   = nan(1, nZones_all);

fprintf('MLPRS excitation: %d steps (%d weeks)  —  %s\n', max_steps, n_weeks, testcase);
fprintf('  LFSR: degree=10  polynomial=x^10+x^7+1  period=%d ticks (%.1f days)\n', ...
        2^10-1, (2^10-1)*clock_period*step_period/86400);
fprintf('  Clock: %d steps = %.0f h per tick\n', clock_period, clock_period*step_period/3600);
fprintf('  Levels: %s\n\n', mat2str(valve_levels,3));

%% Data-collection loop
for i = 1:max_steps

    %% 1. Disturbance forecast (one-step-ahead)
    fcast = handler.get_forecast(forecast_names, step_period, step_period);

    %% 2. Build override struct : direct valve commands, PI bypassed
    u = struct();
    for k = 1:nZones_ctrl
        zn    = ctrl_zones{k};
        v_cmd = valve_levels(valve_idx(k));
        u.(['conHea' zn '_oveActHea_activate']) = 1;
        u.(['conHea' zn '_oveActHea_u'])         = v_cmd;
        u.(['conCoo' zn '_oveTSetCoo_activate']) = 1;
        u.(['conCoo' zn '_oveTSetCoo_u'])         = T_cool_off;
    end
    %  Keep boiler and pump running at all times.
    u.oveTSetPumBoi_activate = 1;
    u.oveTSetPumBoi_u        = T_boi;

    %% 3. Advance simulation  (with retry on solver timeout)
    max_retries = 3;
    res = [];
    for retry = 1:max_retries
        try
            res = handler.advance(u);
            break;
        catch ME
            if retry < max_retries
                warning('[MLPRS] Step %d: advance timed out (%s). Retrying %d/%d in 10 s...', ...
                        i, ME.message, retry, max_retries);
                pause(10);
            else
                % All retries exhausted, try a zero-valve rescue step.
                % A problematic valve combination may have stalled the solver;
                % sending all-zero valves often lets it recover.
                warning('[MLPRS] Step %d: all retries failed. Trying zero-valve rescue step...', i);
                u_rescue = u;
                for k = 1:nZones_ctrl
                    zn = ctrl_zones{k};
                    u_rescue.(['conHea' zn '_oveActHea_u']) = 0;
                end
                pause(5);
                try
                    res = handler.advance(u_rescue);
                    warning('[MLPRS] Step %d: rescue step succeeded. Continuing.', i);
                    break;
                catch
                    % Truly stuck, save partial data and stop.
                    fprintf('[MLPRS] Rescue step also failed. Saving partial data...\n');
                    if height(all_data) > 0
                        mid_crash = floor(height(all_data)/2);
                        if mid_crash > 0
                            writetable(all_data(1:mid_crash,:),     'data/multilevel_identification.csv');
                            writetable(all_data(mid_crash+1:end,:), 'data/multilevel_validation.csv');
                            fprintf('[MLPRS] Partial data saved (id:%d / val:%d rows).\n', ...
                                    mid_crash, height(all_data)-mid_crash);
                        end
                    end
                    rethrow(ME);
                end
            end
        end
    end

    %% 4. Stall detection
    curr_temps = nan(1, nZones_all);
    for k = 1:nZones_all
        if isfield(res, zone_temp_sig{k})
            curr_temps(k) = res.(zone_temp_sig{k});
        end
    end
    if ~any(isnan(prev_temps)) && max(abs(curr_temps(1:nZones_ctrl) - prev_temps(1:nZones_ctrl))) < 1e-4
        warning('[MLPRS] Simulation frozen at step %d — forcing full heat.', i);
        valve_idx(:)   = nLevels;
        clock_count(:) = 0;
    end
    prev_temps = curr_temps;

    %% 5. Log row
    new_row = struct2table(res, 'AsArray', true);
    for k = 1:nZones_ctrl
        zn = ctrl_zones{k};
        new_row.(['valve_' zn]) = valve_levels(valve_idx(k));
        new_row.(['Qdot_'  zn]) = 0;   % placeholder
    end
    for k = 1:nZones_ctrl
        if isfield(res, zone_heat_sig{k})
            new_row.(['Qdot_' ctrl_zones{k}]) = res.(zone_heat_sig{k});
        end
    end
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

    %% 6. Update LFSR state : clock tick every clock_period steps
    for k = 1:nZones_ctrl
        clock_count(k) = clock_count(k) + 1;

        if clock_count(k) >= clock_period
            clock_count(k) = 0;

            %  Read two consecutive bits from the LFSR.
            %  Two bits → one of 4 patterns → 4 amplitude levels.
            %  Using consecutive bits preserves PRBS spectral flatness.
            [b1, lfsr_regs{k}] = lfsr_next(lfsr_regs{k}, LFSR_TAPS);
            [b2, lfsr_regs{k}] = lfsr_next(lfsr_regs{k}, LFSR_TAPS);
            %  Binary word to level index (1-based):
            %    (0,0)→1  (0,1)→2  (1,0)→3  (1,1)→4
            valve_idx(k) = 2*b1 + b2 + 1;
        end

        %  Safety floor: zone too cold → open valve fully until next tick.
        if ~isnan(curr_temps(k)) && curr_temps(k) < T_min
            valve_idx(k)      = nLevels;
            clock_count(k)    = 0;         % reset so we tick again soon
            override_log(i,k) = 1;
        end

        %  Safety ceiling: zone too hot → close valve until next tick.
        %  Prevents temperatures above T_max that destabilise the Modelica solver.
        if ~isnan(curr_temps(k)) && curr_temps(k) > T_max
            valve_idx(k)      = 1;         % level 1 = v=0 (off)
            clock_count(k)    = 0;
            override_log(i,k) = 1;
        end
    end

    %  Checkpoint every 336 steps (3.5 days), preserves data if later crash.
    if mod(i, 336) == 0 && height(all_data) > 0
        mid_cp = floor(height(all_data)/2);
        if mid_cp > 0
            writetable(all_data(1:mid_cp,:),     'data/multilevel_identification.csv');
            writetable(all_data(mid_cp+1:end,:), 'data/multilevel_validation.csv');
            fprintf('[MLPRS] Checkpoint saved at step %d (%d rows).\n', i, height(all_data));
        end
    end

    if mod(i, 48) == 0   % every 12 h
        vstr = arrayfun(@(k) sprintf('%.2f', valve_levels(valve_idx(k))), ...
                        1:nZones_ctrl, 'UniformOutput', false);
        tstr = arrayfun(@(t) sprintf('%.1f', t-273.15), curr_temps, 'UniformOutput', false);
        fprintf('Step %4d/%d | T=[%s]°C | v=[%s] | or=[%s]\n', ...
                i, max_steps, ...
                strjoin(tstr(1:nZones_ctrl),' '), ...
                strjoin(vstr,' '), ...
                sprintf('%.0f%% ', mean(override_log(1:i,:))*100));
    end
end

%% Override summary
fprintf('\n--- Safety-floor override rates (target < 5%%) ---\n');
for k = 1:nZones_ctrl
    fprintf('  %s: %.1f%%\n', ctrl_zones{k}, mean(override_log(:,k))*100);
end

%% Save
if ~exist('data', 'dir'), mkdir('data'); end
mid = floor(height(all_data) / 2);
writetable(all_data(1:mid,:),     'data/multilevel_identification.csv');
writetable(all_data(mid+1:end,:), 'data/multilevel_validation.csv');
fprintf('\nDone.  %d rows saved  (ident: %d | val: %d)\n', ...
        height(all_data), mid, height(all_data)-mid);

handler.stop_test(tid);

%% Local functions
function [bit, reg_out] = lfsr_next(reg, taps)
%LFSR_NEXT  Advance a Fibonacci LFSR one step.
%
%  reg     [1×n]   current shift-register state  (all elements ∈ {0,1})
%  taps    [1×t]   1-based tap positions for XOR feedback
%
%  bit              output bit = reg(end)  (oldest element)
%  reg_out [1×n]   updated register: [new_bit, reg(1:end-1)]
%
%  For x^10 + x^7 + 1 use taps = [3, 10].
%  In this convention reg(end) is the "oldest" (output) bit and
%  reg(1) is where the new feedback bit enters.
%  The feedback = XOR of the tap positions, which for taps=[3,10] gives
%  feedback = reg(3) XOR reg(10), matching the recurrence of x^10+x^7+1.

bit     = reg(end);
feedback = mod(sum(reg(taps)), 2);
reg_out  = [feedback, reg(1:end-1)];
end
