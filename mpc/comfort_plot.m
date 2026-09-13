function fig = comfort_plot(name, log_Tz, log_ref, log_ub, log_Tamb, kpi, ZONE_NAMES, Ts, varargin)
%COMFORT_PLOT  The standard comfort figure, drawn the same way for every controller.
%
%  A 3x2 grid: one panel per controlled zone with its comfort band, plus a last
%  panel holding the uncontrolled hallway and the outdoor temperature. The
%  BOPTEST KPIs of the run go in the title, so a figure always carries the row
%  it belongs to in the KPI tables.
%
% INPUTS
%   name        controller key; known keys map to the label used in the thesis
%   log_Tz      N x nZ zone temperatures [K], uncontrolled zones in the last columns
%   log_ref     N x nU lower setpoint [K] (raw schedule, uncapped)
%   log_ub      N x nU upper setpoint [K] (raw schedule, uncapped)
%   log_Tamb    N x 1  outdoor temperature [K]
%   kpi         BOPTEST KPI struct, shown in the title
%   ZONE_NAMES  cell array of zone labels, length nZ
%   Ts          sample time [s]
%
% OPTIONAL name-value pairs
%   'SaveTo'    file path for PNG export, '' disables
%   'Position'  figure position [x y w h]

    ip = inputParser;
    addParameter(ip, 'SaveTo',   '',                @(x) ischar(x)||isstring(x));
    addParameter(ip, 'Position', [20 20 1300 780],  @(x) isnumeric(x)&&numel(x)==4);
    parse(ip, varargin{:});
    opts = ip.Results;

    ZCOL = [0.12 0.47 0.71;
            0.86 0.35 0.19;
            0.27 0.63 0.27;
            0.58 0.40 0.74;
            0.80 0.47 0.65;
            0.50 0.50 0.50];

    C_LB   = [0.25 0.25 0.25];   % dark grey:  lower-setpoint reference
    C_BAND = [0.78 0.94 0.78];   % light green: comfort band
    C_AMB  = [0.85 0.33 0.10];   % orange:      ambient
    C_HAL  = [0.30 0.60 0.90];   % steel blue:  hallway (uncontrolled)

    nU      = size(log_ref, 2);   % controlled zones
    nZ_log  = size(log_Tz,  2);   % columns present in the temperature log
    n_steps = size(log_Tz,  1);
    t_days  = (1:n_steps).' * Ts / 86400;

    T_C  = log_Tz(:, 1:nU) - 273.15;
    lb_C = log_ref - 273.15;
    ub_C = log_ub  - 273.15;

    label = display_name(name);
    fig = figure('Name', sprintf('%s — temperature', label), 'Color', 'w', ...
                 'Position', opts.Position);

    %% 1. One panel per controlled zone
    for j = 1:min(nU, 5)
        subplot(3, 2, j);
        hold on; grid on; box on;
        fill([t_days; flipud(t_days)], [lb_C(:,j); flipud(ub_C(:,j))], C_BAND, ...
             'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Comfort band');
        plot(t_days, lb_C(:,j), '-', 'Color', C_LB, 'LineWidth', 1.0, ...
             'DisplayName', 'Lower setpoint');
        plot(t_days, T_C(:,j),  '-', 'Color', ZCOL(j,:), 'LineWidth', 1.1, ...
             'DisplayName', strrep(ZONE_NAMES{j},'_','\_'));
        ylabel('Temperature [°C]');
        title(sprintf('\\bf%s', strrep(ZONE_NAMES{j},'_','\_')), 'FontWeight','normal');
        if j > 3, xlabel('Time [days]'); end
        legend('Location','best','FontSize',7,'Box','off');

        % The y range follows the zone temperature, not the band: the wide
        % unoccupied-hour band would otherwise flatten the trace.
        ylim([min(min(T_C(:,j)), min(lb_C(:,j))) - 1, ...
              max(max(T_C(:,j)), max(lb_C(:,j))) + 1]);
    end

    %% 2. Last panel: the hallway, plus the ambient temperature on the right axis
    subplot(3, 2, 6);
    yyaxis left;
    if nZ_log > nU
        hal = strrep(ZONE_NAMES{min(nU+1,end)},'_','\_');
        plot(t_days, log_Tz(:, nU+1) - 273.15, '-', 'Color', C_HAL, ...
             'LineWidth', 1.0, 'DisplayName', hal);
        ylabel(sprintf('T_{%s} [°C]', hal));
        side_title = sprintf('\\bf%s\\rm  (uncontrolled)', hal);
    else
        plot(t_days, T_C(:, end), '-', 'Color', ZCOL(min(nU,6),:), ...
             'LineWidth', 1.0, 'DisplayName', strrep(ZONE_NAMES{nU},'_','\_'));
        ylabel('Temperature [°C]');
        side_title = '';
    end
    yyaxis right;
    plot(t_days, log_Tamb - 273.15, '-.', 'Color', C_AMB, 'LineWidth', 0.8, ...
         'DisplayName','T_{amb}');
    ylabel('Outdoor temperature [°C]');
    xlabel('Time [days]'); grid on; box on;
    title([side_title '  +  ambient'], 'FontWeight','normal');
    legend('Location','best','FontSize',7,'Box','off');

    sgtitle(sprintf('%s   |   BOPTEST KPIs:  %s', label, kpi_str(kpi)), ...
            'FontWeight','bold');

    if ~isempty(opts.SaveTo)
        exportgraphics(fig, opts.SaveTo, 'Resolution', 150);
        fprintf('[comfort_plot] saved → %s\n', opts.SaveTo);
    end
end


%% Helpers
function s = display_name(name)
%DISPLAY_NAME  Map a controller key to its figure label. These strings appear
%  in the headers of the figures reproduced in the thesis, so they are left as
%  they were when those figures were exported.
    map = struct( ...
        'pi_baseline',       'PI baseline', ...
        'mpc_baseline',      'MPC (no bias)', ...
        'ekf_mpc',           'Offset-free MPC (AKF)', ...
        'rl_mpc',            'RL-tuned MPC (training episode)', ...
        'rl_mpc_bias',       'RL-tuned offset-free MPC (training episode)', ...
        'rl_mpc_bench',      'RL-tuned MPC (greedy eval)', ...
        'rl_mpc_bench_bias', 'RL-tuned offset-free MPC (greedy eval)', ...
        'bo_mpc',            'BO-MPC', ...
        'bo_mpc_bias',       'BO-MPC (offset-free)');
    key = matlab.lang.makeValidName(char(name));
    if isfield(map, key)
        s = map.(key);
    else
        s = char(name);
    end
end

function s = kpi_str(kpi)
%KPI_STR  Headline string with the four reported KPIs, in TeX-safe math labels.
    flds   = {'tdis_tot',    'ener_tot', 'cost_tot', 'emis_tot'};
    labels = {'T_{dis,tot}', 'E_{tot}',  'C_{tot}',  'M_{CO_2,tot}'};
    units  = {'K·h/zone',    'kWh/m²',   '$/m²',     'kgCO_2/m²'};
    s = '';
    for k = 1:numel(flds)
        if isfield(kpi, flds{k})
            s = sprintf('%s  %s=%.3f %s ', s, labels{k}, kpi.(flds{k}), units{k});
        end
    end
    s = strtrim(s);
end
