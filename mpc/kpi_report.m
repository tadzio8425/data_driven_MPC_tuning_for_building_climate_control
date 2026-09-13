function metrics = kpi_report(name, kpi, log_Tz, log_Qhea, log_ref, log_ub, ZONE_NAMES, Ts, occ_sp_threshold)
%KPI_REPORT  Standard printout: the official BOPTEST KPIs plus per-zone detail.
%
%  The kpi.* fields are the headline numbers and the ones the thesis tables
%  quote. The per-zone block below them is computed here from the raw schedule,
%  so it does not reproduce kpi.tdis_tot: BOPTEST applies its own occupancy
%  definition and normalises per zone. It is kept because the headline KPI is a
%  single number and the per-zone split is what identifies which zone is
%  responsible for it.
%
% INPUTS
%   name              controller label, e.g. 'pi_baseline'
%   kpi               struct from handler.get_kpi()
%   log_Tz            N x nZ zone temperatures [K], controlled zones first
%   log_Qhea          N x nU heat delivered per controlled zone [W]
%   log_ref           N x nU lower setpoint schedule [K] (raw, uncapped)
%   log_ub            N x nU upper setpoint schedule [K] (raw, uncapped)
%   ZONE_NAMES        cell array of zone labels, controlled zones first
%   Ts                sample time [s]
%   occ_sp_threshold  a step counts as occupied when lb >= threshold [K]
%
% OUTPUT
%   metrics.kpi, .tdis_under, .tdis_over, .E_zone, .tdis_tot_local, .E_tot_local

    nU       = size(log_ref, 2);
    occ_mask = log_ref >= occ_sp_threshold;   % >= so a lb exactly at the threshold counts

    %% 1. Per-zone discomfort and energy
    tdis_under = zeros(nU, 1);
    tdis_over  = zeros(nU, 1);
    E_zone     = zeros(nU, 1);
    for j = 1:nU
        T_j = log_Tz(:, j);
        oj  = occ_mask(:, j);
        tdis_under(j) = sum(max(0, log_ref(:,j) - T_j) .* oj) * Ts/3600;
        tdis_over(j)  = sum(max(0, T_j - log_ub(:,j))  .* oj) * Ts/3600;
        E_zone(j)     = sum(log_Qhea(:, j)) * Ts/3.6e6;
    end
    tdis_tot_local = sum(tdis_under) + sum(tdis_over);
    E_tot_local    = sum(E_zone);

    %% 2. BOPTEST KPI block
    fprintf('\n══════════════════════════════════════════════════════\n');
    fprintf('  Controller: %s\n', name);
    fprintf('══════════════════════════════════════════════════════\n');
    fprintf('  BOPTEST KPIs (the numbers the thesis tables quote):\n');
    flds  = {'tdis_tot','idis_tot','ener_tot','cost_tot','emis_tot','time_rat'};
    units = {'K·h/zone','ppmh/zone','kWh/m²','$/m²','kgCO₂/m²','-'};
    for k = 1:numel(flds)
        if isfield(kpi, flds{k})
            fprintf('    %-10s = %10.4f   [%s]\n', flds{k}, kpi.(flds{k}), units{k});
        end
    end

    %% 3. Per-zone diagnostics
    fprintf('\n  Per-zone detail (occupied steps, raw schedule):\n');
    fprintf('  %-6s  %12s  %12s  %12s  %12s\n', ...
            'Zone', 'Tdis_under', 'Tdis_over', 'Tdis_tot', 'E [kWh]');
    for j = 1:nU
        fprintf('  %-6s  %12.2f  %12.2f  %12.2f  %12.2f\n', ZONE_NAMES{j}, ...
                tdis_under(j), tdis_over(j), tdis_under(j)+tdis_over(j), E_zone(j));
    end
    fprintf('  %-6s  %12.2f  %12.2f  %12.2f  %12.2f\n\n', ...
            'TOTAL', sum(tdis_under), sum(tdis_over), tdis_tot_local, E_tot_local);

    %% 4. Pack the output
    metrics.kpi            = kpi;
    metrics.tdis_under     = tdis_under;
    metrics.tdis_over      = tdis_over;
    metrics.E_zone         = E_zone;
    metrics.tdis_tot_local = tdis_tot_local;
    metrics.E_tot_local    = E_tot_local;
end
