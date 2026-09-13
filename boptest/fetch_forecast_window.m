function [Tamb_w, HGlo_w, Qint_w, lb_w, ub_w] = fetch_forecast_window( ...
    handler, dist_sigs, setp_sigs, ZONE_NAMES, nU, N_p, Ts)
%FETCH_FORECAST_WINDOW  Per-step forecast fetch used in noisy-forecast runs.
%
%  BOPTEST caps the /forecast horizon at 48 h whenever forecast uncertainty is
%  active, so the whole-episode bulk prefetch is not available. This helper
%  fetches only the N_p-step lookahead the QP needs, 6 h by default. The
%  returned arrays keep the shapes of the bulk arrays sliced 1:N_p, so the
%  runners can swap the helper in without touching any downstream indexing:
%    Tamb_w  N_p × 1     HGlo_w  N_p × 1     Qint_w  nU × N_p
%    lb_w    N_p × nU    ub_w    N_p × nU    (raw setpoint schedule, K)
%
%  T_Hal is a measurement, not a forecast, and is read from the advance response.
%  Cost is 2 HTTP calls per MPC step.

    % 1. Disturbance window
    fc_dist = handler.get_forecast(dist_sigs, (N_p - 1) * Ts, Ts);
    Tamb_w  = fc_dist.TDryBul(1:N_p);
    HGlo_w  = fc_dist.HGloHor(1:N_p);

    Qint_w = zeros(nU, N_p);
    for j = 1:nU
        zn = ZONE_NAMES{j};
        qc = fc_dist.(sprintf('InternalGainsCon_%s_', zn))(1:N_p);
        qr = fc_dist.(sprintf('InternalGainsRad_%s_', zn))(1:N_p);
        Qint_w(j, :) = qc(:).' + qr(:).';
    end

    % 2. Schedule window
    fc_setp = handler.get_forecast(setp_sigs, (N_p - 1) * Ts, Ts);
    lb_w = zeros(N_p, nU);
    ub_w = zeros(N_p, nU);
    for j = 1:nU
        zn = ZONE_NAMES{j};
        lb_w(:, j) = fc_setp.(sprintf('LowerSetp_%s_', zn))(1:N_p);
        ub_w(:, j) = fc_setp.(sprintf('UpperSetp_%s_', zn))(1:N_p);
    end
end
