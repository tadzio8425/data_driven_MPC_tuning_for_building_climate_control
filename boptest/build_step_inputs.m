function [Tz_meas, dvec, lb_h] = build_step_inputs(res, lb_mpc, i, p, ...
                                Tamb_full, HGlo_full, Qint_full, n_total, SIG_TZ, nU, nYc)
%BUILD_STEP_INPUTS  Assemble the four inputs an MPC step function expects.
%
%  Shared by the RL and BO runners, so that all of them present the controller
%  with exactly the same measurement, disturbance window and comfort band.
%
%    Tz_meas  measured zone temperatures [K]
%    dvec     disturbance forecast stacked over the horizon, nD per step, with
%             the measured hallway temperature held constant across it
%    lb_h     comfort bound over the horizon: the occupied target during
%             occupancy, the raw setback schedule otherwise
%
%  i indexes the forecast arrays: the absolute step under a bulk prefetch, and 1
%  when the caller refetches a window at every step under forecast uncertainty.

    Tz_meas = zeros(p.nYm, 1);
    for zi = 1:nU, Tz_meas(zi) = res.(SIG_TZ{zi}); end
    T_Hal = res.(SIG_TZ{end});

    dvec = zeros(p.nD * p.N_p, 1);
    for k = 1:p.N_p
        ii = min(i + k - 1, n_total);
        dvec((k-1)*p.nD + (1:p.nD)) = [Tamb_full(ii); HGlo_full(ii);
                                        Qint_full(:,ii); T_Hal];
    end
    dvec(~isfinite(dvec)) = 0;

    lb_h = zeros(nYc * p.N_p, 1);
    w_target = p.w_occ * ones(nU, 1);
    for k = 1:p.N_p
        ii    = min(i + k - 1, size(lb_mpc, 1));
        lb_k  = lb_mpc(ii, :).';
        occ_k = double(lb_k > p.occ_sp_threshold);
        lb_h((k-1)*nYc + (1:nYc)) = occ_k .* w_target + (1 - occ_k) .* lb_k;
    end
end
