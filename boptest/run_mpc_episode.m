function L_total = run_mpc_episode(x, use_bias, BOPTEST_TUNC, BOPTEST_SUNC, start_time_in)
%RUN_MPC_EPISODE  The black-box objective the Bayesian optimisation minimises.
%
%  Deploys a fresh seven-day BOPTEST episode, runs the MPC with the cost
%  parameters bayesopt proposes, and returns the realised closed-loop cost
%
%    J_real = w_E·E_tot [kWh] + w_d·Tdis_under [K·h] + w_o·Tdis_over [K·h]
%
%  accumulated over the week and the five zones, with all three weights at one.
%  Nothing in the search reads the model: the only signal returned is what the
%  plant actually cost.
%
%  INPUTS
%    x             1x1 table from bayesopt with fields
%                    q_Liv q_Ro1 q_Ro2 q_Ro3 q_Bth   occupied tracking weights
%                    r_u                             shared input penalty
%                    w_Liv w_Ro1 w_Ro2 w_Ro3 w_Bth   reference offsets [K]
%    use_bias      bias states active in the KF (default false, see bo_mpc_tune)
%    BOPTEST_TUNC  temperature forecast noise: '' | 'medium' | 'high'
%    BOPTEST_SUNC  the same options for solar irradiation
%    start_time_in scenario start [s], 0 for the week the search is run on
%
%  OUTPUT
%    L_total  realised J_real over the episode, Inf if the BOPTEST run fails, so
%             that a failed evaluation costs the search one point and not the run
%
%  q_unocc and the three slack penalties stay at their defaults from
%  tracking_mpc_params.mat: only the eleven variables above are optimised.

if nargin < 2 || isempty(use_bias),      use_bias      = true; end
if nargin < 3 || isempty(BOPTEST_TUNC),  BOPTEST_TUNC  = '';   end
if nargin < 4 || isempty(BOPTEST_SUNC),  BOPTEST_SUNC  = '';   end
if nargin < 5 || isempty(start_time_in), start_time_in = 0;    end

%% 1. Real-cost weights (identical to rl_mpc_loop)
w_E = 1.0;
w_d = 1.0;
w_o = 1.0;

%% 2. Load the model and the precomputed QP matrices
load('models/rc_params_2C6z.mat', 'ZONE_NAMES', 'Ts', 'nZ');
p   = load('models/tracking_mpc_params.mat');
nU  = p.nU;   nYc = p.nYc;
nZb = nZ;     % 6 BOPTEST zones, 5 controlled plus Hal

clear rl_mpc_step   % reset the persistent parameter cache between episodes

%% 3. Build θ from the BO table row
theta.q_diag   = [[x.q_Liv; x.q_Ro1; x.q_Ro2; x.q_Ro3; x.q_Bth], ...
                  diag(p.Q_tilde(1:nYc, 1:nYc))];   % occupied from BO, unoccupied default
theta.r_u      = [x.r_u; x.r_u];                    % one penalty, occupied and unoccupied
theta.w_offset = [x.w_Liv; x.w_Ro1; x.w_Ro2; x.w_Ro3; x.w_Bth];
theta.g1_over  = p.g1_over;
theta.g1_under = p.g1_under;
theta.g2       = p.g2;

%% 4. Signal names
SIG_TZ  = {'conHeaLiv_reaTZon_y','conHeaRo1_reaTZon_y','conHeaRo2_reaTZon_y', ...
           'conHeaRo3_reaTZon_y','conHeaBth_reaTZon_y','reaTHal_y'};
SIG_QH  = {'reaHeaLiv_y','reaHeaRo1_y','reaHeaRo2_y','reaHeaRo3_y','reaHeaBth_y'};
SIG_CON = {'InternalGainsCon[Liv]','InternalGainsCon[Ro1]','InternalGainsCon[Ro2]', ...
           'InternalGainsCon[Ro3]','InternalGainsCon[Bth]','InternalGainsCon[Hal]'};
SIG_RAD = {'InternalGainsRad[Liv]','InternalGainsRad[Ro1]','InternalGainsRad[Ro2]', ...
           'InternalGainsRad[Ro3]','InternalGainsRad[Bth]','InternalGainsRad[Hal]'};

%% 5. Scenario
start_time = start_time_in;
warm_up    = 7*24*3600;
n_days     = 7;
n_steps    = round(n_days * 86400 / Ts);

%% 6. Run the episode, wrapped so that bayesopt never crashes on a failed one
try
    handler = RequestHandler('http://127.0.0.1', 80);
    tid     = handler.deploy_test('multizone_residential_hydronic');
    handler.set_scenario(start_time, warm_up, BOPTEST_TUNC, BOPTEST_SUNC);
    handler.set_step(Ts);
    res = handler.advance(struct());

    USE_NOISY_FCAST = ~isempty(BOPTEST_TUNC) || ~isempty(BOPTEST_SUNC);
    n_total = n_steps + p.N_p;

    setp_sigs = cell(1, 2*nU);
    for j = 1:nU
        setp_sigs{2*j-1} = sprintf('LowerSetp[%s]', ZONE_NAMES{j});
        setp_sigs{2*j  } = sprintf('UpperSetp[%s]', ZONE_NAMES{j});
    end
    dist_sigs = [{'TDryBul','HGloHor'}, SIG_CON(1:nU), SIG_RAD(1:nU)];

    if USE_NOISY_FCAST
        Tamb_full = zeros(p.N_p, 1);      % refilled each step by fetch_forecast_window
        HGlo_full = zeros(p.N_p, 1);
        Qint_full = zeros(nU, p.N_p);
        lb_mpc    = zeros(p.N_p, nU);
        ub_mpc    = zeros(p.N_p, nU);
    else
        lb_mpc  = zeros(n_total, nU);
        ub_mpc  = zeros(n_total, nU);
        fc_setp = handler.get_forecast(setp_sigs, n_total*Ts, Ts);
        for j = 1:nU
            lb_mpc(:,j) = fc_setp.(sprintf('LowerSetp_%s_', ZONE_NAMES{j}))(1:n_total);
            ub_mpc(:,j) = fc_setp.(sprintf('UpperSetp_%s_', ZONE_NAMES{j}))(1:n_total);
        end

        fc_dist   = handler.get_forecast(dist_sigs, n_total*Ts, Ts);
        Tamb_full = fc_dist.TDryBul(1:n_total);
        HGlo_full = fc_dist.HGloHor(1:n_total);
        Qint_full = zeros(nU, n_total);
        for j = 1:nU
            qc = fc_dist.(regexprep(SIG_CON{j},'[\[\]]','_'))(1:n_total);
            qr = fc_dist.(regexprep(SIG_RAD{j},'[\[\]]','_'))(1:n_total);
            Qint_full(j,:) = qc(:).' + qr(:).';
        end
    end

    Tz_all = zeros(nZb, 1);
    for i = 1:nZb, Tz_all(i) = res.(SIG_TZ{i}); end
    z_aug  = [Tz_all(1:nU); Tz_all(1:nU) + 1.5; zeros(p.nd_bias, 1)];
    u_prev = zeros(nU, 1);

    L_total = 0;

    for i = 1:n_steps

        if USE_NOISY_FCAST
            [Tamb_full, HGlo_full, Qint_w, lb_mpc, ub_mpc] = fetch_forecast_window( ...
                handler, dist_sigs, setp_sigs, ZONE_NAMES, nU, p.N_p, Ts);
            Qint_full(1:nU, :) = Qint_w;
            i_loc = 1;
            n_loc = p.N_p;
        else
            i_loc = i;
            n_loc = n_total;
        end

        [Tz_meas, dvec, lb_h] = build_step_inputs( ...
            res, lb_mpc, i_loc, p, Tamb_full, HGlo_full, Qint_full, n_loc, SIG_TZ, nU, nYc);

        % Greedy MPC step, no learning
        [v_cmd, ~, z_aug_k] = ...
            rl_mpc_step(Tz_meas, dvec, lb_h, z_aug, u_prev, theta, [], use_bias);

        ov = struct('oveTSetPumBoi_activate', 1, 'oveTSetPumBoi_u', p.w_occ);
        for j = 1:nU
            zn = ZONE_NAMES{j};
            ov.(sprintf('conHea%s_oveActHea_activate',  zn)) = 1;
            ov.(sprintf('conHea%s_oveActHea_u',         zn)) = v_cmd(j);
            ov.(sprintf('conHea%s_oveTSetHea_activate', zn)) = 0;
        end
        res_next = handler.advance(ov);

        Q_heat  = zeros(nU, 1);
        Tz_next = zeros(nU, 1);
        for j = 1:nU
            Q_heat(j)  = res_next.(SIG_QH{j});
            Tz_next(j) = res_next.(SIG_TZ{j});
        end

        % Accumulate the realised cost
        lb_now  = lb_mpc(min(i_loc, end), :).';
        ub_now  = ub_mpc(min(i_loc, end), :).';
        occ     = double(lb_now > p.occ_sp_threshold);
        L_total = L_total ...
                + w_E * sum(Q_heat) * Ts/3.6e6 ...
                + w_d * sum(occ .* max(0, lb_now - Tz_next)) * Ts/3600 ...
                + w_o * sum(occ .* max(0, Tz_next - ub_now)) * Ts/3600;

        z_aug  = z_aug_k;
        u_prev = Q_heat;
        res    = res_next;
    end

    handler.stop_test(tid);

catch ME
    warning('[run_mpc_episode] Episode failed: %s', ME.message);
    try
        handler.stop_test(tid);
    catch
    end
    L_total = Inf;
    return
end

fprintf('[BO]  q=[%5.1f %5.1f %5.1f %5.1f %5.1f]  r_u=%6.3g  w=[%+.2f %+.2f %+.2f %+.2f %+.2f]  J_real=%8.3f\n', ...
        x.q_Liv, x.q_Ro1, x.q_Ro2, x.q_Ro3, x.q_Bth, x.r_u, ...
        x.w_Liv, x.w_Ro1, x.w_Ro2, x.w_Ro3, x.w_Bth, L_total);
end
