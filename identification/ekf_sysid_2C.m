%% ekf_sysid_2C.m: Grey-box EKF identification of the 2C zone model
%
% Estimates the RC parameters of the five controlled zones of the BOPTEST
% multizone_residential_hydronic testcase from the multilevel excitation
% experiment in data/multilevel_identification.csv, and validates them on
% data/multilevel_validation.csv.
%
% It also exports every figure of Chapter 3 into figures/.
%
% BEFORE RUNNING:
%   1. startup
% Run ONCE, before mpc/tracking_mpc_params_build.
% Saves  models/rc_params_2C6z.mat
%
%% Physical model (2C per zone, 5 zones)
%
%  Each zone has TWO lumped thermal masses:
%    T1_i    air temperature          (fast dynamics, MEASURED)
%    T2_i    solid mass / slab temp   (slow dynamics, UNMEASURED)
%
%  For each modelled zone i in {Liv, Ro1, Ro2, Ro3, Bth}:
%
%    C1_i * dT1_i/dt = UA_ext_i   * (Tamb  - T1_i)           [exterior envelope]
%                    + UA_toHal_i * (T_Hal - T1_i)           [hallway exchange]
%                    + sum_{j in adj5} UA_ij * (T1_j - T1_i) [direct zone coupling]
%                    + UA_int_i   * (T2_i  - T1_i)           [slab->air exchange]
%                    + Qdot_hea_i + Qsol_i + Qint_i          [heat sources]
%
%    C2_i * dT2_i/dt = UA_int_i * (T1_i - T2_i)             [air->slab exchange]
%
%  where adj5 = direct zone-to-zone walls EXCLUDING Hal:
%    Liv-Ro1,  Liv-Bth,  Ro1-Ro2,  Ro3-Bth   (4 shared walls)
%
%  Key structural properties:
%    * UA_ij = UA_ji     (symmetry, one parameter per shared wall)
%    * UA_toHal_i > 0    (all modelled zones border the hallway)
%    * UA_int_i > 0      (slab always exchanges with its air node)
%    * C1_i, C2_i > 0   (positivity enforced via EKF bounds)
%    * Only T1 (air) is measured . T2 (slab) is hidden, estimated by EKF
%
%% Parameter vector (29 parameters)
%
%  C1       [J/K]    5  air capacitances (fast mass)
%  C2       [J/K]    5  slab capacitances (slow mass)
%  UA_ext   [W/K]    5  exterior envelope conductances
%  UA_int   [W/K]    5  internal air-slab conductances
%  UA_toHal [W/K]    5  conductances to the hallway
%  UA_c     [W/K]    4  direct zone-to-zone coupling conductances
%                   ----
%                   29  total
%
%% EKF augmented state (nz = 39)
%
%    z = [ T1_Liv .. T1_Bth   <- indices  1: 5   air temps      (MEASURED)
%          T2_Liv .. T2_Bth   <- indices  6:10   slab temps     (HIDDEN)
%          C1_Liv .. C1_Bth   <- indices 11:15   air capacitances
%          C2_Liv .. C2_Bth   <- indices 16:20   slab capacitances
%          UAe_Liv .. UAe_Bth <- indices 21:25   exterior conductances
%          UAi_Liv .. UAi_Bth <- indices 26:30   internal air-slab conductances
%          UAh_Liv .. UAh_Bth <- indices 31:35   conductances to Hal
%          UAc(1) .. UAc(4) ] <- indices 36:39   direct zone coupling
%
%  Temperatures: Euler-integrated each step.
%  Parameters:   random-walk  p(k+1) = p(k) + noise.
%
%% Analytical Jacobian
%
%  Defining:
%    UA_sum1_i = UA_ext_i + UA_toHal_i + sum_j UA_ij + UA_int_i
%    f1_i = (UA_ext_i*(Tamb-T1_i) + UA_toHal_i*(T_Hal-T1_i)
%            + sum_j UA_ij*(T1_j-T1_i) + UA_int_i*(T2_i-T1_i) + Q_i) / C1_i
%    f2_i = UA_int_i*(T1_i - T2_i) / C2_i
%
%  Non-zero Jacobian entries for air node row i:
%    dT1_i_next/dT1_i       = 1 - Ts*UA_sum1_i/C1_i
%    dT1_i_next/dT1_j       = Ts*UA_ij/C1_i              (direct adjacent j)
%    dT1_i_next/dT2_i       = Ts*UA_int_i/C1_i
%    dT1_i_next/dC1_i       = -Ts*f1_i/C1_i
%    dT1_i_next/dUA_ext_i   = Ts*(Tamb-T1_i)/C1_i
%    dT1_i_next/dUA_int_i   = Ts*(T2_i-T1_i)/C1_i
%    dT1_i_next/dUA_toHal_i = Ts*(T_Hal-T1_i)/C1_i
%    dT1_i_next/dUA_ij      = Ts*(T1_j-T1_i)/C1_i
%
%  Non-zero Jacobian entries for slab node row i:
%    dT2_i_next/dT1_i     = Ts*UA_int_i/C2_i
%    dT2_i_next/dT2_i     = 1 - Ts*UA_int_i/C2_i
%    dT2_i_next/dC2_i     = -Ts*f2_i/C2_i
%    dT2_i_next/dUA_int_i = Ts*(T1_i-T2_i)/C2_i
%
%  Parameter rows: identity (random walk).
%

clear; clc; close all;

%% Section 1. Building topology
%  Edit this section if the floor plan changes.

% Modelled zones only, Hal is not a state
ZONE_NAMES = {'Liv','Ro1','Ro2','Ro3','Bth'};
nZ = numel(ZONE_NAMES);   % 5

iLiv=1; iRo1=2; iRo2=3; iRo3=4; iBth=5;

% Direct zone-to-zone shared walls (Hal excluded, handled via UA_toHal).
% ORDER MATTERS: defines UAc indices. Do not reorder.
ADJ_PAIRS = [iLiv, iRo1;   % 1  UA_Liv-Ro1
             iLiv, iBth;   % 2  UA_Liv-Bth
             iRo1, iRo2;   % 3  UA_Ro1-Ro2
             iRo3, iBth];  % 4  UA_Ro3-Bth
nP = size(ADJ_PAIRS, 1);   % 4 direct coupling conductances

% All modelled zones have exterior walls
HAS_EXT_WALL = logical([1, 1, 1, 1, 1]);
nE = sum(HAS_EXT_WALL);    % 5
EXT_ZONE_IDX = find(HAS_EXT_WALL);   % [1,2,3,4,5]

% All modelled zones border Hal, one UA_toHal per zone
nH = nZ;   % 5

% Solar gain coefficients: A_win * g_value [m^2]
ALPHA_SOL = [4.0, 1.5, 1.5, 1.0, 0.5] * 0.6;

% BOPTEST signal names, 5 modelled zones only
SIG_TZ  = {'conHeaLiv_reaTZon_y','conHeaRo1_reaTZon_y','conHeaRo2_reaTZon_y', ...
            'conHeaRo3_reaTZon_y','conHeaBth_reaTZon_y'};
SIG_QH  = {'reaHeaLiv_y','reaHeaRo1_y','reaHeaRo2_y','reaHeaRo3_y','reaHeaBth_y'};
SIG_CON = {'fcast_InternalGainsCon_Liv_','fcast_InternalGainsCon_Ro1_', ...
           'fcast_InternalGainsCon_Ro2_','fcast_InternalGainsCon_Ro3_', ...
           'fcast_InternalGainsCon_Bth_'};
SIG_RAD = {'fcast_InternalGainsRad_Liv_','fcast_InternalGainsRad_Ro1_', ...
           'fcast_InternalGainsRad_Ro2_','fcast_InternalGainsRad_Ro3_', ...
           'fcast_InternalGainsRad_Bth_'};
SIG_THAL = 'reaTHal_y';   % Hal temperature (measured)

% State vector index ranges
IDX_T1  = 1                  : nZ;     %  1: 5
IDX_T2  = nZ               + (1:nZ);   %  6:10
IDX_C1  = 2*nZ              + (1:nZ);  % 11:15
IDX_C2  = 3*nZ              + (1:nZ);  % 16:20
IDX_UAe = 4*nZ              + (1:nE);  % 21:25
IDX_UAi = 4*nZ+nE           + (1:nZ);  % 26:30
IDX_UAh = 4*nZ+nE+nZ        + (1:nH);  % 31:35
IDX_UAc = 4*nZ+nE+nZ+nH     + (1:nP);  % 36:39
nz      = 4*nZ + nE + nZ + nH + nP;   %     39  total
%  ASOL is NOT in the EKF state . it is a fixed physical constant identified
%  via a separate batch optimisation (fminsearch) after the EKF converges.

fprintf('2C model: %d modelled zones\n', nZ);
fprintf('Direct shared walls (excl. Hal): %d\n', nP);
fprintf('EKF augmented state dimension: nz = %d\n', nz);
fprintf('Parameters: C1x%d + C2x%d + UAex%d + UAix%d + UAhx%d + UAcx%d = %d total\n\n', ...
        nZ, nZ, nE, nZ, nH, nP, nZ+nZ+nE+nZ+nH+nP);

%% Section 2. Initial parameter guesses

% Air capacitances C1 [J/K]
% Air + light furniture: rho_air*cp*V * 2-3x factor
%   Liv~50m^3, Ro1/Ro2~30m^3, Ro3~20m^3, Bth~15m^3
C1_0 = [500, 300, 420, 100, 80] * 1e3;   % [J/K]

% Slab capacitances C2 [J/K]
% Hydronic concrete floor slab: rho*cp*thickness*A_floor
%   rho~2300 kg/m^3, cp=880 J/(kg·K), ~5cm thick
%   Liv ~16m^2 -> 1.62 MJ/K, bedrooms ~10m^2 -> ~1 MJ/K
C2_0 = [3.058, 1.5, 1.5, 1.2, 0.8] * 1e6;   % [J/K]

% Exterior conductances UA_ext [W/K]
% U_wall * A_ext, insulated wall: 0.3-0.8 W/(m^2·K)
UAe_0 = [40, 12, 15, 15, 8] * 1.0;   % [W/K]

% Internal air-slab conductances UA_int [W/K]
% Convective transfer floor->air: h~5 W/(m^2·K)
%   Liv: 5*16=80 W/K, bedrooms: 5*10=50 W/K
UAi_0 = [80, 50, 50, 45, 35] * 1.0;   % [W/K]

% Conductances to Hal UA_toHal [W/K]
% Each zone shares a wall/door opening with Hal.
UAh_0 = [60, 30, 50, 40, 60] * 1.0;   % [W/K]

% Direct zone-to-zone coupling conductances UA_c [W/K]
% [Liv-Ro1, Liv-Bth, Ro1-Ro2, Ro3-Bth]
UAc_0 = [40, 25, 25, 25] * 1.0;   % [W/K]

theta0 = [C1_0, C2_0, UAe_0, UAi_0, UAh_0, UAc_0];   % 29 parameters

%% Section 3. EKF noise tuning

R_meas = 0.1^2;    % measurement noise [K^2]

Q_T1 = 0.01^2;    % air temperature process noise [K^2/step]
Q_T2 = 0.001^2;   % slab temperature process noise (10x smaller)

Q_C1  = (50)^2;     % [J/K]^2
Q_C2  = (600)^2;    % [J/K]^2
Q_UAe = (0.01)^2;   % [W/K]^2
Q_UAi = (0.01)^2;   % [W/K]^2
Q_UAh = (0.01)^2;   % [W/K]^2
Q_UAc = (0.01)^2;   % [W/K]^2

Q_ekf = diag([ Q_T1 * ones(1,nZ),  ...
               Q_T2 * ones(1,nZ),  ...
               Q_C1 * ones(1,nZ),  ...
               Q_C2 * ones(1,nZ),  ...
               Q_UAe* ones(1,nE),  ...
               Q_UAi* ones(1,nZ),  ...
               Q_UAh* ones(1,nH),  ...
               Q_UAc* ones(1,nP) ]);

%% Section 4. Load data

data_id  = readtable('data/multilevel_identification.csv');   % days  0-14
data_val = readtable('data/multilevel_validation.csv');       % days 14-28
Ts = data_id.time(2) - data_id.time(1);   % 900 s

function D = extract_signals(data, SIG_TZ, SIG_THAL, SIG_QH, SIG_CON, SIG_RAD, ALPHA_SOL)
% EXTRACT_SIGNALS  Pull all zone signals into a struct.
%   D.Tz   [N x nZ]  measured air temperatures [K]  (5 modelled zones)
%   D.THal [N x 1]   measured Hal temperature [K]
%   D.Qhea [N x nZ]  heating power [W]
%   D.Qsol [N x nZ]  solar gain [W]
%   D.Qint [N x nZ]  internal gains [W]
%   D.Tamb [N x 1]   outdoor temperature [K]
    nZ_ = numel(SIG_TZ);
    N   = height(data);
    D.N    = N;
    D.Tamb = data.weatherStation_reaWeaTDryBul_y;
    D.THal = data.(SIG_THAL);
    HGlo   = data.weatherStation_reaWeaHGloHor_y;
    D.Tz   = zeros(N, nZ_);
    D.Qhea = zeros(N, nZ_);
    D.Qsol = zeros(N, nZ_);
    D.Qint = zeros(N, nZ_);
    for i = 1:nZ_
        D.Tz(:,i)   = data.(SIG_TZ{i});
        D.Qhea(:,i) = data.(SIG_QH{i});
        D.Qsol(:,i) = HGlo * ALPHA_SOL(i);
        D.Qint(:,i) = data.(SIG_CON{i}) + data.(SIG_RAD{i});
    end
end

id  = extract_signals(data_id,  SIG_TZ, SIG_THAL, SIG_QH, SIG_CON, SIG_RAD, ALPHA_SOL);
val = extract_signals(data_val, SIG_TZ, SIG_THAL, SIG_QH, SIG_CON, SIG_RAD, ALPHA_SOL);

fprintf('Data: id=%d steps (%.1fd), val=%d steps (%.1fd), Ts=%gs\n\n', ...
        id.N, id.N*Ts/86400, val.N, val.N*Ts/86400, Ts);

% --- BASELINE dataset for the MLPRS excitation-comparison plots (§3.3.1) ---
% Loaded best-effort: if the file is missing the comparison plots are skipped.
try
    data_base_id = readtable('data/baseline_identification.csv');
    base_id     = extract_signals(data_base_id, SIG_TZ, SIG_THAL, SIG_QH, SIG_CON, SIG_RAD, ALPHA_SOL);
    have_base   = true;
    fprintf('Baseline dataset (PI-controlled): %d steps (%.1fd)\n\n', ...
            base_id.N, base_id.N*Ts/86400);
catch
    have_base = false;
    warning('baseline_identification.csv not found — MLPRS excitation-comparison plots will be skipped.');
end

%% Section 5. EKF core functions

function [f1, f2, UA_sum1] = zone_rates_2C(T1, T2, theta, Tamb, T_Hal, Q_net, ...
                                            nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL)
% ZONE_RATES_2C  Compute dT1/dt and dT2/dt for all 5 modelled zones.
%   T_Hal   scalar measured hallway temperature.

    T1    = reshape(T1,    nZ, 1);
    T2    = reshape(T2,    nZ, 1);
    Q_net = reshape(Q_net, nZ, 1);

    % Unpack theta: [C1, C2, UAe, UAi, UAh, UAc]
    C1  = reshape(theta(1:nZ),                      nZ, 1);
    C2  = reshape(theta(nZ       + (1:nZ)),          nZ, 1);
    UAe = reshape(theta(2*nZ     + (1:nE)),          nE, 1);
    UAi = reshape(theta(2*nZ+nE  + (1:nZ)),          nZ, 1);
    UAh = reshape(theta(2*nZ+nE+nZ   + (1:nH)),      nH, 1);
    UAc = reshape(theta(2*nZ+nE+nZ+nH + (1:nP)),     nP, 1);

    % Total conductance draining each air node
    UA_sum1 = UAi + UAh;   % slab coupling + Hal exchange always present
    for i = 1:nZ
        if HAS_EXT_WALL(i)
            e = find(EXT_ZONE_IDX == i);
            UA_sum1(i) = UA_sum1(i) + UAe(e);
        end
    end
    for p = 1:nP
        i = ADJ_PAIRS(p,1);  j = ADJ_PAIRS(p,2);
        UA_sum1(i) = UA_sum1(i) + UAc(p);
        UA_sum1(j) = UA_sum1(j) + UAc(p);
    end

    % Heat injected into each air node
    Q_in1 = UAi .* T2 + UAh * T_Hal;   % slab contribution + Hal exchange
    for i = 1:nZ
        if HAS_EXT_WALL(i)
            e = find(EXT_ZONE_IDX == i);
            Q_in1(i) = Q_in1(i) + UAe(e) * Tamb;
        end
    end
    for p = 1:nP
        i = ADJ_PAIRS(p,1);  j = ADJ_PAIRS(p,2);
        Q_in1(i) = Q_in1(i) + UAc(p) * T1(j);
        Q_in1(j) = Q_in1(j) + UAc(p) * T1(i);
    end

    f1 = (Q_in1 - UA_sum1.*T1 + Q_net) ./ C1;
    f2 = UAi .* (T1 - T2) ./ C2;
end


function F = build_jacobian_2C(z, Tamb, T_Hal, f1, f2, UA_sum1, Ts, ...
                                nz, nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL)
% BUILD_JACOBIAN_2C  Analytical Jacobian F = dz_next/dz.
%   T_Hal is not a state, so no column for it in F.

    T1  = reshape(z(1:nZ),                   nZ, 1);
    T2  = reshape(z(nZ+(1:nZ)),              nZ, 1);
    C1  = reshape(z(2*nZ+(1:nZ)),            nZ, 1);
    C2  = reshape(z(3*nZ+(1:nZ)),            nZ, 1);
    UAi = reshape(z(4*nZ+nE+(1:nZ)),         nZ, 1);
    UAc = reshape(z(4*nZ+nE+nZ+nH+(1:nP)),   nP, 1);

    F = eye(nz);

    for i = 1:nZ
        % Air node row i

        % Self-decay of T1_i
        F(i, i) = 1 - Ts * UA_sum1(i) / C1(i);

        % T2_i -> T1_i  (slab->air)
        F(i, nZ+i) = Ts * UAi(i) / C1(i);

        % C1_i -> T1_i
        F(i, 2*nZ+i) = -Ts * f1(i) / C1(i);

        % UA_ext_i -> T1_i
        if HAS_EXT_WALL(i)
            e = find(EXT_ZONE_IDX == i);
            F(i, 4*nZ+e) = Ts * (Tamb - T1(i)) / C1(i);
        end

        % UA_int_i -> T1_i
        F(i, 4*nZ+nE+i) = Ts * (T2(i) - T1(i)) / C1(i);

        % UA_toHal_i -> T1_i  (T_Hal is a known scalar, not a state)
        F(i, 4*nZ+nE+nZ+i) = Ts * (T_Hal - T1(i)) / C1(i);

        % Slab node row nZ+i

        % T1_i -> T2_i
        F(nZ+i, i) = Ts * UAi(i) / C2(i);

        % Self-decay of T2_i
        F(nZ+i, nZ+i) = 1 - Ts * UAi(i) / C2(i);

        % C2_i -> T2_i
        F(nZ+i, 3*nZ+i) = -Ts * f2(i) / C2(i);

        % UA_int_i -> T2_i
        F(nZ+i, 4*nZ+nE+i) = Ts * (T1(i) - T2(i)) / C2(i);
    end

    % Direct zone coupling
    for p = 1:nP
        i       = ADJ_PAIRS(p,1);  j = ADJ_PAIRS(p,2);
        col_UAc = 4*nZ + nE + nZ + nH + p;

        % T1_j -> T1_i and vice versa
        F(i, j) = Ts * UAc(p) / C1(i);
        F(j, i) = Ts * UAc(p) / C1(j);

        % UA_ij -> T1_i and T1_j
        F(i, col_UAc) = Ts * (T1(j) - T1(i)) / C1(i);
        F(j, col_UAc) = Ts * (T1(i) - T1(j)) / C1(j);
    end

    % Parameter rows: identity (already set by eye(nz))
end

%% Section 6. Run EKF on identification set

C_obs = [eye(nZ), zeros(nZ, nz-nZ)];   % 5 x 39
R_ekf = R_meas * eye(nZ);

T2_init = id.Tz(1,:)';   % slab starts at air temperature

z0 = [ id.Tz(1,:)';   % T1
        T2_init;       % T2
        theta0' ];     % parameters

P0 = diag([ 0.5^2     * ones(1,nZ),  ...   % T1 [K^2]
            1.0^2     * ones(1,nZ),  ...   % T2 [K^2]
            (50e3)^2  * ones(1,nZ),  ...   % C1 [J/K]^2
            (0.5e6)^2 * ones(1,nZ),  ...   % C2 [J/K]^2
            5^2       * ones(1,nE),  ...   % UAe [W/K]^2
            10^2      * ones(1,nZ),  ...   % UAi [W/K]^2
            10^2      * ones(1,nH),  ...   % UAh [W/K]^2
            5^2       * ones(1,nP) ]);     % UAc [W/K]^2

% Physical bounds
C1_air  = [50, 35, 35, 30, 20] * 1231;   % air-only floor [J/K]
UAi_LB  = [40, 25, 25, 20, 15];          % h*A_floor lower bound [W/K]

% Per-zone C2 bounds
%   LB: 0.04m slab × 0.5 safety  UB: 3× physical (covers walls+ceiling mass)
%   Liv~16m², Ro1/Ro2~10m², Ro3~8m², Bth~5m²
C2_LB = [0.8, 0.4, 0.4, 0.3, 0.2]' * 1e6;   % [J/K]
C2_UB = [8.0, 5.0, 5.0, 4.0, 2.5]' * 1e6;   % [J/K]

% Per-zone UAi upper bounds: prevents C2/UAi degeneracy
%   (tau_slab = C2/UAi is what data constrains; without UB both can inflate)
UAi_UB = [240, 150, 150, 80, 75];   % [W/K]

% Per-zone UAh upper bounds: UA_toHal must be capped or it dominates
%   UA_int for Ro3/Bth, making those params degenerate in open-loop.
%   Wall+door geometry: U=1.5 W/m²K × A_wall + 8 W/K door
UAh_UB = [80, 60, 60, 50, 60]';   % [W/K]

LB = [ -inf(2*nZ, 1);
        C1_air(:);           % C1 > air-only floor
        C2_LB;               % C2 > zone-specific slab floor
        0.1  * ones(nE,1);   % UAe > 0.1 W/K
        UAi_LB(:);           % UAi > zone-specific floor
        1.0  * ones(nH,1);   % UAh > 1 W/K
        0.1  * ones(nP,1) ]; % UAc > 0.1 W/K

UB = [  inf(2*nZ, 1);
        2e6    * ones(nZ,1);  % C1 < 2 MJ/K
        C2_UB;                % C2 per-zone ceiling
        150    * ones(nE,1);  % UAe < 150 W/K
        UAi_UB(:);            % UAi per-zone ceiling
        UAh_UB;               % UAh per-zone ceiling (prevents Hal shortcut dominating)
        200    * ones(nP,1) ]; % UAc < 200 W/K

z_log = zeros(nz, id.N);
P_log = zeros(nz, nz, id.N);
z_log(:,1)   = z0;
P_log(:,:,1) = P0;
z_hat = z0;
P_hat = P0;

fprintf('Running EKF (2C) on identification set...\n');
t_start = tic;

for k = 1:id.N-1

    T1_k    = z_hat(IDX_T1);
    T2_k    = z_hat(IDX_T2);
    theta_k = z_hat(2*nZ+1:end);
    Tamb_k  = id.Tamb(k);
    THal_k  = id.THal(k);
    Q_net_k = id.Qhea(k,:)' + id.Qsol(k,:)' + id.Qint(k,:)';

    % PREDICT (RK4)
    [T1_new, T2_new, f1_eff, f2_eff, UA_sum1] = rk4_step( ...
        T1_k, T2_k, theta_k, Tamb_k, THal_k, Q_net_k, Ts, ...
        nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

    z_pred            = z_hat;
    z_pred(IDX_T1)    = max(min(T1_new, 320), 250);
    z_pred(IDX_T2)    = max(min(T2_new, 320), 250);

    F = build_jacobian_2C(z_hat, Tamb_k, THal_k, f1_eff, f2_eff, UA_sum1, Ts, ...
                           nz, nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

    P_pred = F * P_hat * F' + Q_ekf;

    % UPDATE
    y_k    = id.Tz(k+1,:)';
    innov  = y_k - C_obs * z_pred;
    S_inn  = C_obs * P_pred * C_obs' + R_ekf;
    K_gain = P_pred * C_obs' / S_inn;

    z_hat  = z_pred + K_gain * innov;
    P_hat  = (eye(nz) - K_gain * C_obs) * P_pred;

    z_hat  = max(z_hat, LB);
    z_hat  = min(z_hat, UB);

    z_log(:,k+1)   = z_hat;
    P_log(:,:,k+1) = P_hat;
end

fprintf('EKF complete in %.1f s.\n\n', toc(t_start));

%% Section 7. Extract converged estimates

i_conv = floor(0.8*id.N) : id.N;

C1_est  = mean(z_log(IDX_C1,  i_conv), 2)';
C2_est  = mean(z_log(IDX_C2,  i_conv), 2)';
UAe_est = mean(z_log(IDX_UAe, i_conv), 2)';
UAi_est = mean(z_log(IDX_UAi, i_conv), 2)';
UAh_est = mean(z_log(IDX_UAh, i_conv), 2)';
UAc_est = mean(z_log(IDX_UAc, i_conv), 2)';

C1_std  = sqrt(diag(P_log(IDX_C1,  IDX_C1,  end)))';
C2_std  = sqrt(diag(P_log(IDX_C2,  IDX_C2,  end)))';
UAe_std = sqrt(diag(P_log(IDX_UAe, IDX_UAe, end)))';
UAi_std = sqrt(diag(P_log(IDX_UAi, IDX_UAi, end)))';
UAh_std = sqrt(diag(P_log(IDX_UAh, IDX_UAh, end)))';
UAc_std = sqrt(diag(P_log(IDX_UAc, IDX_UAc, end)))';

fprintf('=========================================================\n');
fprintf(' CONVERGED PARAMETER ESTIMATES\n');
fprintf('=========================================================\n\n');

fprintf(' Air capacitances C1 [kJ/K]:\n');
fprintf('   %-6s  %8s +/- %-8s\n','Zone','Est.','1sigma');
for i=1:nZ
    fprintf('   %-6s  %8.2f +/- %.2f\n', ZONE_NAMES{i}, C1_est(i)/1e3, C1_std(i)/1e3);
end

fprintf('\n Slab capacitances C2 [MJ/K]:\n');
fprintf('   %-6s  %8s +/- %-8s\n','Zone','Est.','1sigma');
for i=1:nZ
    fprintf('   %-6s  %8.3f +/- %.3f\n', ZONE_NAMES{i}, C2_est(i)/1e6, C2_std(i)/1e6);
end

fprintf('\n Exterior conductances UA_ext [W/K]:\n');
fprintf('   %-6s  %8s +/- %-8s\n','Zone','Est.','1sigma');
for i=1:nZ
    fprintf('   %-6s  %8.3f +/- %.3f\n', ZONE_NAMES{i}, UAe_est(i), UAe_std(i));
end

fprintf('\n Internal (air-slab) conductances UA_int [W/K]:\n');
fprintf('   %-6s  %8s +/- %-8s\n','Zone','Est.','1sigma');
for i=1:nZ
    fprintf('   %-6s  %8.3f +/- %.3f\n', ZONE_NAMES{i}, UAi_est(i), UAi_std(i));
end

fprintf('\n Hallway conductances UA_toHal [W/K]:\n');
fprintf('   %-6s  %8s +/- %-8s\n','Zone','Est.','1sigma');
for i=1:nZ
    fprintf('   %-6s  %8.3f +/- %.3f\n', ZONE_NAMES{i}, UAh_est(i), UAh_std(i));
end

fprintf('\n Direct coupling conductances UA_ij [W/K]:\n');
fprintf('   %-12s  %8s +/- %-8s\n','Wall','Est.','1sigma');
for p=1:nP
    i=ADJ_PAIRS(p,1); j=ADJ_PAIRS(p,2);
    fprintf('   %-12s  %8.3f +/- %.3f\n', ...
            sprintf('%s-%s',ZONE_NAMES{i},ZONE_NAMES{j}), UAc_est(p), UAc_std(p));
end

fprintf('\n Implied time constants:\n');
for i=1:nZ
    UA_tot = UAi_est(i) + UAe_est(i) + UAh_est(i);
    for p=1:nP
        if ADJ_PAIRS(p,1)==i || ADJ_PAIRS(p,2)==i
            UA_tot = UA_tot + UAc_est(p);
        end
    end
    tau_air  = C1_est(i) / UA_tot / 3600;
    tau_slab = C2_est(i) / UAi_est(i) / 3600;
    fprintf('   %-6s  tau_air=%.2fh   tau_slab=%.1fh\n', ZONE_NAMES{i}, tau_air, tau_slab);
end

%% Section 8. Build state-space model + save
%
%  Linearises the 2C model at the converged EKF parameters and saves the
%  discrete ZOH matrices for use by the MPC params build scripts.
%
%  Topology:
%    x = [T_air(5); T_wall(5)]               nx_ss = 10
%    u = Q_heat(5) [W]                        nU_ss =  5
%    d = [Tamb; HGloHor; Qint(5)]             nD_ss =  7
%    y = T_air(5)      (C_ctrl = C_full)
%
%  T_Hal is included as the last disturbance.  It is measured at every step
%  and held constant over the MPC horizon (valid since T_Hal changes slowly
%  relative to the 15-min sample time).  Excluding it from Bd while keeping
%  UAh in the A self-drain would force d_bias to absorb ~160 K per zone.

nx_ss = 2 * nZ;              %  10: [T_air(5); T_wall(5)]
nU_ss = nZ;                  %   5: Q_heat per zone
nD_ss = 1 + 1 + nZ + 1;     %   8: [Tamb; HGloHor; Qint(5); T_Hal]

% Continuous A (10×10)
Ac_ss = zeros(nx_ss, nx_ss);
for i = 1:nZ
    UA_drain = UAe_est(i) + UAh_est(i) + UAi_est(i);
    for pp = 1:nP
        if ADJ_PAIRS(pp,1)==i || ADJ_PAIRS(pp,2)==i
            UA_drain = UA_drain + UAc_est(pp);
        end
    end
    Ac_ss(i,    i)    = -UA_drain    / C1_est(i);   % air self-drain
    Ac_ss(i,  nZ+i)   =  UAi_est(i) / C1_est(i);   % slab → air
    Ac_ss(nZ+i, i)    =  UAi_est(i) / C2_est(i);   % air → slab
    Ac_ss(nZ+i, nZ+i) = -UAi_est(i) / C2_est(i);   % slab self-drain
end
%  Zone-to-zone air coupling (off-diagonal entries of the air block)
for pp = 1:nP
    a = ADJ_PAIRS(pp,1);  b = ADJ_PAIRS(pp,2);
    Ac_ss(a, b) = Ac_ss(a, b) + UAc_est(pp) / C1_est(a);
    Ac_ss(b, a) = Ac_ss(b, a) + UAc_est(pp) / C1_est(b);
end

% Continuous Bu (10×5) : heat injected into each air node
Buc_ss = zeros(nx_ss, nU_ss);
for i = 1:nZ
    Buc_ss(i, i) = 1 / C1_est(i);
end

% Continuous Bd (10×8) : d = [Tamb; HGloHor; Qint(5); T_Hal]
Bdc_ss = zeros(nx_ss, nD_ss);
for i = 1:nZ
    Bdc_ss(i, 1)      = UAe_est(i)   / C1_est(i);   % Tamb    → air (exterior wall)
    Bdc_ss(i, 2)      = ALPHA_SOL(i) / C1_est(i);   % HGloHor → air (solar gain)
    Bdc_ss(i, 2+i)    = 1            / C1_est(i);   % Qint_i  → air node i
    Bdc_ss(i, nD_ss)  = UAh_est(i)   / C1_est(i);   % T_Hal   → air (hallway conduction)
end

% Output matrices : all 5 air temperatures
%  C_ctrl = C_full since Hal is not a state (no 6th measurement row needed)
C_ctrl_ss = [eye(nZ), zeros(nZ, nZ)];   % 5×10
C_full_ss  = C_ctrl_ss;

% Discrete ZOH via single augmented expm
nIn_ss = nU_ss + nD_ss;   % 13
M_zoh  = [Ac_ss, [Buc_ss, Bdc_ss]; zeros(nIn_ss, nx_ss + nIn_ss)];
E_zoh  = expm(M_zoh * Ts);
Ad     = E_zoh(1:nx_ss,          1:nx_ss);
Bud    = E_zoh(1:nx_ss, nx_ss         + (1:nU_ss));   % 10×5
Bdd    = E_zoh(1:nx_ss, nx_ss+nU_ss   + (1:nD_ss));   % 10×8

%  Diagnostics
eigs_c = eig(Ac_ss);
tau_h  = sort(-1 ./ real(eigs_c(real(eigs_c) < 0)) / 3600, 'ascend');
fprintf('[ss_build] Time constants: %.1f – %.1f h  |  rho(Ad) = %.6f\n', ...
        tau_h(1), tau_h(end), max(abs(eig(Ad))));

% Save
%  nZ is saved as nZ+1 = 6 so that pi_baseline (nU = nZ-1 = 5) and
%  ekf_mpc (nZ_boptest = nZ = 6) work correctly after loading.
%  ZONE_NAMES gets Hal appended for the same reason.
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'models');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

s_rc               = struct();
s_rc.nZ            = nZ + 1;
s_rc.ZONE_NAMES    = [ZONE_NAMES(:)', {'Hal'}];
s_rc.nP            = nP;
s_rc.Ts            = Ts;
s_rc.ADJ_PAIRS     = ADJ_PAIRS;
s_rc.nU_ss         = nU_ss;
s_rc.nD_ss         = nD_ss;
s_rc.nx_ss         = nx_ss;
s_rc.Ac            = Ac_ss;
s_rc.Buc           = Buc_ss;
s_rc.Bdc           = Bdc_ss;
s_rc.Ad            = Ad;
s_rc.Bud           = Bud;
s_rc.Bdd           = Bdd;
s_rc.C_ctrl        = C_ctrl_ss;
s_rc.C_full        = C_full_ss;
s_rc.C1_est        = C1_est;
s_rc.C2_est        = C2_est;
s_rc.UAe_est       = UAe_est;
s_rc.UAi_est       = UAi_est;
s_rc.UAh_est       = UAh_est;
s_rc.UAc_est       = UAc_est;
s_rc.ALPHA_SOL     = ALPHA_SOL;

save(fullfile(out_dir, 'rc_params_2C6z.mat'), '-struct', 's_rc');
fprintf('Saved  →  %s\n\n', fullfile(out_dir, 'rc_params_2C6z.mat'));
clear s_rc Ac_ss Buc_ss Bdc_ss C_ctrl_ss C_full_ss Ad Bud Bdd

%% Section 9. Open-loop validation

function [T1_new, T2_new, f1_eff, f2_eff, UA_sum1] = rk4_step( ...
        T1, T2, theta, Tamb, T_Hal, Q_net, Ts, ...
        nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL)
% RK4_STEP  Classic 4th-order Runge-Kutta step for the 2C thermal model.
%   Parameters are held constant within the step (random-walk in EKF).
%   Intermediate temperatures clamped to [230, 340] K to prevent NaN/Inf
%   if parameters are temporarily unphysical during EKF transients.

    T_LO = 230;  T_HI = 340;
    clamp = @(T) max(min(T, T_HI), T_LO);

    [k1_f1, k1_f2, UA_sum1] = zone_rates_2C(T1, T2, theta, Tamb, T_Hal, Q_net, ...
                                              nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

    T1_k2 = clamp(T1 + (Ts/2)*k1_f1);  T2_k2 = clamp(T2 + (Ts/2)*k1_f2);
    [k2_f1, k2_f2] = zone_rates_2C(T1_k2, T2_k2, theta, Tamb, T_Hal, Q_net, ...
                                     nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

    T1_k3 = clamp(T1 + (Ts/2)*k2_f1);  T2_k3 = clamp(T2 + (Ts/2)*k2_f2);
    [k3_f1, k3_f2] = zone_rates_2C(T1_k3, T2_k3, theta, Tamb, T_Hal, Q_net, ...
                                     nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

    T1_k4 = clamp(T1 + Ts*k3_f1);      T2_k4 = clamp(T2 + Ts*k3_f2);
    [k4_f1, k4_f2] = zone_rates_2C(T1_k4, T2_k4, theta, Tamb, T_Hal, Q_net, ...
                                     nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

    f1_eff = (k1_f1 + 2*k2_f1 + 2*k3_f1 + k4_f1) / 6;
    f2_eff = (k1_f2 + 2*k2_f2 + 2*k3_f2 + k4_f2) / 6;

    T1_new = clamp(T1 + Ts*f1_eff);
    T2_new = clamp(T2 + Ts*f2_eff);
end


function T1_sim = simulate_rc_2C(D, Ts, C1_est, C2_est, UAe_est, UAi_est, UAh_est, UAc_est, ...
                                   nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL)
% SIMULATE_RC_2C  Open-loop RK4 simulation, no KF correction.
%   Uses measured T_Hal at each step.
%   Returns T1_sim [nZ x N].

    theta = [C1_est(:)', C2_est(:)', UAe_est(:)', UAi_est(:)', UAh_est(:)', UAc_est(:)'];
    T1 = D.Tz(1,:)';
    T2 = D.Tz(1,:)';

    T1_sim = zeros(nZ, D.N);
    T1_sim(:,1) = T1;

    for k = 1:D.N-1
        Q_net_k = D.Qhea(k,:)' + D.Qsol(k,:)' + D.Qint(k,:)';
        [T1, T2] = rk4_step(T1, T2, theta, D.Tamb(k), D.THal(k), Q_net_k, Ts, ...
                              nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);
        T1 = max(min(T1, 330), 240);
        T2 = max(min(T2, 330), 240);
        T1_sim(:,k+1) = T1;
    end
end

fprintf('\nRunning open-loop simulation...\n');
T_sim_id  = simulate_rc_2C(id,  Ts, C1_est, C2_est, UAe_est, UAi_est, UAh_est, UAc_est, ...
                             nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);
T_sim_val = simulate_rc_2C(val, Ts, C1_est, C2_est, UAe_est, UAi_est, UAh_est, UAc_est, ...
                             nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);

rmse_fn = @(a,b) sqrt(mean((a-b).^2));
fprintf('\n Open-loop RMSE:\n');
fprintf('   %-6s  %8s  %8s\n','Zone','RMSE-ID','RMSE-Val');
for i=1:nZ
    fprintf('   %-6s  %8.4f  %8.4f\n', ZONE_NAMES{i}, ...
            rmse_fn(id.Tz(:,i),  T_sim_id(i,:)'), ...
            rmse_fn(val.Tz(:,i), T_sim_val(i,:)'));
end

%% Section 10. Plots
%  Every figure of Chapter 3 is exported straight into figures/, so they can be
%  uploaded to Overleaf without being re-saved by hand.

fig_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

t_id  = (0:id.N-1)'  * Ts/3600;
t_val = (0:val.N-1)' * Ts/3600;

% Calendar timestamps for the date axes. The .time column of the CSVs is
% relative to the start of recording, and both campaigns start at t = 0, so the
% identification set covers days 0-14 and the validation set days 14-28. The
% year is arbitrary: BOPTEST uses TMY weather and only the month and day are
% shown on the axes.
T_REF   = datetime(2025,1,1,0,0,0,'Format','dd-MMM');   % day 0 of the simulated year
dt_id   = T_REF + seconds(data_id.time);
dt_val  = T_REF + seconds(data_val.time);
if exist('data_base_id','var')
    dt_base = T_REF + seconds(data_base_id.time);
end

ZONE_COLORS = [0.12 0.47 0.71;
               0.86 0.35 0.19;
               0.27 0.63 0.27;
               0.58 0.40 0.74;
               0.80 0.47 0.65];
c_meas = [0.15 0.15 0.15];
c_slab = [0.70 0.87 0.95];
c_hal  = [0.50 0.50 0.50];

% Fig 0a & 0b: MLPRS excitation comparison (Section 3.3.1)
% These figures support the choice of MLPRS over baseline PI operation for
% grey-box identification: 0a is the qualitative view (input/output traces),
% 0b is the quantitative view (singular values of the input covariance).
if have_base
    % --- Fig 0a: input/output time-trace comparison (Liv as representative) ---
    %  Both panels share the SAME colour scheme so the reader compares the
    %  pattern, not the colours:  heating power = blue (left axis),
    %  zone temperature = orange (right axis).  The dataset is identified
    %  by the panel title only.
    zone_show = 1;   % Liv
    n_show    = min([base_id.N, id.N, round(4 * 86400 / Ts)]);   % first ~4 days
    t_show    = (0:n_show-1)' * Ts/3600;
    dt_base_show = dt_base(1:n_show);
    dt_id_show   = dt_id(1:n_show);
    col_Q     = [0.12 0.47 0.71];   % blue:   heating power
    col_T     = [0.86 0.35 0.19];   % orange: zone temperature
    col_base  = [0.86 0.35 0.19];   % (kept for Fig 0b, distinguishes datasets there)
    col_mlp   = [0.12 0.47 0.71];

    figure('Name','MLPRS excitation — input/output traces','Color','white', ...
           'Position',[100 100 1200 600]);

    subplot(2,1,1);
    yyaxis left;
    plot(dt_base_show, base_id.Qhea(1:n_show, zone_show)/1e3, '-', 'LineWidth',1.0, ...
         'Color',col_Q);
    ylabel('Q_{heat,Liv} [kW]','FontSize',9);
    yyaxis right;
    plot(dt_base_show, base_id.Tz(1:n_show, zone_show)-273.15, '-', 'LineWidth',1.0, ...
         'Color',col_T);
    ylabel('T_{Liv} [\circC]','FontSize',9);
    % Force axis colours to match the line colours so labels and lines agree
    ax = gca;
    ax.YAxis(1).Color = col_Q;
    ax.YAxis(2).Color = col_T;
    xlabel('Date','FontSize',9);
    title('\bfBaseline closed-loop\rm — built-in PI tracking a comfort setpoint', ...
          'FontSize',10,'FontWeight','normal');
    grid on; box on;

    subplot(2,1,2);
    yyaxis left;
    plot(dt_id_show, id.Qhea(1:n_show, zone_show)/1e3, '-', 'LineWidth',1.0, ...
         'Color',col_Q);
    ylabel('Q_{heat,Liv} [kW]','FontSize',9);
    yyaxis right;
    plot(dt_id_show, id.Tz(1:n_show, zone_show)-273.15, '-', 'LineWidth',1.0, ...
         'Color',col_T);
    ylabel('T_{Liv} [\circC]','FontSize',9);
    ax = gca;
    ax.YAxis(1).Color = col_Q;
    ax.YAxis(2).Color = col_T;
    xlabel('Date','FontSize',9);
    title('\bfMLPRS open-loop\rm — multi-level pseudo-random excitation', ...
          'FontSize',10,'FontWeight','normal');
    grid on; box on;
    sgtitle('Excitation comparison: heating power and resulting zone temperature (Liv)', ...
            'FontSize',12,'FontWeight','bold');
    exportgraphics(gcf, fullfile(fig_dir,'excitation_plot.png'), 'Resolution', 200);

    % --- Fig 0b: information matrix singular-value comparison ---
    %  Per-step input vector  u(k) = [Qhea(k,1) ... Qhea(k,5)]^T.
    %  Sample covariance  Sigma_u = (1/N) * sum_k (u(k)-u_mean)(u(k)-u_mean)^T
    %  is proportional to the Fisher information of any parameter that
    %  enters the model linearly through the heating inputs.  Its smallest
    %  eigenvalue measures the worst-excited direction of the input space;
    %  the condition number kappa = lambda_max / lambda_min quantifies how
    %  ill- or well-posed the identification problem is.
    U_base    = base_id.Qhea - mean(base_id.Qhea, 1);
    U_mlp     = id.Qhea      - mean(id.Qhea,      1);
    Sigma_b   = (U_base' * U_base) / size(U_base, 1);
    Sigma_m   = (U_mlp'  * U_mlp ) / size(U_mlp,  1);
    sv_base   = sort(eig(Sigma_b), 'descend');
    sv_mlp    = sort(eig(Sigma_m), 'descend');
    cond_base = sv_base(1) / max(sv_base(end), eps);
    cond_mlp  = sv_mlp(1)  / max(sv_mlp(end),  eps);

    figure('Name','MLPRS excitation — input information matrix','Color','white', ...
           'Position',[140 140 720 480]);
    semilogy(1:nZ, sv_base, 'o-', 'Color',col_base, 'LineWidth',1.5, ...
             'MarkerFaceColor',col_base, 'MarkerSize',7, ...
             'DisplayName', sprintf('Baseline  (\\kappa = %.2g)', cond_base));
    hold on;
    semilogy(1:nZ, sv_mlp,  's-', 'Color',col_mlp,  'LineWidth',1.5, ...
             'MarkerFaceColor',col_mlp,  'MarkerSize',7, ...
             'DisplayName', sprintf('MLPRS  (\\kappa = %.2g)', cond_mlp));
    xlabel('Mode index','FontSize',9);
    ylabel('Eigenvalue of \Sigma_u  [W^2]','FontSize',9);
    title('Input information matrix: sample covariance of zone heating commands', ...
          'FontSize',10,'FontWeight','normal');
    legend('Location','northeast','FontSize',9,'Box','off');
    grid on; box on;
    set(gca,'XTick',1:nZ);
    exportgraphics(gcf, fullfile(fig_dir,'input_information.png'), 'Resolution', 200);

    fprintf('\nMLPRS excitation comparison:\n');
    fprintf('   Baseline  cond(Sigma_u) = %.6g\n', cond_base);
    fprintf('   MLPRS     cond(Sigma_u) = %.6g\n', cond_mlp);
    fprintf('   Improvement factor    = %.2gx\n\n', cond_base/cond_mlp);

end


% Fig 2: RC parameter trajectories, single thesis-ready figure
%   3 rows x 2 cols, all 6 RC parameter groups in one compact figure.
%     Row 1:  phi_air (C1)        |  phi_mass (C2)
%     Row 2:  UA_e (envelope)     |  UA_i (air <-> slab)
%     Row 3:  UA_h (zone <-> Hal) |  UA_c (zone <-> zone)
figure('Name','EKF RC-parameter trajectories','Color','white', ...
       'Position',[40 40 960 960]);

% --- (1,1)  phi_air  (= C1) ---------------------------------------------
subplot(3,2,1);
hold on; box on; grid on;
for i = 1:nZ
    plot(dt_id, z_log(IDX_C1(i),:)/1e3, '-', ...
         'Color', ZONE_COLORS(i,:), 'LineWidth', 1.3);
end
ylabel('$\varphi_{\mathrm{air}}$ [kJ/K]','Interpreter','latex','FontSize',11);
title('Air-node capacitance $\varphi_{\mathrm{air}}$', ...
      'Interpreter','latex','FontSize',11,'FontWeight','normal');
legend(ZONE_NAMES,'Location','best','FontSize',8,'Box','off');
set(gca,'FontSize',9);

% --- (1,2)  phi_mass  (= C2) --------------------------------------------
subplot(3,2,2);
hold on; box on; grid on;
for i = 1:nZ
    plot(dt_id, z_log(IDX_C2(i),:)/1e6, '-', ...
         'Color', ZONE_COLORS(i,:), 'LineWidth', 1.3);
end
ylabel('$\varphi_{\mathrm{mass}}$ [MJ/K]','Interpreter','latex','FontSize',11);
title('Slab capacitance $\varphi_{\mathrm{mass}}$', ...
      'Interpreter','latex','FontSize',11,'FontWeight','normal');
legend(ZONE_NAMES,'Location','best','FontSize',8,'Box','off');
set(gca,'FontSize',9);

% --- (2,1)  UA_e  (envelope) --------------------------------------------
subplot(3,2,3);
hold on; box on; grid on;
for i = 1:nE
    plot(dt_id, z_log(IDX_UAe(i),:), '-', ...
         'Color', ZONE_COLORS(i,:), 'LineWidth', 1.2);
end
ylabel('$U\!A_{e}$ [W/K]','Interpreter','latex','FontSize',11);
title('Envelope conductance $U\!A_{e}$ (zone $\leftrightarrow$ ambient)', ...
      'Interpreter','latex','FontSize',11,'FontWeight','normal');
legend(ZONE_NAMES,'Location','best','FontSize',8,'Box','off');
set(gca,'FontSize',9);

% --- (2,2)  UA_i  (air <-> slab) ----------------------------------------
subplot(3,2,4);
hold on; box on; grid on;
for i = 1:nZ
    plot(dt_id, z_log(IDX_UAi(i),:), '-', ...
         'Color', ZONE_COLORS(i,:), 'LineWidth', 1.2);
end
ylabel('$U\!A_{i}$ [W/K]','Interpreter','latex','FontSize',11);
title('Internal conductance $U\!A_{i}$ (air $\leftrightarrow$ slab)', ...
      'Interpreter','latex','FontSize',11,'FontWeight','normal');
legend(ZONE_NAMES,'Location','best','FontSize',8,'Box','off');
set(gca,'FontSize',9);

% --- (3,1)  UA_h  (zone <-> Hal) ----------------------------------------
subplot(3,2,5);
hold on; box on; grid on;
for i = 1:nH
    plot(dt_id, z_log(IDX_UAh(i),:), '-', ...
         'Color', ZONE_COLORS(i,:), 'LineWidth', 1.2);
end
xlabel('Date','FontSize',10);
ylabel('$U\!A_{h}$ [W/K]','Interpreter','latex','FontSize',11);
title('Hallway conductance $U\!A_{h}$ (zone $\leftrightarrow$ Hal)', ...
      'Interpreter','latex','FontSize',11,'FontWeight','normal');
legend(ZONE_NAMES,'Location','best','FontSize',8,'Box','off');
set(gca,'FontSize',9);

% --- (3,2)  UA_c  (zone <-> zone) ---------------------------------------
subplot(3,2,6);
hold on; box on; grid on;
pair_labels = cell(nP,1);
for p = 1:nP
    i = ADJ_PAIRS(p,1); j = ADJ_PAIRS(p,2);
    col = (ZONE_COLORS(i,:) + ZONE_COLORS(j,:)) / 2;
    plot(dt_id, z_log(IDX_UAc(p),:), '-', ...
         'Color', col, 'LineWidth', 1.2);
    pair_labels{p} = sprintf('%s-%s', ZONE_NAMES{i}, ZONE_NAMES{j});
end
xlabel('Date','FontSize',10);
ylabel('$U\!A_{c}$ [W/K]','Interpreter','latex','FontSize',11);
title('Inter-zone coupling $U\!A_{ij}$ (zone $\leftrightarrow$ zone)', ...
      'Interpreter','latex','FontSize',11,'FontWeight','normal');
legend(pair_labels,'Location','best','FontSize',8,'Box','off');
set(gca,'FontSize',9);

% No overall title: the thesis supplies it through the LaTeX caption.
exportgraphics(gcf, fullfile(fig_dir,'parameter_convergence.png'), 'Resolution', 200);

%% Section 10b. Validation on unseen data (one-step + 6 h horizon)
%  Two complementary checks of the identified 2C model on the validation set:
%    (1) One-step-ahead prediction: feed the measured air temperature back at
%        every step and ask the model only to predict the next sample.  This
%        is the natural diagnostic for an estimator trained on a Markovian
%        update.
%    (2) Multi-step open-loop forecast over the MPC horizon (6 h = 24 steps):
%        from many starting points along the validation set, simulate the
%        model forward for 24 steps and compare against the measured
%        temperatures.  This is the regime the MPC actually exercises.
%
%  Both use the converged parameter estimates (C1_est, ..., UAc_est).
%  Disturbances and heating inputs are taken from the validation CSV.

theta_est = [C1_est(:)', C2_est(:)', UAe_est(:)', UAi_est(:)', UAh_est(:)', UAc_est(:)'];

%% One-step-ahead prediction on the validation set
fprintf('\nRunning one-step-ahead prediction on validation set...\n');
T1_one  = zeros(val.N, nZ);
T2_one  = zeros(val.N, nZ);
T1_one(1,:) = val.Tz(1,:);
T2_one(1,:) = val.Tz(1,:);          % init slab at measured air temperature
for k = 1:val.N-1
    T1k = val.Tz(k,:)';             % MEASURED air temperature at step k
    T2k = T2_one(k,:)';              % propagated slab estimate
    Q_net_k = val.Qhea(k,:)' + val.Qsol(k,:)' + val.Qint(k,:)';
    [T1n, T2n] = rk4_step(T1k, T2k, theta_est, val.Tamb(k), val.THal(k), ...
                           Q_net_k, Ts, ...
                           nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);
    T1_one(k+1,:) = T1n';
    T2_one(k+1,:) = T2n';
end

e_one    = val.Tz - T1_one;          % residuals [val.N x nZ]
rmse_one = sqrt(mean(e_one.^2, 1));
mae_one  = mean(abs(e_one), 1);

% Ljung-style NRMSE fit percentage:
%   fit% = 100 * ( 1 - || y - y_hat || / || y - mean(y) || )
fit_pct_fn = @(y, yhat) 100 * (1 - norm(y - yhat) / norm(y - mean(y)));

fit_one = zeros(1, nZ);
for i = 1:nZ
    fit_one(i) = fit_pct_fn(val.Tz(:,i), T1_one(:,i));
end

fprintf('\n One-step-ahead error on validation set:\n');
fprintf('   %-6s  %10s  %10s  %8s\n','Zone','RMSE [K]','MAE [K]','fit [%%]');
for i=1:nZ
    fprintf('   %-6s  %10.4f  %10.4f  %8.2f\n', ...
            ZONE_NAMES{i}, rmse_one(i), mae_one(i), fit_one(i));
end

%% Multi-step open-loop forecast over the MPC horizon (6 h)
%  Use stride=1 so that we get a forecast landing on every step k,
%  enabling a "6 h-ahead prediction" time series in addition to the
%  RMSE-vs-horizon curve.
H_steps     = 24;                     % 24 x 15 min = 6 h
stride      = 1;
step_starts = 1:stride:val.N - H_steps;
n_starts    = numel(step_starts);

fprintf('\nRunning %d-step (%.1f h) forecasts at %d starting points...\n', ...
        H_steps, H_steps*Ts/3600, n_starts);

T1_h      = zeros(H_steps+1, nZ, n_starts);
T1_meas_h = zeros(H_steps+1, nZ, n_starts);

for s = 1:n_starts
    k0 = step_starts(s);
    T1 = val.Tz(k0,:)';
    T2 = T2_one(k0,:)';
    T1_h(1,:,s)      = T1';
    T1_meas_h(1,:,s) = val.Tz(k0,:);
    for h = 1:H_steps
        kk = k0 + h - 1;
        Q_net_k = val.Qhea(kk,:)' + val.Qsol(kk,:)' + val.Qint(kk,:)';
        [T1, T2] = rk4_step(T1, T2, theta_est, val.Tamb(kk), val.THal(kk), ...
                             Q_net_k, Ts, ...
                             nZ, nE, nH, nP, EXT_ZONE_IDX, ADJ_PAIRS, HAS_EXT_WALL);
        T1_h(h+1,:,s)      = T1';
        T1_meas_h(h+1,:,s) = val.Tz(k0+h,:);
    end
end

err_h  = T1_h - T1_meas_h;            % [H+1 x nZ x n_starts]
rmse_h = sqrt(squeeze(mean(err_h.^2, 3)));   % [H+1 x nZ]

% Build the "6 h-ahead prediction" time series:
%   T1_6h_pred(k,:) is the model output at step k, obtained from a forecast
%   started H_steps steps earlier from the measured air temperature.
T1_6h_pred = nan(val.N, nZ);
for s = 1:n_starts
    k_end = step_starts(s) + H_steps;
    T1_6h_pred(k_end,:) = T1_h(H_steps+1,:,s);
end

% Per-zone metrics on the 6 h-ahead forecast (only where defined)
mask_6h  = ~isnan(T1_6h_pred(:,1));
e_6h     = val.Tz(mask_6h,:) - T1_6h_pred(mask_6h,:);
rmse_6h  = sqrt(mean(e_6h.^2, 1));
mae_6h   = mean(abs(e_6h), 1);
fit_6h   = zeros(1, nZ);
for i = 1:nZ
    fit_6h(i) = fit_pct_fn(val.Tz(mask_6h,i), T1_6h_pred(mask_6h,i));
end

fprintf('\n 6 h-ahead error on validation set:\n');
fprintf('   %-6s  %10s  %10s  %8s\n','Zone','RMSE [K]','MAE [K]','fit [%%]');
for i=1:nZ
    fprintf('   %-6s  %10.4f  %10.4f  %8.2f\n', ...
            ZONE_NAMES{i}, rmse_6h(i), mae_6h(i), fit_6h(i));
end

% Fig V1: one-step-ahead prediction, same layout as the open-loop figs
%    Measured (solid) vs one-step predicted (dashed) per zone, RMSE in legend.
%    Every panel keeps its own date axis to stay self-contained.
figure('Name','Validation — one-step-ahead prediction','Color','white', ...
       'Position',[100 100 1380 760]);
for i = 1:nZ
    subplot(3,2,i);
    plot(dt_val, val.Tz(:,i)-273.15, '-', 'Color', c_meas, 'LineWidth', 0.8, ...
         'DisplayName','Measured');
    hold on;
    plot(dt_val, T1_one(:,i)-273.15, '--', 'Color', ZONE_COLORS(i,:), ...
         'LineWidth', 1.4, ...
         'DisplayName', sprintf('2RC  RMSE=%.3fK', rmse_one(i)));
    ylabel('T_{air} [C]','FontSize',9); grid on; box on;
    title(sprintf('\\bf%s\\rm',ZONE_NAMES{i}),'FontSize',10,'FontWeight','normal');
    xlabel('Date','FontSize',9);
    legend('Location','best','FontSize',8,'Box','off');
    set(gca,'FontSize',9);
end
subplot(3,2,6);
plot(dt_val, val.THal-273.15, '-', 'Color', c_hal, 'LineWidth', 1.0);
ylabel('T_{air} [C]','FontSize',9); xlabel('Date','FontSize',9);
title('\bfHall\rm  (measured)','FontSize',10,'FontWeight','normal');
grid on; box on; set(gca,'FontSize',9);
% No overall title: the thesis supplies it through the LaTeX caption (Fig. 3.5).
exportgraphics(gcf, fullfile(fig_dir,'one_step_ahead.png'), 'Resolution', 200);

% Fig V2: 6 h MPC-horizon RMSE growth (averaged over all starts)
horizon_hours = (0:H_steps) * Ts/3600;

figure('Name','Validation — error growth over MPC horizon','Color','white', ...
       'Position',[120 120 760 460]);
hold on; box on; grid on;
for i = 1:nZ
    plot(horizon_hours, rmse_h(:,i), '-o', ...
         'Color', ZONE_COLORS(i,:), 'LineWidth', 1.3, 'MarkerSize', 4, ...
         'DisplayName', ZONE_NAMES{i});
end
xlabel('Prediction horizon [h]','FontSize',10);
ylabel('RMSE [K]','FontSize',10);
legend('Location','northwest','FontSize',9,'Box','off');
% No overall title: the thesis supplies it through the LaTeX caption (Fig. 3.7).
set(gca,'FontSize',10);
exportgraphics(gcf, fullfile(fig_dir,'rmse_growth.png'), 'Resolution', 200);

% Fig V3: 6 h-ahead prediction, same layout as V1
%    At every step k where a full 6 h forecast exists, plot the predicted
%    air temperature (dashed) against the measured one (solid).  RMSE in
%    the legend; fit percentage is reported in the console only.  Every
%    panel keeps its own date axis to stay self-contained.
figure('Name','Validation — 6 h-ahead prediction','Color','white', ...
       'Position',[140 140 1380 760]);
for i = 1:nZ
    subplot(3,2,i);
    plot(dt_val, val.Tz(:,i)-273.15, '-', 'Color', c_meas, 'LineWidth', 0.8, ...
         'DisplayName','Measured');
    hold on;
    plot(dt_val, T1_6h_pred(:,i)-273.15, '--', 'Color', ZONE_COLORS(i,:), ...
         'LineWidth', 1.4, ...
         'DisplayName', sprintf('2RC  RMSE=%.3fK', rmse_6h(i)));
    ylabel('T_{air} [C]','FontSize',9); grid on; box on;
    title(sprintf('\\bf%s\\rm',ZONE_NAMES{i}),'FontSize',10,'FontWeight','normal');
    xlabel('Date','FontSize',9);
    legend('Location','best','FontSize',8,'Box','off');
    set(gca,'FontSize',9);
end
subplot(3,2,6);
plot(dt_val, val.THal-273.15, '-', 'Color', c_hal, 'LineWidth', 1.0);
ylabel('T_{air} [C]','FontSize',9); xlabel('Date','FontSize',9);
title('\bfHall\rm  (measured)','FontSize',10,'FontWeight','normal');
grid on; box on; set(gca,'FontSize',9);
% No overall title: the thesis supplies it through the LaTeX caption (Fig. 3.6).
exportgraphics(gcf, fullfile(fig_dir,'6_h_ahead.png'), 'Resolution', 200);

fprintf('\nAll figures generated.\n');

%% Section 11. EKF tuning summary

fprintf('\n');
fprintf('=========================================================\n');
fprintf(' TUNING GUIDE  (2C model)\n');
fprintf('=========================================================\n');
fprintf('\n Current noise settings:\n');
fprintf('   Q_T1  = %.0e K^2/step   (air temperature process noise)\n',   Q_T1);
fprintf('   Q_T2  = %.0e K^2/step   (slab temperature process noise)\n',  Q_T2);
fprintf('   Q_C1  = %.0e (J/K)^2   (air capacitance random walk)\n',      Q_C1);
fprintf('   Q_C2  = %.0e (J/K)^2   (slab capacitance random walk)\n',     Q_C2);
fprintf('   Q_UAe = %.0e (W/K)^2   (exterior conductance)\n',             Q_UAe);
fprintf('   Q_UAi = %.0e (W/K)^2   (internal conductance)\n',             Q_UAi);
fprintf('   Q_UAh = %.0e (W/K)^2   (UA_toHal)\n',                         Q_UAh);
fprintf('   Q_UAc = %.0e (W/K)^2   (coupling conductance)\n',             Q_UAc);
fprintf('   R     = %.0e K^2        (measurement noise)\n\n',             R_meas);
fprintf(' Symptom -> Remedy\n\n');
fprintf(' tau_slab > 100h  (unrealistically slow)\n');
fprintf('   -> C2 too large or UAi too small; tighten UAi_LB per zone\n');
fprintf('   -> typical hydronic slab: 10-40h\n\n');
fprintf(' tau_air < 0.5h  (unrealistically fast)\n');
fprintf('   -> UAh or direct UAc too large; check for bound saturation\n\n');
fprintf(' UAh saturating at upper bound\n');
fprintf('   -> raise UB for UAh, or check T_Hal variance in data\n\n');
fprintf(' C1 hitting lower bound\n');
fprintf('   -> C1_air floor may be too high; check volume estimates\n\n');
fprintf(' Open-loop RMSE > 1K\n');
fprintf('   -> expected for PRBS winter data; compare with N4SID baseline\n');
