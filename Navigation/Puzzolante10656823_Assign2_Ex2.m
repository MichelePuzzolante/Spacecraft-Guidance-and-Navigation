% Guidance and Navigation (2024/2025)
% Assignment # 2
% Exercise 2
% Author : Michele Puzzolante
%% REQUEST 1

% Initializzation of the MATLAB enviroment
clear;
close all;
clc;
cspice_kclear();

cspice_furnsh('assignment02.tm');               % Load kernels
addpath('sgp4');                                % Add the path of SGP4

% Set default options for plots
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');

fprintf('\nEXERCISE 2\n');
fprintf('\nREQUEST 1\n');

SAVEPLOTS = false;  % Flag to use when plots needs to be saved as png file

% Initializzation of constatnt and structs 
arcsec2rad = pi / (180*3600);       % Constant for arcseconds to radians conversions
mu = cspice_bodvrd('EARTH','GM',1); % Intialization of global variables
typerun    = 'u';                   % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';                   % afspc approach ('air force space command')
whichconst =  72;                   % WGS72 constants (radius, gravitational parameter)
t0 = '2024-11-18T20:30:00.000';     % Initial time in UTC
tf = '2024-11-18T22:15:00.000';     % Final time in UTC
t0_et = cspice_str2et(t0);          % Initial time in ephemerides time
tf_et = cspice_str2et(tf);          % Final time in ephemerides time
% Set nutation corrections parameters for the 2024-11-18 [rad] 
% taken from https://celestrak.org/SpaceData/EOP-All.txt
ddpsi = -0.115178 * arcsec2rad;       
ddeps = -0.007395 * arcsec2rad; 

% Initiailizzation of the struct for the KOUROU station
KOUROU.name = 'KOUROU';                         % Station name
KOUROU.lat = 5.25144;                           % Latitude                     [deg]
KOUROU.lon = -52.80466;                         % Longitude                    [deg]
KOUROU.alt = -14.67;                            % Altitude                     [m]
KOUROU.sigma_AzEl = 125 * 10^-3 * cspice_rpd;   % Noise on angles measurements [rad]
KOUROU.sigma_r = 0.01;                          % Noise on range measurements  [km]
KOUROU.noiseMatrix = ...                        % Noise matrix
    diag([KOUROU.sigma_AzEl * ones(2, 1); KOUROU.sigma_r]);
KOUROU.minEl = 6 * cspice_rpd;                  % Minimum elevation            [rad]
KOUROU.freq = 60;                               % Measurements frequency       [s]
KOUROU.cost = 30000;                            % Cost per pass                [EUR]

% Initiailizzation of the struct for the TROLL station
TROLL.name = 'TROLL';                          % Station name
TROLL.lat = -72.011977;                        % Latitude                     [deg]
TROLL.lon = 2.536103;                          % Longitude                    [deg]
TROLL.alt =  1298 ;                            % Altitude                     [m]
TROLL.sigma_AzEl = 125 * 10^-3 * cspice_rpd;   % Noise on angles measurements [rad]
TROLL.sigma_r = 0.01;                          % Noise on range measurements  [km]
TROLL.noiseMatrix = ...                        % Noise matrix
    diag([TROLL.sigma_AzEl * ones(2, 1); TROLL.sigma_r]);
TROLL.minEl = 0 * cspice_rpd;                  % Minimum elevation            [rad]
TROLL.freq = 30;                               % Measurements frequency       [s]
TROLL.cost = 35000;                            % Cost per pass                [EUR]

% Initiailizzation of the struct for the SVALBARD station
SVALBARD.name = 'SVALBARD';                       % Station name
SVALBARD.lat = 78.229772;                         % Latitude                     [deg]
SVALBARD.lon = 15.407786;                         % Longitude                    [deg]
SVALBARD.alt = 458;                               % Altitude                     [m]
SVALBARD.sigma_AzEl = 125 * 10^-3 * cspice_rpd;   % Noise on angles measurements [rad]
SVALBARD.sigma_r = 0.01;                          % Noise on range measurements  [km]
SVALBARD.noiseMatrix = ...                        % Noise matrix
    diag([SVALBARD.sigma_AzEl * ones(2, 1); SVALBARD.sigma_r]);
SVALBARD.minEl = 8 * cspice_rpd;                  % Minimum elevation            [rad]
SVALBARD.freq = 60;                               % Measurements frequency       [s]
SVALBARD.cost = 35000;                            % Cost per pass                [EUR]

% read SMOS's TLE
SMOS_TLE = read_3LE('36036', 'tle\36036.3le', whichconst);                                                              
% Get SMOS's TLE epoch
[year_SMOS, mon_SMOS, day_SMOS, hr_SMOS, min_SMOS, sec_SMOS] = invjday(SMOS_TLE.jdsatepoch, SMOS_TLE.jdsatepochf);    
% Convert SMOS's TLE epoch into UTC
SMOS_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year_SMOS, mon_SMOS, day_SMOS, hr_SMOS, min_SMOS, sec_SMOS]);
% Convert SMOS's TLE epoch into eph time
SMOS_epoch_et = cspice_str2et(SMOS_epoch_str);                                                                           

% Get SMOS reference state at a reference epoch of t_ref
% Propagate using SGP4 SMOS' state until t_ref
[SMOS_TLE, rr_TEME_SMOS, vv_TEME_SMOS] = sgp4(SMOS_TLE,  0);                                                             
% Compute centuries from 2000-01-01T00:00:00.00 TDT
ttt = cspice_unitim(SMOS_epoch_et, 'ET', 'TDT')/cspice_jyear()/100;                                                      
% TEME to ECI conversion
[SMOS_rr0_ref, SMOS_vv0_ref] = teme2eci(rr_TEME_SMOS,  vv_TEME_SMOS, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);                 
% Assemle the state at t_ref
SMOS_xx0_ref = [SMOS_rr0_ref; SMOS_vv0_ref];                                                                             

% Computation of the obsculating elemnts at t_ref
SMOS_orb_elem = cspice_oscelt( [rr_TEME_SMOS; vv_TEME_SMOS], SMOS_epoch_et, SMOS_TLE.mu);

% propagation of the state ftom t_ref to t0
[xxf, ~] = propagate2BP(SMOS_xx0_ref, [SMOS_epoch_et t0_et], mu, 'Un');
SMOS_xx0_t0_Kep = xxf(end,:)';

% Computation of the visibility window for KOUROU
[VIS_WIND_KOUROU] = visibilityWindow(KOUROU, SMOS_xx0_t0_Kep, t0_et, tf_et, mu, 'Un');
plotVisbilityWindow(KOUROU, VIS_WIND_KOUROU, [1 2 3]);

% Computation of the visibility window for TROLL
[VIS_WIND_TROLL] = visibilityWindow(TROLL, SMOS_xx0_t0_Kep, t0_et, tf_et, mu, 'Un');
plotVisbilityWindow(TROLL, VIS_WIND_TROLL, [4 5 6]);

% Computation of the visibility window for SVALBARD
[VIS_WIND_SVALBARD] = visibilityWindow(SVALBARD, SMOS_xx0_t0_Kep, t0_et, tf_et, mu, 'Un');
plotVisbilityWindow(SVALBARD, VIS_WIND_SVALBARD, [7 8 9]);

% Display results
fprintf('\nSMOS STATE at t_ref')
fprintf('\nrx\t%.8f\t[km]', SMOS_rr0_ref(1));
fprintf('\nry\t%.8f\t[km]', SMOS_rr0_ref(2));
fprintf('\nrz\t%.8f\t[km]', SMOS_rr0_ref(3));
fprintf('\nvx\t%.8f\t[km/s]', SMOS_vv0_ref(1));
fprintf('\nvy\t%.8f\t[km/s]', SMOS_vv0_ref(2));
fprintf('\nvz\t%.8f\t[km/s]', SMOS_vv0_ref(3));

fprintf('\n\nOBSCUATING ORBITAL ELEMENTS')
fprintf('\na\t\t\t\t%.8f\t[km]', SMOS_orb_elem(1))
fprintf('\ne\t\t\t\t%.8f\t\t[-]', SMOS_orb_elem(2))
fprintf('\ni\t\t\t\t%.4f\t\t\t[deg]', SMOS_orb_elem(3) * cspice_dpr)
fprintf('\nRAAN\t\t\t%.4f\t\t[deg]', SMOS_orb_elem(4) * cspice_dpr)
fprintf('\nArgOfPer\t\t%.4f\t\t\t[deg]', SMOS_orb_elem(5) * cspice_dpr)
fprintf('\nMean Anomaly\t%.4f\t\t[deg]', SMOS_orb_elem(6) * cspice_dpr)

% Save plots
if SAVEPLOTS

    % Saving 1st figure
    figure(1)
    exportgraphics(gcf,'KOUROU_az.png','Resolution',600);

    % Saving 2nd figure
    figure(2)
    exportgraphics(gcf,'KOUROU_el.png','Resolution',600);

    % Saving 3rd figure
    figure(3)
    exportgraphics(gcf,'KOUROU_rng.png','Resolution',600)

    % Saving 4th figure
    figure(4)
    exportgraphics(gcf,'TROLL_az.png','Resolution',600)

    % Saving 5th figure
    figure(5)
    exportgraphics(gcf,'TROLL_el.png','Resolution',600)

    % Saving 6th figure
    figure(6)
    exportgraphics(gcf,'TROLL_rng.png','Resolution',600)

    % Saving 7th figure
    figure(7)
    exportgraphics(gcf,'SVALBARD_az.png','Resolution',600)

    % Saving 8th figure
    figure(8)
    exportgraphics(gcf,'SVALBARD_el.png','Resolution',600)

    % Saving 9th figure
    figure(9)
    exportgraphics(gcf,'SVALBARD_rng.png','Resolution',600)
end

%% REQUEST 2

fprintf('\nREQUEST 2\n');

% Computation of the measurements for the visibility window for KOUROU
[VIS_WIND_KOUROU] = simulatedMeasurements(KOUROU, VIS_WIND_KOUROU, SMOS_epoch_et, SMOS_TLE, ddpsi, ddeps);
plotVisbilityWindow(KOUROU, VIS_WIND_KOUROU.measurement, []);

% Computation of the measurements for the visibility window for TROLL
[VIS_WIND_TROLL] = simulatedMeasurements(TROLL, VIS_WIND_TROLL, SMOS_epoch_et, SMOS_TLE, ddpsi, ddeps);
plotVisbilityWindow(TROLL, VIS_WIND_TROLL.measurement, []);

% Computation of the measurements for the visibility window for SVALBARD
[VIS_WIND_SVALBARD] = simulatedMeasurements(SVALBARD, VIS_WIND_SVALBARD, SMOS_epoch_et, SMOS_TLE, ddpsi, ddeps);
plotVisbilityWindow(SVALBARD, VIS_WIND_SVALBARD.measurement, []);

%% REQUEST 3

fprintf('\nREQUEST 3\n');

% Propagate using SGP4 SMOS' state until t0
[SMOS_TLE, rr_TEME_SMOS, vv_TEME_SMOS] = sgp4(SMOS_TLE,  (t0_et - SMOS_epoch_et)/60.0);                                
% Compute centuries from 2000-01-01T00:00:00.00 TDT
ttt = cspice_unitim(t0_et, 'ET', 'TDT')/cspice_jyear()/100;
% TEME to ECI conversion
[SMOS_rr0_t0, SMOS_vv0_t0] = teme2eci(rr_TEME_SMOS,  vv_TEME_SMOS, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
% Assemble SMOS's the reference state
SMOS_xx0_t0 = [SMOS_rr0_t0; SMOS_vv0_t0];                                                     

% Set options for lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'StepTolerance', 1e-14, 'FunctionTolerance', 1e-14, 'MaxIterations', 200);

% Call lsqnonlin for KOUROU staion only
[NAV_SOL(1).x, NAV_SOL(1).resnorm, NAV_SOL(1).residual ,NAV_SOL(1).exitflag, ~, ~, NAV_SOL(1).jacobian] =...
    lsqnonlin(@(xx) costFunction(xx, mu, 'Un', VIS_WIND_KOUROU, KOUROU), SMOS_xx0_t0, [], [], options);%

% Call lsqnonlin for all the stations without J2
[NAV_SOL(2).x, NAV_SOL(2).resnorm, NAV_SOL(2).residual ,NAV_SOL(2).exitflag, ~, ~, NAV_SOL(2).jacobian] =...
    lsqnonlin(@(xx) costFunction(xx, mu, 'Un', [VIS_WIND_KOUROU; VIS_WIND_TROLL; VIS_WIND_SVALBARD], [KOUROU; TROLL; SVALBARD]), SMOS_xx0_t0, [], [], options);%

% Call lsqnonlin for all the stations with J2
[NAV_SOL(3).x, NAV_SOL(3).resnorm, NAV_SOL(3).residual ,NAV_SOL(3).exitflag, ~, ~, NAV_SOL(3).jacobian] =...
    lsqnonlin(@(xx) costFunction(xx, mu, 'J2', [VIS_WIND_KOUROU; VIS_WIND_TROLL; VIS_WIND_SVALBARD], [KOUROU; TROLL; SVALBARD]), SMOS_xx0_t0, [], [], options);%

% Elaboration of the results
for i = 1 : length(NAV_SOL)

    % Jacobian
    Jac = full(NAV_SOL(i).jacobian);

    % Computation of the covariance matrix
    NAV_SOL(i).P = NAV_SOL(i).resnorm/(length(NAV_SOL(i).residual) - length(SMOS_xx0_t0)).*inv(Jac.'*Jac);

    % Sigma r
    NAV_SOL(i).sigmaP = sqrt(trace(NAV_SOL(i).P(1:3, 1:3)));

    % Sigma v
    NAV_SOL(i).sigmaV = sqrt(trace(NAV_SOL(i).P(4:6, 4:6)));

    % Computation of the covariance matriox for the keplerian elements
    NAV_SOL(i).Pkep = coovarianceMapping(NAV_SOL(i).x, NAV_SOL(i).P, mu);
end

% Display results
fprintf('\nSTATE\t\tREF STATE\t\tKOUROUN only\t\tALL STATIONS\t\tALL STATIONS with J2')
fprintf('\nx[km]\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t\t%.4f', SMOS_xx0_t0(1), NAV_SOL(1).x(1), NAV_SOL(2).x(1), NAV_SOL(3).x(1))
fprintf('\ny[km]\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t\t%.4f', SMOS_xx0_t0(2), NAV_SOL(1).x(2), NAV_SOL(2).x(2), NAV_SOL(3).x(2))
fprintf('\nz[km]\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t\t%.4f', SMOS_xx0_t0(3), NAV_SOL(1).x(3), NAV_SOL(2).x(3), NAV_SOL(3).x(3))
fprintf('\nvx[km/s]\t%.4f\t\t\t%.4f\t\t\t\t%.4f\t\t\t\t%.4f', SMOS_xx0_t0(4), NAV_SOL(1).x(4), NAV_SOL(2).x(4), NAV_SOL(3).x(4))
fprintf('\nvy[km/s]\t%.4f\t\t\t%.4f\t\t\t\t%.4f\t\t\t\t%.4f', SMOS_xx0_t0(5), NAV_SOL(1).x(5), NAV_SOL(2).x(5), NAV_SOL(3).x(5))
fprintf('\nvz[km/s]\t%.4f\t\t\t%.4f\t\t\t\t%.4f\t\t\t\t%.4f', SMOS_xx0_t0(6), NAV_SOL(1).x(6), NAV_SOL(2).x(6), NAV_SOL(3).x(6))
fprintf('\nsigma_r[km]\t\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL(1).sigmaP, NAV_SOL(2).sigmaP, NAV_SOL(3).sigmaP)
fprintf('\nsigma_v[km/s]\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL(1).sigmaV, NAV_SOL(2).sigmaV, NAV_SOL(3).sigmaV)
fprintf('\nsigma_a[km]\t\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL(1).Pkep(1, 1), NAV_SOL(2).Pkep(1, 1), NAV_SOL(3).Pkep(1, 1))
fprintf('\nsigma_i[deg]\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL(1).Pkep(3, 3) * cspice_dpr, NAV_SOL(2).Pkep(3, 3) * cspice_dpr, NAV_SOL(3).Pkep(3, 3) * cspice_dpr)

%% REQUEST 4

fprintf('\nREQUEST 4\n');

% Set options for lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'StepTolerance', 1e-14, 'FunctionTolerance', 1e-14, 'MaxIterations', 200);

% Call lsqnonlin for KOUROU and SVALBARD with J2
[NAV_SOL_REAL(1).x, NAV_SOL_REAL(1).resnorm, NAV_SOL_REAL(1).residual ,NAV_SOL_REAL(1).exitflag, ~, ~, NAV_SOL_REAL(1).jacobian] =...
    lsqnonlin(@(xx) costFunction(xx, mu, 'J2', [VIS_WIND_KOUROU; VIS_WIND_SVALBARD], [KOUROU; SVALBARD]), SMOS_xx0_t0, [], [], options);%

% Call lsqnonlin for KOUROU and TROLL with J2
[NAV_SOL_REAL(2).x, NAV_SOL_REAL(2).resnorm, NAV_SOL_REAL(2).residual ,NAV_SOL_REAL(2).exitflag, ~, ~, NAV_SOL_REAL(2).jacobian] =...
    lsqnonlin(@(xx) costFunction(xx, mu, 'J2', [VIS_WIND_KOUROU; VIS_WIND_TROLL], [KOUROU; TROLL]), SMOS_xx0_t0, [], [], options);%

% Call lsqnonlin for SVALBARD and TROLL with J2
[NAV_SOL_REAL(3).x, NAV_SOL_REAL(3).resnorm, NAV_SOL_REAL(3).residual ,NAV_SOL_REAL(3).exitflag, ~, ~, NAV_SOL_REAL(3).jacobian] =...
    lsqnonlin(@(xx) costFunction(xx, mu, 'J2', [VIS_WIND_SVALBARD; VIS_WIND_TROLL], [SVALBARD; TROLL]), SMOS_xx0_t0, [], [], options);%

% Elaboration of ther results
for i = 1 : length(NAV_SOL_REAL)

    % Jacobian
    Jac = full(NAV_SOL_REAL(i).jacobian);

    % Computation of the covariance matrix
    NAV_SOL_REAL(i).P = NAV_SOL_REAL(i).resnorm/(length(NAV_SOL_REAL(i).residual) - length(SMOS_xx0_t0)).*inv(Jac.'*Jac);

    % Sigma r
    NAV_SOL_REAL(i).sigmaP = sqrt(trace(NAV_SOL_REAL(i).P(1:3, 1:3)));

    % Sigma v
    NAV_SOL_REAL(i).sigmaV = sqrt(trace(NAV_SOL_REAL(i).P(4:6, 4:6)));

    % Computation of the covariance matriox for the keplerian elements
    NAV_SOL_REAL(i).Pkep = coovarianceMapping(NAV_SOL_REAL(i).x, NAV_SOL_REAL(i).P, mu);
end

% Display results
fprintf('\n\nSTATE\t\tREF STATE\t\tKOUROUN-SVAL\t\tKOUROU-TROLL\t\tSVALBARD-TROLL')
fprintf('\nx[km]\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t\t%.4f', SMOS_xx0_t0(1), NAV_SOL_REAL(1).x(1), NAV_SOL_REAL(2).x(1), NAV_SOL_REAL(3).x(1))
fprintf('\ny[km]\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t\t%.4f', SMOS_xx0_t0(2), NAV_SOL_REAL(1).x(2), NAV_SOL_REAL(2).x(2), NAV_SOL_REAL(3).x(2))
fprintf('\nz[km]\t\t%.4f\t\t%.4f\t\t\t%.4f\t\t\t%.4f', SMOS_xx0_t0(3), NAV_SOL_REAL(1).x(3), NAV_SOL_REAL(2).x(3), NAV_SOL_REAL(3).x(3))
fprintf('\nvx[km/s]\t%.4f\t\t\t%.4f\t\t\t\t%.4f\t\t\t\t%.4f', SMOS_xx0_t0(4), NAV_SOL_REAL(1).x(4), NAV_SOL_REAL(2).x(4), NAV_SOL_REAL(3).x(4))
fprintf('\nvy[km/s]\t%.4f\t\t\t%.4f\t\t\t\t%.4f\t\t\t\t%.4f', SMOS_xx0_t0(5), NAV_SOL_REAL(1).x(5), NAV_SOL_REAL(2).x(5), NAV_SOL_REAL(3).x(5))
fprintf('\nvz[km/s]\t%.4f\t\t\t%.4f\t\t\t\t%.4f\t\t\t\t%.4f', SMOS_xx0_t0(6), NAV_SOL_REAL(1).x(6), NAV_SOL_REAL(2).x(6), NAV_SOL_REAL(3).x(6))
fprintf('\nsigma_r[km]\t\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL_REAL(1).sigmaP, NAV_SOL_REAL(2).sigmaP, NAV_SOL_REAL(3).sigmaP)
fprintf('\nsigma_v[km/s]\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL_REAL(1).sigmaV, NAV_SOL_REAL(2).sigmaV, NAV_SOL_REAL(3).sigmaV)
fprintf('\nsigma_a[km]\t\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL_REAL(1).Pkep(1, 1), NAV_SOL_REAL(2).Pkep(1, 1), NAV_SOL_REAL(3).Pkep(1, 1))
fprintf('\nsigma_i[deg]\t---\t\t\t%.4e\t\t\t%.4e\t\t\t%.4e', NAV_SOL_REAL(1).Pkep(3, 3) * cspice_dpr, NAV_SOL_REAL(2).Pkep(3, 3) * cspice_dpr, NAV_SOL_REAL(3).Pkep(3, 3) * cspice_dpr)

%% REQUEST 5

fprintf('\nREQUEST 5\n');

% Computation of the SMOS orbital period
T = 2 * pi * sqrt(SMOS_orb_elem(1)^3 / mu); 

% Defniition of the new final time
nDays = 8;
tf_new_et = t0_et + cspice_convrt(nDays, 'DAYS', 'SECONDS');

% Computation of the visibility window for KOUROU
[VIS_WIND_KOUROU_REAL] = visibilityWindow(KOUROU, SMOS_xx0_t0, t0_et, tf_new_et, mu, 'J2', SMOS_TLE, SMOS_epoch_et, ddpsi, ddeps);
plotVisbilityWindow(KOUROU, VIS_WIND_KOUROU_REAL, [10 11 12]);

% Plot Orbital revolutions
for i = 1 : round(cspice_convrt(nDays, 'DAYS', 'SECONDS')/T)
    figure(11)
    if i == 1
        xline(t0_et + i * T, 'HandleVisibility', 'on', 'DisplayName', 'i-th SMOS orbital revolution');
    else
        xline(t0_et + i * T, 'HandleVisibility', 'off');
    end
end

% Plot days
for i = 1 : nDays
    figure(11)
    if i == 1
        xline(t0_et + cspice_convrt(i, 'DAYS', 'SECONDS') , 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','on', 'DisplayName', 'j-th Day')
    else
        xline(t0_et + cspice_convrt(i, 'DAYS', 'SECONDS') , 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off')
    end
end

% Computation of the visibility window for TROLL
[VIS_WIND_TROLL_REAL] = visibilityWindow(TROLL, SMOS_xx0_t0, t0_et, tf_new_et, mu, 'J2', SMOS_TLE, SMOS_epoch_et, ddpsi, ddeps);
plotVisbilityWindow(TROLL, VIS_WIND_TROLL_REAL, [13 14 15]);

% Plot Orbital revolutions
for i = 1 : round(cspice_convrt(nDays, 'DAYS', 'SECONDS')/T)
    figure(14)
    if i == 1
        xline(t0_et + i * T, 'HandleVisibility', 'on', 'DisplayName', 'i-th SMOS orbital revolution');
    else
        xline(t0_et + i * T, 'HandleVisibility', 'off');
    end
end

% Plot days
for i = 1 : nDays
    figure(14)
    if i == 1
        xline(t0_et + cspice_convrt(i, 'DAYS', 'SECONDS') , 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','on', 'DisplayName', 'j-th Day')
    else
        xline(t0_et + cspice_convrt(i, 'DAYS', 'SECONDS') , 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off')
    end
end

% Computation of the visibility window for SVALBARD
[VIS_WIND_SVALBARD_REAL] = visibilityWindow(SVALBARD, SMOS_xx0_t0, t0_et, tf_new_et, mu, 'J2', SMOS_TLE, SMOS_epoch_et, ddpsi, ddeps);
plotVisbilityWindow(SVALBARD, VIS_WIND_SVALBARD_REAL, [16 17 18]);

% Plot Orbital revolutions
for i = 1 : round(cspice_convrt(nDays, 'DAYS', 'SECONDS')/T)
    figure(17)
    if i == 1
        xline(t0_et + i * T, 'HandleVisibility', 'on', 'DisplayName', 'i-th SMOS orbital revolution');
    else
        xline(t0_et + i * T, 'HandleVisibility', 'off');
    end
end

% Plot days
for i = 1 : nDays
    figure(17)
    if i == 1
        xline(t0_et + cspice_convrt(i, 'DAYS', 'SECONDS') , 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','on', 'DisplayName', 'j-th Day')
    else
        xline(t0_et + cspice_convrt(i, 'DAYS', 'SECONDS') , 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off')
    end
end

fprintf('\n\t\tKOUROU\t\tTROLL\t\tSVALBARD')
fprintf('\nNpoints\t%d\t\t%d\t\t%d',sum(VIS_WIND_KOUROU_REAL.i_visibility), sum(VIS_WIND_TROLL_REAL.i_visibility), sum(VIS_WIND_SVALBARD_REAL.i_visibility))

close 10 12 13 15 16 18

% Save plots
if SAVEPLOTS

    % Saving 11th figure
    figure(11)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    exportgraphics(gcf,'eleK_req5.png','Resolution',600);

    % Saving 14th figure
    figure(14)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    exportgraphics(gcf,'eleT_req5.png','Resolution',600);

    % Saving 17th figure
    figure(17)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    exportgraphics(gcf,'eleS_req5.png','Resolution',600)
end 
%% FUNCTIONS

%-------------------------------------------------------------------------%
% function: TwoBodyProblem.m  
function [dxdt] = TwoBodyProblem(t, xx, mu, type)
%-------------------------------------------------------------------------%
%                                                                         
% 2BP: function dynamics of 2 Body Problem with the possibility to take into account J2 effect                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [dxdt] = TwoBodyProblem(t, xx, mu, 'Un') 
%   [dxdt] = TwoBodyProblem(t, xx, mu, 'J2') 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  t            [1]     time                                       [s]
%  xx           [6x1]   state vector                               [km;km;km;km/s;km/s;km/s]
%  mu           [1]     gravitational constant                     [km^3/s^2]       
%  type         [1]     type of motion with or without J2          [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  dxdt          [6x1]     ODE        [km/s;km/s;km/s;km/s^2;km/s^2;km/s^2]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Get position vector from the state
rr = xx(1:3); 

% Compute the norm of the position
r = norm(rr); 

% Initialize dxdt equal to zero
dxdt = zeros(6, 1); 

% Initialize acceleration due to J2 effect equal to zero
aJ2 = zeros(3,1);   

% J2 constant
J2 = 0.0010826269;  

% Assign the velocity to the first 3 equation of the ODE
dxdt(1:3) = xx(4:6);

if strcmpi(type, 'Un')

    % Unperturbed case
    aJ2 = zeros(3,1);  

elseif strcmpi(type, 'J2')  

    % Perturbed case
    % Get rotation matrix from ECI to ECEF ref. frame
    R_ECI_ECEF = cspice_pxform('J2000', 'ITRF93', t); 

    % Rotate position into ECEF ref. frame
    rr_ECEF = R_ECI_ECEF * rr;                

    % Compute the norm of the position
    r_ECEF = norm(rr_ECEF);           

    % Get the radius of the Earth from spice
    RE = cspice_bodvrd( 'EARTH', 'RADII', 3);         
    RE = RE(1);      

    % Compute the acceleration due to J2 effect
    aJ2_ECEF = 3/2 * mu * J2 * (RE/r_ECEF)^2 *...     
        [rr_ECEF(1)/r_ECEF^3 * (5*(rr_ECEF(3)/r_ECEF)^2 - 1); rr_ECEF(2)/r_ECEF^3 * (5*(rr_ECEF(3)/r_ECEF)^2 - 1); rr_ECEF(3)/r_ECEF^3 * (5*(rr_ECEF(3)/r_ECEF)^2 - 3)];
    
    % Get rotation matrix from ECEF to ECI ref. frame
    R_ECEF_ECI = cspice_pxform('ITRF93', 'J2000', t); 

    % Rotate the accelartion from ECEF to ECI ref. frame
    aJ2 = R_ECEF_ECI * aJ2_ECEF;                      
else

    fprintf('Error while calling TwoBodyProblem.m'); 
end

% Assign the acceleration due to attraction force and J2 effect to the last 3 equation of the ODE
dxdt(4:6) = -mu/r^3 * rr + aJ2; 

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: propagate2BP.m 
function [xx, tt] = propagate2BP(xx0, tspan, mu, type)
%-------------------------------------------------------------------------%
%                                                                         
% propagate2BP: function for propagating the ode for 2 body problem                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xx, tt] = propagate2BP(x0, tspan, mu, 'Un')
%   [xx, tt] = propagate2BP(x0, tspan, mu, 'J2')
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0           [6x1] initial state                               [km;km;km;km/s;km/s;km/s]
%  tspan        [nx1] time vector                                 [s]
%  mu           [1]   gravitational constant                      [km^3/s^2]
%  type         [1]   type of motion with or without J2           [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  xx          [nx6]   state                                      [km;km;km;km/s;km/s;km/s]
%  tt          [1xn]   time                                       [s]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2023: First version
%
%-------------------------------------------------------------------------%

% Set options
options = odeset('reltol', 1e-13, 'abstol', 1e-14);                           

% Perform integration
[tt, xx] = ode113(@(t,x) TwoBodyProblem(t, x, mu, type), tspan, xx0, options);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: antennaMeasurement.m 
function [azimuth, elevation, range] = antennaMeasurement(station, rr_ECI, vv_ECI, et0)
%-------------------------------------------------------------------------%
%                                                                         
% antennaMEasurements: function for getting measurements wrt to a station, 
% in terms of azimuth, elevation and range, by giving the state of the object
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [azimuth, elevation, range] = antennaMeasurement(stationName, rr_ECI, vv_ECI, et0) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  station      [1]   station struct                            [-]
%  rr_ECI       [3x1] position vector                           [km; km; km]
%  vv_ECI       [3x1] velocity vector                           [km/s; km/s; km/s]
%  et0          [1]   time                                      [s]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  azimuth          [1]   azimuth                               [rad]
%  elevation        [1]   elevation                             [rad]
%  range            [1]   range                                 [km]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Define station name
topoFrame = [station.name, '_TOPO'];          

% Transformation from ECI to topocentric frame
ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et0);  

% Get station state in J2000
rv_station_eci = cspice_spkezr(station.name, et0, 'J2000', 'NONE', 'EARTH'); 

% Compute station-satellite vector in ECI
rv_station_sat_eci = [rr_ECI; vv_ECI] - rv_station_eci;                    

% Convert state into topocentric frame
rv_station_sat_topo = ROT_ECI2TOPO * rv_station_sat_eci; 

% Get azimuth, elevation and range by converting rectangular coordinates into latitudinal ones
[range, azimuth,  elevation] = cspice_reclat(rv_station_sat_topo(1:3)); 

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: visibilityWindow.m 
function [VIS_WIND] = visibilityWindow(station, xx0, t0, tf, mu, type, SAT_TLE, SAT_epoch_et, ddpsi, ddeps)
%-------------------------------------------------------------------------%
%                                                                         
% visibilityWindow: function which computes the visibility window of a
% satellite wrt to a ground station
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [VIS_WIND] = visibilityWindow(station, xx0, t0, tf, mu, 'Un') 
%   [VIS_WIND] = visibilityWindow(station, xx0, t0, tf, mu, 'J2')
%   [VIS_WIND] = visibilityWindow(station, xx0, t0, tf, mu, 'SGP4', SAT_TLE, SAT_epoch_et, ddpsi, ddeps)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  station      [-]   struct with station data                  [-]
%  xx0          [6x1] initial state                             [km;km;km;km/s;km/s;kms/s]
%  t0           [1]   initial time                              [s]
%  tf           [1]   final time                                [s]
%  mu           [1]   gravitational constant                    [km^3/s^2]
%  type         [1]   Type of 2BP with or withouth J2           [-]
%  SAT_TLE      [1]   struct with data of TLE                   [-]
%  SAT_epoch_et [1]   epoch of TLE                              [s]
%  ddpsi        [1]   nutation correction parameter             [rad]
%  ddeps        [1]   nutation correction parameter             [rad]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  VIS_WIND     [-]   struct with visbility window data         [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Get the measurement frequency from the struct
delta_t = station.freq;

% Computation of the time grid
tspan = t0 : delta_t : tf;
nPoints = length(tspan); % Get the number of points in the time vector

% Initializzation of variable for the SAT struct
VIS_WIND.tspan = tspan;
VIS_WIND.azimuth = zeros(nPoints, 1);         % Azimuth
VIS_WIND.elevation = zeros(nPoints, 1);       % Elevation
VIS_WIND.range = zeros(nPoints, 1);           % Range
VIS_WIND.i_visibility.KO = false(nPoints, 1); % indices for visibility window
VIS_WIND.rr = zeros(3, nPoints);              % Position
VIS_WIND.vv = zeros(3, nPoints);              % Velocity

if strcmp(type, 'SGP4')

    % SGP4 propagation
    for i = 1 : nPoints
        
        % Time from SAT_epoch_et in minutes
        tsince = (tspan(i) - SAT_epoch_et)/60.0;

        % Propagate using SGP4 SMOS' state until tsince
        [SAT_TLE, rr_TEME_SMOS, vv_TEME_SMOS] = sgp4(SAT_TLE, tsince);                                
        
        % Compute centuries from 2000-01-01T00:00:00.00 TDT
        ttt = cspice_unitim(tsince, 'ET', 'TDT')/cspice_jyear()/100;
        
        % TEME to ECI conversion
        [SMOS_rr, SMOS_vv] = teme2eci(rr_TEME_SMOS,  vv_TEME_SMOS, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
        
        % Assemble SMOS's the reference state
        VIS_WIND.rr(:, i) = SMOS_rr; 
        VIS_WIND.vv(:, i) = SMOS_vv; 
    end
else
    [rv, ~]= propagate2BP(xx0, tspan, mu, type);   % Get the state of the satellite froM the tspan
    % Assign possition and the velocity to the SAT struct
    VIS_WIND.rr = rv(:,1:3)'; % Position
    VIS_WIND.vv = rv(:,4:6)'; % Velocity
end

% Compute the measurement for every state
for i = 1:nPoints
    [VIS_WIND.azimuth(i), VIS_WIND.elevation(i), VIS_WIND.range(i)] = antennaMeasurement(station, VIS_WIND.rr(:, i), VIS_WIND.vv(:, i), tspan(i));
    
end

% Get indeces for the visibility window by cheching if the elevation resect
% the minimum one required by the ground station
VIS_WIND.i_visibility = VIS_WIND.elevation > station.minEl;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: plotVisbilityWindow.m 
function [] = plotVisbilityWindow(station, VIS_WIND, indicesFigure)
%-------------------------------------------------------------------------%
%                                                                         
% plotVisibilityWindow: function which plots azimuth, elevation, range and
% opening and closing dates for the visibility windows of a station
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [] = plotVisbilityWindow(station, VIS_WIND, indicesFigure)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  station      [-]   struct with station data                  [-]
%  VIS_WIND     [-]   struct with visbility window data         [-]
% indicesFigure [3x1] indice of figures for plots               [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Check if the plots need to be made
if ~isempty(indicesFigure)

    % Plot for Azimuth
    figure(indicesFigure(1))
    plot(VIS_WIND.tspan, VIS_WIND.azimuth *cspice_dpr(), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
    hold on
    grid on;
    xlabel('Epoch J2000 [s]', 'FontSize',12); 
    ylabel('Azimuth [deg]', 'FontSize',12);
    xlim([-inf inf]); 
    ylim([-inf inf]);
    
    % Plot for Elevation
    figure(indicesFigure(2))
    plot(VIS_WIND.tspan, VIS_WIND.elevation*cspice_dpr(), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
    hold on
    if ~strcmp(station.name, 'SVALBARD')
        yline(station.minEl*cspice_dpr(), '--', strcat('Minimum elevation', {' '},num2str(station.minEl*cspice_dpr()),'$^o$'),'Color','black','LineWidth', 1.5, 'Interpreter','latex', 'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    else
        yline(station.minEl*cspice_dpr(), '--', strcat('Minimum elevation', {' '},num2str(station.minEl*cspice_dpr()),'$^o$'),'Color','black','LineWidth', 1.5, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    end
    grid on;
    xlabel('Epoch J2000 [s]','FontSize',12); 
    ylabel('Elevation [deg]', 'FontSize',12);
    xlim([-inf inf]); 
    ylim([-inf inf]);
    
    % Plot Range
    figure(indicesFigure(3))
    plot(VIS_WIND.tspan, VIS_WIND.range, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
    hold on
    grid on;
    xlabel('Epoch J2000 [s]', 'FontSize',12); 
    ylabel('Range [km]', 'FontSize',12); 
    xlim([-inf inf]); 
    ylim([-inf inf]);
end

tk = 0;
fprintf('\n\n%s visibility windows', station.name)

% For cycle for printing the time windows for the station and
% updating the plot for azimuth, elevationa and range
for i = 1 : length(VIS_WIND.tspan) - 1

    % Check if the i-th point is in the visibility window
    if VIS_WIND.i_visibility(i)

        tk = tk + 1;

        % Check if the i-th+1 point is in the visibility window
        if VIS_WIND.i_visibility(i + 1)           
        else

            % Print the visibility window
            fprintf('\n The visbility window for %s opens on %s:', station.name, cspice_et2utc(VIS_WIND.tspan(i- tk + 1), 'ISOC', 3))
            fprintf('\n The visbility window closes for %s on %s:', station.name, cspice_et2utc(VIS_WIND.tspan(i), 'ISOC', 3))
            fprintf('\n Its duration is : %d [s], it contains %d points \n', (tk - 1) * station.freq, tk);

            % Plots
            if ~isempty(indicesFigure)

                % Plots update for Azimuth
                figure(indicesFigure(1))
                plot(VIS_WIND.tspan(i- tk + 1 : i), VIS_WIND.azimuth(i- tk + 1 : i)*cspice_dpr(), 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840])
        
                % Plots update for Elevation
                figure(indicesFigure(2))
                plot(VIS_WIND.tspan(i- tk + 1 : i), VIS_WIND.elevation(i- tk + 1 : i)*cspice_dpr(), 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840])
        
                % Plots update for Range
                figure(indicesFigure(3))
                plot(VIS_WIND.tspan(i- tk + 1 : i), VIS_WIND.range(i- tk + 1 : i), 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840])
            end
            tk = 0;
        end
    end
end

% Check if the last point is in the visibility window
if VIS_WIND.i_visibility(end)
    tk = tk + 1;
end

if tk ~= 0
    
    % Print the visibility window
    fprintf('\n The visbility window for %s opens on %s:', station.name, cspice_et2utc(VIS_WIND.tspan(i - tk + 1), 'ISOC', 3))
    fprintf('\n The visbility window closes for %s on %s:', station.name, cspice_et2utc(VIS_WIND.tspan(i + 1), 'ISOC', 3))
    fprintf('\n Its duration is : %d [s], it contains %d points \n', (tk - 1) * station.freq, tk);

    % Plots
    if ~isempty(indicesFigure)

        % Plots update for Azimuth
        figure(indicesFigure(1))
        plot(VIS_WIND.tspan(i- tk + 1 : i + 1), VIS_WIND.azimuth(i- tk + 1 : i + 1)*cspice_dpr(), 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840])

        % Plots update for Elevation
        figure(indicesFigure(2))
        plot(VIS_WIND.tspan(i- tk + 1 : i + 1), VIS_WIND.elevation(i- tk + 1 : i + 1)*cspice_dpr(), 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840])

        % Plots update for Range
        figure(indicesFigure(3))
        plot(VIS_WIND.tspan(i- tk + 1 : i + 1), VIS_WIND.range(i- tk + 1 : i + 1), 'LineWidth', 2.5, 'Color', [0.6350 0.0780 0.1840])
    end
end

% Add legends  to the plots
if ~isempty(indicesFigure)

    % Azimuth
    figure(indicesFigure(1))
    legend('Azimuth', 'Predicted azimuth', 'Location', 'best')
    
    % Elevation
    figure(indicesFigure(2))
    legend('Elevation', 'Predicted elevation', 'Location', 'best')
    
    % Range
    figure(indicesFigure(3))
    legend('Range', 'Predicted range', 'Location', 'best')
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: simulatedMeasurements.m 
function [VIS_WIND] = simulatedMeasurements(station, VIS_WIND, SAT_epoch_et, SAT_TLE, ddpsi, ddeps)
%-------------------------------------------------------------------------%
%                                                                         
% simulatedMeasurements: function which propagate the state of the satellite with SGP4 and by checking
% the visibility window compute the measuremnts associated to the state and
% then add noises 
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [VIS_WIND] = simulatedMeasurements(station, VIS_WIND, starting_epoch_TLE_et, SAT_TLE, ddpsi, ddeps)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  station             [-]   struct with station data                  [-]
%  VIS_WIND            [-]   struct with visbility window data         [-]
%  SAT_epoch_et        [1]   epoch of TLE                              [s]
%  SAT_TLE             [1]   struct with data of TLE                   [-]
%  ddpsi               [1]   nutation correction parameter             [rad]
%  ddeps               [1]   nutation correction parameter             [rad]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  VIS_WIND            [-]   struct with visbility window data         [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Get the number of points in the time vector
npoints = length(VIS_WIND.tspan);

% Assign time vector to the measurements field
VIS_WIND.measurement.tspan = VIS_WIND.tspan;

% Initializzation of temporary variables for azimuth, elevation and range
measured_azimuth = NaN(npoints, 1);   % Azimuth
measured_elevation = NaN(npoints, 1); % Elevation
measured_range = NaN(npoints, 1);     % Range


% Initializzation of variable for the VIS_WIND struct
VIS_WIND.measurement.azimuth = NaN(npoints, 1);          % Azimuth
VIS_WIND.measurement.elevation = NaN(npoints, 1);        % Elevation
VIS_WIND.measurement.range = NaN(npoints, 1);            % Range
VIS_WIND.measurement.i_visibility = false(npoints, 1);   % Indices for visibility window  

VIS_WIND.measurement.rr = NaN(3, npoints);                    % Position
VIS_WIND.measurement.vv = NaN(3, npoints);                    % Velocity

for i = 1:npoints

    % Chech if in that instant of the time the satellite is in the visibility window
    if VIS_WIND.i_visibility(i) 

        % Minutes from TLE epoch
        tsince = (VIS_WIND.tspan(i) - SAT_epoch_et)/60.0;  

        % SGP4 propagation
        [~, rrTEME_SAT, vvTEME_SAT] = sgp4(SAT_TLE,  tsince);  

        % Compute centuries from 2000-01-01T00:00:00.00 TDT
        ttt = cspice_unitim(VIS_WIND.tspan(i), 'ET', 'TDT')/cspice_jyear()/100;  

        % TEME to ECI conversion
        [rrECI_SAT, vvECI_SAT] = teme2eci(rrTEME_SAT, vvTEME_SAT, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);   
        
        VIS_WIND.measurement.rr(:,i) = rrECI_SAT; % Assign position to the struct
        VIS_WIND.measurement.vv(:,i) = vvECI_SAT; % Assign velocity to the struct
        
        % Compute the measurements
        [measured_azimuth(i), measured_elevation(i), measured_range(i)] = antennaMeasurement(station, rrECI_SAT, vvECI_SAT, VIS_WIND.tspan(i));    
       
    end
end

% Adding noises to the measurements
VIS_WIND.measurement.azimuth = mvnrnd(measured_azimuth, station.sigma_AzEl^2 * eye(npoints));       % Azimuth
VIS_WIND.measurement.elevation = mvnrnd(measured_elevation, station.sigma_AzEl^2 * eye(npoints));   % Elevation
VIS_WIND.measurement.range = mvnrnd(measured_range, station.sigma_r^2 * eye(npoints));              % Range

% Get indeces for the visibility window by cheching if the elevation resect
% the minimum one required by the ground station
VIS_WIND.measurement.i_visibility = VIS_WIND.measurement.elevation > station.minEl;


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: costFunction.m 
function [residuals] = costFunction(xx0, mu, typePropagation, VIS_WIND, station)
%-------------------------------------------------------------------------%
%                                                                         
% costFunction: function used inside lsqnonlinto solve the navigation
% problem
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [residuals] = costFunction(xx0, mu, 'Un', VIS_WIND, station)
%   [residuals] = costFunction(xx0, mu, 'J2', VIS_WIND, station)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx0               [6x1] initial state                             [km;km;km;km/s;km/s;kms/s]
%  mu                [1]   gravitational constant                    [km^3/s^2]
%  typePropagation   [1]   Type of 2BP with or withouth J2           [-]
%  VIS_WIND          [-]   struct with visbility window data         [-]
%  station           [-]   struct with station data                  [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  residual          [mx3]   residuals                       [rad, rad, km]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Get number of the stations
nStations = length(station);     

% Initilizzation of variables
nPoints = zeros(nStations, 1); % Number of points present in the time grid for each station    
nMeas = 0;                     % Number of measurements
WeightMatrix = zeros(3 * nStations, 3); % Weigth matrices

for i = 1 : nStations

    % Get the number of points in the time vector
    nPoints(i) = length(VIS_WIND(i).tspan);
    
    % Get number of measurements for the i-th station
    nMeas = nMeas + sum(VIS_WIND(i).measurement.i_visibility);

    % Get the weigth matrix for each station
    WeightMatrix(3*i - 2: 3 * i, :) = diag([ones(1, 2) * 1/station(i).sigma_AzEl, 1/station(i).sigma_r]);
end

% Initializzation of the residuals
residuals = NaN(nMeas, 3);
k = 1;

% Computation of tehe residuals
for i = 1 : nStations

    % Propagate the state over tspan
    [xx, ~] = propagate2BP(xx0, VIS_WIND(i).tspan, mu, typePropagation);
    
    % Extract the matrix
    W = WeightMatrix(3*i - 2: 3 * i, :);
    
    for j = 1 : nPoints(i)
    
        % Chech if in that instant of the time the satellite is in the visibility window
        if VIS_WIND(i).measurement.i_visibility(j) 
    
            % Compute predicted measurements
            [measurement_predicted_azimuth, measurement_predicted_elevation, measurement_predicted_range] = antennaMeasurement(station(i), xx(j,1:3)', xx(j,4:6)' ,VIS_WIND(i).tspan(j));
    
            % Compute the weighted differences between measurements
            diffMeasWeighted = W * [angdiff(measurement_predicted_azimuth, VIS_WIND(i).measurement.azimuth(j));...
                                                        angdiff(measurement_predicted_elevation, VIS_WIND(i).measurement.elevation(j));...
                                                        (measurement_predicted_range - VIS_WIND(i).measurement.range(j))];
            
            % Update the residual
            residuals(k, :) = diffMeasWeighted';
            k  = k + 1;
        end
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: coovarianceMapping.m 
function Pkep = coovarianceMapping(xx, P, mu)
%-------------------------------------------------------------------------%
%                                                                         
% coovarianceMapping: function which maps the covariance from cartesian to
% keplerian using the linearized matrix evaluated computed as the gradient
% of [a, e, i RAAN, ArgOfPer, theta] w.r.t xx
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   Pkep = coovarianceMapping(xx, P, mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx   [6x1] state                                             [km;km;km;km/s;km/s;kms/s]
%  P    [6x6] Covariance Matrix                                 [km^2 km^2/s; km^2/s (km/s)^2] 
%  mu   [1]   gravitational constant                            [km^3/s^2]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  Pkep [6x6] Covariance Matrix of obsculating orbital elements [-] 
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  03/12/2024: First version
%
%-------------------------------------------------------------------------%

% Extract position from the state
rr = xx(1:3);

% Computation of the norm of the position
r = norm(rr);

% Extract components
rx = rr(1); % x - component
ry = rr(2); % y - component
rz = rr(3); % z - component 

% Extract position from the state
vv = xx(4:6);

% Computation of the norm of the velocity
v = norm(vv);

% Extract components
vx = vv(1); % x - component
vy = vv(2); % y - component
vz = vv(3); % z - component

% derivative of semimajoraxis w.r.t xx
dadxx = [...
    (2*rx)/(r^3*((v^2)/mu - 2/r)^2),...
    (2*ry)/(r^3*((v^2)/mu - 2/r)^2),...
    (2*rz)/(r^3*((v^2)/mu - 2/r)^2),...
    (2*vx)/(mu*((v^2)/mu - 2/r)^2),...
    (2*vy)/(mu*((v^2)/mu - 2/r)^2),...
    (2*vz)/(mu*((v^2)/mu - 2/r)^2)];
 
% derivative of eccentricity w.r.t xx
dedxx = [...
    (2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*((rx*ry)/r^3 - (vx*vy)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*((rx*rz)/r^3 - (vx*vz)/mu) + 2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx^2/r^3 + (vy^2 + vz^2)/mu - 1/r))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)),...
    (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*((rx*ry)/r^3 - (vx*vy)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*((ry*rz)/r^3 - (vy*vz)/mu) + 2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry^2/r^3 + (vx^2 + vz^2)/mu - 1/r))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)),... 
    (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*((rx*rz)/r^3 - (vx*vz)/mu) + 2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*((ry*rz)/r^3 - (vy*vz)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rz^2/r^3 + (vx^2 + vy^2)/mu - 1/r))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)),... 
    -((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vy - 2*ry*vx))/mu + (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(ry*vy + rz*vz))/mu + (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rx*vz - 2*rz*vx))/mu)/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)),...
    -((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vx + rz*vz))/mu - (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(2*rx*vy - ry*vx))/mu + (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(ry*vz - 2*rz*vy))/mu)/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)),... 
    ((2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(2*rx*vz - rz*vx))/mu - (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rx*vx + ry*vy))/mu + (2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(2*ry*vz - rz*vy))/mu)/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2))];

% derivative of inclination w.r.t xx
didxx = [...
    ((ry*vz - rz*vy)*(rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy))/((1 - (rx*vy - ry*vx)^2/((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    -((rx*vz - rz*vx)*(rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy))/((1 - (rx*vy - ry*vx)^2/((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    ((rx*vy - ry*vx)*(2*rz*vx^2 - 2*rx*vz*vx + 2*rz*vy^2 - 2*ry*vz*vy))/(2*(1 - (rx*vy - ry*vx)^2/((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    ((ry*vz - rz*vy)*(vz*rx^2 - rz*vx*rx + vz*ry^2 - rz*vy*ry))/((1 - (rx*vy - ry*vx)^2/((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    -((rx*vz - rz*vx)*(vz*rx^2 - rz*vx*rx + vz*ry^2 - rz*vy*ry))/((1 - (rx*vy - ry*vx)^2/((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    ((rx*vy - ry*vx)*(2*vz*rx^2 - 2*rz*vx*rx + 2*vz*ry^2 - 2*rz*vy*ry))/(2*(1 - (rx*vy - ry*vx)^2/((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vy - ry*vx)^2 + (rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2))];
 
% derivative of RAAN w.r.t xx
dRAANdxx = [...
    -(vz*(ry*vz - rz*vy)^2)/((1 - (rx*vz - rz*vx)^2/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),... 
    (vz*(rx*vz - rz*vx)*(ry*vz - rz*vy))/((1 - (rx*vz - rz*vx)^2/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    -(vz*(rx*vy - ry*vx)*(ry*vz - rz*vy))/((1 - (rx*vz - rz*vx)^2/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    (rz*(ry*vz - rz*vy)^2)/((1 - (rx*vz - rz*vx)^2/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    -(rz*(rx*vz - rz*vx)*(ry*vz - rz*vy))/((1 - (rx*vz - rz*vx)^2/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)),...
    (rz*(rx*vy - ry*vx)*(ry*vz - rz*vy))/((1 - (rx*vz - rz*vx)^2/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2))^(1/2)*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2))];

% derivative of ArgOfPer w.r.t xx
dargOfPerdxx = [...
    (r*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)^(1/2)*abs(mu)*(((vz*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx)^2)/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2) - (vz*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) - ((rx*vz - rz*vx)*(rx^2/r^3 + (vy^2 + vz^2)/mu - 1/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (vz*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vz - rz*vx)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2) + ((vx*vy*r^3 - mu*rx*ry)*(ry*vz - rz*vy))/(mu*r^3*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)))/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2) + (((((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2))*(2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*((rx*ry)/r^3 - (vx*vy)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*((rx*rz)/r^3 - (vx*vz)/mu) + 2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx^2/r^3 + (vy^2 + vz^2)/mu - 1/r)))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2))))/(mu^2*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)*(rx^2 + ry^2 + rz^2) - (r*rx^2*vy^2*vz + r*rx^2*vz^3 - mu*rx^2*vz - 2*r*rx*ry*vx*vy*vz - 2*r*rx*rz*vx*vz^2 + mu*rx*rz*vx + r*ry^2*vx^2*vz + r*ry^2*vz^3 - mu*ry^2*vz - 2*r*ry*rz*vy*vz^2 + mu*ry*rz*vy + r*rz^2*vx^2*vz + r*rz^2*vy^2*vz)^2/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2))^(1/2),... 
    (r*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)^(1/2)*abs(mu)*(((vz*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy)^2)/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2) - (vz*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) - ((ry*vz - rz*vy)*(ry^2/r^3 + (vx^2 + vz^2)/mu - 1/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (vz*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2) + ((vx*vy*r^3 - mu*rx*ry)*(rx*vz - rz*vx))/(mu*r^3*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)))/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2) + (((((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2))*(2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*((rx*ry)/r^3 - (vx*vy)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*((ry*rz)/r^3 - (vy*vz)/mu) + 2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry^2/r^3 + (vx^2 + vz^2)/mu - 1/r)))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2))))/(mu^2*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)*(rx^2 + ry^2 + rz^2) - (r*rx^2*vy^2*vz + r*rx^2*vz^3 - mu*rx^2*vz - 2*r*rx*ry*vx*vy*vz - 2*r*rx*rz*vx*vz^2 + mu*rx*rz*vx + r*ry^2*vx^2*vz + r*ry^2*vz^3 - mu*ry^2*vz - 2*r*ry*rz*vy*vz^2 + mu*ry*rz*vy + r*rz^2*vx^2*vz + r*rz^2*vy^2*vz)^2/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2))^(1/2),... 
    (r*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)^(1/2)*abs(mu)*(((vx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (vy*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx)*(2*rz*vx^2 - 2*rx*vz*vx + 2*rz*vy^2 - 2*ry*vz*vy))/(2*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy)*(2*rz*vx^2 - 2*rx*vz*vx + 2*rz*vy^2 - 2*ry*vz*vy))/(2*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)) + ((vx*vz*r^3 - mu*rx*rz)*(rx*vz - rz*vx))/(mu*r^3*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)) + ((vy*vz*r^3 - mu*ry*rz)*(ry*vz - rz*vy))/(mu*r^3*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)))/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2) + (((((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2))*(2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*((rx*rz)/r^3 - (vx*vz)/mu) + 2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*((ry*rz)/r^3 - (vy*vz)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rz^2/r^3 + (vx^2 + vy^2)/mu - 1/r)))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2))))/(mu^2*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)*(rx^2 + ry^2 + rz^2) - (r*rx^2*vy^2*vz + r*rx^2*vz^3 - mu*rx^2*vz - 2*r*rx*ry*vx*vy*vz - 2*r*rx*rz*vx*vz^2 + mu*rx*rz*vx + r*ry^2*vx^2*vz + r*ry^2*vz^3 - mu*ry^2*vz - 2*r*ry*rz*vy*vz^2 + mu*ry*rz*vy + r*rz^2*vx^2*vz + r*rz^2*vy^2*vz)^2/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2))^(1/2),... 
    (r*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)^(1/2)*abs(mu)*(((rz*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + ((rx*vy - 2*ry*vx)*(ry*vz - rz*vy))/(mu*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)) + ((rx*vz - rz*vx)*(ry*vy + rz*vz))/(mu*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)) - (rz*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx)^2)/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2) - (rz*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vz - rz*vx)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2))/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2) - (((((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2))*((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vy - 2*ry*vx))/mu + (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(ry*vy + rz*vz))/mu + (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rx*vz - 2*rz*vx))/mu))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2))))/(mu^2*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)*(rx^2 + ry^2 + rz^2) - (r*rx^2*vy^2*vz + r*rx^2*vz^3 - mu*rx^2*vz - 2*r*rx*ry*vx*vy*vz - 2*r*rx*rz*vx*vz^2 + mu*rx*rz*vx + r*ry^2*vx^2*vz + r*ry^2*vz^3 - mu*ry^2*vz - 2*r*ry*rz*vy*vz^2 + mu*ry*rz*vy + r*rz^2*vx^2*vz + r*rz^2*vy^2*vz)^2/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2))^(1/2),... 
    -(r*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)^(1/2)*abs(mu)*((((2*rx*vy - ry*vx)*(rx*vz - rz*vx))/(mu*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)) - (rz*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) - ((rx*vx + rz*vz)*(ry*vz - rz*vy))/(mu*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)) + (rz*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy)^2)/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2) + (rz*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2))/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2) + (((((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2))*((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vx + rz*vz))/mu - (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(2*rx*vy - ry*vx))/mu + (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(ry*vz - 2*rz*vy))/mu))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2))))/(mu^2*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)*(rx^2 + ry^2 + rz^2) - (r*rx^2*vy^2*vz + r*rx^2*vz^3 - mu*rx^2*vz - 2*r*rx*ry*vx*vy*vz - 2*r*rx*rz*vx*vz^2 + mu*rx*rz*vx + r*ry^2*vx^2*vz + r*ry^2*vz^3 - mu*ry^2*vz - 2*r*ry*rz*vy*vz^2 + mu*ry*rz*vy + r*rz^2*vx^2*vz + r*rz^2*vy^2*vz)^2/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2))^(1/2),... 
    -(r*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)^(1/2)*abs(mu)*(((rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) - (((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx)*(2*vz*rx^2 - 2*rz*vx*rx + 2*vz*ry^2 - 2*rz*vy*ry))/(2*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)) - (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy)*(2*vz*rx^2 - 2*rz*vx*rx + 2*vz*ry^2 - 2*rz*vy*ry))/(2*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(3/2)) + ((rx*vz - rz*vx)*(2*rx*vz - rz*vx))/(mu*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)) + ((ry*vz - rz*vy)*(2*ry*vz - rz*vy))/(mu*((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2)))/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2) - (((((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx*vz - rz*vx))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2) + (((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry*vz - rz*vy))/((rx*vz - rz*vx)^2 + (ry*vz - rz*vy)^2)^(1/2))*((2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(2*rx*vz - rz*vx))/mu - (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rx*vx + ry*vy))/mu + (2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(2*ry*vz - rz*vy))/mu))/(2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2))))/(mu^2*abs(rx^2*vz^2 - 2*rx*rz*vx*vz + ry^2*vz^2 - 2*ry*rz*vy*vz + rz^2*vx^2 + rz^2*vy^2)*(rx^2 + ry^2 + rz^2) - (r*rx^2*vy^2*vz + r*rx^2*vz^3 - mu*rx^2*vz - 2*r*rx*ry*vx*vy*vz - 2*r*rx*rz*vx*vz^2 + mu*rx*rz*vx + r*ry^2*vx^2*vz + r*ry^2*vz^3 - mu*ry^2*vz - 2*r*ry*rz*vy*vz^2 + mu*ry*rz*vy + r*rz^2*vx^2*vz + r*rz^2*vy^2*vz)^2/(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2))^(1/2)];

% derivative of true anomaly w.r.t xx
dthetadxx = [...
    ((- 2*r^3*rx*vy^2 - 2*r^3*rx*vz^2 + 2*vx*r^3*ry*vy + 2*vx*r^3*rz*vz + mu*rx^3 + mu*rx*ry^2 + mu*rx*rz^2)/(mu*r^4*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + (rx*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(r^3*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + ((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*((rx*ry)/r^3 - (vx*vy)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*((rx*rz)/r^3 - (vx*vz)/mu) + 2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(rx^2/r^3 + (vy^2 + vz^2)/mu - 1/r))*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(2*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2)))/(1 - (rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r))^2/(r^2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)))^(1/2),... 
    ((2*vy*r^3*rx*vx - 2*r^3*ry*vx^2 - 2*r^3*ry*vz^2 + 2*vy*r^3*rz*vz + mu*rx^2*ry + mu*ry^3 + mu*ry*rz^2)/(mu*r^4*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + (ry*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(r^3*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + ((2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*((rx*ry)/r^3 - (vx*vy)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*((ry*rz)/r^3 - (vy*vz)/mu) + 2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(ry^2/r^3 + (vx^2 + vz^2)/mu - 1/r))*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(2*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2)))/(1 - (rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r))^2/(r^2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)))^(1/2),... 
    ((2*vz*r^3*rx*vx + 2*vz*r^3*ry*vy - 2*r^3*rz*vx^2 - 2*r^3*rz*vy^2 + mu*rx^2*rz + mu*ry^2*rz + mu*rz^3)/(mu*r^4*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + (rz*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(r^3*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + ((2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*((rx*rz)/r^3 - (vx*vz)/mu) + 2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*((ry*rz)/r^3 - (vy*vz)/mu) + 2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rz^2/r^3 + (vx^2 + vy^2)/mu - 1/r))*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(2*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2)))/(1 - (rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r))^2/(r^2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)))^(1/2),... 
    -((2*(vx*ry^2 - rx*vy*ry + vx*rz^2 - rx*vz*rz))/(mu*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + (((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vy - 2*ry*vx))/mu + (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(ry*vy + rz*vz))/mu + (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rx*vz - 2*rz*vx))/mu)*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(2*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2)))/(1 - (rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r))^2/(r^2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)))^(1/2),... 
    -((2*(vy*rx^2 - ry*vx*rx + vy*rz^2 - ry*vz*rz))/(mu*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) + (((2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(rx*vx + rz*vz))/mu - (2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(2*rx*vy - ry*vx))/mu + (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(ry*vz - 2*rz*vy))/mu)*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(2*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2)))/(1 - (rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r))^2/(r^2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)))^(1/2),... 
    -((2*(vz*rx^2 - rz*vx*rx + vz*ry^2 - rz*vy*ry))/(mu*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(1/2)) - (((2*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)*(2*rx*vz - rz*vx))/mu - (2*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)*(rx*vx + ry*vy))/mu + (2*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)*(2*ry*vz - rz*vy))/mu)*(rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)))/(2*r*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)^(3/2)))/(1 - (rx*((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r) + ry*((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r) + rz*((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r))^2/(r^2*(((rx*vy^2 - ry*vx*vy + rx*vz^2 - rz*vx*vz)/mu - rx/r)^2 + ((ry*vx^2 - rx*vy*vx + ry*vz^2 - rz*vy*vz)/mu - ry/r)^2 + ((rz*vx^2 - rx*vz*vx + rz*vy^2 - ry*vz*vy)/mu - rz/r)^2)))^(1/2)];
 
% Assemble the jacobian
J = [dadxx; dedxx; didxx; dRAANdxx; dargOfPerdxx; dthetadxx];

% Computation of the Covariance matrix associated to the keplerian elements
Pkep = J * P * J';
 
end
%-------------------------------------------------------------------------%