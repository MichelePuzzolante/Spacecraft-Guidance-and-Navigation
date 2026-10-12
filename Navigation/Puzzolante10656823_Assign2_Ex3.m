% Guidance and Navigation (2024/2025)
% Assignment # 2
% Exercise 3
% Author : Michele Puzzolante
%% REQUEST 1

% Initializzation of the MATLAB enviroment
clear;
close all;
clc;
cspice_kclear();

cspice_furnsh('assignment02.tm');               % Load kernels

% Set default options for plots
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');

fprintf('\nEXERCISE 3\n');
fprintf('\nREQUEST 1\n');

SAVEPLOTS = false;  % Flag to use when plots needs to be saved as png file

% Initializzation of constatnts and structs 
mu = cspice_bodvrd('MOON','GM',1); % Intialization of global variables
t0 = '2024-11-18T16:30:00.000';     % Initial time in UTC
tf = '2024-11-18T20:30:00.000';     % Final time in UTC
t0_et = cspice_str2et(t0);          % Initial time in ephemerides time
tf_et = cspice_str2et(tf);          % Final time in ephemerides time
freq = 30;

% Initializzation of the variables for the LUNAR ORBITER
LUNAR_ORBITER.rr0 = [4307.844185282820; -1317.980749248651; 2109.210101634011];   % Initial position   [km]
LUNAR_ORBITER.vv0 = [-0.110997301537882; -0.509392750828585; 0.815198807994189];  % Initial velocity   [km/s]
LUNAR_ORBITER.xx0 = [LUNAR_ORBITER.rr0; LUNAR_ORBITER.vv0];                       % Initial state      [km;km/s]
LUNAR_ORBITER.sigma_r = 100 * 1e-3;                                               % Measurements noise [km] 
LUNAR_ORBITER.P0 = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]);              % Initial covariance

% Initiailizzation of the struct for the LUNAR_LANDER
LUNAR_LANDER.name = 'MOONLANDER';                     % Station name
LUNAR_LANDER.lat = 78.0 * cspice_rpd;                 % Latitude                     [deg]
LUNAR_LANDER.lon = 15.0 * cspice_rpd;                 % Longitude                    [deg]
LUNAR_LANDER.alt = 0.00;                              % Altitude                     [km]
LUNAR_LANDER.minEl = 0 * cspice_rpd;                  % Minimum elevation            [rad]
LUNAR_ORBITER.noiseMatrix = diag(ones(4, 1) * LUNAR_ORBITER.sigma_r^2);

% Initailizzation of the time vector
tspan = t0_et : freq : tf_et;

% Propagation of the LUNAR ORBITER
[xx, tt] = propagate2BP(LUNAR_ORBITER.xx0, tspan, mu);

% Extract the radius of the Moon from the jernel
R_Moon = cspice_bodvrd( 'MOON', 'RADII', 3);
RM = R_Moon(1);

% Computation of the Moon's flatness
f_Moon = (RM - R_Moon(3)) / RM;

% Conversion of the lander coordinates from geotedic to rectangular
rr_MOONLANDER_IAU_MOON = cspice_georec(LUNAR_LANDER.lon, LUNAR_LANDER.lat, LUNAR_LANDER.alt, RM, f_Moon);

% Extract the rotatio  matrix from IAU MOON to J2000
R_IAUMOON2J2000 = cspice_sxform('IAU_MOON', 'J2000', t0_et);

% Rotation of the position of the lander in J2000 reference frame
rrMOONLANDER_J2000 = R_IAUMOON2J2000 * [rr_MOONLANDER_IAU_MOON; 0; 0; 0];

% Initializzation of the visibility window struct
VIS_WIND.estimated.elevation = NaN(length(tspan), 1);

% Estimation of the visibility window
for i = 1 : length(tspan)

    % Extract the rotation matrix from the kernels from J2000 to MOON
    % CENTERED FIXED frame
    R_J20002IAUMOON = cspice_sxform('J2000', 'IAU_MOON', tspan(i));

    % Roatation of the orbiter position into IAU MOON
    xx_IAU_MOON = R_J20002IAUMOON * xx(i, :)';

    % Computation of the relative position between Lander and Orbiter in IAU MOON
    rr_relative_LANDER_ORBITER = xx_IAU_MOON(1:3) - rr_MOONLANDER_IAU_MOON;

    % Extract rotation matrix from IAU MOON to reference frame fixed with
    % the LANDER
    R_IAUMOON2TOPO = cspice_pxform('IAU_MOON', 'MOONLANDER_TOPO', tspan(i));

    % Rotation of the raltive position in the Lander reference frame
    rr_MOONLANDERETOPO = R_IAUMOON2TOPO * rr_relative_LANDER_ORBITER;

    % Computation of the elevation
    [~, ~,  VIS_WIND.estimated.elevation(i)] = cspice_reclat(rr_MOONLANDERETOPO); % Get azimuth, elevation and range by converting rectangular coordinates into latitudinal ones
end

% Check the visibility window
VIS_WIND.estimated.i_visibility = VIS_WIND.estimated.elevation > LUNAR_LANDER.minEl;

% Print results
fprintf('\nLANDER STATE')
fprintf('\nx\t%.8f\t\t[km]', rrMOONLANDER_J2000(1))
fprintf('\ny\t%.8f\t[km]',rrMOONLANDER_J2000(2))
fprintf('\nz\t%.8f\t[km]', rrMOONLANDER_J2000(3))
fprintf('\nvx\t%.8f\t\t[km/s]', rrMOONLANDER_J2000(4))
fprintf('\nvy\t%.8f\t\t[km/s]', rrMOONLANDER_J2000(5))
fprintf('\nvz\t%.8f\t\t[km/s]', rrMOONLANDER_J2000(6))

% Creation of a sphere for simulating Moon
[X,Y,Z] = sphere(40);
X = RM * X;
Y = RM * Y;
Z = RM * Z;

% Plot the Orbiter and the lander
figure(1)
surf(X, Y, Z, 'FaceColor',[0.8, 0.8, 0.8], 'EdgeAlpha', 1)
hold on
plot3(xx(1, 1), xx(1, 2), xx(1, 3), 'LineWidth', 2, 'LineStyle', 'none', 'MarkerSize', 20, 'Marker', '.')
plot3(xx(:, 1), xx(:, 2), xx(:, 3), 'LineWidth', 2)
plot3(xx(end, 1), xx(end, 2), xx(end, 3), 'LineWidth', 2, 'LineStyle', 'none', 'MarkerSize', 20, 'Marker', '.')
plot3(rrMOONLANDER_J2000(1), rrMOONLANDER_J2000(2), rrMOONLANDER_J2000(3), 'LineWidth', 2, 'LineStyle', 'none', 'MarkerSize', 30, 'Marker', '.')
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
legend('Moon', 'Orbiter position at $t_0$', 'Orbiter trajectory', 'Orbiter position at $t_f$', 'Lander position at $t_0$')
axis equal

% Plot the elevation
figure(2)
plot(tspan, VIS_WIND.estimated.elevation * cspice_dpr, 'LineWidth', 2)
hold on
yline(LUNAR_LANDER.minEl*cspice_dpr(), '--', 'Minimum elevation 0$^o$','LineWidth', 1.5, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'Color', ' black');
xlabel('Epoch J2000 [s]');
ylabel('Elevation [deg]');
xlim([-Inf Inf]);
ylim([-Inf Inf])

% Save plots
if SAVEPLOTS

    % Saving 1st figure
    figure(1)
    exportgraphics(gcf,'trajectory.png','Resolution',600);

    % Saving 2nd figure
    figure(2)
    exportgraphics(gcf,'estimated_el.png','Resolution',600);
end
%% REQUEST 2

fprintf('\nREQUEST 2\n');

% Simulation of the measurements
[VIS_WIND] = simulatedMeasurements(VIS_WIND, LUNAR_LANDER, LUNAR_ORBITER, tspan, mu);

% Computation of the measurements errors
error_x = (VIS_WIND.real.rr(1, :) - VIS_WIND.measurements.rr(1, :));    % x - position
error_y = (VIS_WIND.real.rr(2, :) - VIS_WIND.measurements.rr(2, :));    % y - position
error_z = (VIS_WIND.real.rr(3, :) - VIS_WIND.measurements.rr(3, :));    % z - psotion
error_range = (VIS_WIND.real.range - VIS_WIND.measurements.range);      % range

% Plot Orbiter x-position measurements
figure(3)
plot(tspan, VIS_WIND.measurements.rr(1, :), 'LineWidth', 2)
xlabel('Epoch J2000 [s]');
ylabel('$r_x$ [km]', 'FontSize', 12);
xlim([-Inf Inf]);

% Plot Orbiter y-position measurements
figure(4)
plot(tspan, VIS_WIND.measurements.rr(2, :), 'LineWidth', 2)
xlabel('Epoch J2000 [s]');
ylabel('$r_y$ [km]', 'FontSize', 12);
xlim([-Inf Inf]);

% Plot Orbiter z-position measurements
figure(5)
plot(tspan, VIS_WIND.measurements.rr(3, :), 'LineWidth', 2)
xlabel('Epoch J2000 [s]');
ylabel('$r_z$ [km]', 'FontSize', 12);
xlim([-Inf Inf]);

% Plot the range
figure(6)
plot(tspan, VIS_WIND.measurements.range, 'LineWidth', 2)
xlabel('Epoch J2000 [s]');  
ylabel('$\rho$ [km]', 'FontSize', 12);
xlim([-Inf Inf]);

% Plot errors
figure(7)
plot(tspan, error_x, 'LineStyle', 'none', 'Marker','.', 'MarkerSize', 10, 'Color', [0.6350 0.0780 0.1840])
hold on
yline(LUNAR_ORBITER.sigma_r, '--', '$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
yline(-LUNAR_ORBITER.sigma_r, '--', '-$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
xlabel('Epoch J2000 [s]');  
ylabel('$\sigma_rx$ [km]', 'FontSize', 12);

% Plot errors
figure(8)
plot(tspan, error_y, 'LineStyle', 'none', 'Marker','.', 'MarkerSize', 10, 'Color', [0.6350 0.0780 0.1840])
hold on
yline(LUNAR_ORBITER.sigma_r, '--', '$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
yline(-LUNAR_ORBITER.sigma_r, '--', '-$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
xlabel('Epoch J2000 [s]');  
ylabel('$\sigma_ry$ [km]', 'FontSize', 12);

% Plot errors
figure(9)
plot(tspan, error_z, 'LineStyle', 'none', 'Marker','.', 'MarkerSize', 10, 'Color', [0.6350 0.0780 0.1840])
hold on
yline(LUNAR_ORBITER.sigma_r, '--', '$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
yline(-LUNAR_ORBITER.sigma_r, '--', '-$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
xlabel('Epoch J2000 [s]');  
ylabel('$\sigma_rz$ [km]', 'FontSize', 12);

% Plot errors
figure(10)
plot(tspan, error_range, 'LineStyle', 'none', 'Marker','.', 'MarkerSize', 10, 'Color', [0.6350 0.0780 0.1840])
hold on
yline(LUNAR_ORBITER.sigma_r, '--', '$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
yline(-LUNAR_ORBITER.sigma_r, '--', '-$\sigma_p$','LineWidth', 2, 'Interpreter','latex', 'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','top', 'FontSize', 12)
xlabel('Epoch J2000 [s]');  
ylabel('$\sigma_{\rho}$ [km]', 'FontSize', 12);

% Save plots
if SAVEPLOTS

    % Saving 3rd figure
    figure(3)
    exportgraphics(gcf,'rx_meas.png','Resolution',600)

    % Saving 4th figure
    figure(4)
    exportgraphics(gcf,'ry_meas.png','Resolution',600)

    % Saving 5th figure
    figure(5)
    exportgraphics(gcf,'rz_meas.png','Resolution',600)

    % Saving 6th figure
    figure(6)
    exportgraphics(gcf,'rng_meas.png','Resolution',600)

    % Saving 7th figure
    figure(7)
    exportgraphics(gcf,'err_rx_meas.png','Resolution',600)

    % Saving 8th figure
    figure(8)
    exportgraphics(gcf,'err_ry_meas.png','Resolution',600)

    % Saving 9th figure
    figure(9)
    exportgraphics(gcf,'err_rz_meas.png','Resolution',600)

    % Saving 10th
    figure(10)
    exportgraphics(gcf, 'err_rng_meas.png','Resolution', 600)
end

%% REQUEST 3

fprintf('\nREQUEST 3\n');

% Call UKF 
[x_hat_ORB, P_ORB] = UKF_ORBITER(VIS_WIND, LUNAR_ORBITER, tspan, mu, 0.01, 2, 0);

% Compute the error in position and velocity and the 3 sigma for position
% and velocity
[errPos_ORB, sigma3Pos_ORB, errVel_ORB, sigma3Vel_ORB] = Errorand3SigmaComputation(x_hat_ORB, P_ORB, VIS_WIND.real.rr, VIS_WIND.real.vv);

% Plot errrors with 3sigma position
figure(11)
semilogy(tspan, errPos_ORB, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
hold on
semilogy(tspan, sigma3Pos_ORB, 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
grid on
xlabel('Epoch J2000 [s]', 'FontSize',12); 
ylabel('[km]', 'FontSize',12); 
legend('Error', '$3 \sigma$')
xlim([-inf inf]); 
ylim([-inf inf]);

% Plot errrors with 3sigma velocity
figure(12)
semilogy(tspan, errVel_ORB, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
hold on
semilogy(tspan, sigma3Vel_ORB, 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
grid on
xlabel('Epoch J2000 [s]', 'FontSize',12); 
ylabel('[km/s]', 'FontSize',12); 
legend('Error', '$3 \sigma$')
xlim([-inf inf]); ylim([-inf inf]);

% Save plots
if SAVEPLOTS

    % Saving 11th figure
    figure(11)
    exportgraphics(gcf,'sigma_r_req3.png','Resolution',600);

    % Saving 12th figure
    figure(12)
    exportgraphics(gcf,'sigma_v_req3.png','Resolution',600);
end
%% REQUEST 4

fprintf('\nREQUEST 4\n');

% Call UFK
[x_hat_ORB_LAND, P_ORB_LAND] = UKF_ORBITER_LANDER(VIS_WIND, LUNAR_ORBITER, LUNAR_LANDER, tspan, mu, 0.01, 2, 0);

% Extraction of the covariance matrices
for i = 1 : length(tspan)

    % Orbiter state
    P_ORB_2(1:6, 1 +(i-1)*6:6 +(i-1)*6) = P_ORB_LAND(1:6, 1 +(i-1)*8:6 +(i-1)*8);

    % Lander coordinates
    P_LAND(1:2, 1 +(i-1)*2:2 +(i-1)*2) = P_ORB_LAND(7:8, 7 +(i-1)*8:8 +(i-1)*8);
end

% Compute the error in position and velocity and the 3 sigma for position
% and velocity
[errPos_ORB_2, sigma3Pos_ORB_2, errVel_ORB_2, sigma3Vel_ORB_2] = Errorand3SigmaComputation(x_hat_ORB_LAND(1:6, :), P_ORB_2, VIS_WIND.real.rr, VIS_WIND.real.vv);

% Compute the error lander coordinates
[errLONLAT, sigma3LONLAT] = Errorand3SigmaComputationLATLON(x_hat_ORB_LAND(7:8, :), P_LAND, tspan);

% Plot errrors with 3sigma position
figure(13)
semilogy(tspan, errPos_ORB_2, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
hold on
semilogy(tspan, sigma3Pos_ORB_2, 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
grid on
xlabel('Epoch J2000 [s]', 'Interpreter','latex', 'FontSize',12); ylabel('[km]', 'Interpreter','latex', 'FontSize',12); 
legend('Error', '$3 \sigma$')
xlim([-inf inf]); ylim([-inf inf]);
 
% Plot errrors with 3sigma Velocity
figure(14)
semilogy(tspan, errVel_ORB_2, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
hold on
semilogy(tspan, sigma3Vel_ORB_2, 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
grid on
xlabel('Epoch J2000 [s]', 'FontSize',12); 
ylabel('[km/s]',  'FontSize',12); 
legend('Error', '$3 \sigma$')
xlim([-inf inf]); ylim([-inf inf]);

% Plot errrors with 3sigma of lander coordinates
figure(15)
semilogy(tspan, errLONLAT * cspice_dpr, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
hold on
semilogy(tspan, sigma3LONLAT * cspice_dpr, 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
grid on
xlabel('Epoch J2000 [s]', 'FontSize',12); 
ylabel('[deg]', 'Interpreter','latex', 'FontSize',12); 
legend('Error', '$3 \sigma$')
xlim([-inf inf]); ylim([-inf inf]);

% Save plots
if SAVEPLOTS

    % Saving 13th figure
    figure(13)
    exportgraphics(gcf,'sigma_r_req4.png','Resolution',600);

    % Saving 14th figure
    figure(14)
    exportgraphics(gcf,'sigma_v_req4.png','Resolution',600);

    % Saving 15th figure
    figure(15)
    exportgraphics(gcf,'sigma_latlon_req4.png','Resolution',600)
end

%% FUNCTIONS

%-------------------------------------------------------------------------%
% function: TwoBodyProblem.m  
function [dxdt] = TwoBodyProblem(~, xx, mu)
%-------------------------------------------------------------------------%
%                                                                         
% 2BP: function dynamics of 2 Body Problem                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [dxdt] = TwoBodyProblem(t, xx, mu) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  t            [1]     time                       [s]
%  xx           [6x1]   state vector               [km;km;km;km/s;km/s;km/s]
%  mu           [1]     gravitational constant     [km^3/s^2]       
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  dxdt          [6x1]     ODE                     [km/s;km/s;km/s;km/s^2;km/s^2;km/s^2]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  15/12/2024: First version
%
%-------------------------------------------------------------------------%

% Get position vector from the state
rr = xx(1:3); 

% Compute the norm of the position
r = norm(rr); 

% Initialize dxdt equal to zero
dxdt = zeros(6, 1);

% Assign the velocity to the first 3 equation of the ODE
dxdt(1:3) = xx(4:6); 

% Assign the acceleration due to attraction force and J2 effect to the last 3 equation of the ODE
dxdt(4:6) = -mu/r^3 * rr; 

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: propagate2BP.m 
function [xx, tt] = propagate2BP(xx0, tspan, mu)
%-------------------------------------------------------------------------%
%                                                                         
% propagate2BP: function for propagating the ode for 2 body problem                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xx, tt] = propagate2BP(x0, tspan, mu)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx           [6x1] initial state                               [km;km;km;km/s;km/s;km/s]
%  tspan        [nx1] time vector                                 [s]
%  mu           [1]   gravitational constant                      [km^3/s^2]
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
[tt, xx] = ode113(@(t,x) TwoBodyProblem(t, x, mu), tspan, xx0, options);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: orbiterMeasurement.m 
function [rr_ECI, range] = orbiterMeasurement(lander, rr_ECI, vv_ECI, et0)
%-------------------------------------------------------------------------%
%                                                                         
% orbiterMeasurements: function for getting measurements of the orbirter
% position and range w.r.t. the lander 
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [rr_ECI, range] = orbiterMeasurement(lander, rr_ECI, vv_ECI, et0) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  lander       [1]   struct with lander data          [-]
%  rr_ECI       [3x1] position vector                  [km; km; km]
%  vv_ECI       [3x1] velocity vector                  [km/s; km/s; km/s]
%  et0          [1]   time                             [s]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  rr_ECI       [3x1] position                         [km]
%  range        [1]   range                            [km]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  15/12/2024: First version
%
%-------------------------------------------------------------------------%

% Define station name
topoFrame = [lander.name, '_TOPO'];               

% Transformation from ECI to topocentric frame
ROT_ECI2TOPO = cspice_sxform('J2000', topoFrame, et0);  

% Get alnder state in J2000
rv_station_eci = cspice_spkezr(lander.name, et0, 'J2000', 'NONE', 'MOON'); 

 % Compute lander-orbiter vector in ECI
rv_station_sat_eci = [rr_ECI; vv_ECI] - rv_station_eci;                     

% Convert state into topocentric frame
rv_station_sat_topo = ROT_ECI2TOPO * rv_station_sat_eci; % Convert state into topocentric frame

% Computation of the range
range = norm(rv_station_sat_topo(1:3));

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: simulatedMeasurements.m 
function [VIS_WIND] = simulatedMeasurements(VIS_WIND, lander, orbiter, tspan, mu)
%-------------------------------------------------------------------------%
%                                                                         
% simulatedMeasurements: function which propagate the state of the satellite with and by checking
% the visibility window compute the measuremnts associated to the state and
% then add noises 
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [VIS_WIND] = simulatedMeasurements(VIS_WIND, station, orbiter, tspan, mu) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  VIS_WIND     [-]   struct with visbility window data         [-]
%  lander       [1]   struct with lander data                   [-]
%  orbiter      [1]   struct with orbiter data                  [-]
%  tspan        [nx1] time vector                               [s]
%  mu           [1]   gravitational constant                    [km^3/s^2]
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
%  15/12/2024: First version
%
%-------------------------------------------------------------------------%

% Get the number of points in the time vector
npoints = length(tspan);

% Initializzation of temporary variables for state and range
range = NaN(npoints, 1);                 % Range
rr = NaN(3, npoints);                    % Position
vv = NaN(3, npoints);                    % Velocity

% Initializzation of variable for theVIS_WIND struct
VIS_WIND.measurements.range = NaN(1, npoints);                 % Position
VIS_WIND.measurements.rr = NaN(3, npoints);                    % Velocity

% Propagation
[xx, ~] = propagate2BP(orbiter.xx0, tspan, mu);

% Measurements computation
for i = 1 : npoints

    % Chech if in that instant of the time the satellite is in the visibility window
    if VIS_WIND.estimated.i_visibility(i)      
        
        rrECI_SAT = xx(i, 1:3)';
        vvECI_SAT = xx(i, 4:6)';
        vv(:,i) = vvECI_SAT; % Assign velocity to the struct
        
        % Compute the measurements
        [rr(:, i), range(i)] = orbiterMeasurement(lander, rrECI_SAT, vvECI_SAT, tspan(i));    
    end
end

% Save values without noises
VIS_WIND.real.range = range;
VIS_WIND.real.rr = rr;
VIS_WIND.real.vv = vv;

% Adding noises to the measurements
VIS_WIND.measurements.range = mvnrnd(range, orbiter.sigma_r^2);             % Range
VIS_WIND.measurements.rr = mvnrnd(rr, ones(1, npoints) * orbiter.sigma_r^2);% Position

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: UKF_ORBITER.m 
function [x_hat, P] = UKF_ORBITER(VIS_WIND, orbiter, tspan, mu, alpha, beta, k)
%-------------------------------------------------------------------------%
%                                                                         
% UKF_ORBITER: Unscented Kalman Filter for the orbiter state                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [x_hat, P] = UKF_ORBITER(VIS_WIND, orbiter, tspan, mu, alpha, beta, k) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  VIS_WIND     [-]   struct with visbility window data          [-]
%  orbiter      [1]   struct with orbiter data                   [-]
%  tspan        [nx1]  time vector                               [s] 
%  mu           [1]    gravitational constant                    [km^3/s^2]
%  alpha        [1]    parameter for computing weights           [-] 
%  beta         [1]    parameter for computing weights           [-] 
%  k            [1]    parameter for computing weights           [-]    
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  x_hat          [6xn]    mean state vector                     [km; km; km; km/s; km/s; km/s]
%  P              [6xn*6]  covariance matrix vector              [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  16/12/2024: First version
%
%-------------------------------------------------------------------------%

% In the notation used in the function the subscript k refers to t = tk
% while k1 stands for t = tk + 1

% Intializzation of state and covariance matrix at t = tk 
Pk_hat = orbiter.P0(1:6, 1:6);       % Covariance matrix
xk_hat = mvnrnd(orbiter.xx0, Pk_hat);      % Mean state

nMeas = 3;
N = length(xk_hat);                         % Get number of value in the state vector
nPoints = length(tspan);                    % Get number of points in the time vector
nPointsVW = nPoints;                        % Get number of points in the visibility window
c = alpha^2 * (N + k);                      % Compute c value used for weights
Rk = orbiter.noiseMatrix(1:nMeas, 1:nMeas); % Assign noise matrix of the ground station to a variable

% Initializzation of  output variable
x_hat = zeros(N, nPointsVW);    % Mean state
P = zeros(N, N * nPointsVW);    % Covariance matrix

% Mean weigths
W0_m = 1 - N/c;                     % Mean's weight for sigma point n.0
Wi_m = 1/(2 * c);                   % Mean's weight for sigma points

% Covariance weigths
W0_c = (2 - alpha^2 + beta) - N/c;  % Covariance's weight for sigma point n.0
Wi_c = 1/(2 * c);                   % Covariance's weight for sigma points

% Assign to the output the estimated state of t = tk
x_hat(:, 1) = xk_hat;   % Mean state
P(:, 1:N) = Pk_hat;     % Covariance
tk = 1;                 % Index used for UKF 

for j = 1 : nPoints - 1 % Cycle used for the tspan vector

    % Check if the j-th point and %j+1-th point of tspan are in the visibility window 
    if VIS_WIND.estimated.i_visibility(j) && VIS_WIND.estimated.i_visibility(j + 1)   

        xk_hat = x_hat(:, tk);                      % Assign to xk the value of the mean state of t = tk-1
        Pk_hat = P(:, (1+(tk-1)*N: N+(tk-1)*N));    % Assign to xk the value of the mean state of t = tk-1
                                      
        X_k = zeros(N, 2*N + 1);                    % Matrix of sigma points
        X_k(:, 1) = xk_hat; 
        X_k1 = zeros(N, 2*N + 1);                   % Matrix of propagated of sigma points
        Y_k1 = zeros(nMeas, 2*N + 1);               % Matrix of measurement associated to propagated sigma points

        % Creation of sigma points
        for i = 1 : 2*N

            if i <= N
                Pxi = sqrtm(c * Pk_hat);
                dxi = Pxi(:, i);
            else
                Pxi = -sqrtm(c * Pk_hat);
                dxi = Pxi(:, i-N);
            end
        X_k(:, i + 1) = xk_hat + dxi;  % Assemble the matrix of sigma points
        end
    
        % Propagation of sigma points and apply measurement function
        for i = 1 : 2*N+1

            % Propagate the i-th sigma point
            [xx, ~] = propagate2BP(X_k(:, i), [tspan(j) tspan(j + 1)], mu);   

            % Assemble the matrix of propagated sigma points
            X_k1(:, i) = xx(end, :)';                                             

            % Assemble the matrix measurements
            Y_k1(:, i) = X_k1(1:3, i);                                                                           
        end

        % Compute mean of state and measurement using UT
        x_k1_hat_minus = 0;
        y_k1_hat_minus = 0;
        for i = 1 : 2*N + 1
            if i == 1
                x_k1_hat_minus = x_k1_hat_minus + W0_m * X_k1(:, i); 
                y_k1_hat_minus = y_k1_hat_minus + W0_m * Y_k1(:, i); 
            else
                x_k1_hat_minus = x_k1_hat_minus + Wi_m * X_k1(:, i);  
                y_k1_hat_minus = y_k1_hat_minus + Wi_m * Y_k1(:, i);  
            end
        end

        % Compute Covariance, Pee,k Pxy,k
        P_k1_minus = zeros(N, N);
        P_ee_k1 = zeros(nMeas, nMeas);
        P_xy_k1 = zeros(N, nMeas);
        for i = 1 : 2*N + 1
            
            % Create the variable used for the different in the measurements
            diff_yk = Y_k1(:, i) - y_k1_hat_minus;

            if i == 1
                P_k1_minus = P_k1_minus + W0_c * (X_k1(:, i) - x_k1_hat_minus) .* (X_k1(:, i) - x_k1_hat_minus)';
                P_ee_k1 = P_ee_k1 + W0_c * diff_yk .* diff_yk';
                P_xy_k1 = P_xy_k1 + W0_c * (X_k1(:, i) - x_k1_hat_minus) .* diff_yk';
            else
                P_k1_minus = P_k1_minus + Wi_c * (X_k1(:, i) - x_k1_hat_minus) .* (X_k1(:, i) - x_k1_hat_minus)';
                P_ee_k1 = P_ee_k1 + Wi_c * diff_yk .* diff_yk';
                P_xy_k1 = P_xy_k1 + Wi_c * (X_k1(:, i) - x_k1_hat_minus) .* diff_yk';
            end
        end
        
        % Add the noise matrix to the P_ee_k1
        P_ee_k1 = P_ee_k1 + Rk; 
        
        % Update state
        Kk = P_xy_k1 / (P_ee_k1); % Compute Kalman gain
        
        % Get measuremnts at t = tk
        yk = VIS_WIND.measurements.rr(:, j+1);

        % Update mean
        x_hat(:, tk + 1) = x_k1_hat_minus + Kk * (yk - y_k1_hat_minus);

        % Update covariance
        P(:, (1+(tk)*N: N+(tk)*N)) = P_k1_minus - Kk * P_ee_k1 * Kk';
        tk = tk + 1;
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: UKF.m 
function [x_hat, P] = UKF_ORBITER_LANDER(VIS_WIND, orbiter, station, tspan, mu, alpha, beta, k)
%-------------------------------------------------------------------------%
%                                                                         
% UKF_ORBITER_LANDER: function Unscented Kalman Filter for orbiter state and lander coordinates                
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   UKF_ORBITER_LANDER(VIS_WIND, orbiter, station, tspan, mu, alpha, beta, k) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  VIS_WIND     [-]   struct with visbility window data          [-]
%  orbiter      [1]   struct with orbiter data                   [-]
%  lander       [1]   struct with lander data                    [-]
%  tspan        [nx1]  time vector                               [s] 
%  mu           [1]    gravitational constant                    [km^3/s^2]
%  alpha        [1]    parameter for computing weights           [-] 
%  beta         [1]    parameter for computing weights           [-] 
%  k            [1]    parameter for computing weights           [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  x_hat          [6xn]    mean state vector                     [km; km; km; km/s; km/s; km/s]
%  P              [6xn*6]  covariance matrix vector              [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  16/12/2024: First version
%
%-------------------------------------------------------------------------%

% In the notation used in the function the subscript k refers to t = tk
% while k1 stands for t = tk + 1

% Intializzation of state and covariance matrix at t = tk 
Pk_hat = orbiter.P0;       % Covariance matrix
xk_hat = mvnrnd([orbiter.xx0; station.lon; station.lat], Pk_hat);      % Mean state

nMeas = 4;
N = length(xk_hat);                                 % Get number of value in the state vector
nPoints = length(tspan);                            % Get number of points in the time vector
nPointsVW = nPoints;                                % Get number of points in the visibility window
c = alpha^2 * (N + k);                              % Compute c value used for weights
Rk = orbiter.noiseMatrix;                           % Assign noise matrix of the ground station to a variable

% Initializzation of  output variable
x_hat = zeros(N, nPointsVW);    % Mean state
P = zeros(N, N * nPointsVW);    % Covariance matrix

% Mean weigths
W0_m = 1 - N/c;                     % Mean's weight for sigma point n.0
Wi_m = 1/(2 * c);                   % Mean's weight for sigma points

% Covariance weigths
W0_c = (2 - alpha^2 + beta) - N/c;  % Covariance's weight for sigma point n.0
Wi_c = 1/(2 * c);                   % Covariance's weight for sigma points

% Assign to the output the estimated state of t = tk
x_hat(:, 1) = xk_hat;   % Mean state
P(:, 1:N) = Pk_hat;     % Covariance
tk = 1;                 % Index used for UKF 

% Extract Moon radius
R_Moon = cspice_bodvrd( 'MOON', 'RADII', 3);
RM = R_Moon(1);

% Computation of the Moon flatness
f_Moon = (RM - R_Moon(3)) / RM;

% Cycle used for the tspan vector
for j = 1 : nPoints - 1 

    % Check if the j-th point and %j+1-th point of tspan are in the visibility window 
    if VIS_WIND.estimated.i_visibility(j) && VIS_WIND.estimated.i_visibility(j + 1)   

        xk_hat = x_hat(:, tk);                      % Assign to xk the value of the mean state of t = tk-1
        Pk_hat = P(:, (1+(tk-1)*N: N+(tk-1)*N));    % Assign to xk the value of the mean state of t = tk-1
                                      
        X_k = zeros(N, 2*N + 1);                    % Matrix of sigma points
        X_k(:, 1) = xk_hat; 
        X_k1 = zeros(N, 2*N + 1);                   % Matrix of propagated of sigma points
        Y_k1 = zeros(nMeas, 2*N + 1);               % Matrix of measurement associated to propagated sigma points

        % Creation of sigma points
        for i = 1 : 2*N

            if i <= N
                Pxi = (sqrtm(c * Pk_hat));
                dxi = Pxi(:, i);
            else
                Pxi = -(sqrtm(c * Pk_hat));
                dxi = Pxi(:, i-N);
            end
        X_k(:, i + 1) = xk_hat + dxi;  % Assemble the matrix of sigma points
        end
    
        % Propagation of sigma points and apply measurement function
        for i = 1 : 2*N+1
            
            % Propagate the i-th sigma point
            [xx, ~] = propagate2BP(X_k(1:6, i), [tspan(j) tspan(j + 1)], mu);   

            % Assemble the matrix of propagated sigma points
            X_k1(1:6, i) = xx(end, :)';                                         
            
            % Propagation of the lander coordinates
            X_k1(7:8, i) = X_k(7:8, i);

            % Extract the rotation matrix from Moon Centered inertial frame
            % to Moon Centered fixed
            R_J20002IAUMOON = cspice_sxform('J2000', 'IAU_MOON', tspan(j + 1));

            % Conversion of lander coordinates to rectangular
            rr_MOONLANDER_IAU_MOON = cspice_georec(X_k1(7, i), X_k1(8, i), 0, RM, f_Moon);

            % Rotation of the orbiter state into IAU_MOON
            rv_ORBITER_IAU_MOON = R_J20002IAUMOON * X_k1(1:6, i);

            % Computation of the range
            rng = norm(rv_ORBITER_IAU_MOON(1:3) - rr_MOONLANDER_IAU_MOON);
            
            % Assemble the measurements matrix
            Y_k1(:, i) = [X_k1(1:3, i); rng]; 
        end

        % Compute mean of state and measurement using UT
        x_k1_hat_minus = zeros(N, 1);
        y_k1_hat_minus = zeros(4, 1);
        for i = 1 : 2*N + 1

            if i == 1
                x_k1_hat_minus = x_k1_hat_minus + W0_m * X_k1(:, i); 
                y_k1_hat_minus = y_k1_hat_minus + W0_m * Y_k1(:, i); 
            else
                x_k1_hat_minus = x_k1_hat_minus + Wi_m * X_k1(:, i);  
                y_k1_hat_minus = y_k1_hat_minus + Wi_m * Y_k1(:, i);  
            end
        end

        % Compute Covariance, Pee,k Pxy,k
        P_k1_minus = zeros(N, N);
        P_yy_k1 = zeros(nMeas, nMeas);
        P_xy_k1 = zeros(N, nMeas);
        for i = 1 : 2*N + 1

            % Variable used for the difference in the state
            diff_xk1 = X_k1(:, i)  - x_k1_hat_minus;

            % Variable used for the difference in the
            % measurements
            diff_yk1 = Y_k1(:, i) - y_k1_hat_minus;

            if i == 1
                P_k1_minus = P_k1_minus + W0_c * diff_xk1 .* diff_xk1';
                P_yy_k1 = P_yy_k1 + W0_c * diff_yk1 .* diff_yk1';
                P_xy_k1 = P_xy_k1 + W0_c * diff_xk1 .* diff_yk1';
            else
                P_k1_minus = P_k1_minus + Wi_c * diff_xk1 .* diff_xk1';
                P_yy_k1 = P_yy_k1 + Wi_c * diff_yk1 .* diff_yk1';
                P_xy_k1 = P_xy_k1 + Wi_c * diff_xk1 .* diff_yk1';
            end
        end

        % Add the noise matrix to the P_ee_k1
        P_yy_k1 = P_yy_k1 + Rk; 

        % Update state
        Kk = P_xy_k1 / (P_yy_k1); % Compute Kalman gain
        
        % Get measuremnts at t = tk
        yk = [VIS_WIND.measurements.rr(:, j+1); VIS_WIND.measurements.range(j+1)];

        % Update mean
        x_hat(:, tk + 1) = x_k1_hat_minus + Kk * (yk - y_k1_hat_minus);

        % Update covariance
        P(:, (1+(tk)*N: N+(tk)*N)) = P_k1_minus - Kk * P_yy_k1 * Kk';
        tk = tk + 1;
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: Errorand3SigmaComputation.m 
function [err, sigma3] = Errorand3SigmaComputationLATLON(x_hat, P, tspan)
%-------------------------------------------------------------------------%
%                                                                         
% Errorand3SigmaComputationLATLON: function which computes errors between real 
% lander coordinates and estimated one, and computes also the 3sigma associate to the
% angles.
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [err, sigma3] = Errorand3SigmaComputationLATLON(x_hat, P, tspan)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x_hat          [2xn]    mean state vector                     [rad; rad]
%  P              [2xn*2]  covariance matrix vector              [-]
%  tspan          [nx1]  time vector                             [s] 
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  err             [nx1] norm of angles errors                 [rad]
%  sigma3          [nx1] 3 sigma of angles                     [rad]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  17/12/2024: First version
%
%-------------------------------------------------------------------------%

% Initializzation of output variables
err = zeros(size(x_hat, 2), 1);      % Position errors
sigma3 = zeros(size(x_hat, 2), 1);   % 3 sigma position

% Extract Moon radius
R_Moon = cspice_bodvrd( 'MOON', 'RADII', 3);
RM = R_Moon(1);

% Computation of Moon flatness
f_Moon = (RM - R_Moon(3)) / RM;

% Computation
for i = 1 : length(tspan)

    % Extract real coordinates
    rr_LANDER = cspice_spkpos('MOONLANDER', tspan(i), 'IAU_MOON', 'NONE', 'MOON');

    % conversion in geotedic
    [realLON, realLAT, ~] = cspice_recgeo(rr_LANDER, RM, f_Moon);

    % Computation of the error
    err(i) = norm([realLON; realLAT] - x_hat(:, i));

    % extraction of the P
    PLONLAT = P(1:2, (1+(i-1)*2: 2+(i-1)*2));

    % Computation of 3 sigma
    sigma3(i) = 3 * sqrt(trace(PLONLAT));
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: Errorand3SigmaComputation.m 
function [errPos, sigma3Pos, errVel, sigma3Vel] = Errorand3SigmaComputation(x_hat, P, realPos, realVel)
%-------------------------------------------------------------------------%
%                                                                         
% Errorand3SigmaComputation: function which computes errors between real 
% state and estimated one, and computes also the 3sigma associate to the
% position and velocity.
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [errPos, sigma3Pos, errVel, sigma3Vel] = Errorand3SigmaComputation(x_hat, P, realPos, realVel)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x_hat       [6xn]    mean state vector                    [km; km;km;km/s;km/s;km/s]
%  P           [6xn*6]  covariance matrix vector             [-]
%  realPos     [3xn]    matrix with real position            [km; km; km]
%  realVel     [3xn]    matrix with real velocity            [km/s; km/s; km/s]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  errPos      [nx1] norm of position errors                 [km]
%  sigma3Pos   [nx1] 3 sigma of position                     [km]
%  errPos      [nx1] norm of velocity errors                 [km/s]
%  sigma3VEl   [nx1] 3 sigma of velocity                     [km/s]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  17/12/2024: First version
%
%-------------------------------------------------------------------------%

% Initializzation of output variables
errPos = zeros(size(x_hat, 2), 1);      % Position errors
sigma3Pos = zeros(size(x_hat, 2), 1);   % 3 sigma position
errVel = zeros(size(x_hat, 2), 1);      % Velocity errors
sigma3Vel = zeros(size(x_hat, 2), 1);   % 3 sigma velocity

% Computation
for i = 1 : size(x_hat, 2)    

    % Computation of the error relative to the position
    errPos(i) = norm(realPos(:, i) - x_hat(1:3, i));

    % Computation of the error relative to the velocity
    errVel(i) = norm(realVel(:, i) - x_hat(4:6, i));

    % Extract P position
    Ppos = P(1:3, (1+(i-1)*6: 3+(i-1)*6));

    % Extract P velocity
    Pvel = P(4:6, (4+(i-1)*6: 6+(i-1)*6));

    % Computation of 3 sigma position
    sigma3Pos(i) = 3 * sqrt(trace(Ppos));

    % Computation of 3 sigma velocity
    sigma3Vel(i) = 3 * sqrt(trace(Pvel));
end

end
%-------------------------------------------------------------------------%
