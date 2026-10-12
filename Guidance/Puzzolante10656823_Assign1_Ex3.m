% Guidance and Navigation (2024/2025)
% Assignment # 1
% Exercise 3
% Author : Michele Puzzolante
%% REQUEST 1
% Initializzation of the workspace
clear;
close all;
clc;

SAVEPLOTS = false; % Flag to use when plots has to be saved as png file

% Set default options for plots
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');

fprintf('\nREQUEST 1\n')

% Struct which contains the the problem's dimensional data 
DIM_param.hi = 800;                 % [km]
DIM_param.hf = 1000;                % [km]
DIM_param.delta_i = 0.75;           % [deg]
DIM_param.RE = 6378.1366;           % [km]
DIM_param.mu = 398600.435;          % [km^3/s^2]
DIM_param.r0 = 750 + DIM_param.RE;  % [km]
DIM_param.k1 = 1e-5;                % [-]
DIM_param.k2 = 1e-4;                % [-]
DIM_param.m0 = 1000;                % [kg]
DIM_param.T_MAX = 3.000;            % [N]
DIM_param.Isp = 3120;               % [s]
DIM_param.g0 = 9.81 * 10^-3;

% Unit values
DIM_param.DU = 7178.1366;                           % [km]
DIM_param.MU = DIM_param.m0;                        % [kg]
DIM_param.TU = sqrt(DIM_param.DU^3 / DIM_param.mu); % [s]
DIM_param.VU = DIM_param.DU/DIM_param.TU;           % [km/s]
DIM_param.AU = DIM_param.DU/DIM_param.TU^2;         % [km/s^2]
DIM_param.FU = DIM_param.MU * DIM_param.AU * 10^3;  % [N]

% Struct which contains the the problem's adimensional data 
ADIM_param.hi = DIM_param.hi / DIM_param.DU;        % [-]
ADIM_param.hf = DIM_param.hf / DIM_param.DU;        % [-]
ADIM_param.delta_i = 0.75 * cspice_rpd;             % [rad]
ADIM_param.RE = DIM_param.RE / DIM_param.DU;        % [-]
ADIM_param.mu = 1;                                  % [-]
ADIM_param.r0 = DIM_param.r0 / DIM_param.DU;        % [-]
ADIM_param.k1 = DIM_param.k1;                       % [-]
ADIM_param.k2 = DIM_param.k2;                       % [-]
ADIM_param.m0 = DIM_param.m0 / DIM_param.MU;        % [-]
ADIM_param.T_MAX = DIM_param.T_MAX / DIM_param.FU;  % [-]
ADIM_param.Isp = DIM_param.Isp / DIM_param.TU;      % [-]
ADIM_param.g0 = DIM_param.g0 / DIM_param.AU;
% Unit values
ADIM_param.DU = 1;       % [-]
ADIM_param.MU = 1;       % [-]
ADIM_param.TU = 1;       % [-]
ADIM_param.VU = 1;       % [-]
ADIM_param.AU = 1;       % [-]
ADIM_param.FU = 1;       % [-]

% Computation of the initial state
DIM_param.xx0 = [DIM_param.hi + DIM_param.RE; 0; 0;...               % Dimesnional initial position  [km]
            0; sqrt(DIM_param.mu/(DIM_param.hi + DIM_param.RE)); 0]; % Dimensional initial velocity  [km/s]
ADIM_param.xx0 = [DIM_param.xx0(1:3)/DIM_param.DU;...                % Adimesnional initial position [km]
                 DIM_param.xx0(4:6)/DIM_param.VU];                   % Adimensional initial velocity [km/s]

% Computation of the final state
DIM_param.xxf = [DIM_param.hf + DIM_param.RE; 0; 0;...                  % Dimesnional final position  [km]
             0; sqrt(DIM_param.mu/(DIM_param.hf + DIM_param.RE)) * cos(ADIM_param.delta_i); sqrt(DIM_param.mu/(DIM_param.hf + DIM_param.RE)) * sin(ADIM_param.delta_i)]; % Dimensional final velocity  [km/s]
ADIM_param.xxf = [DIM_param.xxf(1:3)/DIM_param.DU;...                   % Adimesnional final position [km]
                  DIM_param.xxf(4:6)/DIM_param.VU];                     % Adimensional final velocity [km/s]

% Display the initial state
fprintf('INITIAL STATE\n')
fprintf('x0 =\t%.6f\t[km]\n', DIM_param.xx0(1)); fprintf('y0 =\t%.6f\t[km]\n', DIM_param.xx0(2)); fprintf('z0 =\t%.6f\t[km]\n', DIM_param.xx0(3));
fprintf('vx0 =\t%.8f\t[km/s]\n', DIM_param.xx0(4)); fprintf('vy0 =\t%.8f\t[km/s]\n', DIM_param.xx0(5)); fprintf('vz0 =\t%.8f\t[km/s]\n', DIM_param.xx0(6));

% Display the final state
fprintf('FINAL STATE\n')
fprintf('xf =\t%.6f\t[km]\n', DIM_param.xxf(1)); fprintf('yf =\t%.6f\t[km]\n', DIM_param.xxf(2)); fprintf('zf =\t%.6f\t[km]\n', DIM_param.xxf(3));
fprintf('vxf =\t%.8f\t[km/s]\n', DIM_param.xxf(4)); fprintf('vyf =\t%.8f\t[km/s]\n', DIM_param.xxf(5)); fprintf('vzf =\t%.8f\t[km/s]\n',DIM_param. xxf(6));

% Function for the spatial density od debries
q = @(r, param) param.k1./ (param.k2 + ((r - param.r0)/param.DU).^2 ); 

% Plot q
rho_vec = DIM_param.RE + linspace(DIM_param.hi - 100, DIM_param.hf + 100, 1000);
figure(1)
plot(rho_vec, q(rho_vec, DIM_param), 'LineWidth', 2, 'Color',  [0.6350 0.0780 0.1840], 'HandleVisibility', 'off')
hold on
plot(norm(DIM_param.xx0(1:3)), q(norm(DIM_param.xx0(1:3)), DIM_param), 'LineStyle', 'none', 'LineWidth', 2, 'Marker','x', 'MarkerSize', 10, 'Color',[0.3010 0.7450 0.9330])
plot(norm(DIM_param.xxf(1:3)), q(norm(DIM_param.xxf(1:3)), DIM_param), 'LineStyle', 'none', 'LineWidth', 2, 'Marker','x', 'MarkerSize', 10, 'Color',[0.4660 0.6740 0.1880])
xlabel('$\rho [km]$',  'FontSize', 15);    
ylabel('$q(\rho)$',  'FontSize', 15);
legend('$\rho_i$', '$\rho_f$', 'FontSize', 12)
xlim([rho_vec(1) rho_vec(end)])

% Saving plots
if SAVEPLOTS

    % save 1st plot
    figure(1)
    exportgraphics(gcf,'q_rho.png','Resolution',600);
end

%% REQUEST 2

fprintf('\nREQUEST 2\n')

fprintf('UNIT VALUES\n')
fprintf('DISTANCE UNIT\t\tDU = %f\t[km]\n', DIM_param.DU)
fprintf('TIME UNIT\t\t\tTU = %f\t\t[s]\n', DIM_param.TU)
fprintf('MASS UNIT\t\t\tMU = %f\t[kg]\n', DIM_param.MU)
fprintf('VELOCITY UNIT\t\tVU = %f\t\t[km/s]\n', DIM_param.VU)
fprintf('ACCELERATION UNIT\tAU = %f\t\t[km/s^2]\n', DIM_param.AU)
fprintf('FORCE UNIT\t\t\tFU = %f\t\t[N]\n', DIM_param.FU)

fprintf('\nADIMENSIONALIZED PARAMETERS\n')
fprintf('ALTITUDE OF DEP. ORB.\t\t\thi\t\t= %f\n', ADIM_param.hi);
fprintf('ALTITUDE OF ARR. ORB.\t\t\thf\t\t= %f\n', ADIM_param.hf);
fprintf('INCLINATION CHANGE\t\t\t\tdelta_i\t= %f\n', ADIM_param.delta_i);
fprintf('EARTH RADIUS\t\t\t\t\tRE\t\t= %f\n', ADIM_param.RE);
fprintf('EARTH GRAVITATIONAL PARAM.\t\tmu\t\t= %f\n', ADIM_param.mu);
fprintf('REF. RADIUS FOR DEBRIES FLUX\trho_0\t= %f\n', ADIM_param.r0);
fprintf('DEBRIES. SPAT.DEN. CONST.1\t\tk1\t\t= %f\n', ADIM_param.k1);
fprintf('DEBRIES. SPAT.DEN. CONST.2\t\tk2\t\t= %f\n', ADIM_param.k2);
fprintf('INITIAL MASS\t\t\t\t\tm_0\t\t= %f\n', ADIM_param.m0);
fprintf('MAX THRUST\t\t\t\t\t\tT_MAX\t= %f\n', ADIM_param.T_MAX);
fprintf('SPECIFIC IMPULSE\t\t\t\tIsp\t\t= %f\n', ADIM_param.Isp);

%% REQUEST 3

fprintf('\nREQUEST 3\n')

%{
PMP stated
%}

%% REQUEST 4

fprintf('\nREQUEST 4\n')

% Initliazzation of variables
lambda_ub = 250;    % Upper bound for costate values
lambda_lb = -250;   % Lower bound for costate values
tf0 = 10 * 2 * pi;  % Initial guess for final time
exitflag = -1;      % exitflag set = 0 to initialize the while
USE_OPTIMAL = true; % flag used to use the optimal initial condition to avoid the random generation of the initial guess for the costate

% Check if the optimal initial condition have to be used
if USE_OPTIMAL

    % Initial condition
    y0 = [-214.9812;...
            10.3659;...
             0.8856;...
           -10.3929;...
          -214.6105;...
          -112.9454;...
             2.5964;...
            64.4801];
end

% Set options for fsolve
optionsFSOLVE = optimoptions('fsolve', 'Algorithm','trust-region-dogleg', 'Display','iter-detailed',...
                             'MaxIterations',500, 'MaxFunctionEvaluations',5000, 'FunctionTolerance', 1e-28,...
                             'SubproblemAlgorithm','factorization', 'SpecifyObjectiveGradient',false,...
                             'StepTolerance', 1e-14, 'UseParallel',true, 'FiniteDifferenceType','central', 'OptimalityTolerance', 1e-8);

% while loop for searching the initial guess
while exitflag <= 0

    % if for optimal initial condition
    if ~USE_OPTIMAL

        % Generation of the initial guess for the costate
        lambda = lambda_lb + (lambda_ub - lambda_lb).* rand(7,1);

        % Assemble the initla guess for fsolve formed by costate and final time
        y0 = [lambda; tf0];
    end

    % Call fsolve
    [yOPT, residualOPT, exitflag, output] = fsolve(@(y) Residuals(y, ADIM_param), y0, optionsFSOLVE);
end

% Extract the oprimal final time and optimal costate initial condition from the solution
tf_OPT = yOPT(end);     % Final time
lambda_OPT = yOPT(1:7); % Costate

% Propagation of the optimal solution
[~, ~, xx_OPT, lambda_vec_OPT, tt_OPT] = TPBVP_propagation([ADIM_param.xx0; ADIM_param.m0; lambda_OPT], [0 tf_OPT], ADIM_param);

% Computation of the error on the final state
error_Position = norm(residualOPT(1:3)) * DIM_param.DU;         % [km]
error_Velocity = norm(residualOPT(4:6)) * DIM_param.VU * 10^3;  % [m/s]

% Computation of the primer vector in NTW reference frame
[~, alpha_NTW] = primerVector(xx_OPT, lambda_vec_OPT, tt_OPT);

% Computation of the Hamiltonian
H = Hamiltonian(xx_OPT, lambda_vec_OPT, tt_OPT, ADIM_param);

% Display results
fprintf('\nOPTIMAL ORBIT RAISING TRANSFER SOLUTION (Tmax = 3.000 N)\n')
fprintf('lambda0_rx\t=\t %.4f\t[-]\n', lambda_OPT(1))
fprintf('lambda0_ry\t=\t %.4f\t\t[-]\n', lambda_OPT(2))
fprintf('lambda0_rz\t=\t %.4f\t\t[-]\n', lambda_OPT(3))
fprintf('lambda0_vx\t=\t %.4f\t\t[-]\n', lambda_OPT(4))
fprintf('lambda0_vy\t=\t %.4f\t[-]\n', lambda_OPT(5))
fprintf('lambda0_vz\t=\t %.4f\t\t[-]\n', lambda_OPT(6))
fprintf('lambda0_m\t=\t %.4f\t\t[-]\n', lambda_OPT(7))
fprintf('tf\t\t\t=\t %.4f\t[min]\n', tf_OPT * DIM_param.TU /60)
fprintf('mf\t\t\t=\t %.4f\t\t[kg]\n', xx_OPT(7, end) * DIM_param.MU)

fprintf('\nFINAL STATE ERROR WITH RESPECT TO TARGET POSITION AND VELOCITY (Tmax = 3.000 N)\n')
fprintf('||r(tf) - rf|| =\t%.4e\t[km]\n', error_Position)
fprintf('||v(tf) - vf|| =\t%.4e\t[m/s]\n', error_Velocity)

% Plot Transfer adimensional
figure(2)
plot3(0, 0, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize', 50, 'Color', [0 0.4470 0.7410])
hold on
plot3(ADIM_param.xx0(1), ADIM_param.xx0(2), ADIM_param.xx0(3), 'LineStyle','none', 'MarkerSize', 25, 'Marker','.', 'Color', [0.4940 0.1840 0.5560]);
plot3(ADIM_param.xxf(1), ADIM_param.xxf(2), ADIM_param.xxf(3), 'LineStyle','none', 'MarkerSize', 25, 'Marker','.', 'Color', [0.6350 0.0780 0.1840]);
plot3(xx_OPT(1,:), xx_OPT(2,:), xx_OPT(3,:), 'LineWidth', 1.2, 'Color', [0.9290 0.6940 0.1250])
hold off
xlabel('x [-]', 'FontSize', 15); 
ylabel('y [-]', 'FontSize', 15); 
zlabel('z [-]',  'FontSize', 15);
legend('Earth''s position', 'Departure', 'Arrival', 'Transfer trajectory', 'FontSize', 12, 'Position', [0.66,0.71,0.33,0.2])

% Plot the primer's vector components
figure(3)
plot(tt_OPT, alpha_NTW(1, :), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840])
hold on
plot(tt_OPT, alpha_NTW(2, :), 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250])
plot(tt_OPT, alpha_NTW(3, :), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
hold off
xlabel('t [-]',  'FontSize', 15);   
ylabel('$\alpha_N$, $\alpha_T$, $\alpha_W$',  'FontSize', 15);
legend('$\alpha_N$', '$\alpha_T$', '$\alpha_W$', 'FontSize', 11, 'Position', [0.1257,0.0976,0.13,0.14])
xlim([0 +Inf])

% Plot the hamiltonian
figure(4)
plot(tt_OPT, H, 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840])
xlabel('t [-]',  'FontSize', 15);   
ylabel('H [-]',  'FontSize', 15);
xlim([0 +Inf])

% Saving plots
if SAVEPLOTS

    % save 2nd plot
    figure(2)
    exportgraphics(gcf,'transfer_T_3N.png','Resolution',600);

    % save 3rd plot
    figure(3)
    exportgraphics(gcf,'primer_vector.png','Resolution',600);

    % save 4th plot
    figure(4)
    exportgraphics(gcf,'Hamiltonian.png','Resolution',600);
end

%% REQUEST 5

fprintf('\nREQUEST 5\n')

% Initializzation of variables
N_points = 14;                                              % Number of points in thrust vector
T_MAX_new = 2.860;                                          % New value of thrust [N]
T_MAX_vec = linspace(DIM_param.T_MAX, T_MAX_new, N_points); % Thrust vector

% Initializzation of the new initial condition with the value found for
% T_MAX = 3 N
y0_new_T_MAX = yOPT;

% For used to exploit numerical continuation used to find the initial guess
% for T_MAX = 2.86 N
for j = 1 : N_points

    % Initializzation of the flag for while
    exitflag_new_T_MAX = 0;

    % Assign the new thrust value to the structs
    DIM_param.T_MAX = T_MAX_vec(j);                     % Dimensional value
    ADIM_param.T_MAX = DIM_param.T_MAX / DIM_param.FU;  % Adimensional value

    % While for finding the initial guess
    while exitflag_new_T_MAX <= 0

        % Call fsolve
        [yOPT_new_T_MAX, residualOPT_new_T_MAX, exitflag_new_T_MAX, output_new_T_MAX] = fsolve(@(y) Residuals(y, ADIM_param), y0_new_T_MAX, optionsFSOLVE);
    end

    % Initializzation of the new initial condition with the value found for
    % the previous value of thrust
    y0_new_T_MAX = yOPT_new_T_MAX;
end

% Extract the oprimal final time and optimal costate initial condition from the solution
tf_OPT_new_T_MAX = yOPT_new_T_MAX(end);     % Final time
lambda_OPT_new_T_MAX = yOPT_new_T_MAX(1:7); % Costate

% Propagation of the optimal solution
[~, ~, xx_OPT_new_T_MAX, lambda_vec_OPT_new_T_MAX, tt_OPT_new_T_MAX] = TPBVP_propagation([ADIM_param.xx0; ADIM_param.m0; lambda_OPT_new_T_MAX], [0 tf_OPT_new_T_MAX], ADIM_param);

% Computation of the error on the final state
error_Position_new_T_MAX = norm(residualOPT_new_T_MAX(1:3)) * DIM_param.DU;         % [km]
error_Velocity_new_T_MAX = norm(residualOPT_new_T_MAX(4:6)) * DIM_param.VU * 10^3;  % [m/s]

% Display results
fprintf('\nOPTIMAL ORBIT RAISING TRANSFER SOLUTION (Tmax = 2.860 N)\n')
fprintf('lambda0_rx\t=\t %.4f\t[-]\n', lambda_OPT_new_T_MAX(1))
fprintf('lambda0_ry\t=\t %.4f\t\t[-]\n', lambda_OPT_new_T_MAX(2))
fprintf('lambda0_rz\t=\t %.4f\t\t[-]\n', lambda_OPT_new_T_MAX(3))
fprintf('lambda0_vx\t=\t %.4f\t\t[-]\n', lambda_OPT_new_T_MAX(4))
fprintf('lambda0_vy\t=\t %.4f\t[-]\n', lambda_OPT_new_T_MAX(5))
fprintf('lambda0_vz\t=\t %.4f\t\t[-]\n', lambda_OPT_new_T_MAX(6))
fprintf('lambda0_m\t=\t %.4f\t\t[-]\n', lambda_OPT_new_T_MAX(7))
fprintf('tf\t\t\t=\t %.4f\t[min]\n', tf_OPT_new_T_MAX * DIM_param.TU /60)
fprintf('mf\t\t\t=\t %.4f\t\t[kg]\n', xx_OPT_new_T_MAX(7, end) * DIM_param.MU)

fprintf('\nFINAL STATE ERROR WITH RESPECT TO TARGET POSITION AND VELOCITY (Tmax = 2.860 N)\n')
fprintf('||r(tf) - rf|| =\t%.4e\t[km]\n', error_Position_new_T_MAX)
fprintf('||v(tf) - vf|| =\t%.4e\t[m/s]\n', error_Velocity_new_T_MAX)

% Plot the transfer in adimensional space
figure(5)
plot3(0, 0, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize', 50, 'Color', [0 0.4470 0.7410])
hold on
plot3(ADIM_param.xx0(1), ADIM_param.xx0(2), ADIM_param.xx0(3), 'LineStyle','none', 'MarkerSize', 25, 'Marker','.', 'Color', [0.4940 0.1840 0.5560]);
plot3(ADIM_param.xxf(1), ADIM_param.xxf(2), ADIM_param.xxf(3), 'LineStyle','none', 'MarkerSize', 25, 'Marker','.', 'Color', [0.6350 0.0780 0.1840]);
plot3(xx_OPT_new_T_MAX(1,:), xx_OPT_new_T_MAX(2,:), xx_OPT_new_T_MAX(3,:), 'LineWidth', 1.2, 'Color', [0.4660 0.6740 0.1880])
hold off
xlabel('x [-]',  'FontSize', 15); 
ylabel('y [-]',  'FontSize', 15); 
zlabel('z [-]',  'FontSize', 15);
legend('Earth''s position', 'Departure', 'Arrival', 'Transfer trajectory', 'FontSize', 12, 'Position', [0.66,0.71,0.33,0.2])

% Saving plots
if SAVEPLOTS

    % save 5th plot
    figure(5)
    exportgraphics(gcf,'transfer_T_2_86N.png','Resolution',600);
end

%% FUNCTIONS

%-------------------------------------------------------------------------%
% function: Residuals.m 
function [R, GR] = Residuals(y, param)
%-------------------------------------------------------------------------%
%                                                                         
% Residuals: Residuals of the problems                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [R, GR] = Residuals(y, param)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y            [8x1] vector of variables                        [-]
%  param        [-]   struct of constant values                  [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  R            [8x1] Residuals                                  [-]
%  GR           [8x8] Gradient's of residuals                    [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  5/11/2024: First version
%
%-------------------------------------------------------------------------%

% Exctract the costate and final time
lambda_rr = y(1:3); % position costate
lambda_vv = y(4:6); % velocity costate
lambda_m  = y(7);   % mass costate
tf = y(8);          % Finale time

% Extract the parameters from the struct
xx0 = param.xx0;       % initial state
m0  = param.m0;        % Initial mass 
xxf_TRUE = param.xxf;  % Final state

% Assemble the initial guess
y0 = [xx0; m0; lambda_rr; lambda_vv; lambda_m];

% Propagation of the dynamics
[xxf, lambda_f, ~, ~, ~] = TPBVP_propagation(y0, [0 tf], param);

% Extract lambda_m at tf from the propagation
lambda_mf = lambda_f(7);

% Computation of the Hamiltonian at final time
Hf = Hamiltonian(xxf, lambda_f, tf, param);

% Residuals' computation
R = [xxf(1:6) - xxf_TRUE;...
    lambda_mf;...
    Hf];

if nargout > 1

    % Gradient's residuals
    GR = [];
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: TPBVP_propagation.m 
function [xxf, lambda_f, xx, lambda, tt] = TPBVP_propagation(y0, tspan, param)
%-------------------------------------------------------------------------%
%                                                                         
% TPBVP_propagation: function for propagate the TPBVP                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xxf, lambda_f, xx, lambda, tt] = TPBVP_propagation(y0, tspan, param)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y0           [14x1] initial state and costate                 [state; costate]
%  tspan        [1x2] initial and fina time                      [-]
%  param        [-]   struct of constant values                  [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  xxf          [7x1] State at tf                                [-]
%  lambdaf      [7x1] Costate at tf                              [-]
%  xx           [7xn] State                                      [-]
%  lambda       [7xn] Costate                                    [-]
%  tt           [nx1] Time vector                                [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  5/11/2024: First version
%
%-------------------------------------------------------------------------%

% Set tolerance
Reltol = 1e-13; % Relative
Abstol = 1e-14; % Absolute

% Set options for ode
optionsODE = odeset('reltol', Reltol, 'abstol', Abstol * ones(14, 1));

% Perform integration
[tt, yy] = ode113(@(t, y) TPBVP_dynamics(t, y, param), tspan, y0, optionsODE);

% Extract state and mass
xx  = yy(:, 1:7)';
xxf = yy(end, 1:7)';

% Extract costate
lambda   = yy(:,8:14)';
lambda_f = yy(end,8:14)';

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: TPBVP_dynamics.m 
function [dydt] = TPBVP_dynamics(~, y, param)
%-------------------------------------------------------------------------%
%                                                                         
% TPBVP_dynamics: Ode for the two points boundary value problem                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [dydt] = TPBVP_dynamics(~, y, param)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  ~            [-]    time                                       [-]
%  y            [17x1] state and costate                          [-]
%  param        [-]    struct of constant values                  [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  dydt         [14x1] Ode                                        [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  5/11/2024: First version
%
%-------------------------------------------------------------------------%

% Extract position, velocity, mass, lambda_r, lambda_v and lambda_m from
% the variable vector
rr = y(1:3);                 % Position
r = norm(rr);                % Position's norm
vv = y(4:6);                 % Velocity
m = y(7);                    % Mass
lambda_rr = y(8:10);         % Lambda r
lambda_vv = y(11:13);        % Lambda v
lambda_v  = norm(lambda_vv); % Lambda v's norm

% Extract parameter from the struct
mu    = param.mu;       % Earth's gravitational parameter
T_MAX = param.T_MAX;    % Thrust
Isp   = param.Isp;      % Specific impulse
g0    = param.g0;       % Gravitational acceleration
r0    = param.r0;       % rho0
k1    = param.k1;       % k1
k2    = param.k2;       % k2
DU    = param.DU;       % Distance unit

% Computation of PRIMER vector
alpha = -lambda_vv/lambda_v;

% Call St switching function
[u, ~] = switchingFunction(y(1:14), param);

% Initializzation of the ode
dydt = zeros(14, 1);

% S/C dynamics
dydt(1:3) = vv;                                  % Velocity
dydt(4:6) = -mu/r^3 * rr + u * T_MAX/m * alpha;  % Acceleration
dydt(7)   = -u * T_MAX/ Isp / g0;                % Variation of mass

% Compute derivative of the debries spatial density w.r.t position
dq_drho = -(2 * k1 * DU^2 * (r - r0)) / ((k2 * DU^2 + (r - r0)^2)^2 * r);

% COSTATE dynamics
dydt(8 : 10)  = -dq_drho * rr - 3 * mu/r^5 * dot(rr, lambda_vv) * rr + mu/r^3 * lambda_vv;  % Lamda_r dot
dydt(11 : 13) = -lambda_rr;                                                                 % Lamda_v dot
dydt(14)      = -u * lambda_v * T_MAX / m^2;                                                % Lamda_m dot

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: switchingFunction.m 
function [u, St] = switchingFunction(y, param)
%-------------------------------------------------------------------------%
%                                                                         
% switchingFunction: Function used for compute the control action and the 
%                    switching function                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [u, St] = switchingFunction(y, param)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y            [14x1] vector of variables                       [-]
%  param        [-]    struct of constant values                 [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  u            [1]    Control action                            [-]
%  St           [1]    Switching function                        [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  5/11/2024: First version
%
%-------------------------------------------------------------------------%

% Exctract costate from the variables
m = y(7);                    % Mass
lambda_vv = y(11:13);        % Lambda v
lambda_v  = norm(lambda_vv); % Lambda v's norm
lambda_m  = y(14);           % Lambda m

% Extract parameters from the struct
Isp = param.Isp;    % Specific impluse
g0  = param.g0;     % Gravitational acceleration

% Compute St 
St = - lambda_v / m * Isp * g0 - lambda_m;

% Thrust throttle policy
if St > 0
    u = 0;
else
    u = 1;
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: Hamiltonian.m 
function [H] = Hamiltonian(xx, lambda, tt, param)
%-------------------------------------------------------------------------%
%                                                                         
% Hamiltonian: function which computes the hamiltonian                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [H] = Hamiltonian(xx, lambda, tt, param)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx           [7xn] State                                      [-]
%  lambda       [7xn] Costate                                    [-]
%  tt           [nx1] Time vector                                [-]
%  param        [-]   struct of constant values                  [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  H            [nx1] Hamiltonian                                [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  5/11/2024: First version
%
%-------------------------------------------------------------------------%

% Get the number of points in tt
nPoints  = length(tt);

% Extract the parameters from the struct
r0    = param.r0;       % rho0
k1    = param.k1;       % k1
k2    = param.k2;       % k2
DU    = param.DU;       % Distance unit

% Initializzation of the Hamiltonian vectors
H  = zeros(nPoints, 1);

% Computation of the Hamiltonian
for j = 1 : nPoints

    % Extract the position, and computation of the norm
    r_j = norm(xx(1:3, j));

    % Computation of the value of q(rho)
    q_j =  (k1 * DU^2)/(k2*DU^2 + (r_j - r0)^2);

    % Computation of the f(xx)
    [dydt] = TPBVP_dynamics(tt(j), [xx(:, j); lambda(:, j)], param);

    % Computation of the i-th value of H
    H(j) = q_j + dot(lambda(:, j), dydt(1:7));
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: primerVector.m 
function [alpha, alpha_NTW] = primerVector(xx, lambda, tt)
%-------------------------------------------------------------------------%
%                                                                         
% primerVector: Function which computes the value of the primer vector in
%               NTW reference frame
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [alpha, alpha_NTW] = primerVector(xx, lambda, tt)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx           [7xn] State                                      [-]
%  lambda       [7xn] Costate                                    [-]
%  tt           [nx1] Time vector                                [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  alpha        [3xn] Primer vector                              [-]
%  alpha_NTW    [3xn] Primer vector  in NTW ref.frame            [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  5/11/2024: First version
%
%-------------------------------------------------------------------------%

% Get the length of tt
nPoints = length(tt);

% Initializzation of the output variables
alpha = zeros(3, nPoints);
alpha_NTW = zeros(3, nPoints);

% Primer's vector  computation
for j = 1 : nPoints

    % Extract lambda_vv from the costate
    lambda_vv = lambda(4:6, j);

    % Extract the position from the state
    rr = xx(1:3, j);

    % Extract the velocity from the state
    vv = xx(4:6, j);

    % Computation of the primer vector
    alpha(:, j) = -lambda_vv / norm(lambda_vv);

    % Computation of the axes of the NTW ref. frame
    T = vv / norm(vv);                          % T - axis
    W = -rr / norm(rr);                         % W - axis
    N = cross(T, W) / norm(cross(T, W));        % N - axis

    % Assemble the rotation matrix
    R = [N, T, W]';

    % Rotation of the primer vector in NTW ref. frame
    alpha_NTW(:, j) = R * alpha(:, j);
end

end
%-------------------------------------------------------------------------%