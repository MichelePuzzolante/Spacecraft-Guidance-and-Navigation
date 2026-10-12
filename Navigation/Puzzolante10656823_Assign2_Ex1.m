% Guidance and Navigation (2024/2025)
% Assignment # 2
% Exercise 1
% Author : Michele Puzzolante
%% REQUEST 1

% Initializzation of the MATLAB enviroment
clear;
close all;
clc;

% Set default options for plots
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');

fprintf('\nEXERCISE 1\n');
fprintf('\nREQUEST 1\n');

SAVEPLOTS = false;  % Flag to use when plots needs to be saved as png file

% Initializzation of the struct PBRFBP_const used in function to pass the
% constant values
PBRFBP_const.mu      = 1.21506683 * 10^-2;  % Earth–Moon mass parameter          [-]
PBRFBP_const.ms      = 3.28900541 * 10^5;   % Scaled masso of the Sun            [-]
PBRFBP_const.rho     = 3.88811143 * 10^2;   % Scaled Sun-(Earth+Moon) distance   [-]
PBRFBP_const.omegaS  = -9.25195985 * 10^-1; % Scaled angular velocity of the Sun [-]
PBRFBP_const.lEM     = 3.84405 * 10^8;      % Earth-Moon distance                [m]
PBRFBP_const.omegaEM = 2.66186135 * 10^-6;  % Earth-Moon angular velocity        [1/s]
PBRFBP_const.RE      = 6378;                % Mean Earth’s radius                [km]
PBRFBP_const.RM      = 1738;                % Mean Moon’s radius                 [km]
PBRFBP_const.hi      = 167;                 % Altitude of departure orbit        [km]
PBRFBP_const.hf      = 100;                 % Altitude of arrival orbit          [km]
PBRFBP_const.DU      = 3.84405000 * 10^5;   % Distance Unit                      [km]
PBRFBP_const.TU      = 4.34811305;          % Time Unit                          [days]
PBRFBP_const.VU      = 1.02323281;          % Velocity Unit                      [km/s]

% Initial state
rr_i = [-0.011965533749906; -0.017025663128129]; % Initial position
vv_i = [10.718855256727338; 0.116502348513671];  % Initial velocity
xx_i = [rr_i; vv_i];

% Initial and final time
t_i = 1.282800225339865; % Initial time
t_f = 9.595124551366348; % Final time

% Covaraince matrix 
P0 = [1.041e-15, 6.026e-17, 5.647e-16, 4.577e-15;...
      6.026e-17, 4.287e-18, 4.312e-17, 1.855e-16;...
      5.647e-16, 4.312e-17, 4.432e-16, 1.455e-15;...
      4.577e-15, 1.855e-16, 1.455e-15, 2.822e-14];

% Initializzation of the variable
nPoints = 5;                         % Number of points
tspan = linspace(t_i, t_f, nPoints); % Time vector
alpha = 1;                           % Constatnt used for the weights in UT 
beta = 2;                            % Constatnt used for the weights in UT
k = 0;                               % Constatnt used for the weights in UT
NMCpoints = 1000;                    % Number of MC samples

% Propagation using LINCOV approach
[x_hat_linCov, P_linCov] = linCov(xx_i, P0, tspan, PBRFBP_const);

% Propagation using UT approach
[x_hat_UT, P_UT] = UT(xx_i, P0, tspan, PBRFBP_const, alpha, beta, k);

% Plot ellipsesfor linCov and UT
figure(1)
hold on

% LinCov
plot(x_hat_linCov(1, end),x_hat_linCov(2, end), 'Marker','.', 'MarkerSize',26, 'Color', [0.9290 0.6940 0.1250], 'LineStyle','none')  % Plot mean
[xLCM, yLCM] = covarianceEllipse(x_hat_linCov(:, end), P_linCov(:, end-3:end), 200, 5.99);                                            % Get covariance ellipse
plot(xLCM,yLCM, 'LineWidth', 1.5, 'LineStyle','-', 'Color', [0.4660 0.6740 0.1880])                                                  % Plot covariance ellipse                                                                                            % Plot covariance ellipse

% UT                       
plot(x_hat_UT(1, end), x_hat_UT(2, end), 'Marker','.', 'MarkerSize',26, 'Color', [0.3010 0.7450 0.9330], 'LineStyle','none')         % Plot mean
[xUTM, yUTM] = covarianceEllipse(x_hat_UT(:, end), P_UT(:, end-3:end), 200, 5.99);                                                    % Get covariance ellipse
plot(xUTM,yUTM, 'LineWidth', 1.5, 'LineStyle','-', 'Color', [0 0.4470 0.7410])                                                       % Plot covariance ellipse                                                                                            % Plot covariance ellipse
xlabel('x [-]')
ylabel('y [-]')
legend('LinCov mean', 'LinCov Confidence ellipse', 'UT mean', 'UT Confidence ellipse')

% Save plots
if SAVEPLOTS

    % Saving 1st figure
    figure(1)
    exportgraphics(gcf,'LinCov_UT_confidence_ellipse.png','Resolution',600);
end

%% REQUEST 2
fprintf('\nREQUEST 2\n');

% Propagation using MC approach
[x_hat_MC, P_MC, MC_Samples] = MonteCarlo(xx_i, P0, tspan, PBRFBP_const, NMCpoints);


% Initializzation of the vectors for the 3sqrt(max(eig(Pr)))
threeSigma_r_linCov = zeros(nPoints, 1);    % LinCov
threeSigma_r_UT = zeros(nPoints, 1);        % UT
threeSigma_r_MC = zeros(nPoints, 1);        % MC

% Initializzation of the vectors for the 3sqrt(max(eig(Pv)))
threeSigma_v_linCov = zeros(nPoints, 1);    % LinConv
threeSigma_v_UT = zeros(nPoints, 1);        % UT
threeSigma_v_MC = zeros(nPoints, 1);        % MC

% Differences between three methods at t_f
error_r_linCov_UT = abs(threeSigma_r_linCov(end) - threeSigma_r_UT(end));
error_r_linCov_MC = abs(threeSigma_r_linCov(end) - threeSigma_r_MC(end));
error_r_MC_UT = abs(threeSigma_r_MC(end) - threeSigma_r_UT(end));
error_v_linCov_UT = abs(threeSigma_v_linCov(end) - threeSigma_v_UT(end));
error_v_linCov_MC = abs(threeSigma_v_linCov(end) - threeSigma_v_MC(end));
error_v_MC_UT = abs(threeSigma_v_MC(end) - threeSigma_v_UT(end));

% Computation  3sqrt(max(eig(Pr))) and 3sqrt(max(eig(Pv)))
for i = 1 : nPoints

    % 3sqrt(max(eig(Pr)))
    threeSigma_r_linCov(i) = 3 * sqrtm( max( eig(P_linCov(1:2, (1+(i-1)*4: 2+(i-1)*4)))));    % LinCov
    threeSigma_r_UT(i)     = 3 * sqrtm( max( eig(P_UT(1:2, (1+(i-1)*4: 2+(i-1)*4)))));        % UT
    threeSigma_r_MC(i)     = 3 * sqrtm( max( eig(P_MC(1:2, (1+(i-1)*4: 2+(i-1)*4)))));        % MC
    
    % 3sqrt(max(eig(Pv)))
    threeSigma_v_linCov(i) = 3 * sqrtm( max( eig(P_linCov(3:4, (3+(i-1)*4: 4+(i-1)*4)))));    % LinCov
    threeSigma_v_UT(i)     = 3 * sqrtm( max( eig(P_UT(3:4, (3+(i-1)*4: 4+(i-1)*4)))));        % UT
    threeSigma_v_MC(i)     = 3 * sqrtm( max( eig(P_MC(3:4, (3+(i-1)*4: 4+(i-1)*4)))));        % MC
end

% Computation of difference in sigma_r
sigma_r_LC_UT = threeSigma_r_linCov(end) - threeSigma_r_UT(end);
sigma_r_MC_LC = threeSigma_r_MC(end) - threeSigma_r_linCov(end);
sigma_r_MC_UT = threeSigma_r_MC(end) - threeSigma_r_UT(end);

% Computation of difference in sigma_v
sigma_v_LC_UT = threeSigma_v_linCov(end) - threeSigma_v_UT(end);
sigma_v_MC_LC = threeSigma_v_MC(end) - threeSigma_v_linCov(end);
sigma_v_MC_UT = threeSigma_v_MC(end) - threeSigma_v_UT(end);

% Display results
fprintf('\n\t\t\tLC-UT\t\t\tMC-LC\t\t\tMC-UT')
fprintf('\nsigma_r\t%.3e\t\t%.3e\t\t%.3e', sigma_r_LC_UT, sigma_r_MC_LC, sigma_r_MC_UT)
fprintf('\nsigma_v\t%.3e\t\t%.3e\t\t%.3e', sigma_v_LC_UT, sigma_v_MC_LC, sigma_v_MC_UT)

% Plot 3sqrt(max(eig(Pr)))
figure(2)
plot(tspan, threeSigma_r_linCov, 'LineStyle','-', 'Marker','o', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5)
hold on
plot(tspan, threeSigma_r_UT, 'LineStyle','-', 'Marker','o', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.5)
plot(tspan, threeSigma_r_MC, 'LineStyle','-', 'Marker','o', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5)
xlabel('t [-]')
xticks(tspan);
ylabel('3 $\sqrt{\max(\lambda(P_r))}$ [-]')
legend('LinCov', 'UT', 'MC', 'Location', 'northwest')

% Plot 3sqrt(max(eig(Pv)))
figure(3)
plot(tspan, threeSigma_v_linCov, 'LineStyle','-', 'Marker','o', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5)
hold on
plot(tspan, threeSigma_v_UT, 'LineStyle','-', 'Marker','o', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.5)
plot(tspan, threeSigma_v_MC, 'LineStyle','-', 'Marker','o', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5)
xlabel('t [-]')
xticks(tspan);
ylabel('3 $\sqrt{\max(\lambda(P_v))}$ [-]')
legend('LinCov', 'UT', 'MC', 'Location', 'northwest')

% PLot ellipsesfor linCov, UT and MC with samples
figure(4)
hold on
% LinCov
plot(x_hat_linCov(1, end),x_hat_linCov(2, end), 'Marker','.', 'MarkerSize',26, 'Color', [0.9290 0.6940 0.1250], 'LineStyle','none')  % Plot mean
plot(xLCM,yLCM, 'LineWidth', 1.5, 'LineStyle','-', 'Color', [0.4660 0.6740 0.1880])                                                  % Plot covariance ellipse
% UT                       
plot(x_hat_UT(1, end), x_hat_UT(2, end), 'Marker','.', 'MarkerSize',26, 'Color',[0.3010 0.7450 0.9330], 'LineStyle','none')         % Plot mean
plot(xUTM,yUTM, 'LineWidth', 1.5, 'LineStyle','-', 'Color', [0 0.4470 0.7410])                                                      % Plot covariance ellipse  

% MC
plot(x_hat_MC(1, end), x_hat_MC(2, end), 'Marker','.', 'MarkerSize',26, 'Color', [0.6350 0.0780 0.1840], 'LineStyle','none')         % Plot mean
[xMCM, yMCM] = covarianceEllipse(x_hat_MC(:, end), P_MC(:, end-3:end), 200, 5.99);                                                    % Get covariance ellipse
plot(xMCM,yMCM, 'LineWidth', 1.5, 'LineStyle','-', 'Color', [0.4940 0.1840 0.5560])                                                  % Plot covariance ellipse  

% Plotting samples 
for i = 1 : NMCpoints
        
    sampleMMC = MC_Samples(1+(i-1)*4:4+(i-1)*4, end);
    plot(sampleMMC(1), sampleMMC(2), 'Marker','.', 'Color',[0.8500 0.3250 0.0980], 'LineStyle','none')
end
xlabel('x [-]')
ylabel('y [-]')
legend('LinCov mean', 'LinCov Confidence ellipse', 'UT mean', 'UT Confidence ellipse', 'MC mean', 'MC Confidence ellipse', 'MC samples')

% Plot the components of the samples at t_f using qqplot
% x
figure(5)
qqplot(MC_Samples(1 : 4 : NMCpoints, end))
title('','HandleVisibility','off')
ylabel('$x(t_f)$ [-]')

% y
figure(6)
qqplot(MC_Samples(2 : 4 : NMCpoints, end))
title('','HandleVisibility','off')
ylabel('$y(t_f)$ [-]')

% vx
figure(7)
qqplot(MC_Samples(3 : 4 : NMCpoints, end))
title('','HandleVisibility','off')
ylabel('$v_x(t_f)$ [-]')

% vy
figure(8)
qqplot(MC_Samples(4 : 4 : NMCpoints, end))
title('','HandleVisibility','off')
ylabel('$v_y(t_f)$ [-]')

% Save plots
if SAVEPLOTS

    % Saving 2nd figure
    figure(2)
    exportgraphics(gcf,'3sigma_r.png','Resolution',600);

    % Saving 3rd figure
    figure(3)
    exportgraphics(gcf,'3sigma_v.png','Resolution',600)

    % Saving 4th figure
    figure(4)
    exportgraphics(gcf,'LinCov_UT_MT_samples_confidence_ellipse.png','Resolution',600)

    % Saving 5th figure
    figure(5)
    exportgraphics(gcf,'x_pos_stdnorm.png','Resolution',600)

    % Saving 6th figure
    figure(6)
    exportgraphics(gcf,'y_pos_stdnorm.png','Resolution',600)

    % Saving 7th figure
    figure(7)
    exportgraphics(gcf,'x_vel_stdnorm.png','Resolution',600)

    % Saving 8th figure
    figure(8)
    exportgraphics(gcf,'y_vel_stdnorm.png','Resolution',600)
end

%% FUNCTION

%-------------------------------------------------------------------------%
% function:PBRFBP_propagation.m 
function [xxf, Phif, tf, xx, tt] = PBRFBP_propagation(xx0, tspan, const, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% PBRFBP_propagation: function for propagating state and the STM using the PBRFBP dynamics                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xxf, Phif, tf, xx, tt] = PBRFBP_propagation(xx0, tspan, const, false)
%   [xxf, Phif, tf, xx, tt] = PBRFBP_propagation(xx0, tspan, const, true)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx0          [4x1] or [20x1] initial state                    [-; -; -; -]
%  t0           [nx1] time vector                                [-]
%  const        [-]   struct of constant values                  [-]
%  flagSTM      [1]   flag used to activate STM propagation      [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  xxf         [4x1]   final state                               [-; -; - ; -]
%  Phif        [4x4]   STM at final time                         [-]
%  tf          [1]     final time                                [-]
%  xx          [nx4]   state from t0 to tf                       [-]
%  tt          [nx1]   time from t0 to tf                        [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  18/10/2024: First version
%
%-------------------------------------------------------------------------%

% Set integartion tolerance
Reltol = 3e-14;     % Relative tolerance
Abstol_xx = 1e-14;   % Absolute tolerance for the state
Abstol_STM = 1e-14; % Absolute tolerance for the STM

% Checking if the integration of the STM dynamics is needed
if flagSTM

    % Set options
    options = odeset('reltol', Reltol, 'abstol', [Abstol_xx * ones(4, 1); Abstol_STM * ones(16,1)]); 

    % Set initial condition for STM
    Phi0 = eye(4);

    % Perform integration
    [tt, xx] = ode113(@(t, x) PBRFBP_dynamics_STM(t, x, const, flagSTM), tspan, [xx0; Phi0(:)], options);
    
    % Extract STM 
    Phif = reshape(xx(end,5:end),4,4);
else

    % Set options
    options = odeset('reltol', Reltol, 'abstol', Abstol_xx * ones(4, 1));

    % Perform integration
    [tt, xx] = ode113(@(t, x) PBRFBP_dynamics_STM(t, x, const, flagSTM), tspan, xx0, options);

    % Set Phif equal to NaN if STM is not needed
    Phif = NaN(4, 4);
end

% Extract final state vector 
xxf = xx(end,1:4)';

% Extract the final time
tf = tt(end);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: PBRFBP_dynamics_STM.m 
function [dxdt] = PBRFBP_dynamics_STM(t, xx, const, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% PBRFBP_dynamics_STM: function dynamics for the Planar Bicircular Restricted Four-Body Problem , with dynamics of STM                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [dxdt] = PBRFBP_dynamics_STM(t, xx, const, false) 
%   [dxdt] = PBRFBP_dynamics_STM(t, xx, const, true)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  t            [1]   integration time                           [-]
%  xx           [4x1] or [20x1] initial state                    [-; -; -; -]
%  const        [-]   struct of constant values                  [-]
%  flagSTM      [1]   flag used to activate STM propagation      [-]              
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  dxdt         [4x1] or [20x1]     ODE                          [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  18/10/2024: First version
%
%-------------------------------------------------------------------------%

% Get the constant values from the struct
mu = const.mu;          % Earth–Moon mass parameter
ms = const.ms;          % scaledmassoftheSun
rho = const.rho;        % scaled Sun–(Earth+Moon) distance
omegaS = const.omegaS;  % scaled angular velocity of the Sun

% Extract position and velocity from the state
x = xx(1);  % x - position
y = xx(2);  % y - position
vx = xx(3); % x - velocity
vy = xx(4); % y - velocity 

% Computation of the distances
r1 = (x + mu)^2 + y^2;                                      % P w.r.t. Earth
r2 = (x + mu - 1)^2 + y^2;                                  % P w.r.t. Moon
r3 = (x - rho*cos(omegaS*t))^2 + (y - rho*sin(omegaS*t))^2;  % P w.r.t. Sun

% Checking if the integration of the STM dynamics is needed
if flagSTM
    
    % Initializzation of the right-hand-side
    dxdt = zeros(20, 1);

    % Computation of matrix A for the dynamic od STM
    A = [0, 0, 1, 0;...
         0, 0, 0, 1;...
         1 - (1 - mu)/r1^(3/2) + (3*(1 - mu)*(x + mu)^2)/r1^(5/2) -...     % dxddot/dx
         mu/r2^(3/2) + (3*mu*(x + mu - 1)^2)/r2^(5/2) -...                 % dxddot/dx
         ms/r3^(3/2) + (3 * ms * (x - rho*cos(omegaS*t))^2)/r3^(5/2), ...  % dxddot/dx
         3*(1 - mu)*(x + mu)*y/r1^(5/2) +...                               % dxddot/dy
         3*mu*y*(x + mu - 1)/r2^(5/2) +...                                 % dxddot/dy
         3*ms*(x - rho*cos(omegaS*t))*(y - rho*sin(omegaS*t))/r3^(5/2),... % dxddot/dy
         0, 2;...
         3*(1 - mu)*(x + mu)*y/r1^(5/2) +...                               % dyddot/dx
         3*mu*y*(x + mu - 1)/r2^(5/2) +...                                 % dyddot/dx
         3*ms*(x - rho*cos(omegaS*t))*(y - rho*sin(omegaS*t))/r3^(5/2),... % dyddot/dx
         1 - (1 - mu)/r1^(3/2) + (3*(1 - mu)*(y)^2)/r1^(5/2) -...          % dyddot/dy
         mu/r2^(3/2) + (3*mu*(y)^2)/r2^(5/2) -...                          % dyddot/dy
         ms/r3^(3/2) + (3 * ms * (y - rho*sin(omegaS*t))^2)/r3^(5/2), ...  % dyddot/dy
         -2, 0];

    % Extract the STM from the state
    Phi = reshape(xx(5:end),4,4); 
    
    % Compute the derivative of the STM
    Phidot = A * Phi ;
    
    % Define dinamyc of STM
    dxdt(5:end) = Phidot(:);
else

    % Initializzation of the right-hand-side
    dxdt = zeros(4, 1);
end

% % Define dinamycs of PBRFBP                                                         
dxdt(1) = vx;   % x - velocity
dxdt(2) = vy;   % y - velocity
dxdt(3) = 2*vy + x - (ms*cos(omegaS*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(omegaS*t)))/(2*((x - rho*cos(omegaS*t))^2 + (y - rho*sin(omegaS*t))^2)^(3/2));
dxdt(4) = y - 2*vx - (ms*sin(omegaS*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(omegaS*t)))/(2*((x - rho*cos(omegaS*t))^2 + (y - rho*sin(omegaS*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: linCov.m 
function [x_hat, P] = linCov(x0_hat, P0, tspan, const)
%-------------------------------------------------------------------------%
%                                                                         
% linCov: function for computing mean and covariance through
%         the linear covariance method               
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [x_hat, P] = linCov(x0_hat, P0, tspan, const)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0_hat       [4x1]  initial mean state vector                       [-; -; -; -]
%  P0           [4x4]  initial covariance matrix                       [-]
%  tspan        [nx1]  time vector                                     [-]
%  const        [-]   struct of constant values                        [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  x_hat          [4xn]    mean state vector                           [-; -; -; -]
%  P              [4xn*4]  covariance matrix vector                    [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/12/2024: First version
%
%-------------------------------------------------------------------------%

nPoints = length(tspan);    % Get number of points in the time vector

% Initializzation of  output variable
x_hat = zeros(4, nPoints);  % Mean state
P = zeros(4, 4 * nPoints);  % Covariance matrix

% Assign the initial value of mean and covariance to the outputs
x_hat(:, 1) = x0_hat;
P(:, 1:4) = P0;

for i = 2 : nPoints

    % Propagation
    [xxf, Phif, ~, ~, ~] = PBRFBP_propagation(x_hat(:, i - 1), [tspan(i - 1)  tspan(i)], const, true);

    % Propagated mean
    x_hat(:, i) = xxf;               

    % Compute the covariance matrix
    P(:, (1+(i-1)*4: 4+(i-1)*4)) = Phif * P(:, (1+(i-2)*4: 4+(i-2)*4)) * Phif';              
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: UT.m 
function [x_hat, P] = UT(x0_hat, P0, tspan, const, alpha, beta, k)
%-------------------------------------------------------------------------%
%                                                                         
% UT: function for computing mean and covariance through
%     the unscented transform               
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [x_hat, P] = UT(x0_hat, P0, tspan, const, alpha, beta, k)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0_hat       [4x1]  initial mean state vector                      [-; -; -; -]
%  P0           [4x4]  initial covariance matrix                      [-]
%  tspan        [nx1]  time vector                                    [-]
%  const        [-]   struct of constant values                       [-]
%  alpha        [1]    parameter for computing weights                [-] 
%  beta         [1]    parameter for computing weights                [-] 
%  k            [1]    parameter for computing weights                [-] 
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  x_hat        [4xn]    mean state vector                            [-; -; -; -]
%  P            [4xn*4]  covariance matrix vector                     [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/12/2024: First version
%
%-------------------------------------------------------------------------%

N = length(x0_hat);         % Get number of value in the state vector
nPoints = length(tspan);    % Get number of points in the time vector

c = alpha^2 * (N + k);      % Compute c value used for weights

% Initializzation of  output variable
x_hat = zeros(4, nPoints);  % Mean state
P = zeros(4, 4 * nPoints);  % Covariance matrix

% Assign the initial value of mean and covariance to the outputs
x_hat(:, 1) = x0_hat;
P(:, 1:4) = P0;

% Creation of sigma points
X0 = x0_hat;                % Sigma point 0
X = zeros(4, 2*N + 1);      % Matrix used for the sigma points
X(:, 1) = X0;

for i = 1:2*N
    if i <= N
        Pxi = sqrtm(c * P0);
        dxi = Pxi(:, i);
    else
        Pxi = -sqrtm(c * P0);
        dxi = Pxi(:, i-N);
    end
X(:, i + 1) = x0_hat + dxi; % Assemble the matrix of sigma points
end


% Propagation simga points
Y = zeros(4 * (2 * N + 1), nPoints); % Initializzation of propagate sigma points matrix
for i = 1 : 2*N+1

    % Propagation
    [~, ~, ~, xx, ~] = PBRFBP_propagation(X(:, i), tspan, const, false);  

    % Assemble the matrix of propagated sigma points
    Y((1+(i-1)*4: 4+(i-1)*4), :) = xx';                                                        
end

% Computation of the mean
W0_m = 1 - N/c;     % Mean's weight for sigma point n.0
Wi_m = 1/(2 * c);   % Mean's weight for sigma points
for j = 2 : nPoints
    for i = 1 : 2*N + 1
        if i == 1
            x_hat(:, j) = x_hat(:, j) + W0_m * Y((1+(i-1)*4: 4+(i-1)*4), j); 
        else
            x_hat(:, j) = x_hat(:, j) + Wi_m * Y((1+(i-1)*4: 4+(i-1)*4), j); 
        end
    end
end

% Computation of the covariance
W0_c = (2 - alpha^2 + beta) - N/c;  % Covariance's weight for sigma point n.0
Wi_c = 1/(2 * c);                   % Covariance's weight for sigma points
for j = 2 : nPoints
    for i = 1 : 2*N + 1
        if i == 1
            P(:, (1+(j-1)*4: 4+(j-1)*4)) = P(:, (1+(j-1)*4: 4+(j-1)*4)) + W0_c * (Y((1+(i-1)*4: 4+(i-1)*4), j) - x_hat(:, j)) .* (Y((1+(i-1)*4: 4+(i-1)*4), j) - x_hat(:, j))'; 
        else
            P(:, (1+(j-1)*4: 4+(j-1)*4)) = P(:, (1+(j-1)*4: 4+(j-1)*4)) + Wi_c * (Y((1+(i-1)*4: 4+(i-1)*4), j) - x_hat(:, j)) .* (Y((1+(i-1)*4: 4+(i-1)*4), j) - x_hat(:, j))'; 
        end
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: MonteCarlo.m 
function [x_hat, P, Y] = MonteCarlo(x0_hat, P0, tspan, const, NMCpoints)
%-------------------------------------------------------------------------%
%                                                                         
% MonteCarlos: function for computing mean and covariance through
% MonteCarlo simulation               
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [x_hat, P, Y] = MonteCarlo(x0_hat, P0, tspan, const, NMCpoints)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0_hat       [4x1]  initial mean state vector                       [-; -; -; -]
%  P0           [4x4]  initial covariance matrix                       [-]
%  tspan        [nx1]  time vector                                     [-] 
%  const        [-]   struct of constant values                       [-]
%  NMCpoints    [1]    number of points for MonteCarlo simulation      [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  x_hat        [4xn]    mean state vector                             [-; -; -; -]
%  P            [4xn*4]  covariance matrix vector                      [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/12/2024: First version
%
%-------------------------------------------------------------------------%

nPoints = length(tspan);    % Get number of points in the time vector

% Initializzation of  output variable
x_hat = zeros(4, nPoints);  % Mean state
P = zeros(4, 4 * nPoints);  % Covariance matrix

% Assign the initial value of mean and covariance to the outputs
x_hat(:, 1) = x0_hat;
P(:, 1:4) = P0;

% Generate points
X = mvnrnd(x0_hat, P0, NMCpoints)';

% Intilaize the outputs
Y = zeros(4 * (2 * NMCpoints + 1), nPoints); % Initializzation of propagate sigma points matrix
for i = 1 : NMCpoints

    % Propagation
    [~, ~, ~, xx, ~] = PBRFBP_propagation(X(:, i), tspan, const, false);  

    % Assemble the matrix of propagated sigma points
    Y((1+(i-1)*4: 4+(i-1)*4), :) = xx';                                 
end

% Compute sample mean
for j = 2 : nPoints
    for i = 1 : NMCpoints
        x_hat(:, j) = x_hat(:, j) + Y((1+(i-1)*4: 4+(i-1)*4), j);
    end
    x_hat(:, j) = x_hat(:, j)/NMCpoints;
end

% Compute sample covariance
for j = 2 : nPoints
    for i = 1 : NMCpoints
         P(:, (1+(j-1)*4: 4+(j-1)*4)) = P(:, (1+(j-1)*4: 4+(j-1)*4)) + (Y((1+(i-1)*4: 4+(i-1)*4), j) - x_hat(:, j)) .* (Y((1+(i-1)*4: 4+(i-1)*4), j) - x_hat(:, j))';
    end
    P(:, (1+(j-1)*4: 4+(j-1)*4)) =  P(:, (1+(j-1)*4: 4+(j-1)*4)) / (NMCpoints - 1);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: covarianceEllipse.m  
function [x, y] = covarianceEllipse(mean, P, n, chi)
%-------------------------------------------------------------------------%
%                                                                         
% covarianceEllipse: function which gives n points of an covariance ellipse 
% with a given mean and covariance matrix
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [x, y] = covarianceEllipse(mean, P, R, n, chi)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  mean       [6x1]         mean vector         [km;km;km;km/s;km/s;km/s]
%  P          [1]           covariance matrix                        [-]  
%  R          [1]           rotation matrix                          [-]
%  n          [1]           nuber of points for x,y coord            [-]  
%  chi        [1]           confidence, tipically 95%                [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  x          [nx1]     x-coordinate                                [km]
%  y          [nx1]     y-coordinate                                [km]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  01/12/2024: First version
%
%-------------------------------------------------------------------------%

P = P(1:2, 1:2); 
lambda = eig(P);      % Compute eigenvalues of P
lambda1 = max(lambda);     % Assign maximum eigenvalue
lambda2 = min(lambda);     % Assign minimum eigenvalue

% Compute rotation angle of the ellipse
if P(1, 2) == 0 && P(1, 1) >= P(2, 2)
    theta = 0;
elseif P(1, 2) == 0 && P(1, 1) < P(2, 2)
    theta = pi/2;
else
    theta = atan2(lambda1 - P(1, 1), P(1, 2));
end
theta = wrapTo2Pi(theta);   % Wrap the angle 

% Compute the x coordinate
x = mean(1) + sqrt(chi * lambda1) * cos(theta) * cos(linspace(0, 2*pi, n)) - sqrt(chi * lambda2) * sin(theta) * sin(linspace(0, 2*pi, n));

% Compute the y coordinate
y = mean(2) + sqrt(chi * lambda1) * sin(theta) * cos(linspace(0, 2*pi, n)) + sqrt(chi * lambda2) * cos(theta) * sin(linspace(0, 2*pi, n));

end
%-------------------------------------------------------------------------%