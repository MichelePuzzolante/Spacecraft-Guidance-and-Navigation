% Guidance and Navigation (2024/2025)
% Assignment # 1
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
mu = 0.012150;     % Earth - Moon mass parameter

% Defnition of the OMEGA function
OMEGA = @(x, y, z) 1/2 * (x^2 + y^2) + (1 - mu)/((x + mu)^2 + y^2 + z^2)^(1/2) + (mu)/((x + mu - 1)^2 + y^2 + z^2)^(1/2) + 1/2 * mu * (1 - mu);

% Definition of the integral of motion
J = @(x, y, z, vx, vy, vz) 2 * OMEGA(x, y, z) - (vx^2 + vy^2 + vz^2);

% Defnition of the derivative of OMEGA w.r.t. x
dOMEGAdx = @(x, y, z) x - (1 - mu) * (x + mu)./((x + mu).^2 + y.^2 + z.^2).^(3/2) - (mu) * (x + mu - 1)./((x + mu - 1).^2 + y.^2 + z.^2).^(3/2);

% Defnition of the derivative of OMEGA w.r.t. y
dOMEGAdy = @(x, y, z) y - (1 - mu) * (y)./((x + mu).^2 + y.^2 + z.^2).^(3/2) - (mu) * (y)./((x + mu - 1).^2 + y.^2 + z.^2).^(3/2);

% Defnition of the derivative of OMEGA w.r.t. z
dOMEGAdz = @(x, y, z) - (1 - mu) * (z)./((x + mu).^2 + y.^2 + z.^2).^(3/2) - (mu) * (z)./((x + mu - 1).^2 + y.^2 + z.^2).^(3/2);

% Definition of the Derivative of J w.r.t x, z, and vy
dJdxzvy = @(x, y, z, vx, vy, vz) [dOMEGAdx(x, y, z), dOMEGAdz(x, y, z), -2*vy];

% Computation of the 5 Lagrangian points

% L1, L2 and L3 whcih lies on the x axis
xL1 = fzero(@(x) dOMEGAdx(x, 0, 0), 0.5);
xL2 = fzero(@(x) dOMEGAdx(x, 0, 0), 1.5);
xL3 = fzero(@(x) dOMEGAdx(x, 0, 0), -1);

% L4 and L5 have y different from 0 and x = 1/2 - mu
xL4L5 = 1/2 - mu;
yL4 = fzero(@(y) dOMEGAdy(xL4L5, y, 0), 1);
yL5 = fzero(@(y) dOMEGAdy(xL4L5, y, 0), -1);

% Computation of the 5 Jacobi constants associated to Lagrangian points
C1 = J(xL1, 0, 0, 0, 0, 0);
C2 = J(xL2, 0, 0, 0, 0, 0);
C3 = J(xL3, 0, 0, 0, 0, 0);
C4 = J(xL4L5, yL4, 0, 0, 0, 0);
C5 = J(xL4L5, yL5, 0, 0, 0, 0);

% Check if the condition on the Jacobi constants is satisfied
if C1 > C2 && C2 > C3 && C3 > C4 && C4 == C5
    fprintf('\nThe condition C1 > C2 > C3 > C4 = C5 is respected\n');
    fprintf('L1:\tx = %.10f,\ty = %.10f,\tz = 0,\tC1 = %.10f\n', xL1, 0, C1)
    fprintf('L2:\tx = %.10f,\ty = %.10f,\tz = 0,\tC2 = %.10f\n', xL2, 0, C2)
    fprintf('L3:\tx = %.10f,\ty = %.10f,\tz = 0,\tC3 = %.10f\n', xL3, 0, C3)
    fprintf('L4:\tx = %.10f,\ty = %.10f,\tz = 0,\tC4 = %.10f\n', xL4L5, yL4, C4)
    fprintf('L5:\tx = %.10f,\ty = %.10f,\tz = 0,\tC5 = %.10f\n', xL4L5, yL5, C5)
else
    fprintf('\nError: The condition C1 > C2 > C3 > C4 = C5 is not respected\n');
end

% vectors of values used to plot dOMEGA/dx and dOMEGA/dy
x1 = linspace(-2, -mu-1e-10, 100);
x2 = linspace(-mu+1e-10, 1-mu-1e-10, 100);
x3 = linspace(1-mu+1e-10, 2, 100);

% Plot dOMGEA/dx
figure(1)
plot(-mu, 0, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', [0 0.4470 0.7410])          % Earth
hold on
plot(1 - mu, 0, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', [0.6350 0.0780 0.1840])  % Moon
plot(x1,  dOMEGAdx(x1, 0, 0), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off')
plot(x2,  dOMEGAdx(x2, 0, 0), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off')
plot(x3,  dOMEGAdx(x3, 0, 0), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off')
plot(xL1, 0, 'LineStyle', 'none', 'Marker','*','MarkerSize', 10, 'Color', 'r', 'LineWidth', 2)
plot(xL2, 0, 'LineStyle', 'none', 'Marker','*','MarkerSize', 10, 'Color', 'g', 'LineWidth', 2)
plot(xL3, 0, 'LineStyle', 'none', 'Marker','*','MarkerSize', 10, 'Color', 'b', 'LineWidth', 2)
xlabel('x [-]', 'FontSize',15);
ylabel('$\frac{\partial \Omega}{\partial x}$ [-]', 'FontSize',15, 'Rotation', 0);
legend('Earth', 'Moon', '$L_1$', '$L_2$', '$L_3$','FontSize',15);
ylim([-10 10])

% Plot dOMEGA/dy
figure(2)
plot(x1, dOMEGAdy(xL4L5, x1, 0), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off')
hold on
plot(x2, dOMEGAdy(xL4L5, x2, 0), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off')
plot(x3, dOMEGAdy(xL4L5, x3, 0), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off')
plot(yL4, 0, 'LineStyle', 'none', 'Marker','x','MarkerSize', 10, 'Color', 'm', 'LineWidth', 2)
plot(yL5, 0, 'LineStyle', 'none', 'Marker','x','MarkerSize', 10, 'Color', 'c', 'LineWidth', 2)
xlabel('y [-]', 'FontSize',15);
ylabel('$\frac{\partial \Omega}{\partial y}$ [-]', 'FontSize',15, 'Rotation', 0);
legend('$L_4$', '$L_5$','FontSize',15);

% Plot the Lagrange points
figure(3)
plot(-mu, 0, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', [0 0.4470 0.7410])          % Earth
hold on
plot(1 - mu, 0, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20, 'Color', [0.6350 0.0780 0.1840])  % Moon
plot(xL1, 0, 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 10, 'Color', 'r', 'LineWidth', 2)        % L1
plot(xL2, 0, 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 10, 'Color', 'g', 'LineWidth', 2)        % L2
plot(xL3, 0, 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 10, 'Color', 'b', 'LineWidth', 2)        % L3
plot(xL4L5, yL4, 'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 10, 'Color', 'm', 'LineWidth', 2)    % L4
plot(xL4L5, yL5, 'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 10, 'Color', 'c', 'LineWidth', 2)    % L5
plot(-mu + cos(0 :0.1:2*pi), sin(0:0.1:2*pi), 'LineStyle','--', 'LineWidth',0.5, 'HandleVisibility','off','Color','black')
line([-mu xL4L5], [0 yL4], 'LineWidth', 0.5, 'Color', 'black', 'LineStyle','--'); 
line([1-mu xL4L5], [0 yL4], 'LineWidth', 0.5, 'Color', 'black', 'LineStyle','--'); 
line([-mu xL4L5], [0 yL5], 'LineWidth', 0.5, 'Color', 'black', 'LineStyle','--'); 
line([1-mu xL4L5], [0 yL5], 'LineWidth', 0.5, 'Color', 'black', 'LineStyle','--');
xlabel('x [-]', 'FontSize',15);
ylabel('y [-]', 'FontSize',15);
legend('Earth', 'Moon', '$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$',...
       'Location','bestoutside', 'FontSize', 15);
hold off
xlim([-1.2 1.2])
axis equal

% Save plots
if SAVEPLOTS

    % Saving 1st figure
    figure(1)
    exportgraphics(gcf,'dOMEGA_dx.png','Resolution',600);

    % Saving 2nd figure
    figure(2)
    exportgraphics(gcf,'dOMEGA_dy.png','Resolution',600);

    % Saving 3rd figure
    figure(3)
    exportgraphics(gcf,'Lagrange_Points.png','Resolution',600);
end

%% REQUEST 2
% The following section has been used to find the periodic orbit passing
% through xx0 with Jacobi constant equal to Cf

fprintf('\nREQUEST 2\n');

% Initializzation fo the initial state
x0 = 1.068792441776;  % x - position
y0 = 0;               % y - position
z0 = 0.071093328515;  % z - position
vx0 = 0;              % x - velocity
vy0 = 0.319422926485; % y - velocity
vz0 = 0;              % z - velocity

% Assemble the initial state
xx0 = [x0; y0; z0; vx0; vy0; vz0];

% Initializzation of the variables
Cf = 3.09;       % Final value of Jacobi constant
error_yf = 1;    % error in y - position
error_vxf = 1;   % error in x - velocity
error_vzf = 1;   % error in z - velocity
error_Cf = 1;    % error for Jacobi constant
tol_R2 = 1e-14;  % tolerance used to stop the correction algorithm
Nmax_iter = 100; % Maximum number of iteration, used to avoid infinte loop
N_iter = 1;      % Variable used to count the number of iterations
t0 = -0.1;       % Initial time of integration

% Compute a first guess for final time
tf = 0;        % Final time of integration, has been set has the time when the orbit cross the x-axis

Phi0 = eye(6); % Initial contition for STM

% Initializzation of the new state vector
xx0_new_R2 = xx0;

% While used to correct the initial state in order to have a periodic orbit
% with C = 3.09
while (error_vxf > tol_R2 || error_vzf > tol_R2 || error_Cf > tol_R2 || error_yf > tol_R2) && N_iter <= Nmax_iter
    
    % Propagation of the initial state corrected
    [xxf, Phif, te, ~, ~] = CRTBP_propagation([xx0_new_R2; reshape(Phi0(:),36,1)], t0, tf, mu, true);

    % Computation of the f(xx0)
    dflow_dt = CRTBP_dynamic_STM(t0, xx0_new_R2, mu, false);

    % Get the component of the STM used to correct the initial state
    STM_correction = [Phif(2,1), Phif(2,3), Phif(2,5), dflow_dt(2);...
                      Phif(4,1), Phif(4,3), Phif(4,5), dflow_dt(4);...
                      Phif(6,1), Phif(6,3), Phif(6,5), dflow_dt(6)];
    
    
    % Computation of the Jacobi constant associated to the corrected
    % initial state
    Cf_Prop = J(xx0_new_R2(1), xx0_new_R2(2), xx0_new_R2(3), xx0_new_R2(4), xx0_new_R2(5), xx0_new_R2(6));
    
    % Assebmble the vector of the errors
    dxxTrue_Prop = [-xxf(2); -xxf(4); -xxf(6); (Cf - Cf_Prop)];
    
    % Assemble the Matrix used for compute the corrections
    CorrectionMatrix = [STM_correction; dJdxzvy(xx0_new_R2(1), xx0_new_R2(2), xx0_new_R2(3), xx0_new_R2(4), xx0_new_R2(5), xx0_new_R2(6)), 0];

    % Compute the corrections
    dxx = CorrectionMatrix\dxxTrue_Prop;
   
    % add the correction to the initial state
    xx0_new_R2 = xx0_new_R2 + [dxx(1); 0; dxx(2); 0; dxx(3); 0];
    
    % Add the correction to the initial time
    t0 = t0 + dxx(4);
   
    % Computation of the errors
    error_yf = abs(xxf(2));       % y - position
    error_vxf = abs(xxf(4));      % x - velocity
    error_vzf = abs(xxf(6));      % z - velocity
    error_Cf = abs(Cf - Cf_Prop); % Jacobi constant
    
    % Update nuber of iterations
    N_iter = N_iter + 1;          
end

% Display the results
fprintf('\nCORRECTED INITIAL STATE OF THE HALO ORBIT WITH C = 3.09\n')
fprintf('x0\t=\t%.10f\t[-]\n',xx0_new_R2(1))
fprintf('y0\t=\t%.10f\t[-]\n',xx0_new_R2(2))
fprintf('z0\t=\t%.10f\t[-]\n',xx0_new_R2(3))
fprintf('vx0\t=\t%.10f\t[-]\n',xx0_new_R2(4))
fprintf('vy0\t=\t%.10f\t[-]\n',xx0_new_R2(5))
fprintf('vz0\t=\t%.10f\t[-]\n',xx0_new_R2(6))
fprintf('T\t=\t%.10f\t[-]\n',abs(2 * (tf-t0)))

% Plot the halo horbit
[~, ~, ~, xx_Halo_forward, ~] = CRTBP_propagation(xx0_new_R2, t0, tf, mu, false);
[~, ~, ~, xx_Halo_backward, ~] = CRTBP_propagation(xx0_new_R2, -t0, -tf, mu, false);
figure(4)
plot3(-mu, 0, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',30, 'Color', [0 0.4470 0.7410])
hold on
plot3(1-mu, 0, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',30, 'Color', [0.6350 0.0780 0.1840])
plot3(xx_Halo_forward(:, 1), xx_Halo_forward(:,2), xx_Halo_forward(:, 3), 'LineWidth',2, 'Color', [0.4940 0.1840 0.5560])
plot3(xx_Halo_backward(:, 1), xx_Halo_backward(:,2), xx_Halo_backward(:, 3), 'LineWidth',2, 'Color', [0.3010 0.7450 0.9330])
xlabel('x [-]', 'FontSize',15);
ylabel('y [-]', 'FontSize',15);
zlabel('z [-]', 'FontSize',15);
legend('Earth', 'Moon', 'Forward', 'Backward','FontSize',12, 'Location', 'northwest');

% Save plots
if SAVEPLOTS

    % Saving 4th figure
    f4 = figure(4);
    exportgraphics(gcf,'Halo_orbit.png','Resolution',600);
end

%% REQUEST 3

fprintf('\nREQUEST 3\n');

% Definition of constant values
N_points = 20;                              % Number of poinyìts for Jacobi constants                      
Cf_end = 3.04;                              % Last value of Jacobi constant
C_vector = linspace(Cf, Cf_end, N_points);  % Vector of Jacobi constants

% Initializzation of the variables
tol_R3 = 1e-12;     % tolerance used to stop the correction algorithm
Nmax_iter = 100;    % Maximum number of iteration, used to avoid infinte loop
t0_new_R3 = t0;

% Initializzation of the new state vector
xx0_new_R3 = xx0_new_R2;

% Plot the position of Earth and Moon
figure(5)
plot3(-mu, 0, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',30, 'Color', [0 0.4470 0.7410])
hold on
plot3(1-mu, 0, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',30, 'Color', [0.6350 0.0780 0.1840])

% For cycle used to perform the numerical continuation for the values of C
for i = 1 : N_points

    % Initializzation of the variables used in the while to avoid numerical
    % issues
    error_yf = 1;  % Error in y - position
    error_vxf = 1; % Error in x - velocity
    error_vzf = 1; % Error in z - velocity
    error_Cf = 1;  % error for Jacobi constant
    
    % Reset number of iterations
    N_iter = 1; 

    % Extarct the current Jacobi constant value
    Cf_Family = C_vector(i);

    % While for the computation of the corrected initial state
    while (error_vxf > tol_R3 || error_vzf > tol_R3 || error_Cf > tol_R3 || error_yf > tol_R3) && N_iter <= Nmax_iter
    
        % Propagation of the initial state corrected
        [xxf, Phif, te, ~, ~] = CRTBP_propagation([xx0_new_R3; reshape(Phi0(:),36,1)], t0_new_R3, tf, mu, true);
    
        % Computation of the f(xx0)
        dflow_dt = CRTBP_dynamic_STM(t0, xx0_new_R3, mu, false);

        % Get the component of the STM used to correct the initial state
        STM_correction = [Phif(2,1), Phif(2,3), Phif(2,5), dflow_dt(2);...
                          Phif(4,1), Phif(4,3), Phif(4,5), dflow_dt(4);...
                          Phif(6,1), Phif(6,3), Phif(6,5), dflow_dt(6)];
        
        % Computation of the Jacobi constant associated to the corrected
        % initial state
        Cf_Prop = J(xx0_new_R3(1), xx0_new_R3(2), xx0_new_R3(3), xx0_new_R3(4), xx0_new_R3(5), xx0_new_R3(6));   
    
        % Assebmble the vector of the errors
        dxxTrue_Prop = [-xxf(2); -xxf(4); -xxf(6); (Cf_Family - Cf_Prop)];
        
        % Assemble the Matrix used for compute the corrections
        CorrectionMatrix = [STM_correction; dJdxzvy(xx0_new_R3(1), xx0_new_R3(2), xx0_new_R3(3), xx0_new_R3(4), xx0_new_R3(5), xx0_new_R3(6)), 0];
        
        % Compute the corrections
        dxx = CorrectionMatrix\dxxTrue_Prop;
        
        % Update the initial state
        xx0_new_R3 = xx0_new_R3 + [dxx(1); 0; dxx(2); 0; dxx(3); 0];
        
        % Update the initial time
        t0_new_R3 = t0_new_R3 + dxx(4);

        % Computation of the errors
        error_yf = abs(xxf(2));              % y - position
        error_vxf = abs(xxf(4));             % x - velocity
        error_vzf = abs(xxf(6));             % z - velocity
        error_Cf = abs(Cf_Family - Cf_Prop); % Jacobi constant

        % Update nuber of iterations
        N_iter = N_iter + 1;
    end
    
    % Propagation of the correted initail state for the i-th value of C
    [~, ~, ~, xx_Family, ~] = CRTBP_propagation(xx0_new_R3, 0, 2*(tf - t0_new_R3), mu, false);

    % Plot the halo orbit
    figure(5)
    plot3(xx_Family(:, 1), xx_Family(:,2), xx_Family(:, 3), 'Color', (Cf - Cf_Family)/(Cf - Cf_end)*[0, 1, 1] + (Cf_Family - Cf_end)/(Cf - Cf_end) * [1, 0, 1], 'LineWidth', 2);
end

% Setting up the plot
figure(5)
xlabel('x [-]', 'FontSize',15);
ylabel('y [-]', 'FontSize',15);
zlabel('z [-]', 'FontSize',15);
legend('Earth', 'Moon','FontSize',12, 'Location', 'northwest');

% Setting the colorbar options
colormap(cool)
colorB = colorbar;                           % Create the colorbar
colorB.Label.String = 'C';                   % Put label
colorB.Label.FontSize = 15;                  % Size of the label
colorB.Label.VerticalAlignment = 'bottom';   % Vertical aligment
colorB.Label.Interpreter = 'latex';          % Label interpreter
colorB.Label.Rotation = 0;                   % Label rotation
colorB.Label.Position(1) = 0.5;              % Label position
colorB.Label.Position(2) = 3.09;             % Label position
colorB.TickLabelInterpreter = 'latex';       % Interpreter for values
colorB.FontSize = 12;                        % Size of values  
clim([Cf_end Cf]);                           % Set limit for the colorbar

% Display the results
fprintf('\nCORRECTED INITIAL STATE OF THE HALO ORBIT WITH C = 3.04\n')
fprintf('x0\t=\t%.10f\t[-]\n',xx0_new_R3(1))
fprintf('y0\t=\t%.10f\t[-]\n',xx0_new_R3(2))
fprintf('z0\t=\t%.10f\t[-]\n',xx0_new_R3(3))
fprintf('vx0\t=\t%.10f\t[-]\n',xx0_new_R3(4))
fprintf('vy0\t=\t%.10f\t[-]\n',xx0_new_R3(5))
fprintf('vz0\t=\t%.10f\t[-]\n',xx0_new_R3(6))
fprintf('T\t=\t%.10f\t[-]\n',abs(2 * t0_new_R3))

% Save plots
if SAVEPLOTS

    % Saving 5th figure
    f5 = figure(5);
    view(-23.4633,34.2745);
    exportgraphics(gcf,'Halo_orbit_family.png','Resolution',600);
end

%% FUNCTIONS

%-------------------------------------------------------------------------%
% function:CRTBP_propagation.m 
function [xxf, Phif, tf, xx, tt] = CRTBP_propagation(xx0, t0, tf, mu, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% CRTBP_propagation: function for propagating state through the CRTBP dynamics equations              
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xxf, Phif, tf, xx, tt] = CRTBP_propagation(xx0, t0, tf, mu, true) 
%   [xxf, Phif, tf, xx, tt] = CRTBP_propagation(xx0, t0, tf, mu, false)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  x0           [42x1] initial state                              [-]
%  t0           [1]    intial time                                [-]
%  tf           [1]    final time                                 [-]
%  mu           [1]    mass parameter                             [-]
%  flagSTM      [1]    flag used to activate event function       [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  xf          [6x1]   final state                               [-]
%  Phif        [6x6]   STM at final time                         [-]
%  tf          [1]     final time                                [-]
%  xx          [nx6]   state from t0 to tf                       [-]
%  tt          [nx1]   time from t0 to tf                        [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  15/10/2024: First version
%
%-------------------------------------------------------------------------%

% Compute integration time
tof = tf - t0;

% Set options
options = odeset('reltol', 1e-12, 'abstol', 1e-12); 

% Perform integration
[tt, xx] = ode78(@(t, x) CRTBP_dynamic_STM(t, x, mu, flagSTM), [0 tof], xx0, options);

% Extract state vector 
xxf = xx(end,1:6)';

if flagSTM

    % Extract STM 
    Phif = reshape(xx(end,7:end),6,6); 
else
    Phif = NaN(6, 6);
end

% Extract the final time
tf = tt(end);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: CRTBP_dynamic_STM.m 
function [dxdt] = CRTBP_dynamic_STM(~, xx, mu, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% CRTBP_dynamic_STM: ode of Circular Restricted 3 Body Problem, with dynamics of STM                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [dxdt] = CRTBP_dynamic_STM(~, xx, mu, true)
%   [dxdt] = CRTBP_dynamic_STM(~, xx, mu, false)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx           [42x1]   state vector                            [-]
%  mu           [1]   mass parameter                             [-]  
%  flagSTM      [1]    flag used to activate event function      [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  dxdt          [42x1]     ODE                                  [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  15/10/2024: First version
%
%-------------------------------------------------------------------------%

% Get components of position  and velocity of P
x = xx(1);  % x - position
y = xx(2);  % y - position
z = xx(3);  % z - position
vx = xx(4); % x - velocity
vy = xx(5); % y - velocity

% Compute position vector of P w.r.t m1
rr1 = [x + mu; y; z]; % vector
r1 = norm(rr1);       % norm 

% Compute position vector of P w.r.t m2
rr2 = [x + mu - 1; y; z]; % vector
r2 = norm(rr2);           % norm

% Check if STM is requested
if flagSTM

    % Initialize dxdt equal to zero
    dxdt = zeros(42, 1);

    % Put PHI in matrix form
    Phi = reshape(xx(7:end),6,6); 
    
    % df/dv
    dfdv = [(mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1,...  
            (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)),...                                                 
            (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2));... 
            (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)),... 
            (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1,...  
            (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2);...
            (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)),...
            (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),...
            (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)];
    
    % Assemble the A matrix = df/dxx
    A = [zeros(3), eye(3); dfdv, [0 2 0; -2 0 0; 0 0 0]];
    
    % Compute the derivative of the STM
    Phidot = A * Phi ;

    % Define dinamyc of STM
    dxdt(7:end) = Phidot(:);
else

    % Initialize dxdt equal to zero
    dxdt = zeros(6, 1);
end

% Define dxdt for dynamics of CR3BP
dxdt(1 : 3) = xx(4 : 6);                                                            % Velocity
dxdt(4 : 6) = 2 * [vy; -vx; 0] + [x; y; 0] - (1 - mu)/r1^3 * rr1 - mu/r2^3 * rr2;   % Acceleration

end
%-------------------------------------------------------------------------%
