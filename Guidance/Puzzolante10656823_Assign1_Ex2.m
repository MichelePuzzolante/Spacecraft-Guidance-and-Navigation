% Guidance and Navigation (2024/2025)
% Assignment # 1
% Exercise 2
% Author : Michele Puzzolante
%% REQUEST 1
% Initializzation of the workspace
clear;
close all;
clc;

% Clear the kernel pool
cspice_kclear();

SAVEPLOTS = false; % Flag to use when plots has to be saved as png file

% Load kernels
cspice_furnsh('ex02.tm');

% Set default options for plots
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');

fprintf('\n EXERCISE 2')
fprintf('\n\n REQUEST 1\n')

% Initializzation of the struct PBRFBP_const used in function to pass the
% constant values
PBRFBP_const.ms      = 3.28900541 * 10^5;   % Scaled masso of the Sun            [-]
PBRFBP_const.rho     = 3.88811143 * 10^2;   % Scaled Sun-(Earth+Moon) distance   [-]
PBRFBP_const.omegaS  = -9.25195985 * 10^-1; % Scaled angular velocity of the Sun [-]
PBRFBP_const.omegaEM = 2.66186135 * 10^-1;  % Earth-Moon angular velocity        [1/s]
PBRFBP_const.lEM     = 3.84405 * 10^8;      % Earth-Moon distance                [m]
PBRFBP_const.hi      = 167;                 % Altitude of departure orbit        [km]
PBRFBP_const.hf      = 100;                 % Altitude of arrival orbit          [km]
PBRFBP_const.DU      = 3.84405000 * 10^5;   % Distance Unit                      [km]
PBRFBP_const.TU      = 4.34256461;          % Time Unit                          [days]
PBRFBP_const.VU      = 1.02454018;          % Velocity Unit                      [km/s]

% Extract the radius of the Earth from the kernel
EARTH_RADII = cspice_bodvrd('EARTH', 'RADII', 3);
PBRFBP_const.RE = EARTH_RADII(1);

% Extract the radius of the Moon from the kernel
MOON_RADII = cspice_bodvrd('MOON', 'RADII', 3);
PBRFBP_const.RM = MOON_RADII(1);

% Computation of the mass parameters using the Gravitational parameter of
% the Earth and the Moon extracted from the kernel
GM_E = cspice_bodvrd('EARTH', 'GM', 1); % Earth gravitational parameter
GM_M = cspice_bodvrd('MOON', 'GM', 1);  % Moon gravitational parameter
mu_EM = GM_M/(GM_E + GM_M);             % mass parameter for the PBRFBP
PBRFBP_const.mu = mu_EM;                % save mass parameter for the PBRFBP in the struct

% Initial parameters to generate a first guess solution
alpha = 0.2 * pi;
beta  = 1.41;
delta = 4;
ti    = 2;

% Computation of the scaled initial and final position and scaled initial
% velocity
ri = (PBRFBP_const.RE + PBRFBP_const.hi)/PBRFBP_const.DU;              % scaled initial position of the initial orbit
v0 = beta * sqrt((1 - mu_EM)/ri);                                      % scaled initial velocity
PBRFBP_const.ri = ri;                                                  % save the initial position in the struct
PBRFBP_const.rf = (PBRFBP_const.RM + PBRFBP_const.hf)/PBRFBP_const.DU; % save the final position in the struct

% Computation of the first guess optimal solution
x0 = ri * cos(alpha) - mu_EM;   % x - position
y0 = ri * sin(alpha);           % y - position
vx0 = -(v0 - ri) * sin(alpha);  % x - velocity
vy0 = (v0 - ri) * cos(alpha);   % y - velocity

% Assemble the first guess
xx0 = [x0; y0; vx0; vy0]; 

% Rotation of the first guess from rotating frame to Earth Inertial Frame
[xx0_ECI] = RotFrame2InFrame(xx0', ti, mu_EM);
% Conversion of the ECI state into dimensional quantities
xx0_ECI = [xx0_ECI(1:2) * PBRFBP_const.DU;...
           xx0_ECI(3:4) * PBRFBP_const.VU];

% Display the outputs
fprintf('\n---------------INITIAL STATE(Rotating Frame)--------------\n')
fprintf('xi\t=\t%.6f\t[-]\n',  x0)
fprintf('yi\t=\t%.6f\t[-]\n',  y0)
fprintf('vxi\t=\t%.6f\t[-]\n', vx0)
fprintf('vyi\t=\t%.6f\t[-]\n', vy0)

fprintf('\n---------------INITIAL STATE(ECI)-------------------------\n')
fprintf('Xi\t= %.6f\t[km]\n', xx0_ECI(1))
fprintf('Yi\t= %.6f\t\t[km]\n', xx0_ECI(2))
fprintf('Vxi\t= %.6f\t[km]\n', xx0_ECI(3))
fprintf('Vyi\t= %.6f\t\t[km]\n', xx0_ECI(4))

fprintf('\nmu\t= %.10f[-]', mu_EM)

% Propagation of the optimal state found with fmincon
[~, ~, ~, xx_xx0, tt_xx0] = PBRFBP_propagation(xx0, [ti, ti+delta], PBRFBP_const, false);
% Rotation of all the states propagated from rotatimg frame to ECI
[xxIF_xx0] = RotFrame2InFrame(xx_xx0, tt_xx0, mu_EM);
% Conversion of the ECI state into dimensional quantities
xxIF_xx0 = [xxIF_xx0(:, 1:2) * PBRFBP_const.DU,...
            xxIF_xx0(:, 3:4) * PBRFBP_const.VU];

% Plot the transfer in the rotating adimensional reference frame
figure(1)
hold on
plot(-mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',20, 'LineWidth', 1, 'Color', [0 0.4470 0.7410])
plot(1 - mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize', 20, 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840])
plot(xx_xx0(:,1), xx_xx0(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
xlabel('x [-]', 'FontSize', 15); 
ylabel('y [-]', 'FontSize', 15);
legend('Departure', 'Arrival','Trajectory','FontSize', 12, 'Position', [0.15,0.13,0.22,0.15])

% Plot the transfer in the dimensional ECI reference frame
figure(2)
plot(xxIF_xx0(1, 1), xxIF_xx0(1, 2), 'LineStyle','none','Marker','.','Color',[0 0.4470 0.7410], 'MarkerSize',20, 'LineWidth',1)
hold on
plot(xxIF_xx0(end, 1), xxIF_xx0(end, 2), 'LineStyle','none','Marker','.','Color',[0.6350 0.0780 0.1840], 'MarkerSize',20,  'LineWidth',1)
plot(PBRFBP_const.lEM * 1e-3 * sin(0 : 0.01 : 2*pi), PBRFBP_const.lEM * 1e-3 * cos(0 : 0.01 : 2*pi),'LineWidth', 2, 'LineStyle','--', 'Color', [0.9290 0.6940 0.1250])
plot(xxIF_xx0(:,1), xxIF_xx0(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
xlabel('x [km]', 'FontSize', 15); 
ylabel('y [km]', 'FontSize', 15);
legend('Departure', 'Arrival', 'Moon''s orbit','Trajectory', 'FontSize', 12, 'Position', [0.64,0.71, 0.25,0.20])

% Saving plots
if SAVEPLOTS

    % save 1st figure
    figure(1)
    exportgraphics(gcf,'transfer_x0_ADIM.png','Resolution',600);

    % Save 2nd figure
    figure(2)
    exportgraphics(gcf,'transfer_x0_DIM.png','Resolution',600);
end

%% REQUEST 2.1
% The section REQUEST 2.1 and REQUEST 2.2 contain the optimizzation for the
% simple shooting method

fprintf('\n REQUEST 2.1 \n')

% Creation of the lower and upper bound for the variables of the
% Optimizzation problem
% Lower bound vector
lb = -inf(6, 1);
% Upper bound vector
ub = +inf(6,1);

% Inequality linear constraint
A_LinIn = [0 0 0 0 1 -1];
b_LinIn = 0;

% Initializzation of the flag for activate or deactivate the gradients computed
% analitically
flag_G = false;

% Assemble the initial guess for the optimizzation problem
y0 = [xx0; ti; ti+delta];

% Set options for fmincon
options_R_21 = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000,...
                            'ConstraintTolerance', 1e-11, 'StepTolerance', 1e-20, 'FunctionTolerance', 1e-20,...
                            'SpecifyObjectiveGradient', flag_G,'SpecifyConstraintGradient',flag_G,...
                            'Algorithm', 'interior-point', 'SubproblemAlgorithm','ldl-factorization',...
                            'HessianApproximation', 'lbfgs', 'EnableFeasibilityMode', true, 'FiniteDifferenceType','forward');

% Start timer
CPU_time_start_SimpleShooting_noG = tic;

% Call fmincon
[yOPT_noG, f_noG, exitflag_noG, output_noG] = fmincon(@(y) ObjectiveFun_SimpleShoot(y , PBRFBP_const, flag_G), y0, A_LinIn, b_LinIn, [], [], lb, ub, @(y) NonLinConstraints_SimpleShoot(y, PBRFBP_const,flag_G), options_R_21);

% Get computational time
CPU_time_fmincon_SingleShooting_noG = toc(CPU_time_start_SimpleShooting_noG);

% Computation of the Constrints value for the optimal solution
[~, error_Ceq_SS_noG] = NonLinConstraints_SimpleShoot(yOPT_noG, PBRFBP_const, false);

% Display the outputs
fprintf('\n---OPTIMAL INITIAL STATE SIMPLE SHOOTING without GRADIENTS (Rotating Frame)---\n')
fprintf('xi\t=\t%.6f\t[-]\n',  yOPT_noG(1))
fprintf('yi\t=\t%.6f\t[-]\n',  yOPT_noG(2))
fprintf('vxi\t=\t%.6f\t[-]\n',  yOPT_noG(3))
fprintf('vyi\t=\t%.6f\t[-]\n',  yOPT_noG(4))
fprintf('ti\t=\t%.3f\t[-]\n',  yOPT_noG(5))
fprintf('tf\t=\t%.3f\t[-]\n',  yOPT_noG(6))
fprintf('dV\t=\t%.6f\t[km/s]\n',  f_noG * PBRFBP_const.VU)
fprintf('dt\t=\t%.3f\t[days]\n',  (yOPT_noG(6) - yOPT_noG(5)) * PBRFBP_const.TU)

fprintf('||psi_i||\t=\t%.4e\n', norm(error_Ceq_SS_noG(1:2)))
fprintf('||psi_f||\t=\t%.4e\n', norm(error_Ceq_SS_noG(3:4)))

% Propagation of the optimal state found with fmincon
[~, ~, ~, xx_noG, tt_noG] = PBRFBP_propagation(yOPT_noG(1:4), [yOPT_noG(5) yOPT_noG(6)], PBRFBP_const, false);
% Rotation of all the states propagated from rotatimg frame to ECI
[xxIF_noG] = RotFrame2InFrame(xx_noG, tt_noG, mu_EM);
% Conversion of the ECI state into dimensional quantities
xxIF_noG = [xxIF_noG(:, 1:2) * PBRFBP_const.DU,...
            xxIF_noG(:, 3:4) * PBRFBP_const.VU];

% Plot the transfer in the rotating adimensional reference frame
figure(3)
hold on
plot(-mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',20, 'LineWidth', 1, 'Color', [0 0.4470 0.7410])
plot(1 - mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize', 20, 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840])
plot(xx_noG(:,1), xx_noG(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
xlabel('x [-]', 'FontSize', 15); 
ylabel('y [-]', 'FontSize', 15);
legend('Departure', 'Arrival','Trajectory','FontSize', 12, 'Position', [0.72,0.8, 0.22,0.15])

% Plot the transfer in the dimensional ECI reference frame
figure(4)
hold on
plot(xxIF_noG(1, 1), xxIF_noG(1, 2), 'LineStyle','none','Marker','.','Color',[0 0.4470 0.7410], 'MarkerSize',20, 'LineWidth',1)
plot(xxIF_noG(end, 1), xxIF_noG(end, 2), 'LineStyle','none','Marker','.','Color',[0.6350 0.0780 0.1840], 'MarkerSize',20,  'LineWidth',1)
plot(PBRFBP_const.lEM * 1e-3 * sin(0 : 0.01 : 2*pi), PBRFBP_const.lEM * 1e-3 * cos(0 : 0.01 : 2*pi),'LineWidth', 2, 'LineStyle','--', 'Color', [0.9290 0.6940 0.1250])
plot(xxIF_noG(:,1), xxIF_noG(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
xlabel('x [km]', 'FontSize', 15); 
ylabel('y [km]', 'FontSize', 15);
legend('Departure', 'Arrival', 'Moon orbit','Trajectory', 'FontSize', 12, 'Position', [0.72,0.8, 0.22,0.15])

% Saving plots
if SAVEPLOTS

    % save 3rd figure
    figure(3)
    exportgraphics(gcf,'transfer_EM_SimpleShooting_No_Gradient_ADIM.png','Resolution',600);

    % Save 4th figure
    figure(4)
    exportgraphics(gcf,'transfer_EM_SimpleShooting_No_Gradient_DIM.png','Resolution',600);
end

%% REQUEST 2.2
fprintf('\n REQUEST 2.2 \n')

% Initializzation of the flag for activate or deactivate the gradients computed
% analitically
flag_G = true;

% Set the tolerance for checking of the gradients
tol_G = 1e-5;

% Set the options for checking of the gradients
optchkG_SimpleShoot = optimoptions('fmincon', 'FiniteDifferenceType','central', 'FiniteDifferenceStepSize', 1e-8);

% Check the gradient of the Objective function for the simple shooting
[chkObjFun_grad_SimpleShoot, errObjF_Grad_SimpleShoot] = checkGradients(@(y) ObjectiveFun_SimpleShoot(y, PBRFBP_const, true), y0, optchkG_SimpleShoot, 'Display','on', 'Tolerance', tol_G);

% Check the gradient of the constraints for the simple shooting
[chkNonLinConstr_grad_SimpleShoot, errNonLinCons_Jac_SimpleShoot] = checkGradients(@(y) NonLinConstraints_SimpleShoot(y, PBRFBP_const, true), y0, optchkG_SimpleShoot, 'IsConstraint',true, 'Display', 'on', 'Tolerance', tol_G);

% Set options for fmincon
options_R_22 = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000,...
                            'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-20, 'FunctionTolerance', 1e-20,...
                            'SpecifyObjectiveGradient', flag_G,'SpecifyConstraintGradient',flag_G,...
                            'Algorithm', 'interior-point', 'SubproblemAlgorithm','ldl-factorization',...
                            'HessianApproximation', 'lbfgs', 'EnableFeasibilityMode',true);

% Start timer
CPU_time_start_SimpleShooting_G = tic;

% Call fmincon
[yOPT_G, f_G, exitflag_G, output_G] = fmincon(@(y) ObjectiveFun_SimpleShoot(y, PBRFBP_const, flag_G), y0, A_LinIn, b_LinIn, [], [], lb, ub, @(y) NonLinConstraints_SimpleShoot(y, PBRFBP_const, flag_G), options_R_22);

% Get computational time
CPU_time_fmincon_SingleShooting_G = toc(CPU_time_start_SimpleShooting_G);

% Computation of the Constrints value for the optimal solution
[~, error_Ceq_SS_G] = NonLinConstraints_SimpleShoot(yOPT_G, PBRFBP_const, false);

% Display the outputs
fprintf('\n---OPTIMAL INITIAL STATE SIMPLE SHOOTING with GRADIENTS (Rotating Frame)---\n')
fprintf('xi\t=\t%.6f\t[-]\n',  yOPT_G(1))
fprintf('yi\t=\t%.6f\t[-]\n',  yOPT_G(2))
fprintf('vxi\t=\t%.6f\t[-]\n',  yOPT_G(3))
fprintf('vyi\t=\t%.6f\t[-]\n',  yOPT_G(4))
fprintf('ti\t=\t%.3f\t[-]\n',  yOPT_G(5))
fprintf('tf\t=\t%.3f\t[-]\n',  yOPT_G(6))
fprintf('dV\t=\t%.6f\t[km/s]\n',  f_G * PBRFBP_const.VU)
fprintf('dt\t=\t%.3f\t[days]\n',  (yOPT_G(6) - yOPT_G(5)) * PBRFBP_const.TU)

fprintf('||psi_i||\t=\t%.4e\n', norm(error_Ceq_SS_G(1:2)))
fprintf('||psi_f||\t=\t%.4e\n', norm(error_Ceq_SS_G(3:4)))

% Propagation of the optimal state found with fmincon
[~, ~, ~, xx_G, tt_G] = PBRFBP_propagation(yOPT_G(1:4), [yOPT_G(5) yOPT_G(6)], PBRFBP_const, false);

% Rotation of all the states propagated from rotatimg frame to ECI
[xxIF_G] = RotFrame2InFrame(xx_G, tt_G, mu_EM);

% Conversion of the ECI state into dimensional quantities
xxIF_G = [xxIF_G(:, 1:2) * PBRFBP_const.DU,...
          xxIF_G(:, 3:4) * PBRFBP_const.VU];

% Plot the transfer in the rotating adimensional reference frame
figure(5)
hold on
plot(-mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',20, 'LineWidth', 1, 'Color', [0 0.4470 0.7410])
plot(1 - mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize', 20, 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840])
plot(xx_G(:,1), xx_G(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
xlabel('x [-]', 'FontSize', 15); 
ylabel('y [-]', 'FontSize', 15);
legend('Departure', 'Arrival','Trajectory', 'FontSize', 12, 'Position', [0.72,0.8, 0.22,0.15]);

% Plot the transfer in the dimensional ECI reference frame
figure(6)
hold on
plot(xxIF_G(1, 1), xxIF_G(1, 2), 'LineStyle','none','Marker','.','Color',[0 0.4470 0.7410], 'MarkerSize',20, 'LineWidth',1)
plot(xxIF_G(end, 1), xxIF_G(end, 2), 'LineStyle','none','Marker','.','Color',[0.6350 0.0780 0.1840], 'MarkerSize',20,  'LineWidth',1)
plot(PBRFBP_const.lEM * 1e-3 * sin(0 : 0.01 : 2*pi), PBRFBP_const.lEM * 1e-3 * cos(0 : 0.01 : 2*pi),'LineWidth', 2, 'LineStyle','--', 'Color', [0.9290 0.6940 0.1250])
plot(xxIF_G(:,1), xxIF_G(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
hold off
xlabel('x [km]', 'FontSize', 15); 
ylabel('y [km]', 'FontSize', 15);
legend( 'Departure', 'Arrival', 'Moon orbit','Trajectory','FontSize', 12, 'Position', [0.72,0.8, 0.22,0.15]);

% Saving plots
if SAVEPLOTS

    % Save 5th figure
    figure(5)
    exportgraphics(gcf,'transfer_EM_SimpleShooting_Gradient_ADIM.png','Resolution',600);
    
    % Save 6th figure
    figure(6)
    exportgraphics(gcf,'transfer_EM_SimpleShooting_Gradient_DIM.png','Resolution',600);
end

%% REQUEST 3
% This section contains the optimizzation for the multiple shooting method

fprintf('\n REQUEST 3 \n')

% Initializzation of the flag for activate or deactivate the gradients computed
% analitically
flag_G = true;

% Set the number of points for used
N = 4;
PBRFBP_const.N = N;

% Update the lower and uooer bound taking nto account the number of points
% Lower bound
lb_MultiShoot = -Inf*ones(4*(N) + 2, 1);
% Upper bound
ub_MultiShoot = +Inf*ones(4*(N)+2, 1);

% Computation of the points used for the initial guess for the multiple shooting
[~, ~, ~, xx0_MultiShoot, tt_multishoot] = PBRFBP_propagation(xx0, linspace(ti, ti+delta, N), PBRFBP_const, false);
% Assemble the initial guess
y0_MultiShoot = [reshape(xx0_MultiShoot', 4*N, 1); ti; ti+delta];

% Set the options for checking of the gradients
optchkG_MultiShoot =optimoptions('fmincon', 'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-8);

% Check the gradient of the Objective function for the multiple shooting
[chkObjFun_grad_MultiShoot, errObjF_Grad_MultiShoot] = checkGradients(@(y) ObjectiveFun_MultiShoot(y, PBRFBP_const, true), y0_MultiShoot, optchkG_MultiShoot, 'Display','on', 'Tolerance', tol_G);

% Check the gradient of the constraints for the multiple shooting
[chkNonLinConstr_grad_MultiShoot, errNonLinCons_Jac_MultiShoot] = checkGradients(@(y) NonLinConstraints_MultiShoot(y, PBRFBP_const, true), y0_MultiShoot, optchkG_MultiShoot, 'IsConstraint',true, 'Display', 'on', 'Tolerance', tol_G);

% Set options for fmincon
options_R3 = optimoptions('fmincon','Display','iter-detailed', 'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-20,...
                          'MaxFunctionEvaluations', 10000,...
                          'SpecifyObjectiveGradient', flag_G,'SpecifyConstraintGradient',flag_G,...
                           'HessianApproximation', 'lbfgs', 'SubproblemAlgorithm','ldl-factorization',...
                           'EnableFeasibilityMode',true, 'MaxIterations',5000,'Algorithm','interior-point', 'FunctionTolerance', 1e-20);

% Start timer
CPU_time_start_MultipleShooting_G = tic;

% Call fmincon
[yOPT_G_mult, f_G_multi, exitflag_G_multi, output_G_multi] = fmincon(@(y) ObjectiveFun_MultiShoot(y, PBRFBP_const, flag_G), y0_MultiShoot, [], [], [], [], lb_MultiShoot, ub_MultiShoot, @(y) NonLinConstraints_MultiShoot(y, PBRFBP_const, flag_G), options_R3);

% Get computational time
CPU_time_fmincon_MultipleShooting_noG = toc(CPU_time_start_MultipleShooting_G);

% Computation of the Constrints value for the optimal solution
[~, error_Ceq_MS_G] = NonLinConstraints_MultiShoot(yOPT_G_mult, PBRFBP_const, false);

% Display the outputs
fprintf('\n---OPTIMAL INITIAL STATE MULTIPLE SHOOTING with GRADIENTS (Rotating Frame)---\n')
fprintf('xi\t=\t%.6f\t[-]\n',  yOPT_G_mult(1))
fprintf('yi\t=\t%.6f\t[-]\n',  yOPT_G_mult(2))
fprintf('vxi\t=\t%.6f\t[-]\n',  yOPT_G_mult(3))
fprintf('vyi\t=\t%.6f\t[-]\n',  yOPT_G_mult(4))
fprintf('ti\t=\t%.3f\t[-]\n',  yOPT_G_mult(end-1))
fprintf('tf\t=\t%.3f\t[-]\n',  yOPT_G_mult(end))
fprintf('dV\t=\t%.6f\t[km/s]\n',  f_G_multi * PBRFBP_const.VU)
fprintf('dt\t=\t%.3f\t[days]\n',  (yOPT_G_mult(end) - yOPT_G_mult(end-1)) * PBRFBP_const.TU)

fprintf('||psi_i||\t=\t%.4e\n', norm(error_Ceq_MS_G(end-3:end-2)))
fprintf('||psi_f||\t=\t%.4e\n', norm(error_Ceq_MS_G(end-1:end)))

% Creation of the time vector
tj_vec = linspace(yOPT_G_mult(end-1), yOPT_G_mult(end), N);

% Rotation of the 1sr and Nst state from rotatimg frame to ECI
[xxIF_G_mult_start] = RotFrame2InFrame(yOPT_G_mult(1:4)', yOPT_G_mult(end-1), mu_EM);
[xxIF_G_mult_end] = RotFrame2InFrame(yOPT_G_mult(4*N - 3:4*N)', yOPT_G_mult(end), mu_EM);

% Conversion of the ECI state into dimensional quantities
xxIF_G_mult_start = [xxIF_G_mult_start(:, 1:2) * PBRFBP_const.DU,...
                     xxIF_G_mult_start(:, 3:4) * PBRFBP_const.VU];
xxIF_G_mult_end = [xxIF_G_mult_end(:, 1:2) * PBRFBP_const.DU,...
                   xxIF_G_mult_end(:, 3:4) * PBRFBP_const.VU];

% Plot the departure and the arrival in the rotating adimensional reference frame
figure(7)
plot(-mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize',20, 'LineWidth', 1, 'Color', [0 0.4470 0.7410])
hold on
plot(1 - mu_EM, 0, 'LineStyle','none', 'Marker','.', 'MarkerSize', 20, 'LineWidth', 1, 'Color', [0.6350 0.0780 0.1840])

% Plot the departure and the arrival in the dimensional reference frame
figure(8)
plot(xxIF_G_mult_start(1, 1), xxIF_G_mult_start(1, 2), 'LineStyle','none','Marker','.','Color',[0 0.4470 0.7410], 'MarkerSize',20, 'LineWidth',1)
hold on
plot(xxIF_G_mult_end(end, 1), xxIF_G_mult_end(end, 2), 'LineStyle','none','Marker','.','Color',[0.6350 0.0780 0.1840], 'MarkerSize',20,  'LineWidth',1)
plot(PBRFBP_const.lEM * 1e-3 * sin(0 : 0.01 : 2*pi), PBRFBP_const.lEM * 1e-3 * cos(0 : 0.01 : 2*pi),'LineWidth', 2, 'LineStyle','--', 'Color', [0.9290 0.6940 0.1250])

% Plots
for j = 1 : N - 1
    
    % Extraction of the j-th state from the solution of fmincon
    xxj = yOPT_G_mult(4*j-3 : 4*j);

    % Proparagation od the j-th state
    [~, ~, ~, xxj_prop, ttj_prop] = PBRFBP_propagation(xxj, [tj_vec(j), tj_vec(j+1)], PBRFBP_const, false);

    % Rotation of the 1sr and Nst state from rotatimg frame to ECI
    [xxIF_G_mult_j] = RotFrame2InFrame(xxj_prop, ttj_prop, mu_EM);

    % Conversion of the ECI state into dimensional quantities
    xxIF_G_mult_j = [xxIF_G_mult_j(:, 1:2) * PBRFBP_const.DU,...
                     xxIF_G_mult_j(:, 3:4) * PBRFBP_const.VU];

    % Set hadle visibility to on only for the 1st mesh point
    if j == 1
        HandleVis = 'on';
    else
        HandleVis = 'off';
    end
    
    % Plot the mesh point and the arc's trajectory in the adimensional
    % reference frame
    figure(7)
    hold on
    plot(xxj_prop(1, 1), xxj_prop(1, 2), 'LineStyle', 'none', 'Marker','*', 'MarkerSize',7, 'LineWidth',1, 'Color','r', 'HandleVisibility',HandleVis);
    plot(xxj_prop(:, 1), xxj_prop(:, 2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'HandleVisibility',HandleVis);
    
    % Plot the mesh point and the arc's trajectory in the dimensional
    % reference frame
    figure(8)
    hold on
    plot(xxIF_G_mult_j(1,1), xxIF_G_mult_j(1,2), 'LineStyle', 'none', 'Marker','*', 'MarkerSize',7, 'LineWidth',1, 'Color', 'r', 'HandleVisibility',HandleVis);
    plot(xxIF_G_mult_j(:,1), xxIF_G_mult_j(:,2), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560], 'HandleVisibility',HandleVis)
end

% Rotation of the 1sr and Nst state from rotatimg frame to ECI
[xxIF_G_mult_N] = RotFrame2InFrame(yOPT_G_mult(4*N-3 : 4*N)', yOPT_G_mult(end), mu_EM);

% Conversion of the ECI state into dimensional quantities
xxIF_G_mult_N = [xxIF_G_mult_N(:, 1:2) * PBRFBP_const.DU,...
                 xxIF_G_mult_N(:, 3:4) * PBRFBP_const.VU];

% Plot the transfer in the rotating adimensional reference frame
figure(7)
plot(yOPT_G_mult(4*N-3),yOPT_G_mult(4*N-2), 'LineStyle', 'none', 'Marker','*', 'MarkerSize',7, 'LineWidth',1, 'Color','r', 'HandleVisibility','off');
hold off
xlabel('x [-]', 'FontSize', 15); 
ylabel('y [-]', 'FontSize', 15);
legend('Departure', 'Arrival','Mesh point', 'Arc trajectory','FontSize', 12, 'Position', [0.72,0.8, 0.22,0.15])

% Plot the transfer in the dimensional ECI reference frame
figure(8)
plot(xxIF_G_mult_N(1), xxIF_G_mult_N(2), 'LineStyle', 'none', 'Marker','*', 'MarkerSize',7, 'LineWidth',1, 'Color', 'r', 'HandleVisibility','off');
hold off
xlabel('x [km]', 'FontSize', 15); 
ylabel('y [km]', 'FontSize', 15);
legend('Departure', 'Arrival','Moon orbit','Mesh point', 'Arc trajectory', 'FontSize', 12, 'Position', [0.72,0.8, 0.22,0.15])

% Saving plots
if SAVEPLOTS

    % Seve 7th plot
    figure(7)
    exportgraphics(gcf,'transfer_EM_MultipleShooting_ADIM.png','Resolution',600);

    % Save 8th plot
    figure(8)
    exportgraphics(gcf,'transfer_EM_MultipleShooting_DIM.png','Resolution',600);
end

%% REQUEST 4
% This section is used to determine an epoch in which the angle between
% MOON-EMB-SUN is equal tot theta_i
fprintf('\n REQUEST 4 \n')

% Rotation of all the states propagated from rotatimg frame to ECI
xxOPT_G_ECI = RotFrame2InFrame(yOPT_G(1:4)', yOPT_G(5), mu_EM);

% Conversion of the ECI state into dimensional quantities
xx0_ECI_OPT = [[xxOPT_G_ECI(1:2)';0] * PBRFBP_const.DU;...
               [xxOPT_G_ECI(3:4)';0] * PBRFBP_const.VU];

% Computation of the initial angle M-EMB-S
theta_i =  PBRFBP_const.omegaS * yOPT_G(5);
% Check if that_i is <0 to wrap it between [0, 2PI]
if yOPT_G(5) > 0            
    theta_i = 2*pi + theta_i;
end

% Defnition of the frame usdd to express the orientation of the rotating
% frame
frame = 'J2000';

% Initializzation of the starting epoch
et_0 = cspice_str2et('2024-Sep-28 00:00:00.000 TDB');

% Set options for fsolve
options_FSOLVE_R1 = optimoptions('fsolve', 'FunctionTolerance', eps , 'Algorithm','levenberg-marquardt', 'Display','iter', 'StepTolerance',1e-20);

% Call fzero
[ti_et, f_theta] = fsolve(@(et) theta_et(et, frame) - theta_i, cspice_str2et('2024-Oct-10 00:00:00.000 TDB'), options_FSOLVE_R1);

% Number of points used for plotting the evolution of theta w.r.t the time
nP = 500;                       

% Vectort fot thetat values for improving numerical computation
theta_vec = zeros(nP, 1);      

% vector of ephemeries time from et_0 up to 40 days after
et_vec = linspace(et_0, et_0 + cspice_convrt(40, 'DAYS', 'SECONDS'), nP); 

% Computation of theta values
for j = 1 : nP

    theta_vec(j) = theta_et(et_vec(j), frame);
end

% Computation of the final time in ephemeries time
tf_et = ti_et + cspice_convrt((yOPT_G(6) - yOPT_G(5)) * PBRFBP_const.TU, 'DAYS', 'SECONDS');

% Definition of the bodies for the N_body propagator
labels_Body = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initializzation of the bodies gravitational constants
Bodies = nbody_init(labels_Body);

% Definition of the center of integration
center = 'Earth';

% Perform propagation
[xx_ECI, tt_ECI] = nBody_propagation(xx0_ECI_OPT, [ti_et tf_et], Bodies, frame, center);

% Extract Moon's position from Spice
MOONrv = cspice_spkezr('MOON',linspace(ti_et, tf_et, 1000), frame, 'none', center);

% Computation of the errors
error_Position_Real_Trajectory = norm(xx_ECI(end, 1:3)' - MOONrv(1:3, end));
error_Velocity_Real_Trajectory = norm(xx_ECI(end, 4:6)' - MOONrv(4:6, end));

% Dispaly the results
fprintf('The date found with a congiguration of theta_i = %f[deg]\nis the %s\n',  theta_i/cspice_rpd, cspice_et2utc(ti_et, 'C', 3));

fprintf('\nDeparture date\t:\t%s',cspice_et2utc(ti_et, 'C', 3))
fprintf('\nArrival date\t:\t%s\n',cspice_et2utc(tf_et, 'C', 3))

fprintf('\n---------------INITIAL STATE(ECI) 3D---------------------\n')
fprintf('Xi\t= %.8f\t[km]\n', xx0_ECI_OPT(1))
fprintf('Yi\t= %.8f\t\t[km]\n', xx0_ECI_OPT(2))
fprintf('Zi\t= %.8f\t\t[km]\n', xx0_ECI_OPT(3))
fprintf('Vxi\t= %.8f\t\t[km/s]\n', xx0_ECI_OPT(4))
fprintf('Vyi\t= %.8f\t\t[km/s]\n', xx0_ECI_OPT(5))
fprintf('Vzi\t= %.8f\t\t[km/s]\n', xx0_ECI_OPT(6))

fprintf('||r(tf) - rM||\t=\t%.4e\n', norm(error_Position_Real_Trajectory))
fprintf('||v(tf) - vM||\t=\t%.4e\n', norm(error_Velocity_Real_Trajectory))

% Plot the evolution of thetat w.r.t the time
figure(9)
plot(et_vec, theta_vec/cspice_rpd, 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'HandleVisibility','off');
hold on
plot(ti_et, theta_et(ti_et, frame)/cspice_rpd, 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', [0.3010 0.7450 0.9330])
xticks(linspace(et_0, et_0 + cspice_convrt(40, 'DAYS', 'SECONDS'), 9))
xlim([et_0 +Inf])
xticklabels( cspice_et2utc( xticks , 'C', 0));
xlabel('Calendar epoch [UTC]', 'FontSize', 15)
ylabel('$\theta(t)$ [deg]', 'FontSize', 15)
legend('Epoch found', 'FontSize',12, 'Position', [0.71,0.9,0.20,0.055])

% Plot the n-body propagation in the ECI reference frame
figure(10)
plot3(0, 0, 0, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 30, 'Color', [0 0.4470 0.7410])
hold on
plot3(MOONrv(1, :), MOONrv(2, :), MOONrv(3, :), 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250], 'LineStyle','--')
plot3(MOONrv(1, end), MOONrv(2, end), MOONrv(3, end), 'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20)
plot3(xx_ECI(:, 1), xx_ECI(:, 2), xx_ECI(:, 3), 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
plot3(xx_ECI(end, 1), xx_ECI(end, 2), xx_ECI(end, 3), 'LineWidth', 2, 'Color', [0.3010 0.7450 0.9330], 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 20)
xlabel('x [km]', 'FontSize', 15)
ylabel('y [km]', 'FontSize', 15)
zlabel('z [km]', 'FontSize', 15)
legend('Deparure',  'Moon orbit', 'Moon position at $t_f$','Trajectory', 'Arrival', 'Position', [0.658, 0.757, 0.34, 0.24], 'FontSize', 12)

% Saving plots
if SAVEPLOTS

    % save 9th plot
    figure(9)
    exportgraphics(gcf,'Theta_vs_et.png','Resolution',600);

    % Save 10th plot
    figure(10)
    view(0,90)
    exportgraphics(gcf,'transfer_EM_3D.png','Resolution',600);

    % Save 10th plot
    figure(10)
    view(62.3, 8.7)
    legend('Position', [0.14, 0.75, 0.34, 0.24])
    exportgraphics(gcf,'transfer_EM_3D_3dview.png','Resolution',600);
end

%% FUNCTIONS

%-------------------------------------------------------------------------%
% function:nBody_propagation.m 
function [xx, tt] = nBody_propagation(xx0, tspan, Bodies, frame, center)
%-------------------------------------------------------------------------%
%                                                                         
% nBody_propagation: function for propagating state with nbody ode                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xx, tt] = nBody_propagation(xx0, tspan, Bodies, frame, center)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xx0          [6x1] initial state   [km; km; km; km/s; km/s; km/s]
%  tspan        [nx1] time vector                                [s]
%  Bodies       [-]   struct with data relative to the bodies    [-]
%  frame        [1]   frame used to retrive state from Spice     [-]
%  frame        [1]   center used to retrive state from Spice    [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:                          
%  xx          [nx6]   state          [km; km; km; km/s; km/s; km/s]
%  tt          [nx1]   time vector                               [-]
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
Abstol_xx = 1e-14;  % Absolute tolerance for the state

% Set options
options = odeset('reltol', Reltol, 'abstol', Abstol_xx * ones(6, 1));

% Perform integration
[tt, xx] = ode113(@(t, x) nbody_shift_rhs(t, x, Bodies, frame, center), tspan, xx0, options);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% nbody_shift_rhs.m 
function [dxdt] = nbody_shift_rhs(t, x, Bodies, frame, center)
%-------------------------------------------------------------------------%
% nbody_shift_rhs Evaluates the right-hand-side of a N-body propagator
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   The integration centre is the specified by the variable center
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center) 
%
%-------------------------------------------------------------------------%
% INPUTS:
%   t      : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
%   x      : [6,1] cartesian state vector wrt desired object
%   bodies : [1,6] cell-array created with function nbody_init
%
% ------------------------------------------------------------------------%
% OUTPUTS:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
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

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the position between center and object
rr_center_obj = x(1:3);
% Computation of the norm
r_center_obj = norm(rr_center_obj);

% Loop over all bodies (except GM0)
for i = 1:length(Bodies)

    % Check if the center coincides with the i-th bodies
    if strcmp(center, Bodies{i, 1}.name)

        % Computation of acceleration due to the center
        dxdt(4:6) = dxdt(4:6) - Bodies{i, 1}.GM * rr_center_obj/r_center_obj^3;
    else

        % Retrieve position  of i-th celestial body wrt center
        pp_i = cspice_spkpos(Bodies{i, 1}.name, t, frame, 'NONE', center);
        
        % Extract object position wrt. i-th celestial body
        dd_i = rr_center_obj - pp_i;
        % Computation of the norm
        d_i = norm(dd_i);

        % Compute the gravitational q to avoid numerical issues between big
        % and small number
        q_i = dot(rr_center_obj, (rr_center_obj - 2 * pp_i)) / dot(pp_i, pp_i);

        % Compute the gravitational f(q) to avoid numerical issues between big
        % and small number
        fq_i = q_i * (3 + 3 * q_i + q_i^2) / (1 + (1 + q_i)^(3/2));

        % Sum up acceleration to right-hand-side
        dxdt(4:6) = dxdt(4:6) - Bodies{i, 1}.GM * 1/d_i^3 * (rr_center_obj + pp_i * fq_i); 
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: nbody_init.m 
function [bodies] = nbody_init(labels)
%-------------------------------------------------------------------------%
% nbody_init Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [bodies] = nbody_init(labels) 
%
%-------------------------------------------------------------------------%
% INPUTS:
%   labels : [1,n] cell-array with object labels
%
%-------------------------------------------------------------------------%
% OUTPUTS:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  18/10/2024: First version
%
%-------------------------------------------------------------------------%

% Initialize cell array bodies
bodies = cell(size(labels, 1), 1);

% Loop over labels
for i = 1:size(labels, 1)

    % Store body label
    bodies{i, 1}.name = labels(i);

    % Store body gravitational constant
    bodies{i, 1}.GM = cspice_bodvrd(labels(i), 'GM', 1);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: theta_et.m 
function [thetaS] = theta_et(et, frame)
%-------------------------------------------------------------------------%
%                                                                         
% theta_et: Function used for computing the SUN-EMB-MOON angle                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [thetaS] = theta_et(et, frame)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  et           [1]   time                                       [s]
%  const        [1]   reference frame                            [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  thetaS       [1]   Angle                                      [rad]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  24/10/2024: First version
%
%-------------------------------------------------------------------------%

% Get from Spice the state of the MOON w.r.t. EMB in the ref frame (ECI)
rv_EMB_to_MOON = cspice_spkezr('MOON', et, frame ,'none', 'EMB');

% Extract the position from the state
rr_EMB_to_MOON = rv_EMB_to_MOON(1:3);

% Extract the velocity from the state
vv_EMB_to_MOON = rv_EMB_to_MOON(4:6);

% Get the position of the SUN w.r.t. EMB in the ref frame (ECI)
rr_EMB_to_SUN = cspice_spkpos('SUN', et, frame ,'none', 'EMB');

% Compute the normal of the obsculating plane in ECI using the position
% and velocity of the Moon
n_OBS_PLANE_ECI = cross(rr_EMB_to_MOON, vv_EMB_to_MOON)/norm(cross(rr_EMB_to_MOON, vv_EMB_to_MOON));

% Projection of the SUN-EMB vector onto the plane
rr_EMB_to_SUN_projected_OBS = rr_EMB_to_SUN - dot(rr_EMB_to_SUN, n_OBS_PLANE_ECI)/norm(n_OBS_PLANE_ECI)^2 * n_OBS_PLANE_ECI;

% Computation of the cosine of the angle
cos_thetaS = dot(rr_EMB_to_SUN_projected_OBS, rr_EMB_to_MOON) / (norm(rr_EMB_to_SUN_projected_OBS) * norm(rr_EMB_to_MOON));

% Check the sign of the cross product between the EMB-MOON and EMB-SUN
% projected to wrap the angle between 0 and 2pi
sign = cross(rr_EMB_to_MOON, rr_EMB_to_SUN_projected_OBS);

% Computation of the value  of the angle
thetaS = acos(cos_thetaS);

% Wrap if the 3rd component of the dot product is lower than zero
if sign(3) < 0
    
    % change the value of thetaS
    thetaS = 2*pi - thetaS; 
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: ObjectiveFun_MultiShoot.m 
function [F, G] = ObjectiveFun_MultiShoot(y, const, ~)
%-------------------------------------------------------------------------%
%                                                                         
% ObjectiveFun_MultiShoot: Objective function for the multiple shooting method                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [F, G] = ObjectiveFun_MultiShoot(y, const, false)
%   [F, G] = ObjectiveFun_MultiShoot(y, const, true)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y            [(4N+2)x1] ovector of variables                  [-; -; -; -; -; -]
%  const        [-]   struct of constant values                  [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  F            [1]        Objective function                    [-]
%  G            [(4N+2)x1] Objective function's gradient         [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  19/10/2024: First version
%
%-------------------------------------------------------------------------%

% Get the constant values from the struct
mu = const.mu;  % Earth–Moon mass parameter
ri = const.ri;  % scaled radius of the circular orbit around the Earth
rf = const.rf;  % scaled radius of the circular orbit around the Moon
N = const.N;    % number of mesh-points

% Extract values from the variable vector
xxi = y(1:4);         % Initial state
xxN = y(4*N-3:4*N);   % Final state  

% Extract initial position and velocity from the initial state
xi = xxi(1);    % x - position
yi = xxi(2);    % y - position
vxi = xxi(3);   % x - velocity
vyi = xxi(4);   % y - velocity

% Extract final position and velocity from the propagated state
xN = xxN(1);    % x - position
yN = xxN(2);    % y - position
vxN = xxN(3);   % x - velocity
vyN = xxN(4);   % y - velocity

% Computation of the deltaV for the departure
vi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2);
dvi = vi - sqrt((1 - mu)/ri);

% Computation of the deltaV for the arrival
vf = sqrt((vxN - yN)^2 + (vyN + xN + mu - 1)^2);
dvf = vf - sqrt((mu)/rf);

% Definition of the objective function
F = dvi + dvf;

% Checking if the gradient is required
if nargout > 1 
    
    % Initializzation of the gradient of F
    G = zeros(4*N+2, 1);

    % Assemble the gradient
    G(1, 1) = (vyi + xi + mu)/vi;           % dF/dxi
    G(2, 1) = -(vxi - yi)/vi;               % dF/dyi
    G(3, 1) = (vxi - yi)/vi;                % dF/dvxi
    G(4, 1) = (vyi + xi + mu)/vi;           % dF/dvyi
    G(4*N - 3, 1) = (vyN + xN + mu - 1)/vf; % dF/dxN
    G(4*N - 2, 1) = -(vxN - yN)/vf;         % dF/dyN
    G(4*N - 1, 1) = (vxN - yN)/vf;          % dF/dvxN
    G(4*N, 1)     = (vyN + xN + mu - 1)/vf; % dF/dvyN
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: NonLinConstraints_MultiShoot.m 
function [C, Ceq, GC, GCeq] = NonLinConstraints_MultiShoot(y, const, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% NonLinConstraints_MultiShoot: Non linear constraints for the multiple shooting method                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [C, Ceq, GC, GCeq] = NonLinConstraints_MultiShoot(y, const, false)
%   [C, Ceq, GC, GCeq] = NonLinConstraints_MultiShoot(y, const, true)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y            [(4N+2)x1] vector of variables                   [-; -; -; -; -; -]
%  const        [-]        struct of constant values             [-]
%  flagSTM      [1]        flag used to activate STM propagation [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  C            [1]   Non linear inequality constraints          [-]
%  Ceq          [4x1] Non linear equality constraints            [-]
%  GC           [6x1] Non linear inequality constraints' gradient[-]
%  GCeq         [6x4] Non linear equality constraints' gradient  [-]
%
% ------------------------------------------------------------------------%
% AUTHOR:
%  Michele Puzzolante
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  19/10/2024: First version
%
%-------------------------------------------------------------------------%

% Get the constant values from the struct
mu = const.mu;  % Earth–Moon mass parameter
ri = const.ri;  % scaled radius of the circular orbit around the Earth
rf = const.rf;  % scaled radius of the circular orbit around the Moon
RE = const.RE;  % Earth radius
RM = const.RM;  % Moon radius
DU = const.DU;  % distance unit
N = const.N;    % number of mesh-points

% Extract t0 and tf from the variable vector
ti = y(end-1);    % Initial time
tN = y(end);      % Final time

% Extract initial state and final state from the variable vector
xxi = y(1:4);         % initial state
xxN = y(4*N - 3:4*N); % final state

% Creation of the time vector
tt = linspace(ti, tN, N);

% Initializzation for the matrix for the N states
xx = zeros(4, N);

% Extraction of the states from the variable vector
for j = 1 : N
    xx(:, j) = y(4*j-3 : 4*j);
end

% Initializzation for the matrix of the propagated states
xxProp = zeros(4, N-1);

% Initializzation for the STM matrix of the propagated states
PhiProp = zeros(4, 4*(N-1));

% Propagation of the states
for j = 1 : N-1

    [xxProp(:, j), PhiProp(:, 4*j-3 : 4*j), ~, ~, ~] = PBRFBP_propagation(xx(:, j), [tt(j) tt(j+1)], const, flagSTM);
end

% Initializzation for the inequality constraints
C = zeros(2*N+1, 1);

% Inequality constraint for ti and tn
C(end) = ti - tN;

% Inequality constraint
for j = 1 : N

    % Extract the position
    xj = xx(1, j);  % x - position
    yj = xx(2, j);  % y - position

    % Computation of the eta_j constraints
    C(2*j-1, 1) = (RE/DU)^2 - (xj + mu)^2 - yj^2;
    C(2*j, 1) = (RM/DU)^2 - (xj + mu - 1)^2 - yj^2;
end


% Initializzation for the equality constraints
Ceq = zeros(4*N, 1);

% Equality constraints
for j = 1 : N-1
    
    % Computation of the zeta_j constraints
    Ceq(4*j-3 : 4*j) = xxProp(:, j) - xx(:, j+1);
end

% Get initial position and velocity
xi = xxi(1);    % x - position
yi = xxi(2);    % y - position
vxi = xxi(3);   % x - velocity
vyi = xxi(4);   % y - velocity

% Get final position and velocity
xN = xxN(1);    % x - position
yN = xxN(2);    % y - position
vxN = xxN(3);   % x - velocity
vyN = xxN(4);   % y - velocity

% Assemble last equality constraints psi(xxi) and psi(xxN)
Ceq(end-3) = (xi + mu)^2 + yi^2 - ri^2;                         % psi(xxi)
Ceq(end-2) = (xi + mu)*(vxi - yi) + yi*(vyi + xi + mu);         % psi(xxi)
Ceq(end-1) = (xN + mu - 1)^2 + yN^2 - rf^2;                     % psi(xxN)
Ceq(end)   = (xN + mu - 1)*(vxN - yN) + yN*(vyN + xN + mu - 1); % psi(xxN)

% Checking if the gradient is required
if nargout > 2
    
    % Initializzation the gradient of the inequality constraints
    GC = zeros(4*N+2, 2*N+1);

    % Assemble the gradient of the inequality constraints    
    for j = 1 : N
        
        % Extract position
        xj = xx(1, j);  % x - position
        yj = xx(2, j);  % y - position

        % Computation of the deta_j/dxxj constraints
        GC(4*j-3 : 4*j,2*j-1, 1) = [-2*(xj + mu); -2*yj; 0; 0];
        GC(4*j-3 : 4*j,2*j, 1) = [-2*(xj + mu - 1); -2*yj; 0; 0];
    end

    % Derivative of the last constraints
    GC(end-1:end, end) = [1; -1];
    
    % Initializzation the gradient of the inequality constraints
    GCeq = zeros(4*N+2, 4*N);
    
    % Assemble the gradient of the equality constraints 
    for j = 1 : N-1
        
        % Extract the STM of the xxj state
        Phif = PhiProp(:, 4*j-3 : 4*j);

        % Computation of the f(xxj)
        f_xxj = PBRFBP_dynamics_STM(tt(j), xx(:, j), const, false);
        
        % Computation of the f(xxj+1)
        f_xxj_next = PBRFBP_dynamics_STM(tt(j+1), xxProp(:, j), const, false);
        
        GCeq(4*j-3 : 4*j, 4*j-3 : 4*j) = Phif';         % d(flow(xx_j,tj,tj+1) - xx_j+1)/dxx_j
        GCeq(4*j+1 : 4*j+4, 4*j-3 : 4*j) = -eye(4);     % d(flow(xx_j,tj,tj+1) - xx_j+1)/dxx_j+1

        GCeq(end-1, 4*j-3 : 4*j) = (-Phif*f_xxj*(N - j)/(N-1))' + f_xxj_next'*(N - j - 1)/(N - 1);     % d(flow(xx_j,tj,tj+1) - xx_j+1)/dt1
        GCeq(end, 4*j-3 : 4*j) = (-Phif*f_xxj*(j - 1)/(N-1))' + f_xxj_next'*(j)/(N - 1);               % d(flow(xx_j,tj,tj+1) - xx_j+1)/dtN
    end

    GCeq(1:4, end-3) = [2*(xi + mu); 2*yi; 0; 0];           % dpsi(xxi)/dxxi
    GCeq(1:4, end-2) = [vxi; vyi; (xi + mu); yi];           % dpsi(xxi)/dxxi
    GCeq(4*N-3:4*N, end-1) = [2*(xN + mu - 1); 2*yN; 0; 0]; % dpsi(xxN)/dxxN
    GCeq(4*N-3:4*N, end) = [vxN; vyN; (xN + mu - 1); yN];   % dpsi(xxN)/dxxN
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: ObjectiveFun_SimpleShoot.m 
function [F, G] = ObjectiveFun_SimpleShoot(y, const, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% ObjectiveFun_SimpleShoot: Objective function for the simple shooting method                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [F, G] = ObjectiveFun_SimpleShoot(y, const, false)
%   [F, G] = ObjectiveFun_SimpleShoot(y, const, true)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y            [6x1] ovector of variables                       [-; -; -; -; -; -]
%  const        [-]   struct of constant values                  [-]
%  flagSTM      [1]   flag used to activate STM propagation      [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  F            [1]   Objective function                         [-]
%  G            [6x1] Objective function's gradient              [-]
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
mu = const.mu;  % Earth–Moon mass parameter
ri = const.ri;  % scaled radius of the circular orbit around the Earth
rf = const.rf;  % scaled radius of the circular orbit around the Moon

% Extract values from the variable vector
xx0 = y(1:4);   % Initial state
t0 = y(5);      % Initial time
tf = y(6);      % Final time

% Extract initial position and velocity from the initial state
xi = xx0(1);    % x - position
yi = xx0(2);    % y - position
vxi = xx0(3);   % x - velocity
vyi = xx0(4);   % y - velocity

% Propagation of the initial state
[xxf, Phif, ~, ~, ~] = PBRFBP_propagation(xx0, [t0 tf], const, flagSTM);

% Extract final position and velocity from the propagated state
xf = xxf(1);    % x - position
yf = xxf(2);    % y - position
vxf = xxf(3);   % x - velocity
vyf = xxf(4);   % y - velocity

% Computation of the deltaV for the departure
vi = sqrt((vxi - yi)^2 + (vyi + xi + mu)^2);
dvi = vi - sqrt((1 - mu)/ri);

% Computation of the deltaV for the arrival
vf = sqrt((vxf - yf)^2 + (vyf + xf + mu - 1)^2);
dvf = vf - sqrt((mu)/rf);

% Definition of the objective function
F = dvi + dvf;

% Checking if the gradient is required
if nargout > 1 
    
    % Get f evaluated at xi and ti
    f_xxi = PBRFBP_dynamics_STM(t0, xx0, const, false);
    
    % Get f evaluated at xf and tf
    f_xxf = PBRFBP_dynamics_STM(tf, xxf, const, false);
    
    % Compute the derivative of the flow w.r.t the initial time
    dflowdti = -Phif * f_xxi;

    % initializzation of the gradient
    G = zeros(6,1);

    % Computation od the derivatives
    dvidxxi = [(vyi + xi + mu)/vi; -(vxi - yi)/vi; (vxi - yi)/vi; (vyi + xi + mu)/vi]; 
    dvfdxxf = [((vyf + xf + mu - 1))/vf,...
                       ( -(vxf - yf))/vf,...
                       ( (vxf - yf))/vf,...
                       ( (vyf + xf + mu - 1))/vf];
    dvfdxxi = (dvfdxxf * Phif);

    % Assemble the gradient
    G(1:4) = dvidxxi + dvfdxxi';
    G(5) = dvfdxxf * dflowdti;
    G(6) = dvfdxxf * f_xxf;
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% function: NonLinConstraints_SimpleShoot.m 
function [C, Ceq, GC, GCeq] = NonLinConstraints_SimpleShoot(y, const, flagSTM)
%-------------------------------------------------------------------------%
%                                                                         
% NonLinConstraints_SimpleShoot: Non linear constraints for the simple shooting method                  
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [C, Ceq, GC, GCeq] = NonLinConstraints_SimpleShoot(y, const, false)
%   [C, Ceq, GC, GCeq] = NonLinConstraints_SimpleShoot(y, const, true)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  y            [6x1] vector of variables                        [-; -; -; -; -; -]
%  const        [-]   struct of constant values                  [-]
%  flagSTM      [1]   flag used to activate STM propagation      [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  C            [1]   Non linear inequality constraints          [-]
%  Ceq          [4x1] Non linear equality constraints            [-]
%  GC           [6x1] Non linear inequality constraints' gradient[-]
%  GCeq         [6x4] Non linear equality constraints' gradient  [-]
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
mu = const.mu;  % Earth–Moon mass parameter
ri = const.ri;  % scaled radius of the circular orbit around the Earth
rf = const.rf;  % scaled radius of the circular orbit around the Moon

% Extract values from the variable vector
xx0 = y(1:4);   % Initial state
t0 = y(5);      % Initial time
tf = y(6);      % Final time

% Extract initial position and velocity from the initial state
xi = xx0(1);    % x - position
yi = xx0(2);    % y - position
vxi = xx0(3);   % x - velocity
vyi = xx0(4);   % y - velocity

% Propagation of the initial state
[xxf, Phif, ~, ~, ~] = PBRFBP_propagation(xx0, [t0 tf], const, flagSTM);

% Extract final position and velocity from the propagated state
xf  = xxf(1);    % x - position
yf  = xxf(2);    % y - position
vxf = xxf(3);    % x - velocity
vyf = xxf(4);    % y - velocity

% Inequality constraint
C = [];

% Equality constraints
Ceq = [(xi + mu)^2 + yi^2 - ri^2;...
       (xi + mu)*(vxi - yi) + yi*(vyi + xi + mu);...
       (xf + mu - 1)^2 + yf^2 - rf^2;...
       (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1)];

% Checking if the gradient is required
if nargout > 2
    
    % Assemble the gradient of the inequality constraints
    GC = [];

    % Get f evaluated at xi and ti
    f_xxi = PBRFBP_dynamics_STM(t0, xx0, const, false);
    
    % Get f evaluated at xf and tf
    f_xxf = PBRFBP_dynamics_STM(tf, xxf, const, false);
    
    % Compute the derivative of the flow w.r.t the initial time
    dflowdti = -Phif * f_xxi;
    
    dC3dxxf = [2*(xf + mu - 1), 2*yf, 0, 0];
    dC4dxxf = [vxf, vyf, (xf + mu - 1), yf];
    % Assemble the gradient of the equality constraints
    GCeq = [2*(xi + mu), 2*yi, 0, 0, 0, 0;...                       % dCeq1/dy
            vxi, vyi, (xi + mu), yi, 0, 0;...                       % dCeq2/dy
            dC3dxxf * Phif, dC3dxxf * dflowdti, dC3dxxf * f_xxf;    % dCeq3/dy
            dC4dxxf * Phif, dC4dxxf * dflowdti, dC4dxxf * f_xxf]';  % dCeq4/dy
end

end
%-------------------------------------------------------------------------%

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
% function: RotFrame2InFrame.m 
function [xxIF] = RotFrame2InFrame(xxRF, tt, mu)
%-------------------------------------------------------------------------%
%                                                                         
% RotFrame2InFrame: function which rotates the state from rotating frame to inertial frame                   
%                                                                         
%-------------------------------------------------------------------------%
% PROTOTYPE:
%   [xxIF] = RotFrame2InFrame(xxRF, tt, mu) 
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  xxRF         [Nx4]   state vector in rotating frame             [-, -, -, -]
%  tt           [Nx1]   mass parameter                             [-]              
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  xxIF         [Nx4]   state vector in inertial frame             [-, -, -, -]
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

% Get the number of the states 
Npoints = length(tt);

% Initializzation of the state in the inertial frame
xxIF = zeros(Npoints, 4);

% For cycle used to rotate the states
for i = 1 : Npoints
    
    % Get position and velocity rotating frame
    x = xxRF(i, 1);  % x - position
    y = xxRF(i, 2);  % y - position
    vx = xxRF(i, 3); % x - velocity
    vy = xxRF(i, 4); % y - velocity

    % Get the time
    t = tt(i);

    % Rotate position and velocity
    xxIF(i, 1) = (x + mu) * cos(t) - y * sin(t);             % x - position
    xxIF(i, 2) = (x + mu) * sin(t) + y * cos(t);             % y - position
    xxIF(i, 3) = (vx - y) * cos(t) - (vy + x + mu) * sin(t); % x - velocity
    xxIF(i, 4) = (vx - y) * sin(t) + (vy + x + mu) * cos(t); % y - velocity
end

end
%-------------------------------------------------------------------------%



