% Spacecraft Guidance and Navigation
% Assignment # 2
% Author: Matteo Portantiolo

%% START
clearvars; close all; clc; format long g; cspice_kclear();

% Setting plot options
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter','latex')
set(0, 'defaultAxesFontSize', 28 ,'defaultAxesFontSizeMode', 'manual');

%% EX 1: Uncertainty propagation

% Initial state
ri = [-0.011965533749906; -0.017025663128129];
vi = [10.718855256727338; 0.116502348513671];
xi = [ri;vi];

% Time span
ti = 1.282800225339865;
tf = 9.595124551366348; 
t_vect = linspace(ti,tf,5);

% Covariance
P0 = [1.041e-15 6.026e-17 5.647e-16 4.577e-15;
    6.026e-17 4.287e-18 4.312e-17 1.855e-16;
    5.647e-16 4.312e-17 4.432e-16 1.455e-15;
    4.577e-15 1.855e-16 1.455e-15 2.822e-14];

% Spice
cspice_furnsh('assignment02.tm');
fprintf('Number of LSK  kernels: %d\n', cspice_ktotal('lsk'));
fprintf('Number of SPK  kernels: %d\n', cspice_ktotal('spk'));
fprintf('Number of PCK  kernels: %d\n', cspice_ktotal('pck'));
fprintf('Number of CK   kernels: %d\n', cspice_ktotal('ck'));
fprintf('Number of TEXT kernels: %d\n', cspice_ktotal('TEXT'));
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

% Gravitational parameter
GM_e = cspice_bodvrd('Earth','GM',1);
GM_m = cspice_bodvrd('Moon','GM',1);
mu = GM_m/(GM_e+GM_m);

% Parameters
m_s = 3.28900541e5;
rho = 3.88811143e2;
omega_s = -9.25195985e-1;


%% EX 1.1.1: LinCov

% Propagation
[xf_LC,Pf_LC,~,~,~,~] = propagatePBR4BP(t_vect,xi,P0,mu,rho,omega_s,m_s);

% Plot ellipse
figure
new_colors = [0 0.4470 0.7410; 
              0.8500 0.3250 0.0980; 
              0.9290 0.6940 0.1250; 
              0.4940 0.1840 0.5560;  
              0.4660 0.6740 0.1880; 
              0.3010 0.7450 0.9330; 
              0.6350 0.0780 0.1840;  
              0.75 0.75 0;           
              0.25 0.25 0.25]; 
set(gca, 'ColorOrder', new_colors);
colors = get(gca, 'ColorOrder');
hold on
grid on
xlabel('x  [-]','Interpreter','latex','FontSize',40)
ylabel('y  [-]','Interpreter','latex','FontSize',40)
plotEllipses(xf_LC,Pf_LC,colors(1,:))
text(xf_LC(1),xf_LC(2)-0.1e-3, '$\mu\left(t_f\right)_{LinCov}$', 'Color', colors(1,:), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',40); 

% Compute principal standard deviations (square root of eigenvalues)
Pr_LC = Pf_LC(1:2,1:2);
sigma_princ_LC = sqrt(eig(Pr_LC));
fprintf('\nLC principal standard deviations:\n');
disp(sigma_princ_LC);


%% EX 1.1.2: Unscented Transform (UT)

% State size
n = length(xi);

% UT parameters
alpha = 1;
beta = 2;

% UT scaling parameters
lambda = n * (alpha^2 - 1);
gamma = sqrt(n + lambda);

% Cholesky decomposition (lower triangular)
try
    L = chol(P0, 'lower');
catch
    error('Initial covariance matrix P0 is not positive definite.');
end

% Weights
B = gamma * L; 
Wm = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];
Wc = [lambda / (n + lambda) + (1 - alpha^2 + beta), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];

% Sigma points generation
sigmas = [xi, xi + B, xi - B]; 

% Propagation
sigmas_f = NaN(n,2*n+1);
for ii = 1:2*n+1
    [sig_f,~,~,~,~,~] = propagatePBR4BP(t_vect,sigmas(:,ii),P0,mu,rho,omega_s,m_s);
    sigmas_f(:,ii) = sig_f;
end
    
% Mean
xf_UT = sigmas_f * Wm';

% Covariance
Pf_UT = zeros(n,n);
for ii = 1:2*n+1
    diff = sigmas_f(:,ii) - xf_UT;
    Pf_UT = Pf_UT + Wc(ii) * (diff * diff');
end

% Plot ellipse
hold on
plotEllipses(xf_UT,Pf_UT,colors(2,:))
text(xf_UT(1),xf_UT(2)-0.1e-3, '$\mu\left(t_f\right)_{UT}$', 'Color', colors(2,:), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'Interpreter', 'latex','FontSize',40); 
legend(["$\overline{r}_{LinCov}$","$P_{LinCov}$","$\overline{r}_{UT}$","$P_{UT}$"], 'Interpreter', 'latex','FontSize',40); 

% Compute principal standard deviations (square root of eigenvalues)
Pr_UT = Pf_UT(1:2,1:2);
sigma_princ_UT = sqrt(eig(Pr_UT));
fprintf('\nUT principal standard deviations:\n');
disp(sigma_princ_UT);


%% EX 1.2.1: MonteCarlo (MC)

% Population
n_pop = 1000;

% Noise matrix
R = mvnrnd(xi,P0,n_pop)';

% Propagation
xx_MC = NaN(n,n_pop);
for ii = 1:n_pop
    [xf_MC,~,~,~,~,~] = propagatePBR4BP(t_vect,R(:,ii),P0,mu,rho,omega_s,m_s);
    xx_MC(:,ii) = xf_MC;
end

% Mean and covariance
mean_MC = sum(xx_MC,2) ./ n_pop;
diff_MC = xx_MC - mean_MC;
P_MC = (diff_MC*diff_MC') ./ (n_pop-1);

% Plot ellipse
hold on
plotEllipses(mean_MC,P_MC,colors(3,:))
plot(xx_MC(1, :),xx_MC(2, :), '.', 'Color', 'k', 'MarkerSize', 3);
text(xx_MC(1,end),xx_MC(2,end)-0.1e-3, '$\mu\left(t_f\right)_{MC}$', 'Color', colors(3,:), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex','FontSize',40); 
legend(["$\overline{r}_{LinCov}$","$P_{LinCov}$","$\overline{r}_{UT}$","$P_{UT}$","$\overline{r}_{MC}$", ...
   "$P_{MC}$","MC samples"],'Interpreter','latex','Location','best','FontSize',40,'Orientation', 'vertical');

% Compute principal standard deviations (square root of eigenvalues)
Pr_MC = P_MC(1:2,1:2);
sigma_princ_MC = sqrt(eig(Pr_MC));
fprintf('\nMC principal standard deviations:\n');
disp(sigma_princ_MC);


%% EX 1.2.2:  Time evolution

% Initialization
mean_LC = zeros(n, length(t_vect)); 
P_LC = cell(1, length(t_vect));
mean_UT = zeros(n, length(t_vect)); 
P_UT = cell(1, length(t_vect));
mean_MC = zeros(n, length(t_vect)); 
P_MC = cell(1, length(t_vect));

% Initial conditions
mu_LC = xi; 
sigma_LC = P0;
mu_UT = xi; 
sigma_UT = P0;
mu_MC = xi; 
sigma_MC = P0;

% Initial values
mean_LC(:, 1) = mu_LC; 
P_LC{1} = sigma_LC;
mean_UT(:, 1) = mu_UT; 
P_UT{1} = sigma_UT;
mean_MC(:, 1) = mu_MC; 
P_MC{1} = sigma_MC;

% Propagation loop (continuous propagation)
for k = 1:length(t_vect)-1

    % Updated time step
    tt = [t_vect(k) t_vect(k+1)];

    % 1) LinCov
    [mu_LC,sigma_LC] = LinCov(tt,mu_LC,sigma_LC,mu,rho,omega_s,m_s);
    
    % Save current step
    mean_LC(:, k+1) = mu_LC;
    P_LC{k+1} = sigma_LC;

    % 2) UT
    [mu_UT,sigma_UT] = UT(tt,mu_UT,sigma_UT,mu,rho,omega_s,m_s);
    
    % Save current step
    mean_UT(:, k+1) = mu_UT;
    P_UT{k+1} = sigma_UT;

    % 3) MC
    [mu_MC,sigma_MC] = MC(tt,mu_MC,sigma_MC,mu,rho,omega_s,m_s);
    
    % Save current step
    mean_MC(:, k+1) = mu_MC;
    P_MC{k+1} = sigma_MC;

end

% semi-major axis of 3\sigma ellipsoid (SMA)
sma_LCr = NaN(1,length(t_vect));
sma_LCv = NaN(1,length(t_vect));
sma_UTr = NaN(1,length(t_vect));
sma_UTv = NaN(1,length(t_vect));
sma_MCr = NaN(1,length(t_vect));
sma_MCv = NaN(1,length(t_vect));

for k = 1:length(t_vect)

    % LinCov
    sma_LCr(k) = 3*sqrt(max(eig(P_LC{k}(1:2,1:2))));
    sma_LCv(k) = 3*sqrt(max(eig(P_LC{k}(3:4,3:4))));

    % UT
    sma_UTr(k) = 3*sqrt(max(eig(P_UT{k}(1:2,1:2))));
    sma_UTv(k) = 3*sqrt(max(eig(P_UT{k}(3:4,3:4))));

    % MC
    sma_MCr(k) = 3*sqrt(max(eig(P_MC{k}(1:2,1:2))));
    sma_MCv(k) = 3*sqrt(max(eig(P_MC{k}(3:4,3:4))));

end

% Plot the position covariance matrix
figure
subplot(3,1,1)
semilogy(t_vect, sma_LCr, 'o--', 'Color', colors(1,:),'LineWidth',3);
ylabel('$3\sqrt{\max(\lambda_{i}(P_r))}$ [-]', 'Interpreter', 'latex',FontSize=30);
title('LinCov','Interpreter','latex','FontSize',40)
grid on

subplot(3,1,2)
semilogy(t_vect, sma_UTr, 'o--', 'Color', colors(2,:),'LineWidth',3);
ylabel('$3\sqrt{\max(\lambda_{i}(P_r))}$ [-]', 'Interpreter', 'latex',FontSize=30);
title('UT','Interpreter','latex','FontSize',40)
grid on

subplot(3,1,3)
semilogy(t_vect, sma_MCr, 'o--', 'Color', colors(3,:),'LineWidth',3);
xlabel('Time [-]', 'Interpreter', 'latex',FontSize=40);
ylabel('$3\sqrt{\max(\lambda_{i}(P_r))}$ [-]', 'Interpreter', 'latex',FontSize=30);
title('MC','Interpreter','latex','FontSize',40)
grid on

% Plot the velocity covariance matrix
figure
subplot(3,1,1)
semilogy(t_vect, sma_LCv, 'o--', 'Color', colors(1,:),'LineWidth',3);
ylabel('$3\sqrt{\max(\lambda_{i}(P_v))}$ [-]', 'Interpreter', 'latex',FontSize=30);
title('LinCov','Interpreter','latex','FontSize',40)
grid on

subplot(3,1,2)
semilogy(t_vect, sma_UTv, 'o--', 'Color', colors(2,:),'LineWidth',3);
ylabel('$3\sqrt{\max(\lambda_{i}(P_v))}$ [-]', 'Interpreter', 'latex',FontSize=30);
title('UT','Interpreter','latex','FontSize',40)
grid on

subplot(3,1,3)
semilogy(t_vect, sma_MCv, 'o--', 'Color', colors(3,:),'LineWidth',3);
xlabel('Time [-]', 'Interpreter', 'latex',FontSize=40);
ylabel('$3\sqrt{\max(\lambda_{i}(P_v))}$ [-]', 'Interpreter', 'latex',FontSize=30);
title('MC','Interpreter','latex','FontSize',40)
grid on

% Uncertainty comparison at final time
comparison_unc_r = NaN(3, length(t_vect));
comparison_unc_v = NaN(3, length(t_vect));
for k=1:length(t_vect)
    comparison_unc_r(:,k) = [sma_LCr(k); sma_UTr(k); sma_MCr(k)];
    comparison_unc_v(:,k) = [sma_LCv(k); sma_UTv(k); sma_MCv(k)];
end
fprintf('\nPosition uncertainty comparison:\n');
disp(comparison_unc_r);
fprintf('\nVelocity uncertainty comparison:\n');
disp(comparison_unc_v);


%% EX 1.2.3: QQ plot

% Position and velocity components at the final time step
x = xx_MC(1,:); 
y = xx_MC(2,:);  
vx = xx_MC(3,:);  
vy = xx_MC(4,:);  

% Standardization
x_std = standardize(x);   
y_std = standardize(y);
vx_std = standardize(vx); 
vy_std = standardize(vy); 

% QQ plot for x
figure
h = qqplot(x_std);  
ylabel('Quantiles of $x$', 'Interpreter', 'latex','FontSize',40)
grid on
title('QQ Plot of Standardized Position $x$ vs Standard Normal', 'Interpreter', 'latex','FontSize',40)

% QQ plot for y
figure
qqplot(y_std)  
ylabel('Quantiles of $y$', 'Interpreter', 'latex','FontSize',40)
grid on
title('QQ Plot of Standardized Position $y$ vs Standard Normal', 'Interpreter', 'latex','FontSize',40)

% QQ plot for vx
figure
qqplot(vx_std)  
ylabel('Quantiles of $v_x$', 'Interpreter', 'latex','FontSize',40)
grid on
title('QQ Plot of Standardized Velocity $v_x$ vs Standard Normal', 'Interpreter', 'latex','FontSize',40)

% QQ plot for vy
figure
qqplot(vy_std) 
ylabel('Quantiles of $v_y$', 'Interpreter', 'latex','FontSize',40)
grid on
title('QQ Plot of Standardized Velocity $v_y$ vs Standard Normal', 'Interpreter', 'latex','FontSize',40)


%% FUNCTIONS

function dxdt = xyPBR4BP_STM(t, xx, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the state derivative and State Transition Matrix (STM) 
%   propagation for the Planar Bicircular Restricted Four-Body Problem 
%   (PBR4BP). The system models the dynamics of a massless particle in the 
%   rotating synodic frame, under the influence of two primaries and a 
%   third perturbing body moving in a circular orbit.
%
% Inputs:
%   t       - Time [non-dimensional]
%   xx      - Extended state vector (20×1):
%               [x; y; vx; vy; Phi(:)]
%             where Phi is the 4×4 State Transition Matrix stored column-wise
%   mu      - Mass parameter, defined as m2 / (m1 + m2) [-]
%   rho     - Distance of the perturbing body from the barycenter [nd]
%   omega_s - Angular velocity of the perturbing body [rad/nd_time]
%   m_s     - Mass of the perturbing body (normalized to total system mass) [-]
%
% Outputs:
%   dxdt - Time derivative of the extended system (20×1):
%             [xdot; ydot; vxdot; vydot; Phidot(:)]
%--------------------------------------------------------------------------

% Extract position and velocity components
x  = xx(1);
y  = xx(2);
vx = xx(3);
vy = xx(4);

% Reshape STM from vector form to 4×4 matrix
Phi = reshape(xx(5:end), 4, 4);

% Partial derivatives of the effective potential (4BP)
dOM4dx = x ...
    - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) ...
    - (m_s*cos(omega_s*t))/rho^2 ...
    + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) ...
    - (m_s*(2*x - 2*rho*cos(omega_s*t))) ...
      /(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2));

dOM4dy = y ...
    - (m_s*sin(omega_s*t))/rho^2 ...
    - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) ...
    - (m_s*(2*y - 2*rho*sin(omega_s*t))) ...
      /(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2)) ...
    + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

% Jacobian of the equations of motion (4×4 matrix)
dfdx = [ 0, 0, 1, 0
         0, 0, 0, 1
         (mu - 1)/((mu + x)^2 + y^2)^(3/2) ...
           - mu/((mu + x - 1)^2 + y^2)^(3/2) ...
           - m_s/((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2) ...
           + (3*m_s*(2*x - 2*rho*cos(omega_s*t))^2) ...
             /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) ...
           - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1, ...
           (3*m_s*(2*x - 2*rho*cos(omega_s*t))*(2*y - 2*rho*sin(omega_s*t))) ...
             /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) ...
           - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), ...
           0, 2
         (3*m_s*(2*x - 2*rho*cos(omega_s*t))*(2*y - 2*rho*sin(omega_s*t))) ...
           /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) ...
           - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), ...
           (mu - 1)/((mu + x)^2 + y^2)^(3/2) ...
           - mu/((mu + x - 1)^2 + y^2)^(3/2) ...
           - m_s/((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2) ...
           + (3*m_s*(2*y - 2*rho*sin(omega_s*t))^2) ...
             /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) ...
           + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, ...
           -2, 0 ];

% STM derivative
Phidot = dfdx * Phi;

% Assemble right-hand side of extended system
dxdt = zeros(20,1);
dxdt(1:2)   = xx(3:4);       % Position derivatives
dxdt(3)     = dOM4dx + 2*vy; % Acceleration in x
dxdt(4)     = dOM4dy - 2*vx; % Acceleration in y
dxdt(5:end) = Phidot(:);     % Flattened STM derivative

end

function [xf, Pf, PHIf, tf, xx, tt] = propagatePBR4BP(t_vect, x0, P0, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the state vector and covariance matrix in the Planar
%   Bicircular Restricted Four-Body Problem (PBR4BP) using the State
%   Transition Matrix (STM). The dynamics include Sun perturbations
%   relative to the Earth–Moon rotating frame.
%
% Inputs:
%   t_vect - Time vector for integration [1×N]
%   x0     - Initial state vector [4×1]
%   P0     - Initial covariance matrix [4×4]
%   mu     - Earth–Moon mass ratio
%   rho    - Sun–Earth/Moon distance (adimensional units)
%   omega_s- Sun angular velocity in rotating frame [rad/s]
%   m_s    - Sun mass parameter
%
% Outputs:
%   xf    - Final state vector at t = t_vect(end) [4×1]
%   Pf    - Final covariance matrix at t = t_vect(end) [4×4]
%   PHIf  - Final State Transition Matrix [4×4]
%   tf    - Final propagation time [scalar]
%   xx    - Full propagated trajectory [N×(4+16)]
%             first 4 columns = state
%             remaining 16    = STM reshaped row-wise
%   tt    - Time vector returned by integrator [N×1]
%--------------------------------------------------------------------------

% --- Initialize STM ------------------------------------------------------
Phi0 = eye(4);

% Append STM to initial state
x0Phi0 = [x0; Phi0(:)];

% --- Perform integration -------------------------------------------------
options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode78(@(t,x) xyPBR4BP_STM(t,x,mu,rho,omega_s,m_s), t_vect, x0Phi0, options_STM);

% --- Extract results -----------------------------------------------------
xf   = xx(end,1:4)';             % Final state
PHIf = reshape(xx(end,5:end),4,4); % Final STM
tf   = tt(end);                  % Final time

% Propagate covariance via STM
Pf = PHIf * P0 * PHIf';

end

function [] = plotEllipses(xx_mean, P, color)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Plots an uncertainty ellipse corresponding to the covariance matrix P
%   around the mean state xx_mean. The ellipse is scaled to 3σ confidence 
%   level and drawn in the position (x,y) plane. The function also marks 
%   the mean state with a colored marker.
%
% Inputs:
%   xx_mean - Mean state vector [4×1]
%               xx_mean(1:2) = mean position [x; y]
%               xx_mean(3:4) = mean velocity (unused here)
%   P       - State covariance matrix [4×4]
%               Only the position covariance submatrix P(1:2,1:2) is used
%   color   - Plot color (MATLAB format: char, RGB triplet, or hex code)
%
% Outputs:
%   None (ellipse and mean are plotted in the current figure)
%--------------------------------------------------------------------------

% Parameters for the ellipse plot
num_points = 100;              % Number of points along ellipse
theta      = linspace(0,2*pi,num_points);
n_sigma    = 3;                % Confidence level (3σ)

% Extract the 2×2 position covariance
P_pos = P(1:2,1:2);

% Eigen-decomposition of covariance (rotation + scaling)
[R,D] = svd(P_pos);
eigenvalues = diag(D);

% Semi-axes of the ellipse
a = sqrt(eigenvalues(1));      % Semi-major axis
b = sqrt(eigenvalues(2));      % Semi-minor axis

% Generate ellipse points in principal-axis coordinates
ellipse_points = [a*cos(theta); b*sin(theta)];

% Rotate and scale ellipse to covariance shape
rotated_ellipse = n_sigma * R * ellipse_points;

% Translate ellipse to the mean position
x_mean = xx_mean(1);
y_mean = xx_mean(2);
translated_ellipse = rotated_ellipse + [x_mean; y_mean];

% --- Plotting ------------------------------------------------------------
% Plot mean position
plot(x_mean, y_mean, 'o', ...
     'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', color, ...
     'MarkerSize', 8);

% Plot uncertainty ellipse
plot(translated_ellipse(1,:), translated_ellipse(2,:), ...
     'Color', color, 'LineWidth', 2);

end

function [mu_LC, sigma_LC] = LinCov(tt, mu_LC, sigma_LC, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the mean state and covariance matrix using the linear
%   covariance (state transition matrix) method in the Planar Bicircular
%   Restricted Four-Body Problem (PBR4BP).
%
% Inputs:
%   tt       - Propagation time [scalar]
%   mu       - Mass ratio Earth–Moon system
%   rho      - Sun–Earth/Moon distance (adimensional units)
%   omega_s  - Sun angular velocity in rotating frame [rad/s]
%   m_s      - Sun mass parameter
%   mu_LC    - Mean state vector [4×1]
%   sigma_LC - Covariance matrix [4×4]
%
% Outputs:
%   mu_LC    - Propagated mean state [4×1]
%   sigma_LC - Propagated covariance matrix [4×4]
%--------------------------------------------------------------------------

% Propagate mean and covariance through STM
[mu_LC, sigma_LC, ~, ~, ~, ~] = propagatePBR4BP(tt, mu_LC, sigma_LC, mu, rho, omega_s, m_s);

end

function [mu_UT, sigma_UT] = UT(tt, mu_UT, sigma_UT, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the mean and covariance of a random variable using the
%   Unscented Transform (UT) in the context of the Planar Bicircular
%   Restricted Four-Body Problem (PBR4BP). The UT approximates the 
%   nonlinear transformation of Gaussian random variables by deterministically 
%   sampling sigma points and propagating them through the dynamics.
%
% Inputs:
%   tt       - Propagation time [scalar]
%   mu_UT    - Mean state vector [n×1]
%   sigma_UT - Covariance matrix [n×n] (must be positive definite)
%   mu       - Earth–Moon mass ratio
%   rho      - Sun–Earth/Moon distance (adimensional units)
%   omega_s  - Sun angular velocity in rotating frame [rad/s]
%   m_s      - Sun mass parameter
%
% Outputs:
%   mu_UT    - Propagated mean state [n×1]
%   sigma_UT - Propagated covariance matrix [n×n]
%--------------------------------------------------------------------------

% UT Parameters
alpha = 1;   % Spread of sigma points (scaling)
beta  = 2;   % Optimal for Gaussian distributions
n     = length(mu_UT);   % State dimension

% Scaling parameters
lambda = n * (alpha^2 - 1);
gamma  = sqrt(n + lambda);

% Cholesky decomposition
try
    L = chol(sigma_UT, 'lower');   % Lower-triangular factor
catch
    error('Initial covariance matrix sigma_UT is not positive definite.');
end

% Weights for mean and covariance 
Wm = [lambda/(n+lambda), repmat(1/(2*(n+lambda)), 1, 2*n)];
Wc = [lambda/(n+lambda) + (1 - alpha^2 + beta), repmat(1/(2*(n+lambda)), 1, 2*n)];

% Sigma point generation 
B      = gamma * L;
sigmas = [mu_UT, mu_UT + B, mu_UT - B];   % [n × (2n+1)]

% Propagation of sigma points
sigmas_f = NaN(n, 2*n+1);
for ii = 1:2*n+1
    [sig_f, ~, ~, ~, ~, ~] = propagatePBR4BP(tt, sigmas(:,ii), sigma_UT, mu, rho, omega_s, m_s);
    sigmas_f(:,ii) = sig_f;
end

% Mean update 
mu_UT = sigmas_f * Wm';

% Covariance update
sigma_UT = zeros(n,n);
for ii = 1:2*n+1
    diff = sigmas_f(:,ii) - mu_UT;
    sigma_UT = sigma_UT + Wc(ii) * (diff * diff');
end

end

function [mu_MC, sigma_MC] = MC(tt, mu_MC, sigma_MC, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates mean and covariance using a Monte Carlo (MC) approach in the
%   Planar Bicircular Restricted Four-Body Problem (PBR4BP). A population of
%   initial samples is drawn from a Gaussian distribution defined by the
%   given mean and covariance, each sample is propagated, and the final mean
%   and covariance are estimated empirically.
%
% Inputs:
%   tt       - Propagation time [scalar]
%   mu_MC    - Mean state vector [n×1]
%   sigma_MC - Covariance matrix [n×n]
%   mu       - Earth–Moon mass ratio
%   rho      - Sun–Earth/Moon distance (adimensional units)
%   omega_s  - Sun angular velocity in rotating frame [rad/s]
%   m_s      - Sun mass parameter
%
% Outputs:
%   mu_MC    - Propagated mean state [n×1]
%   sigma_MC - Propagated covariance matrix [n×n]
%--------------------------------------------------------------------------

% Monte Carlo population 
n     = length(mu_MC);   % State dimension
n_pop = 1000;            % Number of samples
R     = mvnrnd(mu_MC, sigma_MC, n_pop)';   % Random draws [n×n_pop]

% Propagate each sample 
xx_MC = NaN(n, n_pop);
for ii = 1:n_pop
    [xf_MC, ~, ~, ~, ~, ~] = propagatePBR4BP(tt, R(:,ii), sigma_MC, mu, rho, omega_s, m_s);
    xx_MC(:,ii) = xf_MC;
end

% Empirical mean
mu_MC = mean(xx_MC, 2);

% Empirical covariance 
diff_MC  = xx_MC - mu_MC;
sigma_MC = (diff_MC * diff_MC') / (n_pop - 1);

end

function x_std = standardize(x)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Standardizes a vector by removing its mean and scaling it to unit
%   variance. The result has zero mean and unit standard deviation.
%
% Inputs:
%   x     - Input vector of n samples [n×1 or 1×n]
%
% Outputs:
%   x_std - Standardized vector (same size as x), with:
%              mean(x_std) = 0
%              std(x_std)  = 1
%--------------------------------------------------------------------------

% Compute mean and standard deviation
mu    = mean(x);
sigma = std(x);

% Check for zero standard deviation
if sigma == 0
    error('Standard deviation is zero. Cannot standardize the vector.');
end

% Standardization
x_std = (x - mu) / sigma;

end

