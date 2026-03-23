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

%% EX 3: Sequential filters

% Constant for arcseconds to radians conversions
arcsec2rad = pi / (180*3600);

% Spice
cspice_furnsh('assignment02.tm');
addpath('tle')
addpath('sgp4')
fprintf('Number of LSK  kernels: %d\n', cspice_ktotal('lsk'));
fprintf('Number of SPK  kernels: %d\n', cspice_ktotal('spk'));
fprintf('Number of PCK  kernels: %d\n', cspice_ktotal('pck'));
fprintf('Number of CK   kernels: %d\n', cspice_ktotal('ck'));
fprintf('Number of TEXT kernels: %d\n', cspice_ktotal('TEXT'));
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));
mu = cspice_bodvrd('Moon','GM',1);

% Initial conditions for the lunar orbiter
r0 = [4307.844185282820; -1317.980749248651; 2109.210101634011]; % [km]
v0 = [-0.110997301537882; -0.509392750828585; 0.815198807994189]; % [km/s]
x0 = [r0;v0];

% Time span
et0 = cspice_str2et('2024-11-18T16:30:00.000');
etf = cspice_str2et('2024-11-18T20:30:00.000');
freq = 30;
et_vect = et0:freq:etf;
t_plot = datetime(cspice_timout(et_vect,...
    'YYYY-MM-DD HR:MN:SC.###'),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

% Uncertainties
sigma_rho = 100/1000; %Measurements noise [km]
P0 = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]); %Covariance


%% EX 3.1: Visibility windows

% Propagate sc state wrt Moon centered inertial (MCI)
[~,~,sc_mci,~] = propagate(et_vect,x0,mu);
r_sc_mci = sc_mci(:,1:3)';

% Lunar lander
lander.name = 'MOONLANDER';
lander.el_min = 0;
lander.lat = deg2rad(78); % rad
lander.lon = deg2rad(15); % rad
lander.alt = 0/1000;      % km

% Compute flattening
radii = cspice_bodvrd('MOON', 'RADII', 3);
r1 = radii(1); r3 = radii(3);
flat = (r1 - r3) / r1;

% Compute lander pos wrt Moon center (MCMF)
lander_mcmf = cspice_pgrrec('MOON', lander.lon, lander.lat, lander.alt, r1, flat);

% Rotate sc position vector in Moon-fixed frame (MCMF)
mci2MCMF = cspice_pxform( 'J2000', 'IAU_MOON', et_vect );
r_sc_mcmf = NaN(size(r_sc_mci));
for k = 1 : length(et_vect)
    r_sc_mcmf(:,k) = mci2MCMF(:,:,k)*r_sc_mci(:,k);
end

% Sc position wrt lander (via lander coordinates)
Y_pred = pointing_lander(lander, et_vect, r_sc_mci, true);

% Visibility window
[Y_vis, et_vect_new] = visibility(Y_pred, et_vect, lander);
t_vis = datetime(cspice_timout(et_vect_new, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Plot in Moon-fixed frame
[x_sph, y_sph, z_sph] = sphere(30);
x_sph = r1*x_sph;
y_sph = r1*y_sph;
z_sph = r1*z_sph;
figure
surf(x_sph, y_sph, z_sph, 'EdgeColor', 'k','FaceColor',[0.7 0.7 0.7],'EdgeAlpha',0.0001);
hold on
plot3(r_sc_mcmf(1,:), r_sc_mcmf(2,:), r_sc_mcmf(3,:),'LineWidth',4,'Color',"#0072BD") 
plot3(r_sc_mcmf(1,1), r_sc_mcmf(2,1), r_sc_mcmf(3,1),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot3(r_sc_mcmf(1,end), r_sc_mcmf(2,end), r_sc_mcmf(3,end),'diamond','MarkerEdgeColor', "#EDB120",'MarkerFaceColor',"#EDB120",'MarkerSize',20)
axis equal
plot3(lander_mcmf(1), lander_mcmf(2), lander_mcmf(3),'o','MarkerSize',20,'MarkerFaceColor',"#A2142F",'MarkerEdgeColor','k')
xlabel('x [km]', 'Interpreter','latex','FontSize',40)
ylabel('y [km]', 'Interpreter','latex','FontSize',40)
zlabel('z [km]', 'Interpreter','latex','FontSize',40)
text(lander_mcmf(1), lander_mcmf(2), lander_mcmf(3), 'Moon Lander', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
legend(["","Orbiter trajectory","Initial position","Final position",""],'Location','best','Interpreter','Latex',FontSize = 40)
%title('Orbiter and Lander in Moon-fixed frame','Interpreter','latex', 'Interpreter','latex','FontSize',16)

% Plot elevation
figure 
plot(t_vis, Y_vis(2,:), '.', 'MarkerEdgeColor', [0 0.4470 0.7410],'MarkerSize',20)
yline(0, 'k--', 'LineWidth', 4, 'DisplayName','Minimum Elevation')
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
ylim([-10 90])
legend('Elevation','Minimum elevation','Interpreter','latex','FontSize',40)
%title('Elevation profile', 'Interpreter', 'latex','FontSize',16)
grid on


%% EX 3.2: Simulate measurements

% Sc position wrt lander (via kernels)
Y_pred = pointing_lander(lander, et_vect, r_sc_mci, false);

% Visibility window
[Y_vis, et_vect_new] = visibility(Y_pred, et_vect, lander);
t_vis = datetime(cspice_timout(et_vect_new, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Simulate measurements with noise
range = Y_vis(3,:);
mu_pos = [r_sc_mci;range].';
P_pos = (sigma_rho)^2*diag(ones(4,1));
pos_meas = mvnrnd(mu_pos, P_pos).';
r_meas = pos_meas(1:3,:);

% Visibility window after noise
Y_pred_noise = pointing_lander(lander, et_vect, r_meas, false);
[Y_vis_noise, et_vect_noise] = visibility(Y_pred_noise, et_vect, lander); 

% Error vector
err = mu_pos.'-pos_meas;
z_err = reshape(err, 4*size(err, 2), 1);
mu_z = mean(z_err);
sigma_z = std(z_err);
z_err = (z_err-mu_z)/sigma_z;

% Histogram 
nbins = ceil(sqrt(length(z_err)));

% Plot with standard normal pdf
figure
xx = linspace(-3, 3, 1000);
yy = normpdf(xx);
hold on
histogram(z_err, nbins, 'Normalization','pdf')
plot(xx, yy, 'r', 'LineWidth',4)
xline(-3, 'k--', 'LineWidth', 4)
xline(3, 'k--', 'LineWidth', 4)
xlim([-4 4])
xlabel('Values $z_{pos}$ [-]','Interpreter','latex','FontSize',40)
ylabel('pdf [-]','Interpreter','latex','FontSize',40)
legend('Standardized values', 'Standard normal pdf', '$3 \sigma$ boundaries','Interpreter','latex','FontSize',40)
grid on

% Outliers: values out of the 3 sigma
ind_plus = find(z_err>3);
ind_minus = find(z_err<-3);
outliers = 100*(length(ind_minus) + length(ind_plus))/length(z_err); 


%% EX 3.3: Lunar orbiter absolute state

% Perturbed initial state
x0_UKF = mvnrnd(x0, P0(1:6,1:6)).';

% Noise matrix
R = P_pos(1:3,1:3);

% Unscented Kalman Filter (UKF)
[x_hat_sc, P_sc] = UKF(et_vect, x0_UKF, P0(1:6,1:6), r_meas, R, mu, true);

% Errors
err_sc = x_hat_sc - sc_mci';
err_sc_r = vecnorm(err_sc(1:3,:));
err_sc_v = vecnorm(err_sc(4:6,:));

% 3sigma
sigma_sc_r = NaN(size(P_sc,3),1);
sigma_sc_v = NaN(size(P_sc,3),1);
for ii = 1:size(P_sc,3)
    P_sc_r = P_sc(1:3,1:3,ii);
    P_sc_v = P_sc(4:6,4:6,ii);
    sigma_sc_r(ii) = sqrt(max(eig(P_sc_r)));
    sigma_sc_v(ii) = sqrt(max(eig(P_sc_v)));
end

% Plot errors with 3sigma
figure
subplot(1,2,1)
semilogy(t_vis, err_sc_r, 'LineWidth', 4)
hold on
semilogy(t_vis, 3*sigma_sc_r, 'LineWidth', 4)
grid on
ylabel('$\varepsilon_r$, 3$\sigma_r$ [km]','Interpreter','latex','FontSize',40)
legend('Position error','$3\sigma_r$','Interpreter','latex','FontSize',40)

subplot(1,2,2)
semilogy(t_vis, err_sc_v, 'LineWidth', 4)
hold on
semilogy(t_vis, 3*sigma_sc_v, 'LineWidth', 4)
grid on
ylabel('$\varepsilon_v$, 3$\sigma_v$ [km/s]','Interpreter','latex','FontSize',40)
legend('Velocity error','$3\sigma_v$','Interpreter','latex','FontSize',40)


%% EX 3.4: Lunar lander coordinates

% Perturbed coordinates
coord_UKF = mvnrnd([lander.lat;lander.lon], P0(7:8,7:8)).';

% Noise matrix
R = P_pos;

% Unscented Kalman Filter (UKF)
[x_hat_lander, P_lander] = UKF(et_vect, [x0_UKF; coord_UKF], P0, pos_meas, R, mu, false);

% State error
err_lander = NaN(size(x_hat_lander));
err_lander(1:6,:) = x_hat_lander(1:6,:) - sc_mci';
err_sc_r = vecnorm(err_lander(1:3,:));
err_sc_v = vecnorm(err_lander(4:6,:));

% Coordinates error
lat_ex = deg2rad(78.229772);
lon_ex = deg2rad(15.407786);
err_lander(7:8,:) = x_hat_lander(7:8,:) - [lat_ex;lon_ex];
err_lander_lat = abs(err_lander(7,:));
err_lander_lon = abs(err_lander(8,:));

% 3sigma
sigma_sc_r = NaN(size(P_lander,3),1);
sigma_sc_v = NaN(size(P_lander,3),1);
sigma_lander_lat = NaN(size(P_lander,3),1);
sigma_lander_lon = NaN(size(P_lander,3),1);
for ii = 1:size(P_lander,3)
    P_sc_r = P_lander(1:3,1:3,ii);
    P_sc_v = P_lander(4:6,4:6,ii);
    P_lander_lat = P_lander(7,7,ii);
    P_lander_lon = P_lander(8,8,ii);
    sigma_sc_r(ii) = sqrt(max(eig(P_sc_r)));
    sigma_sc_v(ii) = sqrt(max(eig(P_sc_v)));
    sigma_lander_lat(ii) = sqrt(max(eig(P_lander_lat)));
    sigma_lander_lon(ii) = sqrt(max(eig(P_lander_lon)));
end

% Plot errors with 3sigma
figure
subplot(2,2,1)
semilogy(t_vis, err_sc_r, 'LineWidth', 4)
hold on
semilogy(t_vis, 3*sigma_sc_r, 'LineWidth', 4)
grid on
ylabel('$\varepsilon_r$, 3$\sigma_r$ [km]','Interpreter','latex','FontSize',40)
legend('Variable error','Variable $3\sigma$','Interpreter','latex','FontSize',40)
title('Position','Interpreter','latex','FontSize',40)

subplot(2,2,2)
semilogy(t_vis, err_sc_v, 'LineWidth', 4)
hold on
semilogy(t_vis, 3*sigma_sc_v, 'LineWidth', 4)
grid on
ylabel('$\varepsilon_v$, 3$\sigma_v$ [km/s]','Interpreter','latex','FontSize',40)
%legend('Velocity error','$3\sigma_v$','Interpreter','latex','FontSize',40)
title('Velocity','Interpreter','latex','FontSize',40)

subplot(2,2,3)
semilogy(t_vis, err_lander_lat, 'LineWidth', 4)
hold on
semilogy(t_vis, 3*sigma_lander_lat, 'LineWidth', 4)
grid on
ylabel('$\varepsilon_{lat}$, 3$\sigma_{lat}$ [deg]','Interpreter','latex','FontSize',40)
%legend('Latitude error','$3\sigma_{lat}$','Interpreter','latex','FontSize',40)
title('Latitude','Interpreter','latex','FontSize',40)

subplot(2,2,4)
semilogy(t_vis, err_lander_lon, 'LineWidth', 4)
hold on
semilogy(t_vis, 3*sigma_lander_lon, 'LineWidth', 4)
grid on
ylabel('$\varepsilon_{lon}$, 3$\sigma_{lon}$ [deg]','Interpreter','latex','FontSize',40)
%legend('Longitude error','$3\sigma_{lon}$','Interpreter','latex','FontSize',40)
title('Longitude','Interpreter','latex','FontSize',40)


%% FUNCTIONS

function dy = two_body_rhs(~, y, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the right-hand side of the two-body problem equations of
%   motion under Newtonian gravity. The system is expressed in Cartesian
%   coordinates, with the gravitational parameter mu = GM.
%
% Inputs:
%   ~  - Placeholder for time (not used explicitly)
%   y  - State vector [6×1]:
%          y(1:3) = position [x; y; z]
%          y(4:6) = velocity [vx; vy; vz]
%   mu - Gravitational parameter of central body [km^3/s^2]
%
% Outputs:
%   dy - Time derivative of state vector [6×1]:
%          dy(1:3) = velocity components
%          dy(4:6) = acceleration due to gravity
%--------------------------------------------------------------------------

% Position norm
r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);

% Derivatives of state
dy = [ y(4);                 % dx/dt = vx
       y(5);                 % dy/dt = vy
       y(6);                 % dz/dt = vz
      -(mu/r^3) * y(1);      % dvx/dt
      -(mu/r^3) * y(2);      % dvy/dt
      -(mu/r^3) * y(3) ];    % dvz/dt

end

function [xf, tf, xx, tt] = propagate(t_vect, x0, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the two-body dynamics under Newtonian gravity using a
%   numerical ODE solver (ode78). The state includes Cartesian position
%   and velocity components.
%
% Inputs:
%   t_vect - Time vector for integration [1×N]
%   x0     - Initial state vector [6×1]:
%              [r0x; r0y; r0z; v0x; v0y; v0z]
%   mu     - Gravitational parameter of central body [km^3/s^2]
%
% Outputs:
%   xf - Final state vector at t = t_vect(end) [6×1]
%   tf - Final propagation time [scalar]
%   xx - Full propagated trajectory [N×6]
%   tt - Time vector returned by integrator [N×1]
%--------------------------------------------------------------------------

% Integration options (tight tolerances)
options = odeset('reltol', 1e-13, 'abstol', 1e-20);

% Integrate two-body dynamics
[tt, xx] = ode78(@(t,x) two_body_rhs(t,x,mu), t_vect, x0, options);

% Extract final state and final time
xf = xx(end,:)';
tf = tt(end);

end

function [Y_pred] = pointing_lander(lander, et_vec, r_orbiter_in, flag)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the predicted measurements (azimuth, elevation, and range) of
%   an orbiter as observed from a lunar lander, either using the lander’s
%   provided geodetic coordinates or extracting its position from SPICE
%   kernels.
%
% Inputs:
%   lander       - Structure containing lander info. Fields:
%                    .lat   - latitude [rad] (used if flag = true)
%                    .lon   - longitude [rad] (used if flag = true)
%                    .alt   - altitude [km]   (used if flag = true)
%                    .name  - lander name in SPICE kernels (used if flag = false)
%   et_vec       - Ephemeris time vector [1×N] (s past J2000 TDB)
%   r_orbiter_in - Propagated orbiter trajectory in J2000 frame [3×N] (km)
%   flag         - If true, use (lat, lon, alt) guess.
%                  If false, extract lander state from SPICE kernel.
%
% Outputs:
%   Y_pred - Predicted measurements [3×N]:
%              [Azimuth (deg); Elevation (deg); Range (km)]
%--------------------------------------------------------------------------

frame1 = 'J2000';

if flag
    % --- Case 1: Use guessed lander coordinates --------------------------
    lat = lander.lat;
    lon = lander.lon;
    alt = lander.alt;

    % Moon reference radii (requires PCK kernel)
    radii  = cspice_bodvrd('MOON', 'RADII', 3);
    r_eq   = radii(1);   % equatorial radius
    r_pol  = radii(3);   % polar radius
    flat   = (r_eq - r_pol) / r_eq;

    % Lander position in Moon-fixed frame
    r_lander = cspice_pgrrec('MOON', lon, lat, alt, r_eq, flat);

    % Rotation: J2000 → Moon-fixed → TOPO frame
    MJ20002MCMF = cspice_pxform(frame1, 'IAU_MOON', et_vec);
    MCMF2TOPO   = cspice_eul2m(lat - pi, pi - lon, pi/2, 2, 1, 2);

    % Orbiter in TOPO frame
    r_orbiter_topo = NaN(size(r_orbiter_in));
    for k = 1:length(et_vec)
        r_orbiter_mf      = MJ20002MCMF(:,:,k) * r_orbiter_in(:,k);
        r_rel             = r_orbiter_mf - r_lander;
        r_orbiter_topo(:,k) = MCMF2TOPO * r_rel;
    end

else
    % --- Case 2: Use lander SPICE kernel position ------------------------
    lander_name = lander.name;
    frame2 = strcat(lander_name, '_TOPO');

    % Lander position wrt Moon in inertial frame
    r_lander = cspice_spkpos(lander_name, et_vec, frame1, 'NONE', 'Moon');

    % Rotation J2000 → TOPO frame
    R = cspice_pxform(frame1, frame2, et_vec);

    % Orbiter in TOPO frame
    r_rel = r_orbiter_in - r_lander;
    r_orbiter_topo = NaN(size(r_orbiter_in));
    for k = 1:length(et_vec)
        r_orbiter_topo(:,k) = R(:,:,k) * r_rel(:,k);
    end
end

% Measurements 
range = vecnorm(r_orbiter_topo);
El    = asin(r_orbiter_topo(3,:) ./ range);
Az    = wrapToPi(atan2(r_orbiter_topo(2,:) ./ range, ...
                       r_orbiter_topo(1,:) ./ range));

% Convert to degrees
Az = rad2deg(Az);
El = rad2deg(El);

% Predicted measurements
Y_pred = [Az; El; range];

end

function [Y_pred, et] = visibility(Y_pred, et, station)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Filters the predicted measurements by enforcing a minimum elevation
%   angle constraint, returning only the visible portions of the pass.
%
% Inputs:
%   Y_pred  - Predicted measurements [3×N]:
%                [Azimuth (deg); Elevation (deg); Range (km)]
%   et      - Ephemeris time vector [1×N] (s past J2000 TDB)
%   station - Structure containing station info, must include:
%                .el_min : Minimum elevation angle threshold [deg]
%
% Outputs:
%   Y_pred - Filtered predicted measurements [3×M], M ≤ N
%   et     - Filtered ephemeris times corresponding to visible epochs [1×M]
%--------------------------------------------------------------------------

El_min = station.el_min;       % Minimum elevation [deg]
El     = Y_pred(2,:);          % Extract elevation track

% Visibility mask
ind = (El > El_min);

% Apply filter
et     = et(ind);
Y_pred = Y_pred(:, ind);

end

function [Xi, W_m, W_c] = sigma_points(x0, alpha, beta, P0)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the Unscented Transform (UT) sigma points and their associated
%   mean and covariance weights, given an initial Gaussian distribution.
%
% Inputs:
%   x0    - Mean state vector [n×1]
%   alpha - Spreading parameter (typically small, e.g. 1e-3)
%   beta  - Tuning parameter (2 is optimal for Gaussian priors)
%   P0    - Initial covariance matrix [n×n]
%
% Outputs:
%   Xi  - Sigma points matrix [n×(2n+1)]
%           Xi(:,1)     = x0
%           Xi(:,2:n+1) = x0 + gamma*chol(P0)
%           Xi(:,n+2:2n+1) = x0 - gamma*chol(P0)
%   W_m - Weights for mean calculation [1×(2n+1)]
%   W_c - Weights for covariance calculation [1×(2n+1)]
%--------------------------------------------------------------------------

n = length(x0);  % State dimension

% UT scaling parameters
lambda = n * (alpha^2 - 1);
gamma  = sqrt(n + lambda);

% Cholesky decomposition of covariance
try
    L = chol(P0, 'lower');
catch
    error('Initial covariance matrix P0 is not positive definite.');
end

% Weights for mean and covariance
W_m = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2*n)];
W_c = [lambda / (n + lambda) + (1 - alpha^2 + beta), ...
       repmat(1 / (2 * (n + lambda)), 1, 2*n)];

% Sigma points
Xi = [x0, x0 + gamma*L, x0 - gamma*L];

end

function [Xi_k, x_hat_pred, P_pred, Yi_k, y_hat_pred] = predicted_step(t0, tf, Xi, W_m, W_c, mu, flag)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Performs the **prediction step** of the Unscented Kalman Filter (UKF).
%   Each sigma point is propagated through the nonlinear dynamics, and the
%   predicted state mean, covariance, and measurements are computed.
%   Depending on the flag, the measurement model includes only position
%   (orbital tracking) or position + range from a lunar lander.
%
% Inputs:
%   t0   - Initial propagation epoch [s]
%   tf   - Final propagation epoch [s]
%   Xi   - Sigma points at time t0 [n×(2n+1)]
%   W_m  - Weights for mean [1×(2n+1)]
%   W_c  - Weights for covariance [1×(2n+1)]
%   mu   - Gravitational parameter [km^3/s^2]
%   flag - If true:   measurement = spacecraft position only
%          If false:  measurement = spacecraft position + range from lander
%
% Outputs:
%   Xi_k      - Propagated sigma points at time tf [n×(2n+1)]
%   x_hat_pred- Predicted mean state vector [n×1]
%   P_pred    - Predicted covariance matrix [n×n]
%   Yi_k      - Propagated sigma points in measurement space [m×(2n+1)]
%   y_hat_pred- Predicted mean measurement vector [m×1]
%--------------------------------------------------------------------------

if flag
    % --- Case 1: State propagation with position-only measurements -------
    n = size(Xi,1);
    Xi_k       = zeros(size(Xi));
    x_hat_pred = zeros(n,1);
    P_pred     = zeros(n);
    y_hat_pred = zeros(3,1);

    for ii = 1:2*n+1
        % Propagate sigma points
        x0 = Xi(:,ii);
        [xf,~,~,~] = propagate([t0 tf], x0, mu);
        Xi_k(:,ii) = xf;

        % Predicted mean state
        x_hat_pred = x_hat_pred + W_m(ii)*xf;
    end

    for ii = 1:2*n+1
        % Predicted covariance
        diff   = Xi_k(:,ii) - x_hat_pred;
        P_pred = P_pred + W_c(ii) * (diff*diff');
    end

    % Measurement model = position only
    Yi_k = Xi_k(1:3,:);

    for ii = 1:2*n+1
        y_hat_pred = y_hat_pred + W_m(ii)*Yi_k(:,ii);
    end

else
    % --- Case 2: State propagation with lander range measurement ---------
    n = size(Xi,1);
    Xi_k       = zeros(size(Xi));
    x_hat_pred = zeros(n,1);
    P_pred     = zeros(n);
    range      = zeros(1,2*n+1);
    y_hat_pred = zeros(4,1);

    % Moon flattening for conversion to MCMF
    radii = cspice_bodvrd('MOON', 'RADII', 3);
    r_eq  = radii(1); 
    r_pol = radii(3);
    flat  = (r_eq - r_pol) / r_eq;

    for ii = 1:2*n+1
        % Propagate spacecraft sigma point
        x0 = Xi(1:6,ii);
        [xf,~,~,~] = propagate([t0 tf], x0, mu);
        Xi_k(1:6,ii) = xf;

        % Lander position in Moon-fixed frame
        lat     = Xi(7,ii);
        lon     = Xi(8,ii);
        r_lander = cspice_pgrrec('MOON', lon, lat, 0, r_eq, flat);

        % Transform spacecraft to Moon-fixed frame
        r_sc       = xf(1:3);
        mci2MCMF   = cspice_pxform('J2000', 'IAU_MOON', tf);
        r_sc_MCMF  = mci2MCMF * r_sc;

        % Relative position and range
        r_rel       = r_sc_MCMF - r_lander;
        range(1,ii) = norm(r_rel);

        % Augment sigma point (lander lat/lon unchanged)
        Xi_k(7,ii) = lat;
        Xi_k(8,ii) = lon;

        % Predicted mean state
        x_hat_pred = x_hat_pred + W_m(ii)*Xi_k(:,ii);
    end

    for ii = 1:2*n+1
        diff   = Xi_k(:,ii) - x_hat_pred;
        P_pred = P_pred + W_c(ii) * (diff*diff');
    end

    % Measurement model = spacecraft position + range
    Yi_k = [Xi_k(1:3,:); range];

    for ii = 1:2*n+1
        y_hat_pred = y_hat_pred + W_m(ii)*Yi_k(:,ii);
    end
end

end

function [x_hat_up, P_up] = updated_step(Xi_k, x_hat_pred, P_pred, Yi_k, y_hat_pred, R, W_c, y_k)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Performs the **update step** of the Unscented Kalman Filter (UKF).
%   Given the propagated sigma points and predicted measurement statistics,
%   this function updates the state estimate and covariance using the
%   actual measurement y_k.
%
% Inputs:
%   Xi_k       - Propagated sigma points at time k [n×(2n+1)]
%   x_hat_pred - Predicted state mean [n×1]
%   P_pred     - Predicted state covariance [n×n]
%   Yi_k       - Propagated sigma points in measurement space [m×(2n+1)]
%   y_hat_pred - Predicted measurement mean [m×1]
%   R          - Measurement noise covariance [m×m]
%   W_c        - Covariance weights for sigma points [1×(2n+1)]
%   y_k        - Actual measurement at time k [m×1]
%
% Outputs:
%   x_hat_up - Updated state mean [n×1]
%   P_up     - Updated state covariance [n×n]
%--------------------------------------------------------------------------

n      = size(Xi_k,1);   % State dimension
n_meas = size(Yi_k,1);   % Measurement dimension

% Initialize measurement covariance and cross-covariance
P_yy_k = R;              % Innovation covariance
P_xy_k = zeros(n,n_meas);

% Accumulate contributions from sigma points
for ii = 1:2*n+1
    diff_y = Yi_k(:,ii) - y_hat_pred;
    diff_x = Xi_k(:,ii) - x_hat_pred;

    % Measurement covariance
    P_yy_k = P_yy_k + W_c(ii) * (diff_y * diff_y');

    % Cross covariance
    P_xy_k = P_xy_k + W_c(ii) * (diff_x * diff_y');
end

% Kalman gain
K_k = P_xy_k / P_yy_k;

% Updated mean state
x_hat_up = x_hat_pred + K_k * (y_k - y_hat_pred);

% Updated covariance (Joseph form could be used for numerical stability)
P_up = P_pred - K_k * P_yy_k * K_k';

end

function [x_hat, P] = UKF(et_vect, x0, P0, meas, R, mu, flag)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Implements the Unscented Kalman Filter (UKF) for nonlinear state 
%   estimation. The state is propagated using sigma points through the 
%   nonlinear dynamics, and updated with noisy measurements at each epoch.
%
% Inputs:
%   et_vect - Vector of time instants [1×nt]
%   x0      - Initial state estimate [nx×1]
%   P0      - Initial state covariance matrix [nx×nx]
%   meas    - Measurement sequence [m×nt]
%   R       - Measurement noise covariance [m×m]
%   mu      - Gravitational parameter [km^3/s^2]
%   flag    - If true:   measurement = spacecraft position only
%             If false:  measurement = spacecraft position + range from lander
%
% Outputs:
%   x_hat - Filtered state estimates across all time steps [nx×nt]
%   P     - Filtered state covariance matrices [nx×nx×nt]
%--------------------------------------------------------------------------

% Initialization
nx   = length(x0);   % State dimension
nt   = length(et_vect); % Number of epochs
x_hat = NaN(nx, nt);
P     = NaN(nx, nx, nt);

% UT parameters
alpha = 0.01;   % Spread of sigma points (small positive value)
beta  = 2;     % Optimal for Gaussian distributions

% Store initial condition
x_hat(:,1) = x0;
P(:,:,1)   = P0;
t0 = et_vect(1);

% UKF loop over measurements
for k = 2:nt
    tk = et_vect(k);     % Current epoch
    y_k = meas(:,k);     % Current measurement

    % Sigma points generation
    [Xi, W_m, W_c] = sigma_points(x0, alpha, beta, P0);

    % Prediction step (state & measurement propagation)
    [Xi_k, x_hat_pred, P_pred, Yi_k, y_hat_pred] = ...
        predicted_step(t0, tk, Xi, W_m, W_c, mu, flag);

    % Update step (correction with measurement y_k)
    [x_hat_up, P_up] = updated_step(Xi_k, x_hat_pred, P_pred, ...
                                    Yi_k, y_hat_pred, R, W_c, y_k);

    % Store results
    x_hat(:,k) = x_hat_up;
    % Symmetrize covariance to avoid numerical asymmetries
    P(:,:,k)   = (P_up + P_up')/2;  

    % Prepare for next iteration
    t0 = tk;
    x0 = x_hat_up;
    P0 = P(:,:,k);
end

end

