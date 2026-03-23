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

%% EX 2: Batch filter

% Ground stations
KOUROU.name = 'KOUROU';
KOUROU.el_min = 6;
KOUROU.frequency = 60;
KOUROU.cost = 30000;

TROLL.name = 'TROLL';
TROLL.el_min = 0;
TROLL.frequency = 30;
TROLL.cost = 35000;

SVALBARD.name = 'SVALBARD';
SVALBARD.el_min = 8;
SVALBARD.frequency = 60;
SVALBARD.cost = 35000;

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
mu = cspice_bodvrd('Earth','GM',1);

% Inputs
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)


%% EX 2.1: Visibility windows

% Initialize the satrec structure and TLE
[ satrec, longstr1, longstr2 ] = read_3LE( 36036, '36036.3le', whichconst );

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);
spacecraftName = 'SMOS';
fprintf('\nSatellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s\n', sat_epoch_str);

% Evaluate the TLE
[~,rteme,vteme] = sgp4(satrec, 0.0);

% Get the osculating orbital elements
elts = cspice_oscelt( [rteme;vteme], sat_epoch_et, satrec.mu );
fprintf('\n*** Osculating state ***')
fprintf('\nSMA   [km]:  %.5f', elts(1)/(1-elts(2))); % elts(1)=p --> a (SMA)
fprintf('\nECC   [km]:  %.8f', elts(2));
fprintf('\nINC  [deg]: %.5f', elts(3)*cspice_dpr());
fprintf('\nRAAN [deg]: %.5f', elts(4)*cspice_dpr());
fprintf('\nARGP [deg]: %.5f', elts(5)*cspice_dpr());
fprintf('\nM.AN [deg]: %.5f\n', elts(6)*cspice_dpr());

% Nutation correction: coefficients delta-Psi and delta-Epsilon corrections
ddpsi = -0.114761*arcsec2rad; %  [rad]
ddeps = -0.007531*arcsec2rad; %  [rad]

% Precession: centuries from TDT 2000 January 1 00:00:00.000 
et_ref = sat_epoch_et; % Epoch we want to perform our conversion TEME to ECI
ttt = cspice_unitim(et_ref, 'ET', 'TDT')/cspice_jyear()/100;

% Transform TEME to ECI vectors
ateme = [0;0;0];
[reci, veci, aeci] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
fprintf(1,'\nPosition in ECI %14.7f %14.7f %14.7f',reci );
fprintf(1,'\nVelocity in ECI %14.9f %14.9f %14.9f',veci );
fprintf(1,'\nAcceleration in ECI %14.9f %14.9f %14.9f\n',aeci );

% Propagation to initial epoch
xref = [reci;veci]; 
et0 = cspice_str2et('2024-11-18T20:30:00.000');
[x0,~,~,~] = propagate([et_ref et0],xref,mu);

% Final epoch
etf = cspice_str2et('2024-11-18T22:15:00.000');

% Equally time grid vector, depending on the frequency of measurement acquisition
npoints1 = round((etf-et0)/KOUROU.frequency)+1;
et_vect1 = linspace(et0, etf, npoints1);
npoints2 = round((etf-et0)/TROLL.frequency)+1;
et_vect2 = linspace(et0, etf, npoints2);

% Propagation from t0 to tf
[xf1,~,xx1,~] = propagate(et_vect1,x0,mu);
[xf2,~,xx2,~] = propagate(et_vect2,x0,mu);

% Position
r_eci1 = xx1(:,1:3);
r_eci2 = xx2(:,1:3);

% KOUROU station computations
[azimuth_K, elevation_K, range_K] = pointing_antenna(KOUROU, et_vect1, r_eci1);

% Control of minimum elevation for visibility window
[elevation_K,azimuth_K,range_K,et_vect_newK] = visibility(elevation_K,azimuth_K,range_K,et_vect1,KOUROU);
t_visK = datetime(cspice_timout(et_vect_newK, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Plot
figure
subplot(2, 3, 1)
plot(t_visK, azimuth_K, '*', 'Color', [0 0.4470 0.7410],'LineWidth',5)
ylabel('Az [deg]', 'Interpreter', 'latex','FontSize',40)
title(KOUROU.name, 'Interpreter', 'latex','FontSize',40)
grid on

subplot(2, 3, 4)
plot(t_visK, elevation_K, '*', 'Color', [0.8500 0.3250 0.0980],'LineWidth',5)
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
grid on

% TROLL station computations
[azimuth_T, elevation_T, range_T] = pointing_antenna(TROLL, et_vect2, r_eci2);

% Control of minimum elevation for visibility window
[elevation_T,azimuth_T,range_T,et_vect_newT] = visibility(elevation_T,azimuth_T,range_T,et_vect2,TROLL);
t_visT = datetime(cspice_timout(et_vect_newT, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Plot
subplot(2, 3, 2)
plot(t_visT, azimuth_T, '*', 'Color', [0 0.4470 0.7410],'LineWidth',5)
ylabel('Az [deg]', 'Interpreter', 'latex','FontSize',40)
title(TROLL.name, 'Interpreter', 'latex','FontSize',40)
grid on

subplot(2, 3, 5)
plot(t_visT, elevation_T, '*', 'Color', [0.8500 0.3250 0.0980],'LineWidth',5)
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
grid on

% SVALBARD station computations
[azimuth_S, elevation_S, range_S] = pointing_antenna(SVALBARD, et_vect1, r_eci1);

% Control of minimum elevation for visibility window
[elevation_S,azimuth_S,range_S,et_vect_newS] = visibility(elevation_S,azimuth_S,range_S,et_vect1,SVALBARD);
t_visS = datetime(cspice_timout(et_vect_newS, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Plot
subplot(2, 3, 3)
plot(t_visS, azimuth_S, '*', 'Color', [0 0.4470 0.7410],'LineWidth',5)
ylabel('Az [deg]', 'Interpreter', 'latex','FontSize',40)
title(SVALBARD.name, 'Interpreter', 'latex','FontSize',40)
grid on

subplot(2, 3, 6)
plot(t_visS, elevation_S, '*', 'Color', [0.8500 0.3250 0.0980],'LineWidth',5)
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
grid on


%% EX 2.2.1: Simulate measurements

% Acceleration
a_teme = [0;0;0];

% KOUROU station 
reciK = NaN(length(et_vect_newK),3);
veciK = NaN(length(et_vect_newK),3);
for ii = 1:length(et_vect_newK)
    tsince = (et_vect_newK(ii) - et_ref)/60;

    % Propagation via sgp4
    [~, r_teme, v_teme] = sgp4(satrec, tsince);

    % Conversion to ECI
    [reciK(ii,:), veciK(ii,:), ~] = teme2eci( r_teme, v_teme, a_teme, ttt, ddpsi, ddeps);
end
[azimuthK, elevationK, rangeK] = pointing_antenna(KOUROU, et_vect_newK, reciK); 
[elevationK,azimuthK,rangeK,et_vect_newK] = visibility(elevationK,azimuthK,rangeK,et_vect_newK,KOUROU);

% TROLL station 
reciT = NaN(length(et_vect_newT),3);
for ii = 1:length(et_vect_newT)
    tsince = (et_vect_newT(ii) - et_ref)/60;

    % Propagation via sgp4
    [~, r_teme, v_teme] = sgp4(satrec, tsince);

    % Conversion to ECI
    [reciT(ii,:), ~, ~] = teme2eci( r_teme, v_teme, a_teme, ttt, ddpsi, ddeps); 
end
[azimuthT, elevationT, rangeT] = pointing_antenna(TROLL, et_vect_newT, reciT); 
[elevationT,azimuthT,rangeT,et_vect_newT] = visibility(elevationT,azimuthT,rangeT,et_vect_newT,TROLL);

% SVALBARD station 
reciS = NaN(length(et_vect_newS),3);
for ii = 1:length(et_vect_newS)
    tsince = (et_vect_newS(ii) - et_ref)/60;

    % Propagation via sgp4
    [~, r_teme, v_teme] = sgp4(satrec, tsince);

    % Conversion to ECI
    [reciS(ii,:), ~, ~] = teme2eci( r_teme, v_teme, a_teme, ttt, ddpsi, ddeps); 
end
[azimuthS, elevationS, rangeS] = pointing_antenna(SVALBARD, et_vect_newS, reciS);
[elevationS,azimuthS,rangeS,et_vect_newS] = visibility(elevationS,azimuthS,rangeS,et_vect_newS,SVALBARD);


%% EX 2.2.2: Add Noise

% Noise levels
sigma_az = 125 / 1000; 
sigma_el = 125 / 1000; 
sigma_range = 0.01; 

% Covariance matrix
R = diag([sigma_az^2, sigma_az^2, sigma_range^2]); 

% KOUROU station noise 
measK = [azimuthK, elevationK, rangeK];
measK_noisy = mvnrnd(measK,R);

% Control of minimum elevation for visibility window
[measK_noisy,et_vect_newK] = visibility_noisy(measK_noisy,et_vect_newK,KOUROU);
t_vis_kourou_sgp4 = datetime(cspice_timout(et_vect_newK, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
KOUROU.meas_noisy = measK_noisy;
KOUROU.tspan = [et0, et_vect_newK];

% TROLL station noise
measT = [azimuthT, elevationT, rangeT];
measT_noisy = mvnrnd(measT,R);

% Control of minimum elevation for visibility window
[measT_noisy,et_vect_newT] = visibility_noisy(measT_noisy,et_vect_newT,TROLL);
t_vis_troll_sgp4 = datetime(cspice_timout(et_vect_newT, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
TROLL.meas_noisy = measT_noisy;
TROLL.tspan = [et0, et_vect_newT];

% SVALBARD station noise
measS = [azimuthS, elevationS, rangeS];
measS_noisy = mvnrnd(measS,R);

% Control of minimum elevation for visibility window
[measS_noisy,et_vect_newS] = visibility_noisy(measS_noisy,et_vect_newS,SVALBARD);
t_vis_svalbard_sgp4 = datetime(cspice_timout(et_vect_newS, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
SVALBARD.meas_noisy = measS_noisy;
SVALBARD.tspan = [et0, et_vect_newS];

% KOUROU station simulated measurements plot
figure
subplot(3, 3, 1)
plot(t_vis_kourou_sgp4, measK_noisy(:,1), '*', 'Color', [0 0.4470 0.7410],'LineWidth',5)
ylabel('Az [deg]', 'Interpreter', 'latex','FontSize',40)
title(KOUROU.name, 'Interpreter', 'latex','FontSize',40)
grid on

subplot(3, 3, 4)
plot(t_vis_kourou_sgp4, measK_noisy(:,2), '*', 'Color', [0.8500 0.3250 0.0980],'LineWidth',5)
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
grid on

subplot(3, 3, 7)
plot(t_vis_kourou_sgp4, measK_noisy(:,3), '*', 'Color', 'k','LineWidth',5)
ylabel('$\rho$ [km]', 'Interpreter', 'latex','FontSize',40)
grid on

% TROLL station simulated measurements plot
subplot(3, 3, 2)
plot(t_vis_troll_sgp4, measT_noisy(:,1), '*', 'Color', [0 0.4470 0.7410],'LineWidth',5)
ylabel('Az [deg]', 'Interpreter', 'latex','FontSize',40)
title(TROLL.name, 'Interpreter', 'latex','FontSize',40)
grid on

subplot(3, 3, 5)
plot(t_vis_troll_sgp4, measT_noisy(:,2), '*', 'Color', [0.8500 0.3250 0.0980],'LineWidth',5)
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
grid on

subplot(3, 3, 8)
plot(t_vis_troll_sgp4, measT_noisy(:,3), '*', 'Color', 'k','LineWidth',5)
ylabel('$\rho$ [km]', 'Interpreter', 'latex','FontSize',40)
grid on

% SVALBARD station simulated measurements plot
subplot(3, 3, 3)
plot(t_vis_svalbard_sgp4, measS_noisy(:,1), '*', 'Color', [0 0.4470 0.7410],'LineWidth',5)
ylabel('Az [deg]', 'Interpreter', 'latex','FontSize',40)
title(SVALBARD.name, 'Interpreter', 'latex','FontSize',40)
grid on

subplot(3, 3, 6)
plot(t_vis_svalbard_sgp4, measS_noisy(:,2), '*', 'Color', [0.8500 0.3250 0.0980],'LineWidth',5)
ylabel('El [deg]', 'Interpreter', 'latex','FontSize',40)
grid on

subplot(3, 3, 9)
plot(t_vis_svalbard_sgp4, measS_noisy(:,3), '*', 'Color', 'k','LineWidth',5)
ylabel('$\rho$ [km]', 'Interpreter', 'latex','FontSize',40)
grid on


%% EX 2.3.a: Solve the navigation problem with only KOUROU

% Initial guess (x0 from R2BP propagation)
x0_guess = x0;

% Weights matrix
W_m = inv(sqrtm(R));

% Cost function
fun = @(x) costfunction(x, KOUROU, W_m, false);

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[x,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_guess, [], [], options);

% Difference with initial guess
disp('Only Kourou')
disp('Difference with initial guess (x - x_0) =');
disp((x - x0).');

% Difference with reality (sgp4)
[reci_sgp4, veci_sgp4, ~] = propagate_sgp4(et_ref, et0, KOUROU, et_ref, satrec);
x0_sgp4 = [reci_sgp4(:, end); veci_sgp4(:, end)];
disp('Difference with reality (x - x_sgp4) =');
disp((x - x0_sgp4).');

% Errors evaluation
err = x - x0_sgp4; 
err_r = norm(err(1:3)); 
err_v = norm(err(4:6)); 
disp('Position error = ');
disp(err_r); 
disp('Velocity error = ');
disp(err_v); 

% Covariance computation
B = (jacobian.'*jacobian)\eye(size(jacobian, 2));
P_ls = resnorm/(length(residual)-length(x)).*B;
disp('Covariance =');
disp(P_ls);

% Sqrt(tr(P)) computations
Pr = P_ls(1:3, 1:3); 
Pv = P_ls(4:6, 4:6);
sr_tr_Pr = sqrt(trace(Pr));
sr_tr_Pv = sqrt(trace(Pv));
disp('sqrt(tr(Pr)):')
disp(sr_tr_Pr)
disp('sqrt(tr(Pv)):')
disp(sr_tr_Pv)

% Sigma computation
sigma_vect = linear_mapping(x, P_ls, mu);
sigma_a = sigma_vect(1);
sigma_i = rad2deg(sigma_vect(3));
disp('Semi- major axis sigma [km]:')
disp(sigma_a)
disp('Inclination sigma [deg]:')
disp(sigma_i)


%% EX 2.3.b: Solve the navigation problem for all 3 ground stations

% Cost function
stations = [KOUROU; TROLL; SVALBARD];
fun = @(x) costfunction(x, stations, W_m, false);

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[x,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_guess, [], [], options);

% Print difference with initial guess
disp('All 3 stations')
disp('Difference with initial guess (x - x_0) =');
disp((x - x0).');

% Difference with reality (sgp4)
[reci_sgp4, veci_sgp4, ~] = propagate_sgp4(et_ref, et0, KOUROU, et_ref, satrec);
x0_sgp4 = [reci_sgp4(:, end); veci_sgp4(:, end)];
disp('Difference with reality (x - x_sgp4) =');
disp((x - x0_sgp4).');

% Errors evaluation
err = x - x0_sgp4; 
err_r = norm(err(1:3)); 
err_v = norm(err(4:6)); 
disp('Position error = ');
disp(err_r); 
disp('Velocity error = ');
disp(err_v); 

% Covariance computation
B = (jacobian.'*jacobian)\eye(size(jacobian, 2));
P_ls = resnorm/(length(residual)-length(x)).*B;
disp('Covariance =');
disp(P_ls);

% Sqrt(tr(P)) computations
Pr = P_ls(1:3, 1:3); 
Pv = P_ls(4:6, 4:6);
sr_tr_Pr = sqrt(trace(Pr));
sr_tr_Pv = sqrt(trace(Pv));
disp('sqrt(tr(Pr)):')
disp(sr_tr_Pr)
disp('sqrt(tr(Pv)):')
disp(sr_tr_Pv)

% Sigma computation
sigma_vect = linear_mapping(x, P_ls, mu);
sigma_a = sigma_vect(1);
sigma_i = rad2deg(sigma_vect(3));
disp('Semi- major axis sigma [km]:')
disp(sigma_a)
disp('Inclination sigma [deg]:')
disp(sigma_i)


%% EX 2.3.c: Solve the navigation problem for all 3 ground stations with J2 perturbation

% Initial guess (x0 from 2BP propagation with J2)
[x0_J2,~,~,~] = propagate_J2([et_ref et0],xref,mu);
x0_guess_J2 = x0_J2;

% Cost function
stations = [KOUROU; TROLL; SVALBARD];
fun = @(x) costfunction(x, stations, W_m, true);

% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
[x,resnorm,residual,~,~,~,jacobian] = lsqnonlin(fun, x0_guess_J2, [], [], options);

% Print difference with initial guess
disp('All 3 stations + J2')
disp('Difference with initial guess (x - x_0) =');
disp((x - x0).');

% Difference with reality (sgp4)
[reci_sgp4, veci_sgp4, ~] = propagate_sgp4(et_ref, et0, KOUROU, et_ref, satrec);
x0_sgp4 = [reci_sgp4(:, end); veci_sgp4(:, end)];
disp('Difference with reality (x - x_sgp4) =');
disp((x - x0_sgp4).');

% Errors evaluation
err = x - x0_sgp4; 
err_r = norm(err(1:3)); 
err_v = norm(err(4:6)); 
disp('Position error = ');
disp(err_r); 
disp('Velocity error = ');
disp(err_v); 

% Covariance computation
B = (jacobian.'*jacobian)\eye(size(jacobian, 2));
P_ls = resnorm/(length(residual)-length(x)).*B;
disp('Covariance =');
disp(P_ls);

% Sqrt(tr(P)) computations
Pr = P_ls(1:3, 1:3); 
Pv = P_ls(4:6, 4:6);
sr_tr_Pr = sqrt(trace(Pr));
sr_tr_Pv = sqrt(trace(Pv));
disp('sqrt(tr(Pr)):')
disp(sr_tr_Pr)
disp('sqrt(tr(Pv)):')
disp(sr_tr_Pv)

% Sigma computation
sigma_vect = linear_mapping(x, P_ls, mu);
sigma_a = sigma_vect(1);
sigma_i = rad2deg(sigma_vect(3));
disp('Semi- major axis sigma [km]:')
disp(sigma_a)
disp('Inclination sigma [deg]:')
disp(sigma_i)

% Covariance ellipsoid
[ellipsoid_x, ellipsoid_y, ellipsoid_z] = confidence_ellipsoid(P_ls(1:3,1:3), x(1:3), 3);

% Ellipsoid plot 
figure
surf(ellipsoid_x, ellipsoid_y, ellipsoid_z,'FaceColor','b','EdgeColor','k','FaceAlpha',0.3,'EdgeAlpha',0.3)
hold on
plot3(x(1),x(2),x(3),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',20);
plot3(x0_sgp4(1),x0_sgp4(2),x0_sgp4(3),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',20);
axis equal
text(x(1),x(2),x(3), '$\mathbf{\hat{r}}_0$', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',40);
text(x0_sgp4(1),x0_sgp4(2),x0_sgp4(3), '$\mathbf{r}_{0,sgp4}$', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Interpreter', 'latex','FontSize',40);
xlabel('X [km]','Interpreter','latex','FontSize',40)
ylabel('Y [km]','Interpreter','latex','FontSize',40)
zlabel('Z [km]','Interpreter','latex','FontSize',40)
%title('Covariance Ellipsoid','Interpreter','latex','FontSize',16)


%% EX 2.4: Trade-off analysis

% All combinations analysis
station_combinations = ["Kourou", "Troll", "Svalbard", "Kourou + Troll", "Kourou + Svalbard", "Troll + Svalbard", "KTS"];
[cost_K, sigma_a_K, sigma_i_K] = trade_off_analysis(KOUROU, x0_guess_J2, W_m);
[cost_T, sigma_a_T, sigma_i_T] = trade_off_analysis(TROLL, x0_guess_J2, W_m);
[cost_S, sigma_a_S, sigma_i_S] = trade_off_analysis(SVALBARD, x0_guess_J2, W_m);
[cost_KT, sigma_a_KT, sigma_i_KT] = trade_off_analysis([KOUROU;TROLL], x0_guess_J2, W_m);
[cost_KS, sigma_a_KS, sigma_i_KS] = trade_off_analysis([KOUROU;SVALBARD], x0_guess_J2, W_m);
[cost_TS, sigma_a_TS, sigma_i_TS] = trade_off_analysis([TROLL;SVALBARD], x0_guess_J2, W_m);
[cost_KTS, sigma_a_KTS, sigma_i_KTS] = trade_off_analysis([KOUROU;TROLL;SVALBARD], x0_guess_J2, W_m);

% Costs
costs = [cost_K, cost_T, cost_S, cost_KT, cost_KS, cost_TS, cost_KTS];

% Sigma combinations
sigma_a = [sigma_a_K, sigma_a_T, sigma_a_S, sigma_a_KT, sigma_a_KS, sigma_a_TS, sigma_a_KTS];
sigma_i = [sigma_i_K, sigma_i_T, sigma_i_S, sigma_i_KT, sigma_i_KS, sigma_i_TS, sigma_i_KTS];

% Normalization
sma = elts(1)/(1-elts(2));
inc = elts(3)*cspice_dpr();
sigma_a_ad = sigma_a./sma;
sigma_i_ad = sigma_i./inc;

% Sigma comparison plot
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
figure
loglog(sigma_a_ad(1),sigma_i_ad(1),'o','MarkerEdgeColor','k','MarkerFaceColor',colors(1,:),'MarkerSize',20)
hold on
for i=2:length(sigma_a)-1
    loglog(sigma_a_ad(i),sigma_i_ad(i),'o','MarkerEdgeColor','k','MarkerFaceColor',colors(i,:),'MarkerSize',20)
end
grid on
xlabel('$\frac{\sigma_a}{SMA}$','Interpreter','latex','FontSize',40)
ylabel('$\frac{\sigma_i}{INC}$','Interpreter','latex','FontSize',40)
legend('Kourou','Troll','Svalbard','Kourou+Troll','Kourou+Svalbard','Troll+Svalbard','Interpreter','latex','FontSize',40)
%title('Adimensional Standard Deviations Comparison','Interpreter','latex','FontSize',16)

% Trade-off able
results_table = table(station_combinations', costs', sigma_a', sigma_i', ...
    'VariableNames', {'Stations', 'Cost (€)', 'Sigma_a (km)', 'Sigma_i (deg)'});
disp(results_table);


%% EX 2.5: Long-term analysis

% Set final istant 
etf = cspice_str2et('2024-11-20T20:30:00.000');

% Equally time grid vector, depending on the frequency of measurement acquisition
npointsK = round((etf-et0)/KOUROU.frequency)+1;
et_vectK = linspace(et0, etf, npointsK);
npointsT = round((etf-et0)/TROLL.frequency)+1;
et_vectT = linspace(et0, etf, npointsT);
et_vectS = et_vectK;

% Propagation with J2 effect
[~,~,xx_K,~] = propagate_J2(et_vectK,x0_J2,mu);
[~,~,xx_T,~] = propagate_J2(et_vectT,x0_J2,mu);
[~,~,xx_S,~] = propagate_J2(et_vectS,x0_J2,mu);

% Measurements
[azimuthK, elevationK, rangeK] = pointing_antenna(KOUROU, et_vectK, xx_K(:,1:3));
[azimuthT, elevationT, rangeT] = pointing_antenna(TROLL, et_vectT, xx_T(:,1:3));
[azimuthS, elevationS, rangeS] = pointing_antenna(SVALBARD, et_vectS, xx_S(:,1:3));

% Visibility
[elevationK,azimuthK,rangeK,et_vect_newK] = visibility(elevationK,azimuthK,rangeK,et_vectK,KOUROU);
t_vis_K = datetime(cspice_timout(et_vect_newK, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

[elevationT,azimuthT,rangeT,et_vect_newT] = visibility(elevationT,azimuthT,rangeT,et_vectT,TROLL);
t_vis_T = datetime(cspice_timout(et_vect_newT, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

[elevationS,azimuthS,rangeS,et_vect_newS] = visibility(elevationS,azimuthS,rangeS,et_vectS,SVALBARD);
t_vis_S = datetime(cspice_timout(et_vect_newS, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Long-term plot
stations_names = ["Troll", "Kourou", "Svalbard"];
y_coords = [1, 2, 3]; 
figure;
hold on;
scatter(t_vis_K, y_coords(2)*ones(size(t_vis_K)),200, 'o', 'DisplayName', 'Kourou', 'MarkerEdgeColor', colors(1,:), 'MarkerFaceColor', colors(1,:))
scatter(t_vis_T, y_coords(1)*ones(size(t_vis_T)),200, 'o', 'DisplayName', 'Troll', 'MarkerEdgeColor', colors(2,:), 'MarkerFaceColor', colors(2,:))
scatter(t_vis_S, y_coords(3)*ones(size(t_vis_S)),200, 'o', 'DisplayName', 'Svalbard', 'MarkerEdgeColor', colors(3,:), 'MarkerFaceColor', colors(3,:))
yticks(y_coords);
yticklabels(stations_names);
xtickangle(45); 
ytickangle(45);
ylim([0 4])
grid on;
%title('Visibility windows analysis in two days from $t_0$','Interpreter','latex','FontSize',16)

% Ground track of one day 
et0 = cspice_str2et('2024-11-18T00:00:00.000');
etf = cspice_str2et('2024-11-18T23:59:00.000');
n_points = 10000;
et_vect = linspace(et0,etf,n_points);
reci_smos = NaN(3,n_points);
veci_smos = NaN(3,n_points);
for i = 1:n_points

    % Propagation via sgp4
    tsince = (et_vect(i) - et_ref)/60; % minutes from TLE epoch
    [~,rteme_gt,vteme_gt] = sgp4(satrec,  tsince);
    
    % Compute centuries from 2000-01-01T00:00:00.00 TDT
    ttt_gt = cspice_unitim(et_vect(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % Conversion to ECI
    ateme = [0.0;0.0;0.0];
    [reci_smos(:,i), veci_smos(:,i)] = teme2eci(rteme_gt, vteme_gt, ateme,  ttt_gt, ddpsi, ddeps); 

end
xx_smos = [reci_smos;veci_smos];

% IAU conversion
R = cspice_pxform('J2000','IAU_EARTH', et_vect);
r_iau = NaN(size(reci_smos));
for k=1:length(et_vect)
    r_iau(:,k) = R(:,:,k)*reci_smos(:,k);
end

% Coordinates computation
Radii = cspice_bodvrd('EARTH','RADII',3);
flatness = (Radii(1) - Radii(3))/Radii(1);
[lon,lat,~] = cspice_recgeo(r_iau, Radii(1), flatness);
lon = rad2deg(lon);
lat = rad2deg(lat);

% Ground track plot
Earth_image = flip(imread('Earth_NASA.jpg')); 
figure
hold on
image([-180 180],[-90 90], Earth_image)
for k=1:length(lon)
plot(lon(k),lat(k),'.','MarkerEdgeColor',colors(2,:),'MarkerSize',20)
end
plot(-52.80466, 5.25144, 'o', 'MarkerSize',20, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
plot(2.536103, -72.011977, 'o', 'MarkerSize',20, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
plot(15.407786, 78.229772, 'o', 'MarkerSize',20, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')
text(-52.80466, 5.25144, 'KOUROU', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',40)
text(2.536103, -72.011977, 'TROLL', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',40)
text(15.407786, 78.229772, 'SVALBARD', 'Color', 'k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',40)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]', 'Interpreter', 'latex','FontSize',40)
ylabel('Latitude [deg]', 'Interpreter', 'latex','FontSize',40)
%title('SMOS Ground Track', 'Interpreter', 'latex','FontSize',16)
grid on


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

function dy = two_body_rhs_J2(et, y, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the right-hand side of the two-body problem including the
%   J2 perturbation. The gravitational acceleration is computed in the
%   Earth-Centered Inertial (ECI, J2000) frame, with the J2 effect
%   calculated in the Earth-Centered Earth-Fixed (ECEF, ITRF93) frame
%   and then rotated back to ECI.
%
% Inputs:
%   et - Ephemeris time [s past J2000 TDB], used for frame transformation
%   y  - State vector [6×1]:
%           y(1:3) = position vector in ECI [km]
%           y(4:6) = velocity vector in ECI [km/s]
%   mu - Earth’s gravitational parameter [km^3/s^2]
%
% Outputs:
%   dy - Time derivative of state [6×1]:
%           dy(1:3) = velocity in ECI [km/s]
%           dy(4:6) = acceleration in ECI [km/s²]
%--------------------------------------------------------------------------

% Frame transformation: ECI to ECEF 
rotm   = cspice_pxform('J2000', 'ITRF93', et);
r_eci  = y(1:3);
r_ecef = rotm * r_eci;
rr     = norm(r_ecef);

%  J2 perturbation acceleration in ECEF 
R_E = 6378.1366;             % Earth mean equatorial radius [km]
J2  = 0.0010826269;          % J2 zonal harmonic coefficient

a_J2_ecef = (3/2) * mu * J2 * (R_E/rr)^2 / rr^3 .* r_ecef ...
            .* (5 * (r_ecef(3)/rr)^2 - [1;1;3]);

% Rotate J2 acceleration back to ECI
rotm      = cspice_pxform('ITRF93', 'J2000', et);
a_J2_eci  = rotm * a_J2_ecef;

% Two-body central acceleration in ECI 
r = norm(r_eci);
a_two_body = -(mu/r^3) * r_eci;

% Assemble state derivative 
dy = [ y(4:6);                   % Velocity derivatives = current velocity
       a_two_body + a_J2_eci ];  % Acceleration = 2-body + J2 perturbation

end

function [xf, tf, xx, tt] = propagate_J2(t_vect, x0, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the dynamics of a spacecraft under the two-body problem with
%   J2 perturbation using a numerical ODE solver (ode78). The state includes
%   Cartesian position and velocity components in the Earth-Centered Inertial
%   (ECI, J2000) frame.
%
% Inputs:
%   t_vect - Time vector for integration [1×N]
%   x0     - Initial state vector [6×1]:
%              [r0x; r0y; r0z; v0x; v0y; v0z] in ECI [km, km/s]
%   mu     - Earth gravitational parameter [km^3/s^2]
%
% Outputs:
%   xf - Final state vector at t = t_vect(end) [6×1]
%   tf - Final propagation time [scalar]
%   xx - Full propagated trajectory [N×6]
%   tt - Time vector returned by integrator [N×1]
%--------------------------------------------------------------------------
    
% Integration options (tight tolerances)
options = odeset('reltol', 1e-13, 'abstol', 1e-20);

% Integrate dynamics with J2 perturbation
[tt, xx] = ode78(@(t,x) two_body_rhs_J2(t,x,mu), t_vect, x0, options);

% Extract final state and time
xf = xx(end,:)';
tf = tt(end);

end

function [reci_smos, veci_smos, et_vec] = propagate_sgp4(et0, etf, station, sat_epoch_et, satrec)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates a satellite orbit using the SGP4 model within a given
%   ephemeris time window. The output position and velocity are provided
%   in the Earth-Centered Inertial (ECI, J2000) frame after TEME-to-ECI
%   conversion, including precession-nutation corrections.
%
% Inputs:
%   et0          - Initial ephemeris time of the propagation window [s past J2000 TDB]
%   etf          - Final ephemeris time of the propagation window [s past J2000 TDB]
%   station      - Structure containing ground station information
%                    .frequency = sampling frequency [Hz]
%   sat_epoch_et - TLE reference epoch in ephemeris time [s past J2000 TDB]
%   satrec       - Satellite record structure (initialized from TLE, used by sgp4.mex)
%
% Outputs:
%   reci_smos - Propagated position vectors in ECI frame [3×N] (km)
%   veci_smos - Propagated velocity vectors in ECI frame [3×N] (km/s)
%   et_vec    - Ephemeris time vector of propagation instants [1×N] (s past J2000 TDB)
%--------------------------------------------------------------------------

% Precession-nutation corrections (IAU2000A, example date 18-Nov-2024)
arcsec2rad = pi / (180*3600);
ddpsi = -0.114752 * arcsec2rad;   % Δψ correction [rad]
ddeps = -0.007529 * arcsec2rad;   % Δε correction [rad]

% Time vector setup 
f        = station.frequency; 
n_points = round((etf - et0) / f) + 1;
et_vec   = linspace(et0, etf, n_points);

% Preallocate state arrays
reci_smos = NaN(3, n_points);
veci_smos = NaN(3, n_points);

% Propagation loop
for i = 1:n_points
    
    % Time since TLE epoch [minutes]
    tsince = (et_vec(i) - sat_epoch_et) / 60;  
    [~, rteme, vteme] = sgp4(satrec, tsince);
    
    % Julian centuries of TDT since J2000
    ttt = cspice_unitim(et_vec(i), 'ET', 'TDT') / cspice_jyear() / 100;
    
    % TEME → ECI (J2000) conversion (accounting for nutation corrections)
    ateme = [0.0; 0.0; 0.0]; % acceleration (not used in this context)
    [reci_smos(:,i), veci_smos(:,i)] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);
end

end

function [azimuth, elevation, range] = pointing_antenna(station_struct, et_vect, r_eci)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the azimuth, elevation, and range of a satellite as seen from a
%   ground station, given the satellite ECI positions and the ephemeris time
%   vector. Transformations are performed using SPICE kernels from the 
%   inertial J2000 frame to the local TOPO frame of the station.
%
% Inputs:
%   station_struct - Structure containing station information:
%                      .name = station identifier (SPICE name string)
%   et_vect        - Ephemeris time vector [N×1] (s past J2000 TDB)
%   r_eci          - Satellite positions in ECI [N×3] (km)
%
% Outputs:
%   azimuth   - Azimuth angles [N×1] (deg)
%   elevation - Elevation angles [N×1] (deg)
%   range     - Slant ranges from station to satellite [N×1] (km)
%--------------------------------------------------------------------------

% Initialization
azimuth   = NaN(length(et_vect),1);
elevation = NaN(length(et_vect),1); 
range     = NaN(length(et_vect),1);

% Build station TOPO frame name
stationName_TOPO = strcat(station_struct.name, '_TOPO'); 

for ii = 1:length(et_vect)

    % Transform satellite position from ECI to station TOPO frame
    xform_matrix = cspice_pxform('J2000', stationName_TOPO, et_vect(ii));
    r_topo = xform_matrix * r_eci(ii,:)';

    % Station state in TOPO frame
    [state_topo, ~] = cspice_spkezr(station_struct.name, et_vect(ii), stationName_TOPO, 'NONE', 'EARTH');
    r_gs_topo = state_topo(1:3);

    % Relative position
    r_rel_topo = r_topo - r_gs_topo;

    % Range (km)
    range(ii) = norm(r_rel_topo);
    
    % Components
    x = r_rel_topo(1); 
    y = r_rel_topo(2); 
    z = r_rel_topo(3);
    
    % Azimuth (deg)
    azimuth(ii) = rad2deg(wrapToPi(atan2(y, x))); 

    % Elevation (deg)
    elevation(ii) = rad2deg(asin(z / range(ii))); 

end

end

function [elevation, azimuth, range, et_vect_new] = visibility(elevation, azimuth, range, et_vect, station_struct)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Filters the azimuth, elevation, range, and time vectors to return only
%   the instants when the satellite is visible from the ground station,
%   i.e. when the elevation is above the minimum elevation threshold.
%
% Inputs:
%   elevation      - Elevation angles [N×1] (deg)
%   azimuth        - Azimuth angles [N×1] (deg)
%   range          - Slant ranges [N×1] (km)
%   et_vect        - Ephemeris time vector [N×1] (s past J2000 TDB)
%   station_struct - Structure containing station information:
%                       .el_min = minimum elevation for visibility [deg]
%
% Outputs:
%   elevation  - Filtered elevation angles [M×1] (deg), M ≤ N
%   azimuth    - Filtered azimuth angles [M×1] (deg)
%   range      - Filtered slant ranges [M×1] (km)
%   et_vect_new- Filtered ephemeris time vector [M×1] (s past J2000 TDB)
%--------------------------------------------------------------------------

% Visibility mask: elevation above threshold
visibility_mask = elevation > station_struct.el_min;

% Apply filter
elevation   = elevation(visibility_mask);
azimuth     = azimuth(visibility_mask);
range       = range(visibility_mask);
et_vect_new = et_vect(visibility_mask);

end

function [meas_noisy, et_vect_new] = visibility_noisy(meas_noisy, et_vect, station_struct)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Filters noisy measurement data and corresponding ephemeris time vector
%   to return only the instants when the satellite is visible from the
%   ground station, i.e. when the elevation angle is above the minimum
%   elevation threshold.
%
% Inputs:
%   meas_noisy     - Matrix of noisy measurements [N×m]
%                      meas_noisy(:,2) = elevation (deg)
%   et_vect        - Ephemeris time vector [N×1] (s past J2000 TDB)
%   station_struct - Structure containing station information:
%                       .el_min = minimum elevation for visibility [deg]
%
% Outputs:
%   meas_noisy  - Filtered noisy measurements [M×m], M ≤ N
%   et_vect_new - Filtered ephemeris time vector [M×1] (s past J2000 TDB)
%--------------------------------------------------------------------------

% Visibility mask: elevation above threshold
visibility_mask = meas_noisy(:,2) > station_struct.el_min;

% Apply filter
meas_noisy  = meas_noisy(visibility_mask,:);
et_vect_new = et_vect(visibility_mask);

end

function residuals = costfunction(x, stations, W_m, J2)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the measurement residuals for orbit determination by comparing
%   predicted azimuth/elevation/range from propagated states against noisy
%   ground-station measurements. Supports dynamics with or without J2.
%
% Inputs:
%   x        - Initial state vector [6×1] (ECI J2000) [km, km/s]
%   stations - Array of station structs with fields:
%                .name        : SPICE target name (string)
%                .meas_noisy  : measured [az(deg), el(deg), range(km)] [N×3]
%                .tspan       : ephemeris time vector [N×1] (s past J2000 TDB)
%                .el_min      : minimum elevation for visibility [deg]
%   W_m      - 3×3 weighting matrix for [az, el, range]
%   J2       - Logical flag (true = include J2 perturbation, false = 2-body)
%
% Outputs:
%   residuals - Stacked residuals to minimize (3×M):
%                 rows:   [Δaz(rad); Δel(rad); Δrange(km)]
%                 cols:   all measurement epochs across stations
%--------------------------------------------------------------------------

mu = cspice_bodvrd('Earth', 'GM', 1);
options = odeset('Reltol',1e-13,'Abstol',1e-20);

% J2 flag
if J2
    odefun = @(t, x) two_body_rhs_J2(t, x, mu);
else
    odefun = @(t, x) two_body_rhs(t, x, mu);
end

n_stations = length(stations);
n_res = 0;
for kk = 1:n_stations
    n_res = n_res + size(stations(kk).meas_noisy', 2);
end
residuals = NaN(3, n_res);
nprec = 1;

for kk = 1:n_stations
    station = stations(kk);

    % Measurements (as 3×N) and times
    Y     = station.meas_noisy';
    tspan = station.tspan;
    n     = size(Y, 2);

    % Propagate dynamics
    [~, x_prop] = ode78(odefun, tspan, x, options);
    reci_vec    = x_prop(:, 1:3);                 % [N×3] ECI positions

    % Predicted az/el/range at station
    [azimuth, elevation, range] = pointing_antenna(station, tspan, reci_vec);
    Y_pred = [azimuth, elevation, range]';
    Y_pred = Y_pred(:, 2:end);                    % align with central differencing below

    % Build weighted residuals (use epochs 1..N-1)
    res = zeros(size(Y));
    for k = 1:length(tspan)-1
        dAz = angdiff(deg2rad(Y_pred(1, k)), deg2rad(Y(1, k)));  % radians
        dEl = angdiff(deg2rad(Y_pred(2, k)), deg2rad(Y(2, k)));  % radians
        dr  = Y_pred(3, k) - Y(3, k);                            % km
        res(:, k) = W_m * [dAz; dEl; dr];
    end

    % Stack residuals for this station
    residuals(:, nprec:nprec + n - 1) = res;
    nprec = nprec + n;
end

end

function sigma_vec = linear_mapping(xx_car, P_car, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the standard deviations of the Keplerian orbital elements
%   starting from the Cartesian state covariance matrix. The mapping is
%   performed via the Jacobian of the Cartesian-to-Keplerian transformation.
%
% Inputs:
%   xx_car - Cartesian state vector [6×1] [km, km/s]
%              [rx; ry; rz; vx; vy; vz]
%   P_car  - Cartesian covariance matrix [6×6] [km^2, km·km/s, (km/s)^2]
%   mu     - Gravitational parameter [km^3/s^2]
%
% Outputs:
%   sigma_vec - Standard deviations of the Keplerian elements [6×1]
%                 [a; e; i; Ω; ω; ν] (units: [km], [-], [rad], [rad], [rad], [rad])
%--------------------------------------------------------------------------

% Jacobian of Cartesian → Keplerian transformation
J = jacobian_car2kep(xx_car, mu);

% Propagate covariance into Keplerian space
P_kep = J * P_car * J.';

% Extract standard deviations
sigma_vec = sqrt(diag(P_kep));

end

function kep = car2kep(xx, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Converts a Cartesian state vector into Keplerian orbital elements.
%   The returned elements follow the classical (osculating) definition:
%     - a  : semi-major axis [km]
%     - e  : eccentricity [-]
%     - i  : inclination [rad]
%     - Ω  : right ascension of ascending node (RAAN) [rad]
%     - ω  : argument of periapsis [rad]
%     - θ  : true anomaly [rad]
%
% Inputs:
%   xx - Cartesian state vector [6×1] [km, km/s]
%          xx(1:3) = position [rx; ry; rz]
%          xx(4:6) = velocity [vx; vy; vz]
%   mu - Gravitational parameter [km^3/s^2]
%
% Outputs:
%   kep - Vector of Keplerian orbital elements [6×1]:
%            [a; e; i; Ω; ω; θ]
%--------------------------------------------------------------------------

% Split state into position and velocity
rr = xx(1:3);
vv = xx(4:6);

% Ensure column vectors
if size(rr,1) == 1, rr = rr.'; end
if size(vv,1) == 1, vv = vv.'; end

% Magnitudes
r = norm(rr);
v = norm(vv);

% Semi-major axis
a = -mu / (v^2 - 2*mu/r);

% Angular momentum vector
hh = cross(rr, vv);
h  = norm(hh);

% Eccentricity vector
ee = cross(vv, hh)/mu - rr/r;
e  = norm(ee);
if e == 0
    ee = [1; 0; 0]; % default direction if circular
end

% Inclination
i = acos(hh(3)/h);

% Node vector
k = [0; 0; 1];
N = cross(k, hh);
if norm(N) ~= 0
    N = N / norm(N);
else
    % Equatorial orbit: define arbitrary node vector
    N = [1; 0; 0];
    i = eps;
end

% Right ascension of ascending node (RAAN)
Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

% Argument of periapsis
om = acos(dot(N, ee) / e);
if ee(3) < 0
    om = 2*pi - om;
end

% True anomaly
theta = acos(dot(rr, ee) / (r*e));
if dot(rr, vv) < 0
    theta = 2*pi - theta;
end

% Collect elements
kep = [a; e; i; Om; om; theta];

end

function J = jacobian_car2kep(xx, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Estimates the Jacobian of the Keplerian elements with respect to the
%   Cartesian state via a centered finite-difference scheme.
%
% Inputs:
%   xx - Cartesian state vector [6×1]  (km, km/s)
%   mu - Gravitational parameter [km^3/s^2]
%
% Outputs:
%   J  - Jacobian matrix ∂[a,e,i,Ω,ω,θ]/∂[rx,ry,rz,vx,vy,vz]  [6×6]
%--------------------------------------------------------------------------

n = length(xx);
J = NaN(n);

for k = 1:n
    pert = zeros(size(xx));
    ek   = sqrt(eps) * max(1, abs(xx(k)));   % Step size (scaled, centered FD)
    pert(k) = ek;

    xx_plus  = xx + pert;
    xx_minus = xx - pert;

    K_plus  = car2kep(xx_plus,  mu);
    K_minus = car2kep(xx_minus, mu);

    J(:, k) = (K_plus - K_minus) / (2*ek);   % Column k of Jacobian
end

end

function [ellipsoid_x, ellipsoid_y, ellipsoid_z] = confidence_ellipsoid(P, mu, n_sigma)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the 3D confidence ellipsoid corresponding to a Gaussian
%   distribution with mean and covariance. The ellipsoid is scaled to
%   represent the n-sigma confidence region.
%
% Inputs:
%   P       - Covariance matrix [3×3]
%   mu      - Mean vector [3×1]
%   n_sigma - Confidence scaling factor (e.g. 1, 2, 3 for 1σ, 2σ, 3σ)
%
% Outputs:
%   ellipsoid_x - X-coordinates of the ellipsoid surface [31×31]
%   ellipsoid_y - Y-coordinates of the ellipsoid surface [31×31]
%   ellipsoid_z - Z-coordinates of the ellipsoid surface [31×31]
%--------------------------------------------------------------------------

% Generate unit sphere mesh
[xe, ye, ze] = sphere(30);
points = [xe(:), ye(:), ze(:)]';

% Singular Value Decomposition of covariance
[R, D] = svd(P);      % R = rotation, D = eigenvalue diagonal
d = sqrt(diag(D));    % Standard deviations along principal axes

% Transform sphere → ellipsoid
tr_pts = n_sigma * R * diag(d) * points;

% Reshape into grid
X = reshape(tr_pts(1,:), size(xe));
Y = reshape(tr_pts(2,:), size(ye));
Z = reshape(tr_pts(3,:), size(ze));

% Translate to mean
ellipsoid_x = mu(1) + X;
ellipsoid_y = mu(2) + Y;
ellipsoid_z = mu(3) + Z;

end

function [cost, sigma_a, sigma_i] = trade_off_analysis(stations, x0_guess, W_m)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Performs a trade-off analysis for ground-station selection in orbit
%   determination. The function estimates the state covariance using a
%   least-squares fit to station measurements, maps it into Keplerian
%   elements, and returns the standard deviations of SMA and inclination
%   along with the total cost of the station configuration.
%
% Inputs:
%   stations  - Array of station structs, each with fields:
%                 .meas_noisy : noisy measurements [N×3]
%                 .tspan      : ephemeris time vector [N×1]
%                 .el_min     : minimum elevation [deg]
%                 .cost       : station cost [€]
%   x0_guess  - Initial Cartesian state guess [6×1] [km, km/s]
%   W_m       - Weighting matrix for measurement residuals [3×3]
%
% Outputs:
%   cost    - Total cost of selected ground stations [€]
%   sigma_a - Standard deviation of semi-major axis [km]
%   sigma_i - Standard deviation of inclination [deg]
%--------------------------------------------------------------------------

% Earth gravitational parameter
mu = cspice_bodvrd('Earth', 'GM', 1);

% Residual function for nonlinear least squares
fun = @(x) costfunction(x, stations, W_m, true);

% Solve least-squares orbit determination
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');
[x, resnorm, residual, ~, ~, ~, Jac] = lsqnonlin(fun, x0_guess, [], [], options);

% Estimate covariance from Jacobian
B = (Jac.' * Jac) \ eye(size(Jac, 2));
P = resnorm / (length(residual) - length(x)) * B;

% Map Cartesian covariance into Keplerian space
sigma_vec = linear_mapping(x, P, mu);
sigma_a   = sigma_vec(1);           % SMA std deviation [km]
sigma_i   = rad2deg(sigma_vec(3));  % Inclination std deviation [deg]

% Total station cost
cost = 0;
for k = 1:length(stations)
    cost = cost + stations(k).cost;
end

end

