clear all

% Select datafile to load
addpath('simdata2026');
dataname = "da3211_1.mat";

% UNCOMMENT FOR FILE CHOICE
% dataname = uigetfile('.mat');

try
    load(dataname);
catch
    Error('Error loading data file');
end

%%
close all;

% ===== AIRCRAFT DATA =====

Ixx = 15351; % kgm2
Iyy = 22965; % kgm2
Izz = 36220; % kgm2
Ixz = 1908; % kgm2
m = 5000; % kg
b = 13.3250; % m, wingspan
S = 24.9900; % m2, wing area
c = 1.9910; % m, MAC

% ===== UNCERTAINTY MODELLING =====

% IMU bias
lambdax = 0.02; % m/s2
lambday = 0.02; % m/s2
lambdaz = 0.03; % m/s2
lambdap = 0.005; % deg/s
lambdaq = 0.005; % deg/s
lambdar = 0.002; % deg/s

% IMU noise
wx = lambdax * randn(length(t), 1);
wy = lambday * randn(length(t), 1);
wz = lambdaz * randn(length(t), 1);
wp = lambdap * randn(length(t), 1);
wq = lambdaq * randn(length(t), 1);
wr = lambdar * randn(length(t), 1);

% Wind
Wxe = 6; % m/s
Wye = 2; % m/s
Wze = 10; % m/s

% GPS noise
% Position standard deviation
sigmaPos = 1; % m
% Velocity standard deviation
sigmaVel = 0.01; % m/s
% Attitude standard deviation
sigmaAtt = 0.04; % deg

v_gps_pos = sigmaPos * randn(length(t), 3);
v_gps_vel = sigmaVel * randn(length(t), 3);
v_gps_att = sigmaAtt * randn(length(t), 3);

% Other measurments noise

sigmaVTAS = 0.2; % m/s
sigmaAoA = 0.1; % deg
sigmaSideSlip = 0.25; %deg 

v_vtas = sigmaVTAS * randn(length(t), 1);
v_alpha = sigmaAoA * randn(length(t), 1);
v_beta = sigmaSideSlip * randn(length(t), 1);

% ===== COMPUTE AND PLOT THE ACTUAL TRAJECTORY =====

% Velocities in specific directions
xdot = (u_n.*cos(theta)+(v_n.*sin(phi)+w_n.*cos(phi)).*sin(theta)).*cos(psi)-(v_n.*cos(phi)-w_n.*sin(phi)).*sin(psi) + Wxe;
ydot = (u_n.*cos(theta)+(v_n.*sin(phi)+w_n.*cos(phi)).*sin(theta)).*sin(psi)+(v_n.*cos(phi)-w_n.*sin(phi)).*cos(psi) + Wye;
zdot = -u_n.*sin(theta)+(v_n.*sin(phi)+w_n.*cos(phi)).*cos(theta) + Wze;

xdot_interp = @(tq) interp1(t, xdot, tq, 'pchip', 'extrap');
ydot_interp = @(tq) interp1(t, ydot, tq, 'pchip', 'extrap');
zdot_interp = @(tq) interp1(t, zdot, tq, 'pchip', 'extrap');

x0 = 0;
y0 = 0;
z0 = 0;

x = transpose(x0 + arrayfun(@(T) integral(xdot_interp, t(1), T), t));
y = transpose(y0 + arrayfun(@(T) integral(ydot_interp, t(1), T), t));
z = transpose(z0 + arrayfun(@(T) integral(zdot_interp, t(1), T), t));

% Plot the computed trajectories
figure;
plot3(x, y, z);
xlabel('X Position');
ylabel('Y Position');
zlabel('Z Position');
title('3D Trajectory for file ' + dataname, 'Interpreter', 'none');
grid on;

% ===== APPLY UNCERTAINTY =====

% IMU measurments
Axm = Ax + lambdax + wx;
Aym = Ay + lambday + wy;
Azm = Az + lambdaz + wz;
pm = p + lambdap + wp;
qm = q + lambdaq + wq;
rm = r + lambdar + wr;

% GPS measurments
xm = x + v_gps_pos(:, 1);
ym = y + v_gps_pos(:, 2);
zm = z + v_gps_pos(:, 3);
xdotm = xdot + v_gps_vel(:, 1);
ydotm = ydot + v_gps_vel(:, 2);
zdotm = zdot + v_gps_vel(:, 3);
phim = phi + v_gps_att(:, 1);
thetam = theta + v_gps_att(:, 2);
psim = psi + v_gps_att(:, 3);

% Other measurments
vtasm = vtas + v_vtas;
alpham = alpha + v_alpha;
betam = beta + v_beta;

% Plot clean vs noise trajectory comparison
t0 = 1;
t1 = 2;

idx = t>= t0 & t <= t1;

figure;
plot3(xm(idx), ym(idx), zm(idx), 'r--');
hold on;
plot3(x(idx), y(idx), z(idx), 'b-');
xlabel('X Position');
ylabel('Y Position');
zlabel('Z Position');
title(sprintf('Noisy vs Clean Trajectory from t %.2f to %.2f s for file %s', t0, t1, dataname), 'Interpreter', 'none');
legend('Noisy Trajectory', 'Clean Trajectory');
grid on;