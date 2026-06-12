classdef aircraft_dynamics

    properties (Constant)

        % ===== AIRCRAFT DATA =====

        Ixx = 15351; % kgm2
        Iyy = 22965; % kgm2
        Izz = 36220; % kgm2
        Ixz = 1908; % kgm2
        m = 5000; % kg
        b = 13.3250; % m, wingspan
        S = 24.9900; % m2, wing area
        c = 1.9910; % m, MAC

        % ===== TRUE IMU SENSOR BIASES AND DEVIATIONS =====
        lambda_x_true = 0.02;                      % m/s^2
        lambda_y_true = 0.02;                      % m/s^2
        lambda_z_true = 0.03;                      % m/s^2
        lambda_p_true = 0.005 * (pi / 180);        % rad/s
        lambda_q_true = 0.005 * (pi / 180);        % rad/s
        lambda_r_true = 0.002 * (pi / 180);        % rad/s

        % ===== OTHER DEVIATIONS =====
        % Position standard deviation
        sigmaPos = 1; % m
        % Velocity standard deviation
        sigmaVel = 0.01; % m/s
        % Attitude standard deviation
        sigmaAtt = 0.04; % deg
        sigmaVTAS = 0.2; % m/s
        sigmaAoA = 0.1; % deg
        sigmaSideSlip = 0.25; %deg 

    end

    methods (Static)

        function [X_true, U_meas, Z_meas] = preprocess_flight_data(data, W_true)
            % PREPROCESS_FLIGHT_DATA Prepares all vectors for the EKF loop
            % INPUTS:
            %   data   : The structure loaded from the .mat file (e.g. load('da3211_1.mat'))
            %   W_true : 3x1 vector of true wind components [Wxe; Wye; Wze]
            % OUTPUTS:
            %   X_true : (N x 18) Matrix of clean, true states over time
            %   U_meas : (N x 6)  Matrix of noisy IMU inputs fed to the EKF
            %   Z_meas : (N x 12) Matrix of noisy GPS/Airdata measurements
            
            N = length(data.t);
            dt = data.t(2) - data.t(1);
            
            % --- 1. Generate True Trajectory (Part 1, Q1) ---
            % Calculate navigation-frame velocity derivatives using true data
            xdot = (data.u_n.*cos(data.theta) + (data.v_n.*sin(data.phi) + data.w_n.*cos(data.phi)).*sin(data.theta)).*cos(data.psi) - (data.v_n.*cos(data.phi) - data.w_n.*sin(data.phi)).*sin(data.psi) + W_true(1);
            ydot = (data.u_n.*cos(data.theta) + (data.v_n.*sin(data.phi) + data.w_n.*cos(data.phi)).*sin(data.theta)).*sin(data.psi) + (data.v_n.*cos(data.phi) - data.w_n.*sin(data.phi)).*cos(data.psi) + W_true(2);
            zdot = -data.u_n.*sin(data.theta) + (data.v_n.*sin(data.phi) + data.w_n.*cos(data.phi)).*cos(data.theta) + W_true(3);
            
            % Integrate to find true positions
            x_pos = cumtrapz(data.t, xdot);
            y_pos = cumtrapz(data.t, ydot);
            z_pos = cumtrapz(data.t, zdot);
            
            % Assemble the True State Matrix (N x 18)
            X_true = [x_pos, y_pos, z_pos, ...                     % 1:3 Positions
                      data.u_n, data.v_n, data.w_n, ...            % 4:6 Body Velocities
                      data.phi, data.theta, data.psi, ...          % 7:9 Attitudes
                      zeros(N, 3), ...                             % 10:12 True Accel Biases (0)
                      zeros(N, 3), ...                             % 13:15 True Gyro Biases (0)
                      repmat(W_true', N, 1)];                       % 16:18 True Wind States
                  
            % --- 2. Generate Noisy IMU Inputs (Part 1, Q2) ---
            
            % Generate random IMU noise using your existing class function
            imu_noise = aircraft_dynamics.generate_imu_noise(N);
            
            % Corrupt clean IMU data with biases and noise
            Ax_m = data.Ax + aircraft_dynamics.lambda_x_true + imu_noise(:, 1);
            Ay_m = data.Ay + aircraft_dynamics.lambda_y_true + imu_noise(:, 2);
            Az_m = data.Az + aircraft_dynamics.lambda_z_true + imu_noise(:, 3);
            
            p_m  = data.p  + aircraft_dynamics.lambda_p_true + imu_noise(:, 4);
            q_m  = data.q  + aircraft_dynamics.lambda_q_true + imu_noise(:, 5);
            r_m  = data.r  + aircraft_dynamics.lambda_r_true + imu_noise(:, 6);
            
            U_meas = [Ax_m, Ay_m, Az_m, p_m, q_m, r_m];
            
            % --- 3. Generate Noisy Sensor Measurements (Part 1, Q2) ---
            % Standard deviations from assignment specifications
            v_pos = aircraft_dynamics.sigmaPos * randn(N, 3);                 % GPS Position Noise
            v_vel = aircraft_dynamics.sigmaVel * randn(N, 3);                 % GPS Velocity Noise (0.01 m/s)
            v_att = (aircraft_dynamics.sigmaAtt * (pi/180)) * randn(N, 3);    % GPS Attitude Noise (0.04 deg -> rad)
            v_v   = aircraft_dynamics.sigmaVTAS * randn(N, 1);                  % Airdata V Noise (0.2 m/s)
            v_alp = (aircraft_dynamics.AoA * (pi/180)) * randn(N, 1);     % Airdata Alpha Noise (0.1 deg -> rad)
            v_bet = (aircraft_dynamics.SideSlip * (pi/180)) * randn(N, 1);    % Airdata Beta Noise (0.25 deg -> rad)
            
            % Ground speed components (Clean tracking derivatives)
            u_GS_m = xdot + v_vel(:, 1); [cite: 112]
            v_GS_m = ydot + v_vel(:, 2); [cite: 114]
            w_GS_m = zdot + v_vel(:, 3); [cite: 116]
            
            % Assemble Noisy Measurement Matrix (N x 12) [cite: 90]
            Z_meas = [x_pos + v_pos(:, 1), y_pos + v_pos(:, 2), z_pos + v_pos(:, 3), ... % 1:3 GPS Pos
                      u_GS_m, v_GS_m, w_GS_m, ...                                        % 4:6 GPS Vel
                      data.phi + v_att(:, 1), data.theta + v_att(:, 2), data.psi + v_att(:, 3), ... % 7:9 GPS Att
                      data.vtas + v_v, data.alpha + v_alp, data.beta + v_bet];           % 10:12 Airdata
        end

        function xdot = state_propagation(x, u_m)
            % x = [x,y,z,u,v,w,phi,theta,psi,
            % lambdax, lambday, lambdaz, lambdap, lambdaq, lambdar,
            % Wx, Wy, Wz].T
            % u = [Axm, Aym, Azm, pm, qm, rm].T
            % z = [xm, ym, zm, um, vm, wm,
            % phim, thetam, psim, vtasm, alpham, betam].T
        
            % 1. Unpack states
            x_pos = x(1); y_pos = x(2); z_pos = x(3);
            u = x(4); v = x(5); w = x(6);
            phi = x(7); theta = x(8); psi = x(9);
            lambda_Ax = x(10); lambda_Ay = x(11); lambda_Az = x(12);
            lambda_p  = x(13); lambda_q  = x(14); lambda_r  = x(15);
            Wxe = x(16); Wye = x(17); Wze = x(18);
            
            % 2. Bias correction
            Ax = u_m(1) - lambda_Ax;
            Ay = u_m(2) - lambda_Ay;
            Az = u_m(3) - lambda_Az;
            p  = u_m(4) - lambda_p;
            q  = u_m(5) - lambda_q;
            r  = u_m(6) - lambda_r;
        
            g = 9.81; % Gravitational acceleration
            
            % 3. Initialize state derivative vector
            xdot = zeros(18, 1);
            
            % Kinematic translational equations (NED positions) 
            xdot(1) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*cos(psi) - (v*cos(phi) - w*sin(phi))*sin(psi) + Wxe;
            xdot(2) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*sin(psi) + (v*cos(phi) - w*sin(phi))*cos(psi) + Wye;
            xdot(3) = -u*sin(theta) + (v*sin(phi) + w*cos(phi))*cos(theta) + Wze;
            
            % Body velocity equations [cite: 74, 75, 76]
            xdot(4) = Ax - g*sin(theta) + r*v - q*w;
            xdot(5) = Ay + g*cos(theta)*sin(phi) + p*w - r*u;
            xdot(6) = Az + g*cos(theta)*cos(phi) + q*u - p*v;
            
            % Kinematic rotational equations (Attitude) [cite: 76, 77, 78]
            xdot(7) = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
            xdot(8) = q*cos(phi) - r*sin(phi);
            xdot(9) = q*(sin(phi)/cos(theta)) + r*(cos(phi)/cos(theta));
            
            % Biases and Wind are modeled as constants (derivatives are zero) [cite: 163, 165]
            xdot(10:12) = 0; % Accelerometer biases
            xdot(13:15) = 0; % Gyro biases
            xdot(16:18) = 0; % Wind components
        end

        function z_pred = measurement_model(x)
            % x: 18x1 estimated state vector from your EKF
            % z_pred: 12x1 predicted measurement vector h(x)
            
            % 1. Unpack required states for readability
            x_pos = x(1);
            y_pos = x(2);
            z_pos = x(3);
            u     = x(4);
            v     = x(5);
            w     = x(6);
            phi   = x(7);
            theta = x(8);
            psi   = x(9);
            Wxe   = x(16);
            Wye   = x(17);
            Wze   = x(18);
            
            % 2. Pre-calculate trigonometric terms for computational speed
            s_phi = sin(phi);   c_phi = cos(phi);
            s_the = sin(theta); c_the = cos(theta);
            s_psi = sin(psi);   c_psi = cos(psi);
            
            % 3. Initialize the 12x1 measurement prediction vector
            z_pred = zeros(12, 1);
            
            % --- GPS Position Observables ---
            z_pred(1) = x_pos; % x_GPS
            z_pred(2) = y_pos; % y_GPS
            z_pred(3) = z_pos; % z_GPS
            
            % --- GPS Velocity Observables (Transformation + Wind) ---
            % Transforms body velocities to navigation frame and adds constant wind
            z_pred(4) = (u*c_the + (v*s_phi + w*c_phi)*s_the)*c_psi - (v*c_phi - w*s_phi)*s_psi + Wxe; % u_GS
            z_pred(5) = (u*c_the + (v*s_phi + w*c_phi)*s_the)*s_psi + (v*c_phi - w*s_phi)*c_psi + Wye; % v_GS
            z_pred(6) = -u*s_the + (v*s_phi + w*c_phi)*c_the + Wze;                                    % w_GS
            
            % --- GPS Attitude Observables ---
            z_pred(7) = phi;   % phi_GPS
            z_pred(8) = theta; % theta_GPS
            z_pred(9) = psi;   % psi_GPS
            
            % --- Airdata Sensor Observables ---
            z_pred(10) = sqrt(u^2 + v^2 + w^2);         % True Airspeed (V)
            z_pred(11) = atan2(w, u);                    % Angle of Attack (alpha)
            z_pred(12) = atan2(v, sqrt(u^2 + w^2));      % Side slip angle (beta)
        end

        function G = get_noise(x)
            % x = [x,y,z,u,v,w,phi,theta,psi,Wxe,Wye,Wze,
            % lambdax, lambday, lambdaz, lambdap, lambdaq, lambdar].T
            % imu_noise = [wx, wy, wz, wp, wq, wr].T
        
            % 1. Unpack states
            u = x(4); v = x(5); w = x(6);
            phi = x(7); theta = x(8);
        
            % 2. Matrix G definition
            G = zeros(18,6);
        
            % Row 4: u_dot noise coefficients (-wx + w*wq - v*wr)
            G(4, :) = [-1,  0,  0,  0,  w, -v];
            
            % Row 5: v_dot noise coefficients (-wy - w*wp + u*wr)
            G(5, :) = [ 0, -1,  0, -w,  0,  u];
            
            % Row 6: w_dot noise coefficients (-wz + v*wp - u*wq)
            G(6, :) = [ 0,  0, -1,  v, -u,  0];
            
            % Row 7: phi_dot noise coefficients
            G(7, :) = [ 0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta)];
            
            % Row 8: theta_dot noise coefficients
            G(8, :) = [ 0,  0,  0,  0, -cos(phi),  sin(phi)];
            
            % Row 9: psi_dot noise coefficients
            G(9, :) = [ 0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta)];
            
            % Rows 1:3 (positions), 10:15 (biases), and 16:18 (wind) remain 0 
            % because process noise does not directly drive their differential equations.
        end
        
        function imu_noise = generate_imu_noise(N_samples)
            % GENERATE_IMU_NOISE Generates a [samples x 6] matrix of white Gaussian 
            % noise for IMU accelerometers (m/s^2) and gyroscopes (rad/s).
            
            % IMU noise standard deviations (matching bias levels)
            sigma_wx = 0.02; % m/s^2
            sigma_wy = 0.02; % m/s^2
            sigma_wz = 0.03; % m/s^2
            
            % Convert gyro noise standard deviations from deg/s to rad/s
            deg2rad  = pi / 180;
            sigma_wp = 0.005 * deg2rad; % rad/s
            sigma_wq = 0.005 * deg2rad; % rad/s
            sigma_wr = 0.002 * deg2rad; % rad/s
            
            % Generate independent Gaussian noise channels
            wx = sigma_wx * randn(N_samples, 1);
            wy = sigma_wy * randn(N_samples, 1);
            wz = sigma_wz * randn(N_samples, 1);
            wp = sigma_wp * randn(N_samples, 1);
            wq = sigma_wq * randn(N_samples, 1);
            wr = sigma_wr * randn(N_samples, 1);
            
            % Package into a clean matrix column-wise (with terminating semicolon)
            imu_noise = [wx, wy, wz, wp, wq, wr];
        end
    end
end

