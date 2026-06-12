function xdot = aircraft_dynamics(x, u_m)
    % x = [x,y,z,u,v,w,phi,theta,psi,Wxe,Wye,Wze,
    % lambdax, lambday, lambdaz, lambdap, lambdaq, lambdar].T
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

function delta_xdot = aircraft_dynamics_noise(x, imu_noise)
    % x = [x,y,z,u,v,w,phi,theta,psi,Wxe,Wye,Wze,
    % lambdax, lambday, lambdaz, lambdap, lambdaq, lambdar].T
    % imu_noise = [wx, wy, wz, wp, wq, wr].T

    % 1. Unpack states
    x_pos = x(1); y_pos = x(2); z_pos = x(3);
    u = x(4); v = x(5); w = x(6);
    phi = x(7); theta = x(8); psi = x(9);
    lambda_Ax = x(10); lambda_Ay = x(11); lambda_Az = x(12);
    lambda_p  = x(13); lambda_q  = x(14); lambda_r  = x(15);
    Wxe = x(16); Wye = x(17); Wze = x(18);

    % 2. Unpack uncertainties
    wx = imu_noise(1); wy=imu_noise(2); wz=imu_noise(3);
    wp = imu_noise(4); wq=imu_noise(5); wr=imu_noise(6);

    delta_xdot = zeros(18,1);

    % Velocity derivatives
    delta_xdot(4) = -wx - v*wr + w*wq;
    delta_xdot(5) = -wy - w*wp + u*wr;
    delta_xdot(6) = -wz - u*wq + v*wp;

    % Angle derivatives
    delta_xdot(7) = -wp - wq*sin(phi)*tan(theta) - wr*cos(phi)*tan(theta);
    delta_xdot(8) = -wq*cos(phi) + wr*sin(phi);
    delta_xdot(9) = -wq*sin(phi)/cos(theta) - wr*cos(phi)/cos(theta);
end

function imu_noise = generate_imu_noise(samples)
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
    wx = sigma_wx * randn(samples, 1);
    wy = sigma_wy * randn(samples, 1);
    wz = sigma_wz * randn(samples, 1);
    wp = sigma_wp * randn(samples, 1);
    wq = sigma_wq * randn(samples, 1);
    wr = sigma_wr * randn(samples, 1);
    
    % Package into a clean matrix column-wise (with terminating semicolon)
    imu_noise = [wx, wy, wz, wp, wq, wr];
end

