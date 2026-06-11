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
    
    % Kinematic translational equations (NED positions) [cite: 72, 73]
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