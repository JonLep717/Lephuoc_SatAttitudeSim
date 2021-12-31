close all;
clear all;

%% This is the main file to run. You can verify questions 3, 4, and 5 
%% (detumbling, fixed attitude, and open-loop disturbance torques)

%% See Proj3_CrossTrackSlewVerification.m for verification of Q1 and Q2

[date, pos_ECI, vel_ECI] = readCSV('rv_vec_HST.csv');   %% Read in r and v vectors

C_lvlh_eci = A_LVLH(pos_ECI(1,:).', vel_ECI(1,:).');    %% Define initial LVLH frame

q0_lvlh = shepperd_q(C_lvlh_eci);   %% Create corresponding quaternion from inertial to orbital


Jc = diag([77217,77217,25000]);     %% Define moment of inertia of spacecraft

semi_maj = 6920300.;               
mu = 3.9860e14;
n = sqrt(mu)/(semi_maj)^1.5;
T = 2*pi/n;
w_oi = [0; -n; 0];                  %% Defines angular velocity of orbital frame

t = linspace(0, 86400,1441);
t_state_vec = transpose(t);
tspan = [0 86400];
length(tspan);
init_att = input("Input 0 for initial attitude of [0;0;0;1], Input 1 for initial attitude of [.5;.5;.5;.5], Input 2 to input your own initial body-orbit attitude: ");
if init_att == 0
    q0_body = [0;0;0;1];
elseif init_att == 1
    q0_body = [.5;.5;.5;.5];
elseif init_att == 2
    phi_init = input("Input initial attitude phi [deg]: ");
    theta_init = input("Input initial attitude theta [deg]: ");
    psi_init = input("Input initial attitude psi [deg]: ");
    A_init_body = EulerMatrix(deg2rad(phi_init), deg2rad(theta_init), deg2rad(psi_init));
    q0_body = shepperd_q(A_init_body);
end
mode = 0;
cont = input("Input 0 for open-loop, no controller, Input 1 for closed-loop controller: ");  %% cont == 0 for no controller, cont == 1 for controller
q_target = [0;0;0;0]; % initialize target attitude quaternion
w_0 = [0;0;0];      % initialize initial angular rate
if cont == 1
    mode = input("Input 0 for fixed attitude mode (no initial angular velocity), Input 1 for detumble mode, Input 2 to test slewing mode: ");
    if mode == 0 || 2
        w_0 = [0;0;0];
        if mode == 2
            phi_targ = input("Input target attitude phi [deg]: ");
            theta_targ = input("Input target attitude theta [deg]: ");
            psi_targ = input("Input target attitude psi [deg]: ");
            A_targ_body = EulerMatrix(deg2rad(phi_targ), deg2rad(theta_targ), deg2rad(psi_targ));
            q_target = shepperd_q(A_targ_body);
        end
    elseif mode == 1
        w_0 = [.1;.1;.1];
    end
end
options = odeset('RelTol', 1e-10,'AbsTol', 1e-10);

sol = ode45(@(t,q) eq_motion(t, q, q0_body, q_target, w_oi, Jc, pos_ECI,vel_ECI, t_state_vec, cont, mode), tspan, [q0_lvlh;q0_body;w_0], options);

y = deval(sol,t);

figure
plot(t/3600, y(1:4,:))
title("Orbit-Inertial quaternion")
legend('q1', 'q2', 'q3', 'q4')
xlabel("Time [hours]")

figure
plot(t/3600, y(5:8,:))
title("Body-Orbit Quaternion")
xlabel("Time [hours]")
legend('q1', 'q2', 'q3', 'q4')

figure
for i=1:size(y,2)
    [phi, theta, psi] = euler_angles_321(y(5:8,i));
    angles(1:3,i) = [rad2deg(phi), rad2deg(theta), rad2deg(psi)];
end
plot(t/3600, angles)
title("Body-Orbit Angles [degrees]")
legend(["phi","theta","psi"])
xlabel("Time [hours]")

figure
plot(t/3600, y(9:11,:))
title("Body-Orbit Angular Velocity [deg/s]")
legend(["w1", "w2", "w3"])
xlabel("Time [hours]")

function [dydt, Lc] = eq_motion(t, y, q0_body, q_targ, w_oi, J_c,pos_ECI, vel_ECI, t_state, cont, mode)
    x_pos = pos_ECI(:,1);
    y_pos = pos_ECI(:,2);
    z_pos = pos_ECI(:,3);
    x_vel = vel_ECI(:,1);
    y_vel = vel_ECI(:,2);
    z_vel = vel_ECI(:,3);

    x_pos_int = interp1(t_state, x_pos,t);
    y_pos_int = interp1(t_state, y_pos,t);
    z_pos_int = interp1(t_state, z_pos,t);

    x_vel_int = interp1(t_state, x_vel,t);
    y_vel_int = interp1(t_state, y_vel,t);
    z_vel_int = interp1(t_state, z_vel,t);
    
    pos_ECI_int = [x_pos_int, y_pos_int, z_pos_int].';
    vel_ECI_int = [x_vel_int, y_vel_int, z_vel_int].';
    
    q1 = y(1:4);
    q2 = y(5:8);

    w = y(9:11);
    A = [q1(4), -q1(3),  q1(2);
         q1(3),  q1(4), -q1(1);
         -q1(2), q1(1),  q1(4);
         -q1(1), -q1(2), -q1(3)];
    A2 = [q2(4), -q2(3),  q2(2);
         q2(3),  q2(4), -q2(1);
         -q2(2), q2(1),  q2(4);
         -q2(1), -q2(2), -q2(3)];
    dqdt = 0.5*A*w_oi;
    dq2dt = 0.5*A2*w;
    A_bo = A_q(q2);
    A_ob = A_bo.';

    w_bi = w + A_bo*w_oi;

    if cont==1
        constant_att = EulerMatrix(0,0,0);

        if mode == 1                %% detumble the spacecraft from initial ang. velocity
            Kx = .35;
            Kx_d = 32;
        
            Ky = .35;
            Ky_d = 35;
        
            Kz = .35;
            Kz_d = 32;
            q_c = shepperd_q(constant_att);
        elseif mode == 0            %% hold a constant attitude (body=initial attitude)
            Kx = 850;
            Kx_d = 10000;
        
            Ky = 350;
            Ky_d = 10000;
        
            Kz = 650;
            Kz_d = 10000;
            q_c = q0_body;
        elseif mode == 2            %% Slewing mode, closed loop, just for fun
            Kx = 1;
            Kx_d = 2000;
        
            Ky = 1;
            Ky_d = 2000;
        
            Kz = 1;
            Kz_d = 2000;
            %[phi0, theta0, psi0] = euler_angles_321(q0_body);
            %A_30_deg_CT = EulerMatrix(phi0, theta0 - deg2rad(32.4), psi0);
            %q_c = shepperd_q(A_30_deg_CT);
            q_c = q_targ;
        end

        q_e = [q_c(4), q_c(3), -q_c(2), q_c(1);
            -q_c(3), q_c(4), q_c(1), q_c(2);
            q_c(2), -q_c(1), q_c(4), q_c(3);
            -q_c(1), -q_c(2), -q_c(3), q_c(4)] * [-q2(1);-q2(2);-q2(3);q2(4)];
        
        Tx = -2*Kx*q_e(1)*q_e(4) + Kx_d*w(1);
        Ty = -2*Ky*q_e(2)*q_e(4) + Ky_d*w(2);
        Tz = -2*Kz*q_e(3)*q_e(4) + Kz_d*w(3);
        Lc = [Tx;Ty;Tz];

        if norm(Lc) > 1
            norm(Lc)    %% Return the value of the torque if it exceeds 1Nm
        end
    end

    Ld = DisturbanceTorques(q2, pos_ECI_int, vel_ECI_int);

    if cont==1
        dwdt = J_c\(Ld - Lc - cross(w_bi, J_c*w_bi));
        dydt = [dqdt; dq2dt; dwdt];
    else
        dwdt = J_c\(Ld - cross(w_bi, J_c*w_bi));
        dydt = [dqdt; dq2dt; dwdt];
    end
end

function A = A_q(q)
q_skew = [0 -q(3) q(2);
          q(3), 0, -q(1);
          -q(2), q(1), 0];
A = (q(4)^2-norm(q(1:3))^2)*eye(3) - 2*q(4)*q_skew + 2*q(1:3)*q(1:3)';
end

function q = shepperd_q(A)
    % Shepperd quaternion
    trC  = trace(A);
    q_0 = sqrt((1+trC)/4);
    q_1 = sqrt((1+2*A(1,1)-trC)/4);
    q_2 = sqrt((1+2*A(2,2)-trC)/4);
    q_3 = sqrt((1+2*A(3,3)-trC)/4);
    [max_q, max_idx] = max([q_0, q_1, q_2, q_3]);
    if max_idx == 1
        q_1 = (A(2,3) - A(3,2))/4/max_q;
        q_2 = (A(3,1) - A(1,3))/4/max_q;
        q_3 = (A(1,2) - A(2,1))/4/max_q;
    elseif max_idx == 2
        q_0 = (A(2,3) - A(3,2))/4/max_q;
        q_2 = (A(1,2) + A(2,1))/4/max_q;
        q_3 = (A(3,1) + A(1,3))/4/max_q;
    elseif max_idx == 3
        q_0 = (A(3,1) - A(1,3))/4/max_q;
        q_1 = (A(1,2) + A(2,1))/4/max_q;
        q_3 = (A(2,3) + A(3,2))/4/max_q;
    else
        q_0 = (A(1,2) - A(2,1))/4/max_q;
        q_1 = (A(3,1) + A(1,3))/4/max_q;
        q_2 = (A(2,3) + A(3,2))/4/max_q;
    end
    q = [q_1; q_2; q_3; q_0];
end

function [phi, theta, psi] = euler_angles_321(q)
A = A_q(q);
phi = atan2(A(2,3),A(3,3));
theta = -asin(A(1,3));
psi = atan2(A(1,2),A(1,1));
end