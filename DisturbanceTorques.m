function L_dist = DisturbanceTorques(quat, pos_ECI, vel_ECI)
J_c = diag([77217,77217,25000]);
C_D = 2.5;
rho = 3e-13;
plate_areas = [23,18,18];
plate_normals = [[1;0;0],[cos(30*pi/180);sin(30*pi/180);0],[cos(30*pi/180);sin(30*pi/180);0]];
plate_positions = [[1;0;0], [0;0.6;0], [0;-0.6;0]];
plate_absorpt = [0.1,0.2,0.2];
plate_spec = [0.8,0.75, 0.75];
plate_diffus = [0.1,0.05, 0.05];

%%pos_ECI = [-6811.15;-449.726;1097.175];
%%vel_ECI = [0.997668;-6.71024;3.419294];
L_gg = GravGrad(J_c,quat,pos_ECI, vel_ECI);
D = [20;20;20];
L_mag = MagTorq(D, quat, pos_ECI, vel_ECI);
L_drag = AeroTorque(pos_ECI, vel_ECI, rho, C_D, quat, plate_areas, plate_normals, plate_positions);
L_SRP = SolarTorque(pos_ECI, vel_ECI, quat, plate_areas, plate_normals, plate_positions, plate_absorpt, plate_spec, plate_diffus);
L_dist = L_gg + L_mag + L_drag + L_SRP;
end
%% Gravity Gradient Torque

function L_grav = GravGrad(J, quat, position,velocity)
mu = 398600.441;     % Standard Gravity Parameter
%mu = 3.9860e14;
R_c = norm(position);
%R_C = R_c*(10^3)
%phi = Eul(1);
%theta = Eul(2);
%psi = Eul(3);
%gg_x = (J(3)-J(2))*(((cos(theta))^2)*cos(psi)*sin(psi));
%gg_y = (J(3)-J(1))*((cos(theta)*sin(theta)*cos(psi)));
%gg_z = (J(1)-J(2))*((cos(theta)*sin(theta)*sin(psi)));
%L_grav = (3*mu/(R_c^3))*[gg_x;gg_y;gg_z];

A_body_LVLH = A_q(quat);
A_LVLH_ECI = A_LVLH(position,velocity);
A_body_ECI = A_body_LVLH*A_LVLH_ECI;

R_c_body = A_body_ECI * position;
L_grav = (3*mu/(R_c^5))*R_c_body;
L_grav = cross(L_grav, (J*R_c_body));
end
%% Magnetic Torque

function L_mag = MagTorq(D, quat, position,velocity)
year = 2021;
month = 11;
day = 12;
hour = 6;
minute = 0;
second = 0;
A_ECEF_ECI = ECEF_ECI_Matrix(year,month,day,hour,minute,second);
A_ECI_ECEF = transpose(A_ECEF_ECI);
%position = 1000*position;
a = 6371.2;   % magnetic spherical reference radius [m]
g1_0 = -29554.63e-9;     % magnetic coefficient [T]
g1_1 = -1669.05e-9;      % magnetic coefficient [T]
h1_1 = 5077.99e-9;       % magnetic coefficient [T]
m_vec = (a^3)*[g1_1;h1_1;g1_0];
position_ECEF = A_ECEF_ECI*position;

A_body_LVLH = A_q(quat);
A_LVLH_ECI = A_LVLH(position,velocity);

A_LVLH_ECEF = A_LVLH_ECI*A_ECI_ECEF;
A_body_ECEF = A_body_LVLH*A_LVLH_ECEF;
B_r = (1/((norm(position_ECEF))^5)) * ((3*dot(m_vec,position_ECEF))*position_ECEF - ((norm(position_ECEF)^2)*m_vec));
L_mag = cross(D, B_r);
L_mag = (A_body_ECEF*L_mag);
end
%% Aerodynamic Drag Torque

function L_aero = AeroTorque(pos_ECI, vel_ECI, rho, C_D, quat, plate_area, plate_norm, plate_pos)
A_body_LVLH = A_q(quat);
A_LVLH_ECI = A_LVLH(pos_ECI,vel_ECI);

A_body_ECI = A_body_LVLH * A_LVLH_ECI;

w_earth = [0;0;7.2921159e-5];

v_rel_ECI = vel_ECI + cross(w_earth,pos_ECI);
v_rel_body = A_body_ECI*v_rel_ECI;
v_rel_body_norm = norm(v_rel_body);

L_aero = [0;0;0];
for i = 1:length(plate_area)
    cos_theta_i = (dot(plate_norm(:,i),v_rel_body))/v_rel_body_norm;
    F_aero_i = -(1/2)*rho*C_D*v_rel_body_norm*v_rel_body*plate_area(:,i)*max(cos_theta_i,0);
    cross_r_f_i = cross(plate_pos(:,i), F_aero_i);
    L_aero = L_aero + cross_r_f_i;
end
L_aero = L_aero*(10^3);
end
%% Solar Radiation Torque

function L_sol = SolarTorque(pos_ECI, vel_ECI, quat, plate_area, plate_norm, plate_pos, plate_absorpt, plate_spec, plate_diffus)
year = 2021;
month = 11;
day = 12;
hour = 6;
minute = 0;
second = 0;

A_body_LVLH = A_q(quat);
A_LVLH_ECI = A_LVLH(pos_ECI,vel_ECI);

A_body_ECI = A_body_LVLH * A_LVLH_ECI;

JD = julian_date(year,month,day,hour,minute,second);
T_UT1 = (JD-2451545)/36525;

M_dot = 357.5277233 + 35999.05034*T_UT1;

phi_dot = 280.46 + 36000.771*T_UT1;
eps = 23.439291 - 0.0130042*T_UT1;
phi_eclip = phi_dot + 1.914666471*sin(deg2rad(M_dot)) + 0.019994643*sin(2*deg2rad(M_dot));
r_hat_cross_dot = [cos(deg2rad(phi_eclip));cos(deg2rad(eps))*sin(deg2rad(phi_eclip));sin(deg2rad(eps))*sin(deg2rad(phi_eclip))];
r_cross_dot = 1.000140612 - 0.016708617*cos(deg2rad(M_dot)) - 0.000139589*cos(2*deg2rad(M_dot));

r_I = 10^3 * pos_ECI;
r_I_AU = (1/6.68459e12)*r_I;

r_cross_dot_vec = r_cross_dot*r_hat_cross_dot;
r_sat_dot = r_cross_dot_vec - r_I_AU;
r_sat_dot_norm = norm(r_sat_dot);
r_sat_dot_unit = (1/r_sat_dot_norm)*r_sat_dot;

s_vec = A_body_ECI * r_sat_dot_unit;

Phi = 1361;      %% Solar Constant
c = 3e8;         %% Speed of Light
P = Phi/(c*r_sat_dot_norm^2);
L_sol = zeros(3,1);
for i = 1:length(plate_area)
    cos_theta_SRP = (transpose(plate_norm(:,i)))*s_vec;
    F_SRP_i = -P*plate_area(i)*(2*((plate_diffus(:,i)/3) + (plate_spec(:,i)*cos_theta_SRP))*plate_norm(:,i) + (1-plate_spec(:,i))*s_vec)*max(cos_theta_SRP,0);
    cross_r_F_i = cross(plate_pos(:,i), F_SRP_i);
    L_sol = L_sol + cross_r_F_i;
end
end