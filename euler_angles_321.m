function [phi, theta, psi] = euler_angles_321(q)
A = A_q(q);
phi = atan2(A(2,3),A(3,3));
theta = -asin(A(1,3));
psi = atan2(A(1,2),A(1,1));
end