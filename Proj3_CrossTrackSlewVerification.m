close all;
clear all;

J1 = 77217;
J2 = 77217;
J3 = 25000;

h_dot = 0.05;
f = @(t,y) [y(2); h_dot/J2];

t01 = 0;
y01 = [0;0];

T1 = 15*60;
T2 = 30*60;
T3 = 60*60;
[ts1,ys1] = ode45(f,[t01,T1],y01);
t02 = T1;
y02 = [ys1(length(ys1(:,1)),1); ys1(length(ys1(:,2)),2)];

h_dot = .05;
f = @(t,y) [y(2); -h_dot/J2];
[ts2,ys2] = ode45(f,[t02,T2],y02);
t03 = T2;
y03 = [ys2(length(ys2(:,1)), 1); ys2(length(ys2(:,2)),2)];
h_dot = 0;
f = @(t,y) [y(2); h_dot/J2];
[ts3,ys3] = ode45(f,[t03,T3], y03);
plot(ts1/60,ys1(:,1)*(180/pi),'g')
hold on
plot(ts2/60,ys2(:,1)*(180/pi),'y')
plot(ts3/60, ys3(:,1)*(180/pi), 'r')
title('Slewing Angle vs. Time')
xlabel('Time [minutes]')
ylabel('Angle [degrees]')
legend("Angular Acceleration Period", "Angular Deceleration Period", "Settling Period")
legend("Location", "Southeast")
grid on
hold off

final_angle = rad2deg(ys3(length(ys2(:,1)),1))
ss_error = (abs(30 - final_angle))