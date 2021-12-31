function [date, pos_ECI, vel_ECI] = readCSV(filename)
filename = 'rv_vec_HST.csv';
T = readtable(filename);
date = T.Time_UTCG_(:,:);
x_pos = T.x_km_(:,:);
y_pos = T.y_km_(:,:);
z_pos = T.z_km_(:,:);
x_vel = T.vx_km_sec_(:,:);
y_vel = T.vy_km_sec_(:,:);
z_vel = T.vz_km_sec_(:,:);

Y = zeros(length(date),1);
M = zeros(length(date),1);
D = zeros(length(date),1);
hr = zeros(length(date),1);
min = zeros(length(date),1);
sec = zeros(length(date),1);
for i = 1:length(date)
    t = datetime(date{i}, 'InputFormat', 'dd MMM uuuu HH:mm:ss.SSS');
    Y(i) = year(t);
    M(i) = month(t, 'monthofyear');
    D(i) = day(t, 'dayofmonth');
    hr(i) = hour(t);
    min(i) = minute(t);
    sec(i) = second(t);
end

pos_ECI = [x_pos,y_pos,z_pos];
vel_ECI = [x_vel,y_vel,z_vel];
end