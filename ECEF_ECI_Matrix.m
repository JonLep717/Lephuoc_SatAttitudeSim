function A_ECEF_ECI = ECEF_ECI_Matrix(Y,M,D,h,m,s)
JD_0 = julian_date(Y,M,D,0,0,0);

GMST_rad = GMST(JD_0,Y,M,D,h,m,s);

A_ECEF_ECI = [cos(GMST_rad), cos(pi/2)-GMST_rad, 0;
    cos(pi/2)+GMST_rad, cos(GMST_rad),0;
    0,0,1];
end
%%
function GMST = GMST(JD, Y, M, D, h, m, s)
T_0 = (JD - 2451545)/36525;
GMST = 24110.54841 + (8640184.812866*T_0) + (0.093104*T_0^2) - ((6.2*(10^(-6)))*(T_0^3)) + 1.002737909350795*((3600*h)+(60*m)+s);
if (GMST > 86400)
    while (GMST >86400)
        GMST = GMST-86400;
    end

elseif (GMST < 0)
    while (GMST < 0)
        GMST = GMST + 86400;
    end
end
GMST = GMST/240;
GMST = GMST * (pi/180);
end