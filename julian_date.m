function jd = julian_date(Y,M,D,h,m,s)
    jd = 1721013.5 + 367*Y - floor(7/4*(Y+floor((M+9)/12))) + floor(275*M/9) + D + (60*h+m+s/60)/1440;
end