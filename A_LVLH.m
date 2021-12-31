function A_orb_ECI = A_LVLH(r_I, v_I)
x = [1;0;0];
y = [0;1;0];
z = [0;0;1];
r_I_norm = norm(r_I);
r_I_unit = r_I/r_I_norm;
o_3 = -r_I_unit;
v_I_unit = v_I / norm(v_I);
o_2 = -(cross(r_I_unit,v_I_unit));
o_1 = cross(o_2,o_3);

A_orb_ECI = [dot(o_1,x), dot(o_1,y), dot(o_1,z);
    dot(o_2,x), dot(o_2,y), dot(o_2, z);
    dot(o_3,x), dot(o_3,y), dot(o_3,z);];
end