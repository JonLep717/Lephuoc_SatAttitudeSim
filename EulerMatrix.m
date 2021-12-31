function A_body_ECI = EulerMatrix(phi, theta, psi)
A_BI0 = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)];
A_BI1 = [cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)];
A_BI2 = [-sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)];

A_body_ECI = [A_BI0;A_BI1;A_BI2];
end