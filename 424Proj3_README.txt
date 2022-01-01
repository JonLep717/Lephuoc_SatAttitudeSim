To verify Q1 and Q2 (30 degree cross-track slewing), run the Proj3_CrossTrackSlewVerification.m file (HOLDOVER FILE FROM THE ASSIGNMENT ITSELF, NOT IMPORTANT)
To verify Q3, Q4, and Q5, run the ProjectAssignment3.m file (THE IMPORTANT FILE)

A_LVLH.m - calculates Direction Cosine Matrix (DCM) from Earth-Centered Inertial (ECI) frame to Local-Vertical Local-Horizontal (LVLH) frame
A_q.m - calculates DCM corresponding to a given quaternion
DisturbanceTorques.m - calculates total disturbance torques given a quaternion, a position vector in the ECI frame, and a velocity vector in the ECI frame
ECEF_ECI_Matrix - calculates DCM from ECI to ECEF frames
EulerMatrix - calculates 3-2-1 Euler transformation
julian_date - calculates Julian date for use in calculating Earth-Centered Earth-Fixed (ECEF) frame
rv_vec_HST - CSV containing position and velocity vectors in a given timespan (DO NOT CHANGE)
readCSV - processing rv_vec_HST.csv

Sorry for the mess!
