function [X] = COE2RV(COE, GM)
%{
convert state of spacecraft from COE to ECI position and veloctiy
components

Inputs: 
    COE : Matrix of classical orbital elements
    GM : Earth's gravitational parameter [km3/s2] 

Outputs: 
    X : State vector of spacecraft in ECI

Method: 
    
    Calculate r = a*(1-e^2) / (1 + e cosTA) 
    Calculate h = sqrt(mu*a*(1-e^2))
    Rotation matrix from perifocal frame to ECI  : 
        [x y z] = ROTA * (r*cosTA r*sinTA 0)
        [xdot ydot zdot] = ROTA * (-(GM/h)*sinTA (GM/h)*(cosTA+e) 0)
%}

%Classical Orbital paramters
a = COE(1);
ECC = COE(2);
I = COE(3);
RAAN = COE(4);
argP = COE(5);
TA = COE(6);

%Compute position vector in perifocal frame 
pos = a * (1-ECC^2) / (1 + ECC * cos(TA));
r_P = [pos * cos(TA); pos * sin(TA); 0];

%Compute specific angular momentum 
h = sqrt(GM * a * (1 - ECC^2));

%Comput velocity vector (rdot) in perifocal frame
v_P = [-(GM/h)*sin(TA); (GM/h)*(cos(TA) + ECC); 0];

%Rotation matrices
R3_RAAN = [  cos(RAAN)      sin(RAAN)      0;
           -sin(RAAN)     cos(RAAN)      0;
           0              0              1];
R1_I =    [  1        0            0;
             0        cos(I)       sin(I);
             0        -sin(I)      cos(I)];

R3_argP = [  cos(argP)      sin(argP)    0;
           -sin(argP)     cos(argP)      0;
           0              0              1];

ECItoP = R3_argP*R1_I*R3_RAAN; %3-1-3 rotation to go from Perifocal frame to ECI
PtoECI = ECItoP';                 

r_ECI = PtoECI * r_P;
v_ECI = PtoECI * v_P;
X = [r_ECI; v_ECI];