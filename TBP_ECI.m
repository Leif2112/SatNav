function [Xdot_ECI] = TBP_ECI(t, X, GM)

%{ 

Function to output system of first order ODEs to solve equations of motion
in ECI frame

Inputs:
    t : time [s]
    X : State Vector of spacecraft in ECI 
    GM : Earth's gravitational parameter [km3/s2]

Output: 
    Xdot_ECI = dXdt time derivative of cartesian coordinates 

Method:
    rewrite TBP as system of 6 first-order ODEs Xdot = [...]
    Use ODE45 to integrate the equations of motion over time --> X(t)
    
%}

%Define state vector / extract components
x1 = X(1);
x2 = X(2);
x3 = X(3);
x4 = X(4);
x5 = X(5);
x6 = X(6);

r = sqrt(x1^2 + x2^2 + x3^2);

%TBP as system of first-order ODEs
x1dot = x4;
x2dot = x5;
x3dot = x6;
x4dot = -(GM/r^3) * x1;
x5dot = -(GM/r^3) * x2;
x6dot = -(GM/r^3) * x3;

Xdot_ECI = [x1dot, x2dot, x3dot, x4dot, x5dot, x6dot]';



