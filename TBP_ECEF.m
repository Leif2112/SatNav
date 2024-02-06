function [Xdot_ECEF] = TBP_ECEF(t, X, GM)

%{ 

Function to output system of first order ODEs to solve equations of motion
in ECEF frame

Inputs:
    t : time [s]
    X : State Vector of spacecraft in ECEF frame 
    GM : Earth's gravitational parameter [km3/s2]

Output: 
    Xdot_ECEF = dXdt time derivative of cartesian coordinates 

Method:
    rewrite TBP as system of 6 first-order ODEs Xdot = [...]
    Use ODE45 to integrate the equations of motion over time --> X(t)
    
%}

r = norm(X(1:3));

%Angular velocity of the ECEF frame as seen from ECI frame
we = 7.29212351699038e-5;
w = [0 0 we]';

%TBP as system of first-order ODEs
xdot = X(4:6);
vdot = -(GM/r^3) * X(1:3) - cross(2*w, (X(4:6))) - cross(w, cross(w, X(1:3)));

Xdot_ECEF = [xdot; vdot];



