function [Xdot_Att] = AttitudeDynamics(t, X, Im)

%{ 

Function to integrate attitude dynamics of satellite

Inputs:
    t : time [s]
    X : Attitude state vector of spacecraft 
    I : Inertia matrix of spacecraft body

Output: 
    Xdot_Att = attitude dynamics as function of time
    
%}

B = [0          sin(X(3))               cos(X(3)); 
     0          cos(X(2))*cos(X(3))     -cos(X(2))*sin(X(3)); 
     cos(X(2))  sin(X(2))*sin(X(3))     sin(X(2))*cos(X(3))];

AngVel = X(4:6);

EulAng_dot = (1 / cos(X(2))) * B * AngVel;
EulEq1 = -(cross(AngVel, Im*AngVel));
EulEq = EulEq1' / Im;

Xdot_Att = [EulAng_dot; EulEq'];