clearvars; close all; clc; format longG;

%Physical Constants
GM = 398600.4418;                                                        %Earth's gravitational parameter [km^3/s^2]
Re = 6378.137;                                                           %Earth's equatorial radius [km]
Te = 86164;                                                              %Earth's rotational period [s]
we = 2*pi / Te;                                                          %Angular velocity of spacecraft [rad/s]
Im = [2500 0       0; 
      0    5000    0;
      0    0       6500];                                                %Inertia matrix of spacecraft

AngVel = [-3.092e-4; 6.6161e-4; 7.4606e-4];

tol = 10e-10;
ode_opt = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%Classical Orbital Elements
a = 7151.16;                                                              %Semi-major axis [km]
ECC = 0.0008;                                                             %Eccentricity [~]
I = deg2rad(98.39);                                                       %Inclination [rad]
RAAN = deg2rad(10);                                                       %RAAN [rad]
argP = deg2rad(233);                                                      %Argument of periapsis [rad]
MA0 = deg2rad(127);

%Solve Kepler's equation 
[E0, ~] = Kepler(ECC, MA0, tol);                                          %Eccentric Anomaly in [rad]

%Compute the intial True Anomaly of the spacecraft 
TA0 = 2*atan2(sqrt(1+ECC)*tan(E0/2), sqrt(1-ECC));                        %True Anomaly in [rad]
fprintf('Initial True Anomaly TA0: %.4f rads\n', TA0);

%Check location along orbit
if pi/2 < TA0 && TA0 < 1.5*pi
    fprintf('The spacecraft is near its apoapsis\n');
else 
    fprintf('The spacecraft is near its periapsis\n');
end

%Compute orbital period of spacecraft
n = sqrt(GM/a^3);                                                         %Mean  motion [rad/s]
P = 2*pi/n;                                                               %Orbit period [s]

fprintf('Orbit Period P: %.4f (s)\n', P);

%Create time vector
t0 = 0;
t = linspace(t0, P + t0, 1000);                                            %Time vector [s]
MAt = mod(MA0 + n * (t - t0), 2*pi);                                       %Mean anomaly as a funciton of time [rad]

E_matrix = zeros(1, length(t));                                            %Initialising E_matrix

%Propgate state of satellite : MA(t), E(t), TA(t) using COE
for tt = 1:length(t)

    %Solve Kepler's equation
    [E,~] = Kepler(ECC, MAt(tt), tol);                                     %dEdt eccentric anomaly as a function of time
    E_matrix(tt) = E;                                                      %store value in matrix for plotting
    
    TA(tt) = mod(2*atan2(sqrt(1+ECC)*tan(E/2), sqrt(1-ECC)), 2*pi);        %Compute True Anomaly as a function of time

    %Store COE & convert to Cartesian coords / COE2RV
    COE(:, tt) = [a, ECC, I, RAAN, argP, TA(tt)];                          %Classical Orbital Elements matrix
    X(:,tt) = COE2RV(COE(:,tt), GM);                                       %State vector from classical orbital elements to RV propagation
end 


%Integrate equation of motion 
[~, Xout] = ode113(@TBP_ECI, t, X(1:6, 1), ode_opt, GM);                   %Integrate the equation of motion of the ECI TBP using ODE45
Xout = Xout';                                                              %Transpose of the output
v_Xout = sqrt(Xout(4, :).^2 + Xout(5, :).^2 + Xout(6, :).^2);              %Velocity magnitude from Xout integration 
r_Xout = sqrt(Xout(1, :).^2 + Xout(2, :).^2 + Xout(3, :).^2);              %Position magnitude from Xout integration

diff = abs((X(1:6, :) - Xout(1:6, :)));                                    %Position and velocity discrepency from COE2RV propagation
rdiff = sqrt(diff(1, :).^2 + diff(2, :).^2 + diff(3, :).^2);               %Position error magnitude 
vdiff = sqrt(diff(4, :).^2 + diff(5, :).^2 + diff(6, :).^2);               %Velocity error magnitude 

%Compute specific energy of the spacecraft at every time stamp
for i = 1:length(t)
    sp_e(i) = 0.5 * v_Xout(i)^2 - GM / r_Xout(i);                          %Specific energy of the spacecraft, should be approx. constant and can be approximated to -GM / 2*a
    sp_e2(i) = -GM / (2*a);
end 

w = [0 0 we]';
FI = eye(3);

X_ini_ECEF = [FI * X(1:3, 1); FI * X(4:6, 1)-cross(w, X(1:3, 1))];         %Spacecraft initial conditions in ECEF frame

%Integrate the equations of motion of the satellite in ECEF frame
[~, X_ECEF] = ode113(@TBP_ECEF, 10*t, X_ini_ECEF, ode_opt, GM);
X_ECEF = X_ECEF';

%Initial position and velocity in the ECI frame 1h into orbit
r6 = [6768.27 870.90 2153.59]';
v6 = [-2.0519 -1.4150 7.0323]';

%Initial Euler Angles
alpha = deg2rad(30);
beta = deg2rad(20);
gamma = deg2rad(10);

%Rotation matrices 1-2-1 sequence
R1_alpha = [  1             0               0;
              0             cos(alpha)      sin(alpha);
              0             -sin(alpha)     cos(alpha)];

R2_beta =  [  cos(beta)     0               -sin(beta);
              0             1               0;
              sin(beta)    0               cos(beta)];

R1_gamma = [  1             0               0;
              0             cos(gamma )     sin(gamma);
              0             -sin(gamma )    cos(gamma)];;

%DCM from Orbital to Body frame
R_BO = R1_gamma * R2_beta * R1_alpha;

%RSW Orbital frame axes
O1 = r6 / norm(r6);
h6 = cross(r6, v6);
O3 = h6 / norm(h6);
O2 = cross(O3, O1);

%DCM from Inertial to Orbital frame
R_OI = [O1(1) O1(2) O1(3);
        O2(1) O2(2) O2(3);
        O3(1) O3(2) O3(3)];

%DCM from Inertial to Body frame 
R_BI = R_BO *  R_OI;

EulAng_BO = acos((trace(R_BI) - 1) / 2);                            %Determine Euler angles to rotate from Orbital to Body frame

%Euler's principle axes
EulAX_1 = (R_BI(2,3) - R_BI(3,2)) / (2 * sin(EulAng_BO));           
EulAX_2 = (R_BI(3,1) - R_BI(1,3)) / (2 * sin(EulAng_BO));
EulAX_3 = (R_BI(1,2) - R_BI(2,1)) / (2 * sin(EulAng_BO));
EulAx = [EulAX_1 EulAX_2 EulAX_3]';

%Quaternion representation
q123 = EulAx * sin(EulAng_BO / 2);                                  %Determinne Euler parameters / quaternions 
q4 = cos(EulAng_BO/2);                                          
q = [q123; q4];
scale = norm(q);                                                    %nomalised value of quaternion
q = q/scale;                                                        %normalise quaternion values

%DCM from Inertial to Body frame from quaternions
q_BI = [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2 2*(q(1)*q(2) + q(4)*q(3)) 2*(q(1)*q(3) - q(4)*q(2));
       2*(q(1)*q(2) - q(4)*q(3)) -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2 2*(q(2)*q(3) + q(4)*q(1));
       2*(q(1)*q(3) + q(4)*q(2)) 2*(q(2)*q(3) - q(4)*q(1)) -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];

%Euler Angles from Inertial to body frame
psi = atan2(q_BI(1,2), q_BI(1,1));
theta = asin(-q_BI(1,3));
phi = atan2(q_BI(2,3), q_BI(3,3));
EulAng_BI = [psi; theta; phi];

%Initial attitude state
X_ini_Att = [EulAng_BI; AngVel];

t2 = linspace(0, 3600, 360);                                        
[~, X_Att] = ode113(@AttitudeDynamics, t2, X_ini_Att, ode_opt, Im);  %Integrate attitude of spacecraft / propagate in time

X_Att = X_Att';

%Compute Angular momentum and rotational kinetic energy of spacecraft
for ii = 1:length(t2)
    H(:,ii) = Im*X_Att(4:6, ii);                                      %Angular momentum 
    Hmag(:,ii) = norm(H(:,ii));                                       %Magnitude of angular momentum
    T(:,ii) = 0.5 * (X_Att(4:6, ii)' * H(:,ii));                      %Rotational kinentic energy
    Tmag(:,ii) = norm(T(:,ii));                                       %magnitude of kinetic energy
end

%% Plots
%########## Specific Energy Plot ######################################################
hfig = figure;
ax = gca;
hold on;
sp_diff = abs(sp_e - sp_e2);
yyaxis right
plot(t, sp_e, 'LineWidth', 1.5, 'DisplayName', 'Specific energy', 'Color', "#78DCE8");
set(gca, 'YTick', [])
ax.YColor = 'black';
ax.YColor = 'black';

yyaxis left
plot(t, sp_e2, 'LineWidth', 1.5, 'LineStyle', '-','DisplayName', 'Specific energy', 'Color', "#FF6188");
xlabel('time $t$ $(s)$');
ylabel('specific energy $\zeta$ $(Jkg^{-1})$', 'Color', "#3ce0d5");
ax.YColor = 'black';
ax.YColor = 'black';
xlim([min(t), max(t)]);
lgd = legend('Location', 'best');
lgd.EdgeColor = 'black';
xlim([min(t), max(t)]);

legend({'$\zeta_{1}=-\frac{\mu }{2a}$', '${\zeta_{2}}=\frac{\nu ^{2}}{2}-\frac{\mu }{r}$'});
grid on;

picturewidth = 17.6;
hw_ratio = .75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

% print(hfig, 'specific energy plot', '-dpng', '-painters')


%########## VELOCITY AND POSITION ERROR PLOT ######################################################
hfig = figure;
ax = gca;

yyaxis left
plot(t, vdiff, 'Color', "#FF6188", 'LineWidth', 1, 'DisplayName', '$\epsilon_v$');
xlabel('time $t$ $(s)$')
ylabel('velocity error $\epsilon_v$ $(ms^{-1})$')
ax.YColor = 'black';

yyaxis right;
plot(t, rdiff, 'Color', "#78DCE8", 'LineWidth', 1, 'DisplayName', '$\epsilon_r$');
ylabel('position error $\epsilon_r$ $(m)$', 'Color', "#e62740");
ax.YColor = 'black';

lgd = legend('Location', 'best');
lgd.EdgeColor = 'black';
xlim([min(t), max(t)]);
grid on;

picturewidth = 17.6;
hw_ratio = 0.75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'velocity and error plot', '-dpng', '-painters')

%########## ANOMALY VS TIME PLOT ######################################################
hfig = figure;
yyaxis left;
hold on;
plot(t, E_matrix, 'LineWidth', 1.5, 'LineStyle', '-', 'DisplayName', 'EA', 'Color', "#FF6188");
plot(t, TA, 'LineWidth', 1.5, 'LineStyle', '-',  'DisplayName', 'TA', 'Color', "#78DCE8");
plot(t, MAt, 'LineWidth', 1.5, 'LineStyle', '-', 'DisplayName', 'MA', 'Color', "#AB9DF2");
xlabel('time $t$ $(s)$');
ylabel('angle $(rad)$', 'Color', 'k');
lgd = legend('Location', 'best');
xlim([min(t), max(t)]);
ax = gca;
ax.YColor = 'black';
grid on;
yyaxis right;
Ediff = E_matrix - MAt;
TAdiff = TA - MAt;
plot(t, Ediff, 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'EA deviation from MA', 'Color', "#FF6188");
plot(t, TAdiff, 'LineWidth', 1.5, 'LineStyle', '-.', 'DisplayName', 'TA deviation from MA', 'Color', "#78DCE8");
ylim([-5e-3 5e-3]);
ylabel('deviation angle $(rad)$', 'Color', 'k')
ax.YColor = 'black';

picturewidth = 17.6;
hw_ratio = 0.75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
lgd.FontSize = 10;
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'Anomaly plot', '-dpng', '-painters')

%########## ORBIT PLOTS ###############################################################

% create orbit plot
[Xe,Ye,Ze] = sphere(50);

% Create figure
hfig = figure;

% Plot Earth
mesh(Re*Xe, Re*Ye, Re*Ze, 'EdgeColor', 'k'); hold on; axis equal;
xlabel('$I_{ECI} (km)$');
ylabel('$J_{ECI} (km)$');
zlabel('$K_{ECI} (km)$');

% Plot ECI / ECEF axes
quiver3(0, 0, 0, 1, 0, 0, 1e4, 'k', 'LineWidth', 2);    % I/X-axis
quiver3(0, 0, 0, 0, 1, 0, 1e4, 'k', 'LineWidth', 2);    % J/Y-axis
quiver3(0, 0, 0, 0, 0, 1, 1e4, 'k', 'LineWidth', 2);    % K/Z-axis
text(1.22e4, 0, 0, '$I$', 'Color', 'k')
text(0, 1.1e4, 0, '$J$', 'Color', 'k')
text(0, 0, 1.1e4, '$K$', 'Color', 'k')

% Plot Trajectory in ECI frame
plot3(X(1,:), X(2,:), X(3,:),'LineWidth', 2, 'Color', "#FF6188");
plot3(X(1,1), X(2,1), X(3,1), 'ok', 'MarkerFaceColor', "#A9DC76"); %start indicator
plot3(X(1,end), X(2,end), X(3,end), 'ok', 'MarkerFaceColor', "#A9DC76");%end indicator
view(130, 30);

% Calculate eccentricity vector, specific angular momentum vector, and complete the triad
r = X(1:3, 1);
v = X(4:6, 1);
h = cross(r, v);
e = cross(v, h)/GM - r/norm(r);

ie = e/norm(e);
ih = h/norm(h);
ip = cross(ih, ie)/norm(cross(ih, ie));

quiver3(0, 0, 0, ie(1), ie(2), ie(3), 1e4, 'Color', "#78DCE8", 'LineWidth', 2);
quiver3(0, 0, 0, ip(1), ip(2), ip(3), 1e4, 'Color', "#78DCE8", 'LineWidth', 2);
quiver3(0, 0, 0, ih(1), ih(2), ih(3), 1e4, 'Color', "#78DCE8", 'LineWidth', 2);
text(1e4*ie(1), 1e4*ie(2), 1.2e4*ie(3), '$i_{e}$', 'FontSize', 20, 'Interpreter', 'tex', 'Color', "#78DCE8")
text(1e4*ip(1), 1e4*ip(2), 1.1e4*ip(3), '$i_{p}$', 'FontSize', 20, 'Interpreter', 'tex', 'Color', "#78DCE8")
text(1e4*ih(1), 1.15e4*ih(2), 1e4*ih(3), '$i_{h}$', 'FontSize', 20, 'Interpreter', 'tex', 'Color', "#78DCE8")

% Plot Trajectory in ECEF frame
plot3(X_ECEF(1,:), X_ECEF(2,:), X_ECEF(3,:), 'LineWidth', 2, 'Color', "#78DCE8");
plot3(X_ECEF(1,1), X_ECEF(2,1), X_ECEF(3,1), 'ok', 'MarkerFaceColor', "#A9DC76");
plot3(X_ECEF(1,end), X_ECEF(2,end), X_ECEF(3,end), 'ok', 'MarkerFaceColor', "#FF6188");
grid off 

picturewidth = 17.6;
hw_ratio = .9;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'ECI plot', '-dpng', '-vector')

%######################### POLHODE PLOT ############################
%create orbit plot
[H1,H2,H3] = sphere(50);

% Create figure
hfig = figure;

% Plot Earth
mesh(Hmag(1)*H1, Hmag(1)*H2, Hmag(1)*H3, 'EdgeColor', 'k'); hold on; axis equal;
axis off

%plot momentum vector 
plot3(H(1,:), H(2,:), H(3,:),'LineWidth', 5, 'Color', "#FF6188");

% Plot H axes
quiver3(0, 0, 0, 1, 0, 0, 10, 'k', 'LineWidth', 2);    
quiver3(0, 0, 0, 0, 1, 0, 10, 'k', 'LineWidth', 2);   
quiver3(0, 0, 0, 0, 0, 1, 10, 'k', 'LineWidth', 2);    
text(12, 1, 0, '$H_{1}$', 'Color', 'k')
text(0, 10.25, 0, '$H_{2}$', 'Color', 'k')
text(0, 0, 10.25, '$H_{3}$', 'Color', 'k')
grid off 
view(135, 20);

picturewidth = 17.6;
hw_ratio = 1;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

[h1, h2, h3] = ellipsoid(0, 0, 0, sqrt(2*Im(1,1)*T(1)), sqrt(2*Im(2,2)*T(1)), sqrt(2*Im(3,3)*T(1)), 50);
surf(h1, h2, h3, 'LineStyle', 'none', 'FaceColor', "#78DCE8"); hold on; axis equal;

% print(hfig, 'PoleHhode plot', '-dpng', '-painters')

%################ EULER ANGLES PLOT ##########
hfig = figure;
ax = gca;
hold on

plot(t2, X_Att(1,:), 'Color', "#FF6188", 'LineWidth', 1.5, 'DisplayName', '$Yaw$ $\psi$');
plot(t2, X_Att(2,:), 'Color', "#78DCE8", 'LineWidth', 1.5, 'DisplayName', '$Pitch$ $\theta$');
plot(t2, X_Att(3,:), 'Color', "#AB9DF2", 'LineWidth', 1.5, 'DisplayName', '$Roll$ $\phi$');
xlabel('time $t$ $(s)$')
ylabel('Euler angles $\Psi$ $(rad)$')
ax.YColor = 'black';

lgd = legend('Location', 'best');
lgd.EdgeColor = 'black';
xlim([min(t2), max(t2)]);
grid on;

picturewidth = 17.6;
hw_ratio = 0.75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
lgd.FontSize = 15;
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'Euler Ang plot', '-dpng', '-painters')

%################ ANGULAR VELOCITY PLOT ##########
hfig = figure;
ax = gca;
hold on

plot(t2, X_Att(4,:), 'Color', "#FF6188", 'LineWidth', 1.5, 'DisplayName', '$\omega_{\psi}$');
plot(t2, X_Att(5,:), 'Color', "#78DCE8", 'LineWidth', 1.5, 'DisplayName', '$\omega_{\theta}$');
plot(t2, X_Att(6,:), 'Color', "#AB9DF2", 'LineWidth', 1.5, 'DisplayName', '$\omega_{\phi}$');
xlabel('time $t$ $(s)$')
ylabel('angular velocities $\omega$ $(rad\,s^{-1})$')
ax.YColor = 'black';

lgd = legend('Location', 'best');
lgd.EdgeColor = 'black';
xlim([min(t2), max(t2)]);
grid on;

picturewidth = 17.6;
hw_ratio = 0.75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
lgd.FontSize = 15;
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'Ang velocities plot', '-dpng', '-painters')

%################ ANGULAR MOMENTUM AND ROTATIONAL KINETIC ENERGY PLOT ##########
hfig = figure;
ax = gca;
hold on

plot(t2, Hmag, 'Color', "#FF6188", 'LineWidth', 1.5, 'DisplayName', '$H$');
xlabel('time $t$ $(s)$')
ylabel('angular momentum $H$ $(kg\,m^{2}\,s^{-1})$')
ax.YColor = 'black';

lgd = legend('Location', 'best');
lgd.EdgeColor = 'black';
xlim([min(t2), max(t2)]);
grid on;

picturewidth = 17.6;
hw_ratio = .75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
lgd.FontSize = 15;
ax.YAxis.FontSize = 10;
ax.YLabel.FontSize = 20;
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'Ang momentum plot', '-dpng', '-painters')

hfig = figure;
ax = gca;
hold on

plot(t2, Tmag, 'Color', "#78DCE8", 'LineWidth', 1.5, 'DisplayName', '$T$');
xlabel('time $t$ $(s)$')
ylabel('kinetic energy $T$ $(J)$')
ax.YColor = 'black';

lgd = legend('Location', 'best');
lgd.EdgeColor = 'black';
xlim([min(t2), max(t2)]);
grid on;

picturewidth = 17.6;
hw_ratio = .75;
set(findall(hfig, '-property', 'FontSize'), 'FontSize', 20)
lgd.FontSize = 15;
ax.YAxis.FontSize = 10;
ax.YLabel.FontSize = 20;
set(findall(hfig, '-property', 'Box'), 'Box', 'off')
set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig, 'Position');
set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])
% print(hfig, 'Rot kin En plot', '-dpng', '-painters')

