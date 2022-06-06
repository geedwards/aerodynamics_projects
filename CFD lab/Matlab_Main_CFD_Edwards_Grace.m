%% CFD Lab ASEN 3111
% Author: Grace Edwards
% Submitted 4/28/2020
% Utilizing code from Sarah Swenson, namely the PLLT calculation on line 28
%
% Calculate and plot various statistics about the Tempest using some basic
% geometry constants and airfoil data. Compares this to the experimental
% data for the wing.
clear all;
close all;
clc;
% Constants
AR = 16.5; 
b = 3.22; % m
% Got these from Selig p. 91 - I eyeballed the graph for the lift slope
a0 = .125; % deg^-1
% Assume the same cord length and airfoil at the root and tip
c = 0.23; % m
aero = -2; % Same zero lift AoA for each
% Enter experimental data
alpha = [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16]; % degrees
C_L = [-0.32438, -0.21503, -0.10881, 0.010503, 0.12155, 0.24163, 0.34336, 0.45256, 0.56037, 0.66625, 0.76942, 0.86923, 0.96386, 1.0441, 1.0743, 1.0807, 1.0379, 1.034, 1.0156, 0.97946];
C_D = [0.044251, 0.033783, 0.028627, 0.025864, 0.024643, 0.025099, 0.025635, 0.02766, 0.030677, 0.034855, 0.040403, 0.04759, 0.057108, 0.070132, 0.090921, 0.11193, 0.13254, 0.15465, 0.20959, 0.25668];
% Calculate using Pradtl Lifting Line code
N = 100;
for i=1:20
    geo = alpha(i);
    [C_L_th(i), C_D_th(i), e_th(i)] = PLLT(b, a0, a0, c, c, aero, aero, geo, geo, N);
end
% Plot coefficients of lift side by side
plot(alpha, C_L);
hold on;
plot(alpha, C_L_th);
% Curve fit the linear part of C_L
% I think the linear part stops at about X = 8 degrees
coeffs = polyfit(alpha(1:14), C_L(1:14), 1);
alpha_Fit = linspace(-5, 8, 1000);
CL_Fit = polyval(coeffs, alpha_Fit);
hold on;
plot(alpha_Fit, CL_Fit);
legend('CFD', 'PLLT', 'Best Fit');
xlabel("\alpha (degrees)"); ylabel("C_L");
title("Coefficient of Lift vs. \alpha");
% Plot the coefficients of drag vs. alpha
figure();
plot(alpha, C_D);
hold on;
plot(alpha, C_D_th);
xlabel("\alpha (degrees)"); ylabel("C_D");
title("Coefficient of Drag vs. \alpha");
legend('CFD', 'PLLT');
LoverD = C_L./C_D;
LoverD_th = C_L_th./C_D_th;
% Graph the L/D for experimental and computational
% The one using PLLT looks very weird because it only counts induced drag
figure();
plot(alpha, LoverD);
hold on;
plot(alpha, LoverD_th);
legend('CFD', 'PLLT');
xlabel("\alpha (degrees)"); ylabel("L/D");
title("Lift over Drag vs. \alpha");
% Graph the drag polar from CFD data
figure()
plot(C_L, C_D);
title('Drag Polar');
xlabel('C_L'); ylabel('C_D');
hold on;
% Use data in table 1: Cruise speed
speed = 18; % m/s
rho = 1024; % kg/m^3 according to altitude tables in Anderson
W = 4.5; % kg
L = W;% Lift = weight at cruise
S = 0.63; % m^2, total wing area
C_L_test = L/(.5*rho*speed^2*S);
[closest,closestIndex] = min(abs(C_L - C_L_test));
C_D_test = C_D(closestIndex);
% Estimate oswald efficiency factor
e = 1.78*(1-0.045*AR^0.68)-0.64;
% 3. Plot equations 1.1 and 1.2
K = 1/(pi*AR*e);
C_D_min = min(C_D_th);
low_ind = find(C_D_th==C_D_min);
C_L_minDrag = C_L_th(low_ind(1));
CD1 = C_D_min + K*(C_L_th - C_L_minDrag).^2;
CD_0 = C_D_th(6); % When the AoA is zero. I hard coded this in there
CD2 = CD_0 + C_L_th.^2./(pi*e*AR);
plot(CD1);
hold on;
plot(CD2);
legend('Experimental', 'Equation 1', 'Equation 2');
% Now plot just the experimental because the rest is less accurate
figure()
plot(C_L_test, C_D_test, 'o');
hold on;
plot(C_L, C_D);
title('CFD Drag Polar');
xlabel('C_L'); ylabel('C_D');
% Find the maximum speed
% Solve for when power required equals power available
p_avl = 50*26.1; % current*voltage=watts
% Only use the parts of CL that are above zero, only that makes sense
P_req = sqrt(2.*W.^3.*C_D.^2./(rho.*S.*abs(C_L).^3));
T_req = W./LoverD;
V_inf = P_req./T_req;
max_speed = max(V_inf);
% Plot speed vs. power and find where they intersect
figure()
plot(V_inf, P_req);
hold on;
plot(linspace(min(V_inf), max(V_inf), 20), p_avl*ones(1, 20));
xlabel('Speed (m/s)'); ylabel('Power (W)');
legend('Power Required', 'Power Available');
title('Power vs. Speed');