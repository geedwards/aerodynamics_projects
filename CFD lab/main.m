%%ASEN 3111 - Computational Assignment 4 - Main
% Computes the span efficiency factor, coefficient of lift, and induced
% coefficient of drag for a specific wing using the Prandtl Lifting Line
% Theory. The wing is a NACA 2412 at the root and a 0012 at the tip. Bcause
% both of these airfoils are symmetric, the even terms of PLLT will drop
% out. An asymmetric airfoil would have to use a modified version of the
% function below.
% 
% Author: Grace Edwards
% 
clear all;
close all;
clc;

b = 100; % ft
a0_t = 0.1; % from 0012 airfoil
a0_r = 0.2; % From the 2412 airfoil
c_t = 5; % feet
c_r = 15; % feet
aero_t = 0; % radians
aero_r = -2; % degrees--from the 2412 airfoil
geo_t = 0; % Degrees. Not that it matters, it'd be zero anyway
geo_r = 5; % Degrees aoa at the root
N = 100;
[e, C_L, C_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
% Multpily by the area and dynamic pressure to get the lift and drag in
% pounds
V_inf = 150/3600*5280; % ft/s
rho = 0.002378; % slugs/ft^3
area = b*(c_r-c_t)/2;
qinf = .5*rho*V_inf^2;
L = C_L*qinf*area;
Di = C_Di*qinf*area;




