function [cL,cDi,e] = PLLT(b,a0_r,a0_t,c_r,c_t,aero_r,aero_t, geo_r, geo_t,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------Prandtl Lifting Line Theory Code ------------------%
% Code to return the coefficient of lift, coefficient of drag, and        %
% efficiency factor of a finite wing.                                     %
% Assumes wing can be approximated as trapezoidal and that chord and angle%
% of attack varies linearly between the root and the tip.                 %
%                                                                         %
% INPUTS: b    - wing span                                                %
%         a0_r - inifinite wing lift slope at the root (airfoil data)     %
%         a0_t - inifinite wing lift slope at the tip (airfoil data)      %
%         c_r  - chord length at the root                                 %
%         c_t  - chord length at the tip                                  %
%         aero_r - aerodynamic twist at root                              %
%         aero_t - aerodynamic twist at tip                               %
%         geo_r  - geometric twist at root                                %
%         geo_t  - geometric twist at tip                                 %
%         N      - number of terms in sine series                         %
% OUTPUS: cL  - coefficient of lift                                       %
%         cDi - induced coefficient of drag (drag due to lift)            %
%         e   - efficiency factor                                         %
%                                                                         %
% AUTHORS: S. Swenson, 04/27/2020                                         %
% Work based off of notes from ASEN 3111 by Professor B. Argrow.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assemble the wing parameters along the span of the wing
S = (c_t+c_r)*b/2; % Full Wing Area
AR = b^2/S; % Full Wing Aspect Ratio
for i = 1:N
    % Theta: location along wing
    theta(i) = i*pi/(2*N);
    % Wing Taper at theta(i)
    c(i) = c_r-(c_r-c_t)*cos(theta(i));
    % Variable cross-sectional lift slope at theta(i)
    a0(i) = a0_r-(a0_r-a0_t)*cos(theta(i));
    % Aerodynamic twist at theta(i) - airfoil dependent
    aero(i) = aero_r-(aero_r-aero_t)*cos(theta(i));
    % Geometric twist at theta(i)
    geo(i) = geo_r-(geo_r-geo_t)*cos(theta(i));
end

%% Solve for Coefficients
% For each theta(i), must calculate the N odd coefficients to solve the
% matrix equation x*A = B
% B = alpha(theta)-alpha L=0(theta)
% A = set of A coefficients for each theta
% x = matrix multiplying A to get B
x = zeros(N,N);
for i = 1:N
    B(i) = geo(i)-aero(i);
    for j = 1:N
        x(i,j) = (4*b)/(a0(i)*c(i))*sin((2*j-1)*theta(i))+...
            (2*j-1)*(sin((2*j-1)*theta(i))/sin(theta(i)));
    end
end

%Solve the matrix equation xA = B for A
A = x\B';
% Can check wheter this is performed correctly with
test = x*A-B;
% should be matrix of zeros.

%% Calculate lift and drag coefficients
% lift coefficent
cL = A(1)*pi*AR;

% Induced drag factor
delta = 0;
for i = 2:N
    delta = delta+(2*i-1)*(A(i)/A(1))^2;
end

% efficiency factor
e = 1/(1+delta);

% Induced drag coefficient
cDi = cL.^2/(pi*e*AR);
end

    