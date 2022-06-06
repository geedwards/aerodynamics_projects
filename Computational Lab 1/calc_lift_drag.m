function [L_traps, L_simps, D_traps, D_simps] = calc_lift_drag(num_panels)
theta = linspace(0, 2*pi, num_panels); % Radians
R = 0.5; % meters
v_inf = 30; % m/s^2
rho_inf = 1.225; % kg/m^3
p_inf = 101.3*10^3; %Pa
Gamma = 2*pi*R*v_inf;
q_inf = .5*rho_inf*v_inf^2;
% Analytical solution of the coefficient of pressure
C_p = 1 - (4.*sin(theta).^2 + 2.*Gamma.*sin(theta)./(pi.*R.*v_inf) + (Gamma./(2.*pi.*R.*v_inf)).^2);
% Convert this to pressure
p = C_p.* q_inf + p_inf;
% Integrate for lift and drag via two methods
L_traps = trap_int(p.*cos(theta), theta)
L_simps = simpson(p.*cos(theta), theta)
D_traps = trap_int(p.*sin(theta), theta)
D_simps = simpson(p.*sin(theta), theta)
end