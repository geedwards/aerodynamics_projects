%% ASEN 3111 - Computational Assignment 02 - Main
% Plots the streamlines, equipotential lines, and pressure contours of an
% airfoil using thin airfoil theory and superposing N vortices to model the
% flow. Plots the error of the calculation vs. the number of panels used to
% calculate it. Varies the cord length, angle of attack, and freestream
% velocity in order to compare the plot outputs of each variation, and see
% how the variable affects the output. This will take about <time> to run,
% and will output 11 plots in total. 
%
% Author: Grace Edwards
% Collaborators: Keith Poletti
% Date: 2/27/2020
clear all;
close all;
clc;
c = 2; % meters
alpha = 12; % degrees
V_inf = 68; % m/s
p_inf = 101.3*10^3; % Pa
rho_inf = 1.225; %kg/m^3
% Test out using the parameters given in the lab document, as above, with
% 200 panels. Want both plots and stream function here
fprintf('Calculating first airfoil flow...\n');
Plot_Airfoil_Flow(c, alpha, V_inf, p_inf, rho_inf, 200, true, true);
%% Error study
% Set a very high number of panels as the best possible approximation, and
% run the function
fprintf('Calculating the flow with maximum panels...\n');
N_max = 5000;
[str_max, eq_max, c_p_max, vel_max] = Plot_Airfoil_Flow(c, alpha, V_inf, p_inf, rho_inf, N_max, false, false);
% Create empty lists in which to populate the error in both pressure and
% velocity
diff_pressure = [];
diff_velocity = [];
% Create the points of N panels where we want to calculate
n = linspace(10, 2000, (2000-10)/5+1);
% Loop through the points and calculate the pressure and velocity
fprintf('Comparing to lower amounts of panels...\n');
for i=1:numel(n)
    % Don't want plots or stream function calculated
    [str, eq, c_p, vel] = Plot_Airfoil_Flow(c, alpha, V_inf, p_inf, rho_inf, n(i), false, false);
    % Append to each list the RMS difference between this iteration and the
    % ideal, maximum calculation made above
    diff_pressure(i) = sqrt(sum(sum(abs(c_p-c_p_max).^2))./n(i));
    diff_velocity(i) = sqrt(sum(sum(abs(vel-vel_max).^2))./n(i));
    if mod(n(i), 100) == 0
        fprintf('Calculation %i reached...\n', n(i));
    end
end
% Plot the differences in pressure vs. N
figure(4)
plot(n, diff_pressure);
xlabel('Number of panels'); ylabel('RMS');
title('Error of the Coefficient of Pressure');
% Plot the differences in velocity vs. N
figure(5)
plot(n, diff_velocity);
xlabel('Number of panels'); ylabel('RMS');
title('Error of the Velocity');
%% Varying cord length, AoA, and V_inf and comparing them
% For each of these, it was chosen that 500 panels was enough to get the
% error very low, but not too much that the computation time was too high.

% First, cord length. Go from one meter less to one meter greater then the
% cord length we started with. In this lab, this will be 1 to 3. 
cord_length = linspace(c - 1, c + 1, 6);
% Loop through the difference cord lengths
fprintf('Calculating the effects of cord length variation...\n');
for j = 1:numel(cord_length)
    % Calculate the plots for each. We want the stream function, but not
    % the plots. 
    [curStr, curEq, curCp, curVel] = Plot_Airfoil_Flow(cord_length(j), alpha, V_inf, p_inf, rho_inf, 500, false, true);
    % Make a grid of x and y based on the maximum cord length 
    [x,y]=meshgrid(linspace(-cord_length(end),2*cord_length(end),100),linspace(-2*cord_length(end),2*cord_length(end),100));
    % Add a subplot to the streamline window. Make sure the axes are the
    % same, and that the rectangle drawn matches the calculation.
    figure(6)
    sgtitle('Streamlines with Varying Cord Length');
    subplot(2, 3, j), contour(x, y, curStr, 20);
    hold on
    subplot(2, 3, j), rectangle('Position', [0 0 cord_length(j) 0.01]);
    xlim([-(c-1) 2*(c+1)]), ylim([-2*(c+1) 2*(c+1)]);
    xlabel('x (m)'); ylabel('y (m)');
    title(strcat('c = ', num2str(cord_length(j)), ' m'));
    % Now add a subplot to the equipotential window
    figure(7)
    sgtitle('Equipotential Lines with Varying Cord Length');
    subplot(2, 3, j), contour(x, y, curEq, 20);
    hold on
    subplot(2, 3, j), rectangle('Position', [0 0 cord_length(j) 0.01]);
    xlim([-(c-1) 2*(c+1)]), ylim([-2*(c+1) 2*(c+1)]);
    xlabel('x (m)'); ylabel('y (m)');
    title(strcat('c = ', num2str(cord_length(j)), ' m'));
end
% Cord length doesn't make much difference to the shape of the streamlines.
% However the disturbance extends further into the flow the longer the cord
% is. 
% Next, vary the angle of attack
% Go twelve degrees below, and twelve above. That might be a little high.
aoa = linspace(alpha - 12, alpha + 12, 6);
% Loop through the six angles of attack
fprintf('Calculating the effects of angle of attack variation...\n');
for k = 1:numel(aoa)
    % Calculate the flow for each of them
    [curStr, curEq, curCp, curVel] = Plot_Airfoil_Flow(c, aoa(k), V_inf, p_inf, rho_inf, 500, false, true);
    [x,y]=meshgrid(linspace(-c,2*c,100),linspace(-2*c,2*c,100));
    % Add to the streamline window using appropriate labelling
    figure(8)
    sgtitle('Streamlines with Varying Angle of Attack');
    subplot(2, 3, k), contour(x, y, curStr, 20);
    hold on
    subplot(2, 3, k), rectangle('Position', [0 0 c 0.01]);
    xlim([-c 2*c]), ylim([-2*c 2*c]);
    xlabel('x (m)'); ylabel('y (m)');
    title(strcat('\alpha = ', num2str(aoa(k)), ' °'));
    % Add to the equipotential window
    figure(9)
    sgtitle('Equipotential Lines with Varying Angle of Attack');
    subplot(2, 3, k), contour(x, y, curEq, 20);
    hold on
    subplot(2, 3, k), rectangle('Position', [0 0 c 0.01]);
    xlim([-c 2*c]), ylim([-2*c 2*c]);
    xlabel('x (m)'); ylabel('y (m)');
    title(strcat('\alpha = ', num2str(aoa(k)), ' °'));
end
% This is where the thin airfoil assumptions begin to break down. There are
% some weird things happening at the aoa=24 plots, and the discontinuity at
% the leading edge is emphasized the higher the aoa.
% Finally, make the freestream velocity comparison
% I chose to go fifty under and fifty over the initial V_inf, I don't know
% why
fprintf('Calculating the effects of freestream velocity variation...\n');
fs_vel = linspace(V_inf - 50, V_inf + 50, 6);
% Loop through the six calculations to be made
for m = 1:numel(fs_vel)
    [curStr, curEq, curCp, curVel] = Plot_Airfoil_Flow(c, alpha, fs_vel(m), p_inf, rho_inf, 500, false, true);
    [x,y]=meshgrid(linspace(-c,2*c,100),linspace(-2*c,2*c,100));
    % Add the stream function to the subplots
    figure(10)
    sgtitle('Streamlines with Varying Freestream Velocity');
    subplot(2, 3, m), contour(x, y, curStr, 20);
    hold on
    subplot(2, 3, m), rectangle('Position', [0 0 c 0.01]);
    xlim([-c 2*c]), ylim([-2*c 2*c]);
    xlabel('x (m)'); ylabel('y (m)');
    title(strcat('V_i_n_f = ', num2str(fs_vel(m)), ' m/s'));
    % Add the equipotential function to the subplots
    figure(11)
    sgtitle('Equipotential Lines with Varying Freestream Velocity');
    subplot(2, 3, m), contour(x, y, curEq, 20);
    hold on
    subplot(2, 3, m), rectangle('Position', [0 0 c 0.01]);
    xlim([-c 2*c]), ylim([-2*c 2*c]);
    xlabel('x (m)'); ylabel('y (m)');
    title(strcat('V_i_n_f = ', num2str(fs_vel(m)), ' m/s'));
end
fprintf('Done!\n');
% It appears as though the streamlines become smoother when the velocity is
% higher. 
