function [eval] = trap_int(equation, inputs)
    %% trap_int
    % Integrate the given equation from start to end using the trapezoidal
    % rule. Assuming the inputs and equation variables are vectors of the 
    % same length, inputs containing x values, equation containing f(x).
    % Also assuming even spacing.
    % Author: Grace Edwards
    % Collaborators:
    % Date:
    %
    % Set up some variables necessary to calculate the integral
    step_size = inputs(2) - inputs(1);
    sum = 0;
    % Loop through all the elements in the array
    for i = 2:length(inputs)
        sum = sum + step_size*(equation(i-1) + equation(i));
    end
    eval = sum;
end