function [eval] = simpson(equation, inputs)
    %% simpson
    % Integrate the given equation from start to end using Simpson's rule.
    % Author: Grace Edwards
    % Collaborators:
    % Date:
    %
    % Set up some variables to keep track of things
    step_size = inputs(2) - inputs(1);
    sum = 0;
    % Loop through all the elements in the array
    for i = 2:length(inputs)
        sum = sum + step_size/3*(4*equation(i-1) + equation(i));
    end
    eval = sum;
    
end