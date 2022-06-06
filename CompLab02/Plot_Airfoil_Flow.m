function [StreamFunction, EquiFunction, c_p, vel] = Plot_Airfoil_Flow(c, alpha, V_inf, p_inf, rho_inf, N, plots, calc_strf)
    %Plot_Airfoil_Flow plot streamlines, equipotential lines, and pressure
    %contours of a given thin airfoil
    % Assumes inviscid, incompressible flow as well as a thin airfoil with
    % minimal camber and low angle of attack. 
    % Outputs the following:
    % StreamFunction: Values of psi as a function of x and y.
    % EquiFunction: Values of phi as a function of x and y.
    % c_p: Values of the pressure contour as a function of x and y. 
    % vel: Values of velocity as a function of x and y.
    % Requires inputs for the following:
    % c: the cord length of the airfoil to be modelled
    % alpha: The angle of attack to be modelled. Should be fairly small to
    % uphold the assumptions of thin airfoil theory.
    % V_inf: The freestream velocity of the flow. Should be subsonic, to
    % keep the incompressible assumption valid.
    % p_inf: The pressure far upstream from the airfoil.
    % rho_inf: The density of the air far upstream from the airfoil. Not
    % that it would be different close to the airfoil, we are assuming
    % incompressible flow after all.
    % N: the number of vortices that should be used in the approximation.
    % plots: A boolean telling whether or not to plot the data after
    % calculating it. In the error calculation, this function is run
    % through so many times that you wouldn't want so many plots popping
    % up, so you can disable the plotting feature at leisure. 
    % calc_strf: another boolean, this one asking whether you want to
    % calculate the stream function. It's very computationally expensive,
    % so if you're not explicitly going to use it then you should leave it
    % out.
    % Author: Grace Edwards
    % Collaborators: Keith Poletti
    % Date: 2/27/2020
    %% Finally on to the actual code
    % Set up the variables for the grid of points to plot on
    xGamma = linspace(c/N, c-c/N, N);
    yGamma = zeros(numel(xGamma), 1);
    xmin = -c;
    xmax = 2*c;
    ymin = -2*c;
    ymax = 2*c;
    nx = 100;
    ny = 100;
    % The circulation function of the vortices as a function of x/c
    lilGamma = 2.*deg2rad(alpha).*V_inf.*sqrt((1-xGamma./c)./(xGamma./c));
    % Multiply by delta x to get the circulation
    Gamma = lilGamma.*c/N;
    % Set up the grid, then some anonymous functions to calculate r and
    % theta between two points
    [x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
    radius= @(x,y,x1,y1) sqrt((x-x1).^2+(y-y1).^2);
    theta = @(x, y, x1, y1) atan2(y-y1, x-x1);
    %% Creating the stream function with superposition
    % Start with the equation for uniform flow
    Psi_Uniform = V_inf*(y*cos(deg2rad(alpha))- x*sin(deg2rad(alpha)));
    % Only proceed if it was specifically asked for in the arguments
    if calc_strf == true
        StreamFunction = Psi_Uniform;
        % Loop through the vortices
        for j = 1:numel(Gamma)
            % Add to the stream function the streamlines from the vortex in
            % question
            StreamFunction = StreamFunction + Gamma(j)./(2.*pi).*log(radius(x,y,xGamma(j),yGamma(j)));
        end
    else
        % If it wasn't asked for, return nothing. This way it'll be easy to
        % tell if the variable was allocated incorrectly. 
        StreamFunction = [];
    end
    %% Equipotentials with superposition
    % Start with the uniform flow equation
    EquiFunction = V_inf*(x*cos(deg2rad(alpha)) + y*sin(deg2rad(alpha)));
    % Loop through and add the equations for the vortices, as seen before
    for k = 1:numel(Gamma)
        EquiFunction = EquiFunction - Gamma(k)/(2*pi).*theta(x,y,xGamma(k), yGamma(k));
    end
    
    % Pressure contour
    vel = gradient(EquiFunction);
    c_p = 1- (vel/V_inf).^2;
    
    %% Plot streamfunction at levels
    % Only complete this if it was asked for in the arguments. 
    if plots == true
        % Plot the streamlines, a rectangle to show where the airfoil would
        % be, with 20 contour lines. I think that looks good.
        figure(1)
        rectangle('Position', [0 0 c 0.01]);
        hold on
        contour(x,y,StreamFunction, 20)
        xlabel('x (m)'); ylabel('y (m)');
        title('Streamlines on a Thin Airfoil');
        %Plot the equipotential lines
        figure(2)
        rectangle('Position', [0 0 c 0.01]);
        hold on
        contour(x,y,EquiFunction, 20)
        xlabel('x (m)'); ylabel('y (m)');
        title('Equipotential on a Thin Airfoil');
        % Plot the pressure contours
        figure(3)
        rectangle('Position', [0 0 c 0.01]);
        hold on
        contour(x,y,c_p, 20)
        xlabel('x (m)'); ylabel('y (m)');
        title('Pressure Contour on a Thin Airfoil');
    end

