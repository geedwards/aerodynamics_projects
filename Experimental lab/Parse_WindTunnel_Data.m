%% ASEN 3111 MATLAB Code to Process the Wind Tunnel Data from the Pilot Lab Wind Tunnel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 3111 Aerodynamics Spring Semester 2020
% To run:
%   - Create a folder on your computer that contains just the data from the
%   wind tunnel, downloaded from the Google drive folder. Make sure all of
%   the data is in .csv format.
%   - Save code in the same folder as the downloaded Wind Tunnel .csv files
%   - When you run it will ask you what folder the data is in
% Output: 
%   - .mat structure of the data consolidated. You can save that structure
%     or run the processing script from your main code.
%
% Author: Sara Swenson
% Date: 02/25/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Load Data from Current Folder
% Specify folder that contains the data:
FDir = uigetdir(pwd,'Select the main folder that contains ALL of the wind tunnel data');
MFolder = dir(fullfile(FDir));
dirFlags = [MFolder.isdir];
subFolders = MFolder(dirFlags);
% Check to see how data was saved to folder (same as Google drive or all in
% one folder
% Grab all the .csv files from that folder:
if any(strcmp({subFolders.name},'Airfoil')|strcmp({subFolders.name},'Cylinder'))
    files = dir(fullfile(FDir,'*/*.csv'));
else
    files = dir(fullfile(FDir,'*.csv'));
end

% Create new MATLAB structure to contain mean and standard deviation of all
% data
WTData = {};
% Create a flag for bad data
bad = 0;
% Loop through each data file and store properties for configuration in table
for i = 1:length(files)
    % Get data file path
    fpath = fullfile(files(i).folder,files(i).name);
    % Load data file described by variable fpath into an array
    data = csvread(fpath,1,0);
    
    % Parse file name and place in data structure
    parse_fn = split(files(i).name,'_');
    WTData(i).name = files(i).name;
    WTData(i).model = parse_fn{2};
    WTData(i).speed = str2double(parse_fn{3}(2:end));
    WTData(i).x_location = str2double(parse_fn{4}(2:end));
    WTData(i).repetition = str2double(parse_fn{5}(2:end-4));
    
    % Figure out the indices of where the ELD probe is at new y location by
    % looking at the difference in location between data points and
    % thresholding at 0.1mm
    step_i =[1;  find(diff(data(:,end))>0.2)+1; length(data)+1];
    
    % Check for bad data. ie: didn't sample at high enough rate or didn't
    % take enough samples
    if (length(step_i)>10 & length(step_i)<21) | length(step_i)>=22
        % If more/less than specified number data points were taken
        bad = bad+1;
        WTData(i).bad_data = 1;
    elseif length(data)<10000
        % If not enough data was taken
        bad = bad+1;
        WTData(i).bad_data = 2;
    elseif length(step_i)<=10
        bad = bad+1;
        WTData(i).bad_data = 2;
    else 
        WTData(i).bad_data = 0;
    end
 

    for j = 1:length(step_i)-1
        WTData(i).x_mm_mean(j) = mean(data(step_i(j):step_i(j+1)-1,27));
        WTData(i).x_mm_std(j) = std(data(step_i(j):step_i(j+1)-1,27));
        WTData(i).y_mm_mean(j) = mean(data(step_i(j):step_i(j+1)-1,28));
        WTData(i).y_mm_std(j) = std(data(step_i(j):step_i(j+1)-1,28));
        WTData(i).T_K_mean(j) = mean(data(step_i(j):step_i(j+1)-1,2));
        WTData(i).T_K_std(j) = std(data(step_i(j):step_i(j+1)-1,2));
        WTData(i).rhoatm_kgm3_mean(j) = mean(data(step_i(j):step_i(j+1)-1,3));
        WTData(i).rhoatm_kgm3_std(j) = std(data(step_i(j):step_i(j+1)-1,3));
        WTData(i).airspeed_ms_mean(j) = mean(data(step_i(j):step_i(j+1)-1,4));
        WTData(i).airspeed_ms_std(j) = std(data(step_i(j):step_i(j+1)-1,4));
        WTData(i).Patm_Pa_mean(j) = mean(data(step_i(j):step_i(j+1)-1,1));
        WTData(i).Patm_Pa_std(j) = std(data(step_i(j):step_i(j+1)-1,1));
        WTData(i).Pwindtunnel_Pa_mean(j) = mean(data(step_i(j):step_i(j+1)-1,5));
        WTData(i).Pwindtunnel_Pa_std(j) = std(data(step_i(j):step_i(j+1)-1,5));
        WTData(i).Pwake_Pa_mean(j) = mean(data(step_i(j):step_i(j+1)-1,6));
        WTData(i).Pwake_Pa_std(j) = std(data(step_i(j):step_i(j+1)-1,6));
    end
    
        
    WTData(i).dP = WTData(i).Pwindtunnel_Pa_mean- WTData(i).Pwake_Pa_mean(j);
    if mean(WTData(i).rhoatm_kgm3_mean)<0.5
        WTData(i).bad_data = 2;
        bad = bad+1;
    elseif min((WTData(i).dP./max(WTData(i).dP)))>0.5
        WTData(i).bad_data = 2;
        bad = bad+1;
    end
    
    
end

save('Parsed_WindTunnel_Data.mat','WTData')