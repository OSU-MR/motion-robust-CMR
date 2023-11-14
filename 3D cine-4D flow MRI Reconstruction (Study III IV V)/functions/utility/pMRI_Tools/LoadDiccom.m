function [y,info] = LoadDiccom( directory_name )
% function [y,info] = LoadDiccom( diccom_name )
%
% LoadDiccom loads a diccom image series with and returns the images series
% as a 3D matrix
%
% Inputs:
%   directory_name:    A string containing the name of the folder
%                       containing the diccom images
%
% Outputs:
%   y:      The image series stored in a 3D matrix where the thrid
%           dimension has size equal to the length of the image sequence
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  9/10/2014
%**************************************************************************

direct = dir(directory_name);
cur_dir = pwd;
cd(directory_name);

% Main loop for all the frames
for ind = 3:length(direct);
    
    % Get dicominfo
    info = dicominfo(direct(ind).name);

    % Read dicom file and store it in a 3D array
    y(:,:,info.InstanceNumber) = double(dicomread(direct(ind).name));
end
cd(cur_dir);
end