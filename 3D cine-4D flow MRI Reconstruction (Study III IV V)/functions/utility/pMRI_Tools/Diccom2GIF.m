function Diccom2GIF( diccom_name, file_name, frame_rate, clip_factor)
%function diccom2GIF( diccom_name, file_name, frame_rate, clip_factor)
%
% diccom2GIF creates a GIF from a series of diccom images
%
%
% Inputs:
%   diccom_name:    A string containing the name of the first diccom image
%                   in a series of images.
%
%   file_name:      A string for containing the file name that the GIF will
%                   be saved as.
%
%   frame_rate:     The frame_rate of the GIF in frames per second
%
%   clip_factor:    The amount to clip on the high end of the image.  e.g. 
%                   a clip factor of 2 will plot from 0-1/2*max(x(:))
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  7/14/2014
%**************************************************************************

% Get some file information
file_prefix = diccom_name(1:end-3);
start_num = str2num(diccom_name(end-2:end));
info = dicominfo(diccom_name);              % Get Diccom Info
num_frames = info.CardiacNumberOfImages;    % Number of Frames

% Main loop for all the frames
for ind = 1:num_frames
    % Get the filename
    diccom_name = [file_prefix,num2str(start_num)];
    
    % Get dicominfo
    info = dicominfo(diccom_name);
    
    % Read dicom file and store it in a 3D array
    x(:,:,info.InstanceNumber) = double(dicomread(diccom_name));
    
    % Increment through the files
    start_num = start_num+1;
end

% Create the GIF
create_GIF(x,file_name,frame_rate,clip_factor);
end

