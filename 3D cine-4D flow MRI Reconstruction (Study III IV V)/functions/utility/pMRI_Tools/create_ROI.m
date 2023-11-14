function [ mask,position] = create_ROI( x, type, position_in )
%function [ mask ] = create_ROI( x, varargin )
%
% create_ROI creates a a time varying mask by cycling through a dynamic
% image sequence and having the user highlight a regeon of interest in each
% frame
%
%
% Inputs:
%   x:      A dynamic image sequence.  The first two dimension represent
%           the spatial dementions and the third time.
%
% Optional Inputs:
%   ROI_type:   An optional input string to determine the type of ROI to 
%               draw.  Options are 'ellipse', 'polygon', and 'rectangle'
%
% Output
%   mask:   A binary mask that is the same size as the input x.  zeros
%           representing the locations outside the ROI, ones locations 
%           inside the ROI.
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  7/22/2014
%**************************************************************************

% Check for the optional input
if nargin ==1
    ROI_type = 'ellipse';
    position_in = [];
elseif nargin ==2;
    ROI_type = type;
    position_in = [];
elseif nargin ==3;
    ROI_type = type;
else
    error('Incorrect number of inputs')
end

% Initialize mask
mask = zeros(size(x));

switch lower(ROI_type)
    case {'ellipse'}
        for frame_ind = 1:size(x,3)
            
            % Create Figure for the current frame
            figure(1)
            clf
            class(x)
            imagesc((x(:,:,frame_ind)));
            colormap(gray)
            % Create an im class and wait for the user input
            if frame_ind ==1
                fprintf('Zoom to area of Interest and press ENTER to continue\n')
                pause
                Xl = xlim;
                Yl = ylim;
            else
                xlim(Xl);
                ylim(Yl);
            end
            fprintf('Draw ROI and press ENTER to continue\n')
            get_roi = imellipse;
            if frame_ind ==1
                if ~isempty(position_in)
                    get_roi.setPosition(position_in);
                end
                pause
                position = get_roi.getPosition;
                mask(:,:,frame_ind) = get_roi.createMask;
            else
                get_roi.setPosition(position);
                pause
                position = get_roi.getPosition;
                mask(:,:,frame_ind) = get_roi.createMask;
            end
%             pause
%             position = get_roi.getPosition + 5;
%             mask(:,:,frame_ind) = get_roi.createMask;
            
            clc;
        end
        
    case {'polygon'}
        for frame_ind = 1:size(x,3)
            
            % Create Figure for the current frame
            figure(1)
            clf
            imagesc((x(:,:,frame_ind)));
            colormap(gray)
            
            % Create an im class and wait for the user input
            fprintf('Zoom to area of Interest and press ENTER to continue')
            pause
            fprintf('Draw ROI and press ENTER to continue')
            get_roi = impoly;
            pause
            mask(:,:,frame_ind) = get_roi.createMask;
            clc;
        end
        
    case {'rectangle'}
        for frame_ind = 1:size(x,3)
            
            % Create Figure for the current frame
            figure(1)
            clf
            imagesc((x(:,:,frame_ind)));
            colormap(gray)
            
            % Create an im class and wait for the user input
            fprintf('Zoom to area of Interest and press ENTER to continue')
            pause
            fprintf('Draw ROI and press ENTER to continue')
            get_roi = imrect;
            pause
            mask(:,:,frame_ind) = get_roi.createMask;
            clc
        end
        
    otherwise
        error('Unknown ROI type specified')
end