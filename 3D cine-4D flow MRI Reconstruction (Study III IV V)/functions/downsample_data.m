function [ y ] = downsample_data( DATA, patterns )
%DOWNSAMPLE_DATA Summary of this function goes here
%   Detailed explanation goes here
% Find the locations of the data points

% 3D Data
if ndims(patterns)==3
    
%     DATA = bsxfun(@times,DATA,permute(patterns,[1,2,4,3]));
%     y = DATA(abs(DATA)>0);
    for frame_ind = 1:size(patterns,3)
        mask_pattern{frame_ind} = find(patterns(:,:,frame_ind)>0); % find the nonzero components in the mask
    end

    % Downsample
    y = [];
    for frame_ind = 1:size(patterns,3)
        for coil_ind = 1:size(DATA,3)
           tmp = DATA(:,:,coil_ind,frame_ind);
           tmp = tmp(mask_pattern{frame_ind});
           y = [y;tmp];
        end
    end
% 4D data    
elseif ndims(patterns)==4
    DATA = bsxfun(@times,DATA,permute(patterns,[1,2,3,5,4]));
    patterns = repmat(permute(patterns,[1,2,3,5,4]),[1,1,1,size(DATA,4),1]);
    y = DATA(patterns>0);
end
