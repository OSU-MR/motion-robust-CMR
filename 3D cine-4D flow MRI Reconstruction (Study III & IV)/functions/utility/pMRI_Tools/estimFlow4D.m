function [ flow ] = estimFlow4D( data, mask )
%ESTIMFLOW4D Estimate the flow parameters from 4D flow data using a single
%slice of the data volume
%function stats = BlandAltman( s1, s2, new_fig,)
%
%       Creates a Bland-Altman plot for  two dataset s1 and s2.
% 
% Required Inputs:
%   data:   A structure with fields data.thetaX, data.thetaY, and
%           data.thetaZ, where each field is a 4D array of velocities 
% 
%   mask:   A 4D mask the same size as the velocity maps
%
% Outputs:
%   flow:  output object.
%
%                 sv_z: stroke volume in the Z direction
%                 sv_y: stroke volume in the X direction
%                 sv_x: stroke volume in the Z direction
%         theta_z_mean: mean velocity in the ROI in the Z direction
%          theta_z_max: max velocity in the ROI in the Z direction
%         theta_y_mean: mean velocity in the ROI in the Y direction
%          theta_y_max: max velocity in the ROI in the Y direction
%         theta_x_mean: mean velocity in the ROI in the X direction
%          theta_x_max: max velocity in the ROI in the X direction
%     theta_speed_mean: mean speed in the ROI
%      theta_speed_max: max speed in the ROI
% 
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        adamrich12@gmail.com
% Last update:  4/10/2017
%**************************************************************************
    theta_z_masked = bsxfun(@times,mask,data.thetaZ);
    theta_y_masked = bsxfun(@times,mask,data.thetaY);
    theta_x_masked = bsxfun(@times,mask,data.thetaX);
    
    theta_speed = sqrt( data.thetaX.^2 + data.thetaY.^2 + data.thetaZ.^2 );
    theta_speed_average = theta_speed;
    
    % Average over a small region in the phase
    kernel = [0,1,0;1,1,1;0,1,0];
    kernel2 = [0,0.5,0;0.5,1,0.5;0,0.5,0];
    kernel = kernel/sum(kernel(:));
    kernel2 = kernel2/sum(kernel2(:));
    
    for slice = 1:size(theta_speed,3)
        for phase = 1:size(theta_speed,4)
            theta_speed_average(:,:,slice,phase) = conv2(theta_speed_average(:,:,slice,phase),kernel,'same');
        end
    end
    
    for slice = 1:size(theta_speed,3)
        for phase = 1:size(theta_speed,4)
            theta_speed_average2(:,:,slice,phase) = conv2(theta_speed_average(:,:,slice,phase),kernel2,'same');
        end
    end

    theta_speed = bsxfun(@times,theta_speed,mask);
    theta_speed_average = bsxfun(@times,theta_speed_average,mask);
    theta_speed_average2 = bsxfun(@times,theta_speed_average2,mask);
    
    % sum the flow to esitmate stroke volume
    flow.sv_z = sum(theta_z_masked(:));
    flow.sv_y = sum(theta_y_masked(:));
    flow.sv_x = sum(theta_x_masked(:));

    for ind = 1:size(data.thetaZ,4)
        tmp = theta_z_masked(:,:,:,ind);
        tmp = tmp(abs(tmp)>0);
        flow.theta_z_mean(ind) = mean(tmp);
        flow.theta_z_sum(ind) = sum(tmp);
        [~,max_ind] = max(abs(tmp));
        if ~isempty(max_ind)
            flow.theta_z_max(ind) = tmp(max_ind);
        end
        
        tmp = theta_y_masked(:,:,:,ind);
        tmp = tmp(abs(tmp)>0);
        flow.theta_y_mean(ind) = mean(tmp);
        flow.theta_y_sum(ind) = sum(tmp);
        [~,max_ind] = max(abs(tmp));
        if ~isempty(max_ind)
            flow.theta_y_max(ind) = tmp(max_ind);
        end
        
        tmp = theta_x_masked(:,:,:,ind);
        tmp = tmp(abs(tmp)>0);
        flow.theta_x_mean(ind) = mean(tmp);
        flow.theta_x_sum(ind) = sum(tmp);
        [~,max_ind] = max(abs(tmp));
        if ~isempty(max_ind)
            flow.theta_x_max(ind) = tmp(max_ind);
        end
        
        tmp = theta_speed(:,:,:,ind);
        tmp = tmp(abs(tmp)>0);
        flow.theta_speed_mean(ind) = mean(tmp);
        [~,max_ind] = max(abs(tmp));
        if ~isempty(max_ind)
            flow.theta_speed_max(ind) = tmp(max_ind);
        end
        
        tmp = theta_speed_average(:,:,:,ind);
        tmp = tmp(abs(tmp)>0);
%         flow.theta_speed_max_avg(ind) = mean(tmp);
        [~,max_ind] = max(abs(tmp));
        if ~isempty(max_ind)
            flow.theta_speed_max_avg(ind) = tmp(max_ind);
        end
        
        tmp = theta_speed_average2(:,:,:,ind);
        tmp = tmp(abs(tmp)>0);
%         flow.theta_speed_max_avg2(ind) = mean(tmp);
        [~,max_ind] = max(abs(tmp));
        if ~isempty(max_ind)
            flow.theta_speed_max_avg2(ind) = tmp(max_ind);
        end
    end

end

