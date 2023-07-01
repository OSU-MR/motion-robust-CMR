function [v] = EstimFlowParam(mask,velocity)
%ESTIMFLOWPARAM Summary of this function goes here
%   Detailed explanation goes here

% Create a function handle for bsxfun
point_by_point = @(x,y) x.*y;

% Pull out data from a ROI
v_mask = bsxfun(point_by_point,velocity,mask);
for ind = 1:size(velocity,3)
   tmp1 = v_mask(:,:,ind);
   % Filter image to average peak velocity over small range
%    filt = ones([3,3])/9;
   filt = [0,1,0;1,1,1;0,1,0]/5;
   tmp_conv = conv2(tmp1,filt,'same');
   tmp = tmp1(find(tmp1~=0));
   [~,peak_ind] = max(abs(tmp));
   [~,peak_ind_conv] = max(abs(tmp_conv(:)));
   if ~isempty(peak_ind)
        v.peak_avg(ind) = tmp_conv(peak_ind_conv);
        v.peak(ind) = tmp(peak_ind);
        v.mean(ind) = mean(tmp);
        v.flow(ind) = sum(tmp);
        if ind>1
            v.sum(ind) = v.sum(ind-1) + sum(tmp);
        else
            v.sum(ind) = sum(tmp);
        end
        v.var(ind) = var(tmp);
   else
        v.peak_avg(ind) = 0;
        v.peak(ind) = 0;
        v.mean(ind) = 0;
        v.flow(ind) = 0;
        v.sum(ind) = 0;
        v.var(ind) = 0;
   end

   v.flow(ind) = sum(tmp);
end

v.total_flow = sum(v.flow);
% v_lsq.total_flow = sum(v_lsq.flow);

end

