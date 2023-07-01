function [ velocity_norm ] = norm_velocity( velocity)
%NORM_VELOCITY Normalizes the input velocity to be between -pi:pi
%   Detailed explanation goes here

velocity(find(velocity < -pi)) = -pi;
velocity(find(velocity > pi)) = pi;
velocity_norm = (velocity + pi)/(2*pi);

% range = max(velocity(:))-min(velocity(:))
% velocity_norm = velocity/range;
% velocity_norm = abs(velocity_norm+1-max(velocity_norm(:)));

end

