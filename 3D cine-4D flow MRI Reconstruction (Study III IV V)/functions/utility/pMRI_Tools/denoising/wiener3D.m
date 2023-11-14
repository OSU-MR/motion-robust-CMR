function [ f ] = wiener3D( g, noise )
%WIENER3D Summary of this function goes here
%   Detailed explanation goes here

% % Dror Barran like filter
filt_xy = [0, 1, 1, 1, 0;
           1, 1, 2, 1, 1;
           1, 2, 3, 2, 1;
           1, 1, 2, 1, 1;
           0, 1, 1, 1, 0];
% filt_t = [1, 2, 3, 2, 1];
filt_t = [1, 1, 1, 1, 1];
% filt_t = [1, 1, 1];
% filt_t = 1;
filt = repmat(filt_xy,[1,1,length(filt_t)]);

for y_ind = 1:size(filt,2)
    for x_ind = 1:size(filt,1)
        tmp = filt(x_ind,y_ind,:);
        tmp = tmp(:);
        filt(x_ind,y_ind,:) = tmp.*filt_t';
    end
end
filt = filt/sum(filt(:));

% % ones filter
% N = 3;
% filt_size = [N,N,N];
% filt = ones(filt_size)/prod(filt_size);

% Estimate the local mean of f.
% localMean = convn(g,filt,'same');
localMean = 0;
% Estimate of the local variance of f.
localVar = abs(convn(g.^2,filt,'same')) - abs(localMean).^2;
% localVar = (convn(g.^2,filt,'same')) - (localMean).^2;
% Estimate the noise power if necessary.
if (isempty(noise))
  noise = mean2(localVar);
end

% Compute result
% f = localMean + (max(0, localVar - noise) ./ ...
%           max(localVar, noise)) .* (g - localMean);
%
% Computation is split up to minimize use of memory
% for temp arrays.
f = g - localMean;
g = localVar - noise; 
g = max(g, 0);
localVar = max(localVar, noise);
f = f ./ localVar;
f = f .* g;
f = f + localMean;

end

