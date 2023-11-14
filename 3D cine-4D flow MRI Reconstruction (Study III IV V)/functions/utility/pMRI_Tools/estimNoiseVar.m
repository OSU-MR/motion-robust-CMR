function [ wvar ] = estimNoiseVar( noise_ch,sampB)
%ESTIMNOISEVAR Summary of this function goes here
%   Detailed explanation goes here
    
% Find the sampling pattern direction
tmp = squeeze(sampB(:,1,1));
if mean(tmp) == 0 || mean(tmp) ==1
    dim = 2; % Case when columns represent the phase encode direction of the data
else
    dim = 1; % Case when rows represent the phase encode direction of the data
end
clear tmp

% Estimate Noise from the outermost lines of k-space
noise = [];
if dim == 1;
    samp_col = sampB(:,1,:);
    for ind = 1:size(sampB,3)
        tmp = find(samp_col(:,:,ind)==1);
        min_start = min(tmp);
        min_end = min(size(sampB,2)-tmp);
        if min_start<=min_end
            outer_row(ind) = min_start;
        elseif min_start>min_end
            outer_row(ind) = max(tmp);
        end
        noise = [noise,noise_ch(outer_row(ind),:,ind)];
    end
elseif dim == 2;
    samp_row = sampB(1,:,:);
    for ind = 1:size(sampB,3)
        tmp = find(samp_row(:,:,ind)==1);
        min_start = min(tmp);
        min_end = min(size(sampB,2)-tmp);
        if min_start<=min_end
            outer_col(ind) = min_start;
        elseif min_start>min_end
            outer_col(ind) = max(tmp);
        end
        noise = [noise;noise_ch(:,outer_col(ind),ind)];
    end
end
wvar = var(noise);

end


