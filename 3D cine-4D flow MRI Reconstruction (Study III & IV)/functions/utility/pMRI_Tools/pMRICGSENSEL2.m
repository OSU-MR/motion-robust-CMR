function [xhat,maps] = pMRICGSENSEL2(kdata,samp,wvar,lam,maps)

uniform_var = 0;
precision = 'double';
compute = 'mat';

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,4,3]));

if nargin<5
    %% Estimate sensitivity maps
    weights_b = repmat(sum(samp,3),[1,1,size(kdata,3)]);
    weights_b(find(weights_b==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(kdata,4)./weights_b);
    [x0, maps] = WalshCoilCombine(time_av_b,3);
    maps = repmat(maps,[1,1,1,size(samp,3)]);
    x0 = repmat(x0,[1,1,size(samp,3)]);
else
    %% Estimate sensitivity maps
    weights_b = repmat(sum(samp,3),[1,1,size(kdata,3)]);
    weights_b(find(weights_b==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av_b = ifft2_shift(sum(kdata,4)./weights_b);
    [x0, ~] = WalshCoilCombine(time_av_b,3);
%     maps = repmat(maps,[1,1,1,size(samp,3)]);
    x0 = repmat(x0,[1,1,size(samp,3)]);
end

%% fftshift and downsample data
kdata = fftshift(fftshift(kdata,1),2);
samp = fftshift(fftshift(samp,1),2);
x0 = fftshift(fftshift(x0,1),2);
maps = fftshift(fftshift(maps,1),2);
kdata = downsample_data(kdata,samp);

% Normalize The columns of A to be unit norm
R = numel(samp)/length(find(samp ==1));
maps = maps*sqrt(R);
kdata = kdata*sqrt(R);
wvar = wvar*R;

use_gpu = 1;
if use_gpu
   kdata = gpuArray(single(kdata));
   samp = gpuArray(single(samp));
   x0 = gpuArray(single(x0));
   maps = gpuArray(single(maps));
   compute = 'gpu';
   precision = 'double';
end

%% Create Operators =======================================================
pMRI = pMRI_Op_2D_t(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute);


%% Reconstruct
[xhat] = cgSENSEL2(x0(:), kdata, pMRI,lam);
xhat = reshape(xhat,size(samp));
xhat = fftshift(fftshift(xhat,1),2);
xhat = gather(xhat);


end
