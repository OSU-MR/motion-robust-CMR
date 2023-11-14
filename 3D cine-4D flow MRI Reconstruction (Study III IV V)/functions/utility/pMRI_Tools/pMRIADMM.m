function [xhat,maps,zhat] = pMRIADMM(kdata,samp,wvar,lambda,maps)

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

% Create nd-DWT Linear Transform
w_prop.wname = {'db1','db1','db1'};
w_prop.level = 1;
W = ndDWTLinTrans(w_prop,size(samp),'uniform_var', uniform_var,'compute',...
    compute,'precision',precision);

%% Create a Vector of lambda values
% lambda = setLambda(size(samp),lambda);
k = 10;
lambda_band = zeros(8,1);
lambda_band(1) = 0.001*lambda;
lambda_band(2) = lambda;
lambda_band(3) = lambda;
lambda_band(4) = lambda;
lambda_band(5) = k*lambda;
lambda_band(6) = k*lambda;
lambda_band(7) = k*lambda;
lambda_band(8) = k*lambda;

lambda_vec = zeros([size(samp),8]);
for ind = 1:8
   lambda_vec(:,:,:,ind) = lambda_band(ind);
end
lambda_vec = lambda_vec(:);

% initialize with FISTA
[xhat,zhat] = admmpMRI(x0(:),kdata,pMRI,lambda_vec,1,1);
xhat = reshape(xhat,size(samp));
xhat = fftshift(fftshift(xhat,1),2);
xhat = gather(xhat);

zhat = reshape(zhat,size(samp));
zhat = fftshift(fftshift(zhat,1),2);
zhat = gather(zhat);


end
