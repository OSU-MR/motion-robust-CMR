function [xhat,maps] = pMRIL14D(kdata,samp,wvar,lambda,maps)

uniform_var = 0;
precision = 'single';
compute = 'mat';

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,3,5,4]));

if nargin<5
    avg_image = sum(kdata,5);
    avg_pattern = sum(samp,4);
    avg_pattern(avg_pattern==0) = inf;
    avg_image = bsxfun(@rdivide,avg_image,avg_pattern);
    avg_image = ifft3_shift(avg_image);
    
    % Estimate sensitivity maps
    p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
    p.reEst = 0; % Res-estimating sensitivities
    p.fil = 3;
    [maps,x0] = WalshCoilCombine3D(avg_image,p);
    x0 = repmat(x0,[1,1,1,size(samp,4)]);
else
    avg_image = sum(kdata,5);
    avg_pattern = sum(samp,4);
    avg_pattern(avg_pattern==0) = inf;
    avg_image = bsxfun(@rdivide,avg_image,avg_pattern);
    avg_image = ifft3_shift(avg_image);
    
    % Estimate sensitivity maps
    p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
    p.reEst = 0; % Res-estimating sensitivities
    p.fil = 3;
    [maps_notused,x0] = WalshCoilCombine3D(avg_image,p);
    x0 = repmat(x0,[1,1,1,size(samp,4)]);
end

%% fftshift and downsample data
kdata =fftshift(fftshift(fftshift(kdata,1),2),3);
samp = fftshift(fftshift(fftshift(samp,1),2),3);
x0 = fftshift(fftshift(fftshift(x0,1),2),3);
maps = fftshift(fftshift(fftshift(maps,1),2),3);
kdata = downsample_data(kdata,samp);

% % Normalize The columns of A to be unit norm
% R = numel(samp)/length(find(samp ==1));
% maps = maps*sqrt(R);
% kdata = kdata*sqrt(R);
% wvar = wvar*R;

use_gpu = 1;
if use_gpu
   kdata = gpuArray(single(kdata));
   samp = gpuArray((samp));
   x0 = gpuArray(single(x0));
   maps = gpuArray(single(maps));
   compute = 'gpu';
   precision = 'single';
end

%% Create Operators =======================================================
pMRI = pMRI_Op_3D_t(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute);

options.lip = [];
options.nit = 40;
options.lam = lambda;
options.tol = 1e-5;
options.verbose = 1;

% Run GAMP
% xhat = gampEst(inputEst,EstimOut,Op_b,GAMPopt);
t = tic;
xhat = fistaEst_pmri4D(x0(:),kdata,pMRI,options);
display(sprintf('Elapsed Time = %0.2f minutes',toc(t)/60))
xhat = gather(xhat);
xhat = reshape(xhat,size(samp));
xhat = ifftshift(ifftshift(ifftshift(xhat,1),2),3);

end
