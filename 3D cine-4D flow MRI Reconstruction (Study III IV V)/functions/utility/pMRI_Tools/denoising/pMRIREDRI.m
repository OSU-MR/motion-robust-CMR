function [xhat,maps,zhat] = pMRIREDRI(kdata,samp,wvar,lamRED,lamW,denoise,maps)

uniform_var = 0;
precision = 'double';
compute = 'mat';

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,4,3]));

if nargin<7
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
pMRI = pMRI_Op_2D_t_RI2(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute);

switch lower(denoise)
    case ('nd_wavelets')
        % Create nd-DWT Linear Transform
        w_prop.wname = {'db1','db1','db1'};
        w_prop.level = 1;
        W = ndDWTLinTrans(w_prop,size(samp),'uniform_var', uniform_var,'compute',...
                          compute,'precision',precision);
        denoiser = @(x) shrinkWave3DRI2(x,size(samp),W,lamW);
    case ('nd_wavelets_snipe')
        % Create nd-DWT Linear Transform
        w_prop.wname = {'db1','db1','db1'};
        w_prop.level = 1;
        W = ndDWTLinTrans(w_prop,size(samp),'uniform_var', uniform_var,'compute',...
                          compute,'precision',precision);
        denoiser = @(x) shrinkWaveSNIPE3D(x,size(samp),W,lamW);
    case ('nd_wavelets_wiener')
        % Create nd-DWT Linear Transform
        w_prop.wname = {'db1','db1','db1'};
        w_prop.level = 1;
        W = ndDWTLinTrans(w_prop,size(samp),'uniform_var', 1,'compute',...
                          compute,'precision',precision);
        denoiser = @(x) shrinkWaveWiener3D(x,size(samp),W,lamW);
%         denoiser = @(x,lamW) shrinkWaveWiener3D(x,size(samp),W,lamW);
    case ('median')
        denoiser = @(x) shrinkMEDIAN3D(x,size(samp));
    case ('wiener')
        denoiser = @(x) shrinkWiener3D(x,size(samp),lamW);
    case ('dwt')
        denoiser = @(x) shrinkDWT3D(x,size(samp),lamW);
    case ('vbm4d')
        denoiser = @(x) shrinkVBM4D(x,size(samp),lamW);
    case ('dwt_wiener')
        denoiser = @(x) shrinkDWT3DWiener(x,size(samp),lamW);
end

%% Reconstruct

% initialize with FISTA
% options.lip = [];
% options.nit = 40;
% options.lam = 0.001;
% options.tol = 1e-5;
% options.verbose = 1;
% x0 = fistaEst_pmri(x0(:),kdata,pMRI,options);
% x0 = x0(:);

% [x0] = cgSENSE(x0(:), kdata, pMRI);
% x0 = x0(:);

% [x0] = cgSENSEL2(x0(:), kdata, pMRI,0.005);
% x0 = x0(:);

x0 = [x0(:);zeros(size(x0(:)))];

[xhat,zhat, ~] = admmREDRI(x0(:), kdata, pMRI, lamRED, denoiser);
xhat = xhat(1:length(xhat)/2) + 1j*xhat(length(xhat)/2+1:end);
xhat = reshape(xhat,size(samp));
xhat = fftshift(fftshift(xhat,1),2);
xhat = gather(xhat);

zhat = zhat(1:length(zhat)/2) + 1j*zhat(length(zhat)/2+1:end);
zhat = reshape(zhat,size(samp));
zhat = fftshift(fftshift(zhat,1),2);
zhat = gather(zhat);


end
