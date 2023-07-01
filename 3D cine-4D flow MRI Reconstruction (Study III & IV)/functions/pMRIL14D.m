function [xhat,hist] = pMRIL14D(kdata,samp,weights,opt,maps)

uniform_var = 0;
precision = 'single';
compute = 'mat';

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,3,5,4]));

if nargin<5 %if sensitivity maps are not provided in function input
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
    x0 = repmat(x0,[1,1,1,size(samp,4)]); %x0 represents the initial image
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
    x0 = repmat(x0,[1,1,1,size(samp,4)]); %x0 represents the initial image
end

%% fftshift and downsample data
% Using fftshift to bring the center in the middle
kdata =fftshift(fftshift(fftshift(kdata,1),2),3);
samp = fftshift(fftshift(fftshift(samp,1),2),3);
weights = fftshift(fftshift(fftshift(weights,1),2),3);
% apply respiratory weights
kdata = bsxfun(@times,kdata,permute(weights,[1,2,3,5,4]));
x0 = fftshift(fftshift(fftshift(x0,1),2),3);
maps = fftshift(fftshift(fftshift(maps,1),2),3);
kdata = downsample_data(kdata,samp);

% if using GPU, transfering arrays to GPU
if opt.use_gpu
   kdata = gpuArray(single(kdata));
   samp = gpuArray((samp));
   weights = gpuArray(single(weights));
   x0 = gpuArray(single(x0));
   maps = gpuArray(single(maps));
   compute = 'gpu'; 
   precision = 'single';
end

%% Create Operators =======================================================
pMRI = pMRI_Op_3D_t(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute,'weights',weights);

%% Extract parameters

options.use_gpu = opt.use_gpu;
options.transform = opt.transform;
options.sfolder=opt.sfolder;
options.oIter=opt.oIter;
options.iIter=opt.iIter;
options.gStp =opt.gStp;
options.vrb = opt.vrb;
options.readout = opt.readout;
options.lam_l2l1=opt.lam_cs;
options.mu_l2l1=opt.mu_cs;
options.lam1_l2l1g = opt.lam1_core;
options.lam2_l2l1g =opt.lam2_core;
options.mu1_l2l1g =  opt.mu1_core;
options.mu2_l2l1g =  opt.mu2_core;


t = tic; % Starting reconstruction timer
if strcmp(opt.recon,'cs')
disp('CS Estimation Started')
% xhat is reconstructed complex time varying 3D image
% hist is the history of objective
[xhat,hist] = cs(x0(:),kdata,pMRI,options);

elseif strcmp(opt.recon,'core')
disp('CORe Estimation Started')

% xhat is reconstructed complex time varying 3D image
% hist is the values of objective function along iterations 
[xhat,hist] = core(x0(:),kdata,pMRI,options);

else
    error('Wrong recon method.......')
end

fprintf('Elapsed Time = %0.2f minutes\n',toc(t)/60);
xhat = gather(xhat);
xhat = reshape(xhat,[cat(2,size(samp),1)]);
xhat = ifftshift(ifftshift(ifftshift(xhat,1),2),3);

end
