function [xhat,maps,GAMPout,param] = pMRIEMGM(kdata,samp,wvar,params_in,fudge_in,maps_in)

uniform_var = 0;
precision = 'single';
compute = 'mat';
params = [];
fudge = 1;

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,4,3]));

%% Estimate sensitivity maps
weights_b = repmat(sum(samp,3),[1,1,size(kdata,3)]);
weights_b(find(weights_b==0)) = Inf;

% Normalize k_data by weights and take an inverse Fourier Transform
time_av_b = ifft2_shift(sum(kdata,4)./weights_b);
[x0, maps] = WalshCoilCombine(time_av_b,3);
maps = repmat(maps,[1,1,1,size(samp,3)]);
x0 = repmat(x0,[1,1,size(samp,3)]);

if nargin <4
    
elseif nargin == 4
    params = params_in;
elseif nargin == 5
    params = params_in;
    fudge = fudge_in;
else
    params = params_in;
    fudge = fudge_in;
    maps = maps_in;
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
   precision = 'single';
end

%% Create Operators =======================================================
pMRI = pMRI_Op_2D_t(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute);

% Create nd-DWT Linear Transform
w_prop.wname = {'db1','db1','db1'};
w_prop.level = 1;
W = ndDWTLinTrans(w_prop,size(samp),'uniform_var', uniform_var,'compute',...
    compute,'precision',precision);

% Concatonate All The Linear Transform Operator Together
Op_b = LinTransConcat({pMRI;W},[1,1],precision,compute);

% Create Estimation classes
inputEst = NullEstimIn(0,1);
MeasEstimOut = CAwgnEstimOut(kdata,wvar,0);

%% initialize with FISTA
options.lip = [];
options.nit = 40;
options.lam = .1;
options.tol = 1e-5;
options.verbose = 1;

% initialize with FISTA
x0_FISTA = fistaEst_pmri(x0(:),kdata,pMRI,options);
x0_FISTA = reshape(x0_FISTA,[size(x0)]);
% x0 = fftshift(fftshift(x0_FISTA,1),2);
x0 = x0_FISTA;
GAMPopt.xhat0 = x0(:);

% Set GM estimate
if isempty(params)
    [ lambda, ~, ~, sigma1, sigma2,gm]  = estBGParams(x0);
    sigma1 = sigma1*8;
    simga2 = sigma2*8;
else
    lambda = zeros([size(samp),8]);
    for ind = 1:8
       lambda(:,:,:,ind) = params.lambda(ind);
    end
    lambda = lambda(:);

    sigma1 = zeros([size(samp),8]);;
    for ind = 1:8
       sigma1(:,:,:,ind) = params.sigma1(ind);
    end
    sigma1 = sigma1(:);


    sigma2 = zeros([size(samp),8]);
    for ind = 1:8
       sigma2(:,:,:,ind) = params.sigma2(ind);
    end
    sigma2 = sigma2(:);
end

if use_gpu
    lambda = gpuArray(lambda);
    sigma1 = gpuArray(sigma1);
    sigma2 = gpuArray(sigma2);
end

% warning('scaling sigma esitmates')
% sigma1 = 8*sigma1;
% sigma2 = 8*sigma2;
AnaEstimOut  = CGaussMixEstimOut(zeros(8*numel(x0),1),sigma1,sigma2,lambda);
EstimOut = EstimOutConcat({MeasEstimOut;AnaEstimOut},[pMRI.M,W.M],precision,compute);

sigma1_v = reshape(sigma1,[size(samp),8]);
sigma2_v = reshape(sigma2,[size(samp),8]);
lambda_v = reshape(lambda,[size(samp),8]);

sigma1_v = squeeze(sigma1_v(1,1,1,1:8));
sigma2_v = squeeze(sigma2_v(1,1,1,1:8));
lambda_v = squeeze(lambda_v(1,1,1,1:8));

fprintf(sprintf('sigma1\t\t sigma2\t\t lambda\n'))
for ii = 1:8
    fprintf(sprintf('%f\t%f\t%f\n',sigma1_v(ii),sigma2_v(ii),lambda_v(ii)));
end

% Gamp Options
GAMPopt = GampOpt();
GAMPopt.nit = 50;
GAMPopt.stepWindow = 1;
GAMPopt.step = 0.1;
GAMPopt.verbose = 1;
GAMPopt.legacyOut = 0;
GAMPopt.xhat0 = x0(:);

% fudge = sqrt(8); % fudge factor
for ind = 1%:10
    fprintf(sprintf('Iteration = %d\n',ind));
    % Run GAMP
    xhat = gampEst(inputEst,EstimOut,Op_b,GAMPopt);
    
    % Update Prior on wavelet coefficients
%     [ ~,~,lambda] = CGMUpdate_local(fudge*xhat.zhat,fudge^2*xhat.zvar,size(samp),sigma1,sigma2,lambda);
%     warning('not updating parameters');
    [ sigma1,sigma2,lambda] = CGMUpdate(fudge*xhat.zhat,fudge^2*xhat.zvar,size(samp),sigma1,sigma2,lambda);
    if use_gpu
        lambda = gpuArray(lambda);
        sigma1 = gpuArray(sigma1);
        sigma2 = gpuArray(sigma2);
    end

    AnaEstimOut = CGaussMixEstimOut(zeros(8*numel(samp),1),sigma1,sigma2,lambda);
%     lambda = lambda_local;
    EstimOut = EstimOutConcat({MeasEstimOut;AnaEstimOut},[pMRI.M,W.M],precision,compute);
    
    lambda(lambda<1e-4) = 1e-4;
    sigma1(sigma1<1e-5) = 1e-5;
    sigma2(sigma2<1e-5) = 1e-5;

    sigma1(sigma1>1e3) = 1e3;
    sigma2(sigma2>1e3) = 1e3;

    sigma1_v = reshape(sigma1,[size(samp),8]);
    sigma2_v = reshape(sigma2,[size(samp),8]);
    lambda_v = reshape(lambda,[size(samp),8]);

    sigma1_v = squeeze(sigma1_v(1,1,1,1:8));
    sigma2_v = squeeze(sigma2_v(1,1,1,1:8));
    lambda_v = squeeze(lambda_v(1,1,1,1:8));

    fprintf(sprintf('sigma1\t\t sigma2\t\t lambda\n'))
    for ii = 1:8
        fprintf(sprintf('%f\t%f\t%f\n',sigma1_v(ii),sigma2_v(ii),lambda_v(ii)));
    end

    % Warm Start xb		
    GAMPopt.xhat0 = xhat.xhat;		
    GAMPopt.xvar0 = xhat.xvar;		
    GAMPopt.shat0 = xhat.shat;		
    GAMPopt.svar0 = xhat.svar;		
    GAMPopt.xhatPrev0 = xhat.xhatPrev;		
    GAMPopt.scaleFac = xhat.scaleFac;		
    GAMPopt.step = min(max(xhat.step,0.05),xhat.stepMax);
    
    GAMPopt.nit = 10;
end
GAMPout = xhat;
xhat = reshape(xhat.xhat,size(samp));
xhat = fftshift(fftshift(xhat,1),2);
xhat = gather(xhat);
param.sigma1 = sigma1;
param.sigma2 = sigma2;
param.lambda = lambda;
   
end
