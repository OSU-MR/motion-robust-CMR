%{
Reference: 
The algorithms implemented below corresponds to research
currently submitted for MRM Journal publication and available as a 
preprint on arXiv.

Preprint Details:
Title: "Motion-robust free-running cardiovascular MRI"
arXiv: " arXiv:2308.02088v2 [eess.IV]" (Link: https:arxiv.org/abs/2308.02088v2)



The preprint on arXiv provides an accessible overview of the research 
and may be cited for more detailed information pertaining to the 
algorithmic approach used herein.

How to cite:
Arshad SM, Potter LC, Chen C, Liu Y, Chandrasekaran P, Crabtree C, Han Y,
Ahmad R (2023). Motion-robust free-running cardiovascular MRI. 
arXiv preprint arXiv:2308.02088. 

% ========================================================================%
CMR LAB (https://u.osu.edu/ahmad/)
CMR LAB Github: https://github.com/orgs/OSU-MR

The Ohio State University
Written by:
Syed Murtaza Arshad (arshad.32@osu.edu)
(https://github.com/syedmurtazaarshad)

Rizwan Ahmad, PhD (ahmad.46@osu.edu)
%}

%%
%=========================================================================%
%{
About the code:
Main Script for Dynamic phantom study (Study II): using 
1. Compressed Sensing 'CS' (Lustig et al. 2008)
2. Robust Regression 'RR' (Nikolova et al. 2004).
3. Sparse Outliers 'SO' (Dong et al. 2012) .
4.Compressive recovery with Outlier Rejection 'CORe' (Arshad et al. 2023) 
All algorithms are based on ADMM/Split Bregman implementation
%}
%=========================================================================%
%% Workspace Initialization
clear;
% clc;
close all;

%% Load subfolders
restoredefaultpath
addpath(genpath('./functions/'));
addpath(genpath('./methods/'));
%% Parameters

n=5; % no. of realizations

%% Sampling Patten
% GRO Parameters
param_gro.PE   = 256;  % Size of of phase encoding (PE) grid
param_gro.FR   = 26;   % Numbe of frames
param_gro.n    = 10;   % Number of samples (readouts) per frame
param_gro.M    = param_gro.FR*param_gro.n; % Total number of samples
param_gro.E    = 1;    % Number of encoding, E=1 for cine, E=2 for flow (phase-contrast MRI)
param_gro.tau  = 1;    % Extent of shift between frames, tau = 1 or 2: golden ratio shift, tau>2: tiny golden ratio shift
param_gro.s    = 10;  % s>=1. larger values means higher sampling density in the middle (default: 2.2)
param_gro.alph = 5;    % alph>1. larger alpha means sharper transition from high-density to low-density regions (default: 3)
param_gro.PF   = 0;    % for partial fourier; discards PF samples from one side (default: 0)
param_gro.dsp  = 0;    % Display figures: 0 no, 1 yes
% PEInd: indexes of randomly sampled readout in each frame
% sampling: gives the sampling pattern
 [PEInd, sampling] = gro_fun(param_gro);   

%% Phantom parameters
param.nx  = 256;  % Pixels in x-dimension
param.ny  = 256;  % pixels in y-dimension
param.nt  = param_gro.M ;   % Nuber of frames
param.cy  = 5;    % Number of temporal cycles
param.th  = 1*pi/8; % A spatially invariant phase of the phantom to make phantom complex valued
% Multivariate Gaussian used for k-space filtering to avoid aliasing of the phantom
param.sig = [param.nx/4, param.ny/4]; % Std (in k-space) of the filter
param.row = 0;                        % Correlation coefficient of 
param.ctr = [0,0];                    % Center of multivariate Gaussian
param.nor = 3;                        % Type of normalization of Gaussian
param.disp= 0;


p.sig   = 5e-4; % noise std
p.out   = 0.00; % percentage of outliers
p.oIter = 250;  % outer iterations
p.iIter = 4;   % inner iterations
p.gStp  = 0.5; % gradient step size
p.acs   = 12; % size of the fully sampled ACS region
p.R     = 2;  % accleration rate;
p.sd    = 1; % seed for random number generation
p.vrb   = 50; % verbosite level; every vrb-th iteration spits out the objective function
p.N = [param.nx,param.ny]; % image size

%% for CS
p.lam_cs = 2.*1.2e-1*[1e-1, 1, 1, sqrt(2)];% regularization parameter for four wavelet bands
p.mu_cs =   2.*8.8e-1; % lagrange multiplier hyperparameter

%% for RR
p.lam_rr   = 6.3e-1*[1e-2, 1, 1, sqrt(2)];  % regularization parameter for four wavelet bands
p.mu1_rr   = 6e-2; % lagrange multiplier hyperparameter
p.mu2_rr   = 9.1e-1;  % lagrange multiplier hyperparameter

%% for SO
p.lam1_so = 2.*2.3e-2*[1e-2, 1, 1, sqrt(2)];% regularization parameter for four wavelet bands
p.lam2_so = 2.7e-2; % lagrange multiplier hyperparameter
p.mu_so =   2.*8.5e-2; % lagrange multiplier hyperparameter

%% for CORe

p.mu1_core  = 2.*5e-1; % lagrange multiplier hyperparameter
p.mu2_core  = 5e-1; % lagrange multiplier hyperparameter
p.lam1_core = 2.*2.3e-2*[1e-2, 1, 1, sqrt(2)]; % regularization parameter for four wavelet bands
p.lam2_core = 2.1e-1; % lagrange multiplier hyperparameter

%% NMSE & SSIM Array initializations

nmse_rr_arr=zeros(n,1);
nmse_cs_arr=zeros(n,1);
nmse_so_arr=zeros(n,1);
nmse_core_arr=zeros(n,1);

ssim_rr_arr=zeros(n,1);
ssim_cs_arr=zeros(n,1);
ssim_so_arr=zeros(n,1);
ssim_core_arr=zeros(n,1);

%% Loop Start
for k=1:n

fprintf('\nRealization N=%d\n',k);

%% Define the temporal shape of the Respiratory Motion
seed=43;
rng(seed+k);
ti=rand*2*pi;
tcos = cos(linspace(0+ti,2*pi*param.cy+ti,param.nt)); %defining a cos function with random starting phase
transition=0.3; % Transition defines the time spent in each respiratory state
tfun=zeros(size(tcos));
tfun(tcos>transition)=tcos(tcos>transition).^0.001; % Above transition there is stable expiratory state
tfun(tcos<=transition)=sign(tcos(tcos<=transition)).*abs(tcos(tcos<=transition)).^0.05; % Below transition there is inspiratory phase with slight movement
% plot tfun to observe temporal shape
%% Building phantom in k-space
disp("creating digital phantom K space...")
kdata = exp(1j*param.th)*ephantom(tfun,param); %  exp(1j*param.th) to make it complex
x0 = ifft2r(kdata); % image domain

%% True Image
true_ind=find(tcos>0.999); % finding a reference frame from expiratory phase
x=ifft2r(squeeze(kdata(:,:,true_ind(1))/sqrt(numel(kdata(:,:,true_ind(1)))))); % x is reference image
p.x=x; % storing reference image in parameters

%% Sampling problem
binned=zeros(size(PEInd));
weight=zeros(size(PEInd));
binned(tfun>-0.9)=PEInd(tfun>-0.9); %binned contains indexes of the sampled readouts  
weight(tfun>-0.9)=tfun(tfun>-0.9);  %weight represents the motion state of the sampled readout i.e. ~1:represnets Expiratory State, ~-1:represents Inspiratory State, rest is in between transition  
rep_weights=repmat(weight,1,256)'; %repeating the weights across the readouts as the whole readout has the same weight
binning_pattern=logical(binned); %binned sampling pattern
acceleration=param.ny/numel(unique(binned(binned>0))); %acceleration of the data after binning

kdata_ds=zeros(size(kdata(:,:,1)));
weighting_pattern=zeros(size(kdata(:,:,1)));
% Binning process, filling k-space from the data collected from moving phantom, ideally all the
% binned data should come from expiratory phase, but the problem is designed in a way that ~10% of sampled data comes from inspriatory phase 
for i=1:numel(binned)
    if binned(i)>0
        kdata_ds(:,binned(i))=kdata(:,binned(i),i);
        weighting_pattern(:,binned(i))=rep_weights(binned(i),i);
    end
end

kdata_ds=kdata_ds/sqrt(numel(kdata_ds));%normalizing the sampled k-space data
samp=logical(abs(kdata_ds)); % sampling pattern of sampled k-space data

sampInd = find(samp~=0); % vectorized sampled indexes
p.sampInd=sampInd;

kdata_outliers=kdata_ds(weighting_pattern<0.8 & weighting_pattern~=0); %extracting the outliers from k-space data to analyze the efficacy of outlier rejection

%% Define forward model and it's adjoint
p.A  = @(x) funA (x, sampInd, p.N);
p.At = @(x) funAt(x, sampInd, p.N);

%% Figures
if(k==1) %Displaying representative figures for realization 1

figure;
yyaxis right
plot(binned,'.','MarkerSize', 12);
ylim([0 param.ny])
hold on;

yyaxis left
plot(tfun);
yline(-0.9,'-');
hold off;
title("Realization="+k+"   Acceleration="+acceleration);
ylim([-1 1])
legend({'Motion Cycles','Binning Cutoff for motion','Binned Phase Encodings'},'Location','best')


figure;
subplot(221)
imagesc(abs(kdata(:,:,1)).^.2);
 title("Realization="+k+"    Orginal Kspace");
% 
subplot(222)
imagesc(abs(kdata_ds(:,:,1)).^.2);
title("Realization="+k+"     Downsampled Binned Kspace");
% 
subplot(223)
imagesc(abs(samp(:,:,1)).^.2);
title("Realization="+k+"    Sampling Pattern");
% 
subplot(224)
imagesc(weighting_pattern(:,:,1));
title("Realization="+k+" Weighting Pattern ");
subtitle("1~Expiratory State (Good Data), -1~Inspiratory Sate (Outliers), 0=Not sampled")
colorbar;

end

%% Observed Siganl

y0=kdata_ds(sampInd); %vectorized sampled data
y0=y0./max(abs(y0(:)));
x=x./max(abs(y0(:)));
y=y0 + p.sig*(randn(size(y0)) + 1j*randn(size(y0))); %measured k-space data with added Gaussian noise
y = y/p.sig; % making noise variance 1
x  = x/p.sig; %scaling image
ymax = max(abs(y(:))); 
yscl=75; % scaling factor, you can choose any scaling factor to get objective values in certain ranges
y = yscl*y/ymax; %scaling measured data
x = x/ymax; %scaling image
ymax=max(abs(y(:)));
xmax=max(abs(x(:)));
%% Define NDDWT and it's adjoint
p.W = harr_nddwt_2D({'db1','db1'},size(x),'pres_l2_norm','false');
p.M = [p.N, 4]; 

%% CS
disp('CS Reconstruction Started')
%x_cs: reconstructed image
x_cs=cs(y,p);
nmse_cs_arr(k)=20*log10( norm(x(:)-x_cs(:))/norm(x(:)));
ssim_cs_arr(k)=ssim(real(x_cs),real(x));
%% RR 
disp('RR Reconstruction Started')
%x_rr: reconstructed image
x_rr = rr(y,p);
nmse_rr_arr(k)=20*log10( norm(x(:)-x_rr(:))/norm(x(:)) );
ssim_rr_arr(k)=ssim(real(x_rr),real(x));
%% SO
disp('SO Reconstruction Started')
%x_so: reconstructed image
%v_so: rejected outliers
[x_so,v_so]=so(y,p);
nmse_so_arr(k)=20*log10( norm(x(:)-x_so(:))/norm(x(:)) );
ssim_so_arr(k)=ssim(real(x_so),real(x));

%% CORe
disp('CORe Reconstruction Started')
%x_core: reconstructed image
%v_core: rejected outliers
[x_core,v_core]=core(y,p);
nmse_core_arr(k)= 20*log10( norm(x(:)-x_core(:))/norm(x(:)) );
ssim_core_arr(k)=ssim(real(x_core),real(x));



%% 

formatSpec=['\nRealization: %d' ...
    '\nNMSE:\tcs=%.2f\trr=%.2f\tso=%.2f\tcore=%.2f' ...
    '\nSSIM:\tcs=%.2f\trr=%.2f\tso=%.2f\tcore=%.2f'];
fprintf(formatSpec,k...
    ,nmse_cs_arr(k),nmse_rr_arr(k),nmse_so_arr(k),nmse_core_arr(k)...
    ,ssim_cs_arr(k),ssim_rr_arr(k),ssim_so_arr(k),ssim_core_arr(k));

end

%% Mean SSIM and NMSE a

ssim_rr_mean=mean(ssim_rr_arr);
ssim_cs_mean=mean(ssim_cs_arr);
ssim_so_mean=mean(ssim_so_arr);
ssim_core_mean=mean(ssim_core_arr);

nmse_rr_mean=mean(nmse_rr_arr);
nmse_cs_mean=mean(nmse_cs_arr);
nmse_so_mean=mean(nmse_so_arr);
nmse_core_mean=mean(nmse_core_arr);



formatSpec=['\nOverall Results: %d' ...
    '\nNMSE:\tcs=%.2f\trr=%.2f\tso=%.2f\tcore=%.2f' ...
    '\nSSIM:\tcs=%.2f\trr=%.2f\tso=%.2f\tcore=%.2f'];
fprintf(formatSpec,k,...
    nmse_cs_mean,nmse_rr_mean,nmse_so_mean,nmse_core_mean...
    ,ssim_cs_mean,ssim_rr_mean,ssim_so_mean,ssim_core_mean);


%% Saving data
save(['study2_coread',num2str(n),'_',datestr8601]);


%% Plotting


fig2=figure; subplot(2,3,1); imagesc(real(x)); axis('image'); colormap(gray); title("Realizations# ="+n);subtitle("Target")
        subplot(2,3,2); imagesc(real(p.At(y))); axis('image'); colormap(gray); title('Zero-filled');
         subplot(2,3,2); imagesc(real(p.At(y))); axis('image'); colormap(gray); title('Zero-filled');
        subplot(2,3,3); imagesc(real(x_cs)); axis('image'); colormap(gray); title(['cs mean nmse: ',num2str(nmse_cs_mean,3),'+-',num2str(std(nmse_cs_arr),3)]);
        subplot(2,3,4); imagesc(real(x_so)); axis('image'); colormap(gray); title(['so mean nmse: ',num2str(nmse_so_mean,3),'+-',num2str(std(nmse_so_arr),3)]);
        subplot(2,3,5); imagesc(real(x_rr)); axis('image'); colormap(gray); title(['rr mean nmse: ',num2str(nmse_rr_mean,3),'+-',num2str(std(nmse_rr_arr),3)]);
        subplot(2,3,6); imagesc(real(x_core)); axis('image'); colormap(gray); title(['core mean nmse: ',num2str(nmse_core_mean,3),'+-',num2str(std(nmse_core_arr),3)]);

savefig(fig2,fullfile(['study2nmse_real',num2str(n),'_',datestr8601,'.fig']));


fig3=figure; subplot(2,3,1); imagesc(real(x)); axis('image'); colormap(gray); title("Realizations# ="+n);subtitle("Target")
        subplot(2,3,2); imagesc(real(p.At(y))); axis('image'); colormap(gray); title('Zero-filled');
        subplot(2,3,3); imagesc(real(x_cs)); axis('image'); colormap(gray); title(['cs mean nmse: ',num2str(ssim_cs_mean,3),'+-',num2str(std(ssim_cs_arr),3)]);
        subplot(2,3,4); imagesc(real(x_so)); axis('image'); colormap(gray); title(['so mean ssim: ',num2str(ssim_so_mean,3),'+-',num2str(std(ssim_so_arr),3)]);
        subplot(2,3,5); imagesc(real(x_rr)); axis('image'); colormap(gray); title(['rr mean ssim: ',num2str(ssim_rr_mean,3),'+-',num2str(std(ssim_rr_arr),3)]);
        subplot(2,3,6); imagesc(real(x_core)); axis('image'); colormap(gray); title(['core mean ssim: ',num2str(ssim_core_mean,3),'+-',num2str(std(ssim_core_arr),3)]);

savefig(fig3,fullfile(['study2ssim_real',num2str(n),'_',datestr8601,'.fig']));
