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
Main Script for Static phantom study (Study I): Reconstruction of undersampled noisy
 Shepp-Logan Phantom using 
1. Compressed Sensing 'CS' (Lustig et al. 2008)
2. Robust Regression 'RR' (Nikolova et al. 2004).
3. Sparse Outliers 'SO' (Dong et al. 2012) .
4.Compressive recovery with Outlier Rejection 'CORe' (Arshad et al. 2023) 
All algorithms are based on ADMM/Split Bregman implementation
%}
%=========================================================================%
%% Workspace Initialization

clear;
clc;
close all;

%% Load subfolders
addpath(genpath('./functions/'));
addpath(genpath('./methods/'));
%% Parameters
n_real=50; %no. of realizations
p.N = 1*[128, 128]; % image size
acs = 12; % size of acs region (width of fully sampled k-space at center)
R = 2.4; % Acceleration rate
ymax = 10; % normalization factor for k-space
p.vrb = 100; % print iterations details after every 'vrb' iterations 

% CS parameters
p.mu_cs  = 2e-1; % langrange multiplier hyperparameter
p.lam_cs = 2e-2*[1e-2, ones(1,3)]; % 4 regularization parameters for regularization
p.oIter_cs = 500; % outer iterations
p.iIter_cs = 4; % inner iterations
p.gStp_cs = 0.8; % gradient step size

% RR parameters
p.mu1_rr  = 1e-1; % langrange multiplier hyperparameter
p.mu2_rr  = 1; % langrange multiplier hyperparameter
p.lam_rr = 5e-1*[1e-2, ones(1,3)]; % 4 regularization parameters for regularization
p.oIter_rr = 500; % outer iterations
p.iIter_rr = 4;% inner iterations
p.gStp_rr = 0.8; % gradient step size

% SO parameters
p.mu_so  = 4e-2; % langrange multiplier hyperparameter
p.lam1_so = 2e-2*[1e-2, ones(1,3)]; % 4 regularization parameters for regularization 
p.lam2_so = 2e-2; %langrange multiplier hyperparameter
p.oIter_so = 500; % outer iterations
p.iIter_so = 4;% inner iterations
p.gStp_so = 0.8; % gradient step size

% CORe parameters
p.mu1_core  = 4e-2; % langrange multiplier hyperparameter
p.mu2_core  = 2e-1; % langrange multiplier hyperparameter
p.lam1_core = 2e-2*[1e-2, ones(1,3)]; % 4 regularization parameters for regularization
p.lam2_core = 2e-1; %langrange multiplier hyperparameter
p.oIter_core = 500; % outer iterations
p.iIter_core = 4;% inner iterations
p.gStp1_core = 0.8; % gradient step 1 size 
p.gStp2_core = p.gStp1_core; % gradient step 2 size

%% Initializations

seed = 24; % seed for random number generation
rng(seed);
fo_arr=2*randi(10,n_real,1); % Percentage of outliers in each realization
no_arr=zeros(n_real,1); % Array initialized for number of outlier readouts in each realizatin
sig2_arr=randi(10,n_real,1); % scale of std of outlier noise in each realization
ind_arr=strings(n_real,1); % Indices of outlier readouts for each realization

% Initializing zero arrays for NMSE and SSIM values for all realizations
nmse_cs=zeros(n_real,1);
nmse_rr=zeros(n_real,1);
nmse_so=zeros(n_real,1);
nmse_core=zeros(n_real,1);

ssim_cs=zeros(n_real,1);
ssim_rr=zeros(n_real,1);
ssim_so=zeros(n_real,1);
ssim_core=zeros(n_real,1);


%% Generate complex phantom
xo = phantom('modified shepp-logan', p.N(1));
% You can replace phantom with the image of your own choice, do change image size p.N accordingly
%% NDDWT Operator
% Defining 2D NDDWT transform operator
p.W = harr_nddwt_2D({'db1','db1'},size(xo),'pres_l2_norm','false');
p.M = [p.N, 4]; % Size of wavelet bands

%% Realization Loop

for iter=1:n_real
    disp("Realization: "+iter);
%% User defined parameters
fo  = fo_arr(iter)/100; % fraction of outlier readouts
sig1 = 1e-3; % std of Gaussian noise
sig2 = sig2_arr(iter)*sig1; % std of outlier noise


%% Sampling pattern
p.s = zeros(size(xo)); % initializing sampling pattern
p.s(:,randsample(p.N(2),round(p.N(2)/R))) = 1; % downsampling the data according to acceleration
p.s(:,end/2+1-acs/2:end/2+1+acs/2) = 1; % setting acs region to 1

%% Define operators
p.A  = @(xo) funA (xo, p.s); % Defining forward operator
p.At  = @(xo) funAt (xo, p.s); % Defining adjoint operator


%% Simulating data
x = xo + 2*(rand-0.5)*1j*xo; % making image complex
y = p.A(x); % downsampled vectorized data
ymax = max(abs(y(:)));
y = y/ymax; % normalizing y
x = x/ymax; % normalizing x with same factor
noise=(sig1*(randn(size(y)) + 1j*randn(size(y)))).*p.s; % Gaussian noise to be added
yn = y + noise; % adding Gaussian noise to data
yn = yn/sig1; % making noise variance 1
x  = x/sig1; %scaling image
yscl = 10/(max(abs(yn(:)))); % scaling factor 
yn = yscl*yn; %scaling measured data
x = yscl*x; %scaling image
xmax = max(abs(x(:)));
ymax = max(abs(yn(:)));



% Add outliers
no_arr(iter) = round(fo*sum(p.s(1,:)~=0)); % number of outliers
ind = randsample(find(p.s(1,:)),no_arr(iter)); % indices of outliers
ind_arr(iter)=regexprep(num2str(ind),'\s+',','); % storing the array of indices of outliers
outliers=zeros(size(yn)); % intialzing outliers
outliers(:,ind)=sig2*ymax*(randn(p.N(2),no_arr(iter)) + 1j*randn(p.N(2),no_arr(iter))); % creating outlier noise
yn = yn + outliers; % adding outliers to data

%% Zero-filled recon
x_zf = p.At(yn); % Defining zero-filled image

%% cs
disp('======CS Recon Started======')
x_cs = cs(yn,p); % cs recon method
nmse_cs(iter)=20*log10( norm(x(:)-x_cs(:))/norm(x(:))); % calculating nmse
ssim_cs(iter)=ssim(real(x_cs),real(x)); % calculating ssim
disp(nmse_cs);
disp(ssim_cs);
%% rr
disp('======RR Recon Started======')
x_rr = rr(yn,p); % rr recon method
nmse_rr(iter)=20*log10( norm(x(:)-x_rr(:))/norm(x(:)) ); % calculating nmse
ssim_rr(iter)=ssim(real(x_rr),real(x)); % calculating ssim
disp(nmse_rr);
disp(ssim_rr);
%% so
disp('======SO Recon Started======')
[x_so, yv] = so(yn,p); % so recon method
nmse_so(iter)=20*log10( norm(x(:)-x_so(:))/norm(x(:)) ); % calculating nmse
ssim_so(iter)=ssim(real(x_so),real(x)); % calculating ssim
disp(nmse_so);
disp(ssim_so);
%% core
disp('======CORe Recon Started======')
[x_core, yg] = core(yn,p); % core recon method
nmse_core(iter)= 20*log10( norm(x(:)-x_core(:))/norm(x(:)) ); % calculating nmse
ssim_core(iter)=ssim(real(x_core),real(x)); % calculating ssim
disp(nmse_core);
disp(ssim_core);
end

%% Average Values
% Calculating mean and std of N realizations
mean_nmse_cs=mean(nmse_cs);
std_nmse_cs=std(nmse_cs);
mean_nmse_rr=mean(nmse_rr);
std_nmse_rr=std(nmse_rr);
mean_nmse_so=mean(nmse_so);
std_nmse_so=std(nmse_so);
mean_nmse_core=mean(nmse_core);
std_nmse_core=std(nmse_core);

mean_ssim_cs=mean(ssim_cs);
std_ssim_cs=std(ssim_cs);
mean_ssim_rr=mean(ssim_rr);
std_ssim_rr=std(ssim_rr);
mean_ssim_so=mean(ssim_so);
std_ssim_so=std(ssim_so);
mean_ssim_core=mean(ssim_core);
std_ssim_core=std(ssim_core);

%% Plotting
k=6; % Number of tiles in figure
fig=figure; 
subplot(2,k,1); imagesc(abs(x),[0,xmax]); axis('image','off'); 
title("Reference Image");

subplot(2,k,2); imagesc(abs(x_zf),[0,xmax]); axis('image','off'); title('zf'); colormap('gray'); 
title(['ZF, nmse: ' num2str( 10*log10( norm(x(:)-x_zf(:))/norm(x(:)) ),3) ])

subplot(2,k,3); imagesc(abs(x_cs),[0,xmax]); axis('image','off'); title('CS'); colormap('gray');
title(['CS, avg nmse: ' num2str((mean_nmse_cs),3) '+- std:' num2str((std_nmse_cs),3) ])

subplot(2,k,4); imagesc(abs(x_rr),[0,xmax]); axis('image','off'); title('RR'); colormap('gray');
title(['RR, avg nmse: ' num2str((mean_nmse_rr),3) '+- std:' num2str((std_nmse_rr),3) ])

subplot(2,k,5); imagesc(abs(x_so),[0,xmax]); axis('image','off'); title('SO'); colormap('gray');
title(['SO, avg nmse: ' num2str((mean_nmse_so),3) '+- std:' num2str((std_nmse_so),3) ])

subplot(2,k,6); imagesc(abs(x_core),[0,xmax]); axis('image','off'); title('CORe'); colormap('gray');
title(['CORe, avg nmse: ' num2str((mean_nmse_core),3) '+- std:' num2str((std_nmse_core),3) ])

subplot(2,k,8); imagesc(5*abs(x_zf - x),[0,xmax]); axis('image','off'); colormap('gray');
title(['Error Map, nmse: ' num2str( ssim(real(x_zf),real(x)),3) ])

subplot(2,k,9); imagesc(5*abs(x_cs - x),[0,xmax]); axis('image','off'); colormap('gray');
title(['Error Map, avg ssim: ' num2str((mean_ssim_cs),3) '+- std:' num2str((std_ssim_cs),3) ])

subplot(2,k,10); imagesc(5*abs(x_rr - x),[0,xmax]); axis('image','off'); colormap('gray');
title(['Error Map, avg ssim: ' num2str((mean_ssim_rr),3) '+- std:' num2str((std_ssim_rr),3) ])

subplot(2,k,11); imagesc(5*abs(x_so - x),[0,xmax]); axis('image','off'); colormap('gray');
title(['Error Map, avg ssim: ' num2str((mean_ssim_so),3) '+- std:' num2str((std_ssim_so),3) ])

subplot(2,k,12); imagesc(5*abs(x_core - x),[0,xmax]); axis('image','off'); colormap('gray');
title(['Error Map, avg ssim: ' num2str((mean_ssim_core),3) '+- std:' num2str((std_ssim_core),3) ])

savefig(fig,fullfile(['study1_real',num2str(n_real),'_',datestr8601,'.fig']));


%% Saving data
save(['study1_real',num2str(n_real),'_',datestr8601]);


