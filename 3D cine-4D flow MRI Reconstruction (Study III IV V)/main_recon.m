%{
Github Repository: https://github.com/OSU-MR/motion-robust-CMR

Reference: 
The algorithms implemented below corresponds to research
article published in Magnetic Resonance in Medicine (MRM) journal 
(https://doi.org/10.1002/mrm.30123).

Publication Details:
Title: "Motion-robust free-running volumetric cardiovascular MRI"
Authors: Syed M. Arshad,  Lee C. Potter,  Chong Chen,  Yingmin Liu,  Preethi Chandrasekaran, 
Christopher Crabtree,  Matthew S. Tong,  Orlando P. Simonetti,  Yuchi Han,  Rizwan Ahmad

The article provides an accessible overview of the research 
and may be cited for more detailed information pertaining to the 
algorithmic approach used herein.

How to cite:
Arshad SM, Potter LC, Chen C, et al. Motion-robust free-running volumetric cardiovascular MRI. 
Magn Reson Med. 2024; 1-15. doi: 10.1002/mrm.30123

% ========================================================================%

CMR LAB (https://u.osu.edu/ahmad/)
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
Main Script for 3D Cine & 4D Flow Reconstruction using 
1. Compressive recovery with Outlier Rejection 'CORe' (Arshad et al. 2024) 
2. Compressed Sensing 'CS' (Lustig et al. 2008)
 Both algorithms are based on ADMM/Split Bregman implementation
%=========================================================================%
%}

%% Workspace Initialization

clear; 
close all
%Adding paths for relevant functions and recon methods
addpath(genpath('./functions/'));
addpath(genpath('./recon_methods/'));

%Select the methods you want to use, you can mention one or both
%Note: Method names are case sensitive, all small alphabets
%"cs", "core"
recon=["cs","core"]; 

prompt='Name the Reconstruction Dataset\n';
saveName=input(prompt,'s');

prompt='Is the data 4D flow (1) or 3D cine (0)?\n';
flow=input(prompt);

prompt='Is the data Rest (1) or Exercise (0)?\n';
rest=input(prompt);

opt.flow=flow;
% Add your own path where the data is
[fname, pathname] = uigetfile('./', ...
    'Select data file');   % data file (.mat)


%*************************************************
% add folder path to save results
sfolder = ['./',saveName,'_',datestr8601,'/'];               

%*************************************************
%===================================================================================================
% Data input should be type struct() with fields:
%   For 4D Flow
%   kb: complex k-space array (background encoding) - [kx X ky X kz x Coil X Time/Phases]
%   kx: ... (vel encoding 1)
%   ky: ... (vel encoding 2)
%   kz: ... (vel encoding 3)
%   For 3D cine
%   kb: complex k-space array (background encoding) - [kx X ky X kz x Coil X Time/Phases]
% If following available (otherwise estimate it from k-space data)
%   sampB,X,Y,Z: sampling pattern (optional) - [S1 X S2 X S3 X Time/Phases]
%   weightsB,X,Y,Z: data weights (optional) - [S1 X S2 X S3 X Time/Phases]
%       - If you don't need weights, set to sampling pattern

%% Regularization Parameters 
% 16 regularization params for 16 wavelet bands (CS)
% lambda= c.* [Wav. band 1, Wav. band 2, ....., Wav. band 16]
% First try adjusting scalar 'c' for optimizing regularization
% Then if necessary try changing the ratios of wavelet bands

if(flow) % For 4D Flow
    if(rest)  % For rest datasets
        %For CS
        opt.lam_cs = 2e-4*[1e-2, 1 1,1,1,1,1,1,5,5,5,5,5,5,5,5]; 
        %For CORe
        opt.lam1_core = 2e-4*[1e-2,1 1,1,1,1,1,1, 5,5,5,5,5,5,5,5];  
        opt.lam2_core = 7.5e-2; 
    else  % For exercise datasets (using 2x stronger regularization)
        %For CS
        opt.lam_cs = 4e-4*[1e-2, 1 1,1,1,1,1,1,5,5,5,5,5,5,5,5]; 
        %For CORe
        opt.lam1_core = 4e-4*[1e-2,1 1,1,1,1,1,1, 5,5,5,5,5,5,5,5];  
        opt.lam2_core = 7.5e-2;
    end

else % For 3D Cine
    if(rest)  % For rest datasets
        %For CS
        opt.lam_cs = 5e-4*[1e-2, 1 1,1,1,1,1,1,5,5,5,5,5,5,5,5]; 
        %For CORe
        opt.lam1_core = 5e-4*[1e-2,1 1,1,1,1,1,1, 5,5,5,5,5,5,5,5];  
        opt.lam2_core = 7.5e-2; 
    else  % For exercise datasets
        %For CS
        opt.lam_cs = 7e-4*[1e-2, 1 1,1,1,1,1,1,5,5,5,5,5,5,5,5]; 
        %For CORe
        opt.lam1_core = 7e-4*[1e-2,1 1,1,1,1,1,1, 5,5,5,5,5,5,5,5];  
        opt.lam2_core = 7.5e-2;
    end
end


opt.mu_cs =   5e-1; %Lagrange multiplier
opt.mu1_core =   5e-1; %Lagrange multiplier
opt.mu2_core =   5e-1; %Lagrange multiplier

%% Options % parameters
%formatting of numbers
formatSpec = '%.3g';
opt.coil = 12;          % Number of reduced coils
opt.use_gpu = 1;        % gpu (1) or cpu (0)
opt.nit = 50;           % number of outer iterations
opt.spar = 'jtv';
opt.transform = 'harr'; %harr for harrwavelet transform and tv for tv transform (NN)
opt.sfolder=sfolder;

opt.oIter = 50;     % outer iterations
opt.iIter = 4;      % inner iterations
opt.gStp  = 1e-1;   % gradient step size
opt.vrb   = 5;      %p rint iterations details after every 'vrb' iterations 




%% Import motion resolved (binned/sorted data arrays)
%--------------------------------------------------------------------------------------------------

data = importdata(fullfile(pathname, fname));

if(opt.flow)
    % If 4D flow k-space
    % following concatentes the 5D (readout,phase,slice,coil,frame) 
    % v encoded data in 6 th dim data in this dimension
    kdata = cat(6, data.kb, data.kx, data.ky, data.kz); 
else
    %If 3D Cine k-space
    kdata=data.kb;
end

if isfield(data,'sampB') && ~isempty(data.sampB)
    % if sampling pattern available
    % %4 dimensional sampling pattern (x,y,z,frame,venc) being concatenated over 5th dim
    samp = cat(5, data.sampB, data.sampX, data.sampY, data.sampZ); 
else
    % Estimating sampling pattern from the data
    % squeeze removes the coil dimension & logical converts them into 0 
    % or 1 hence we get sampling pattern
    samp = logical(squeeze(abs(kdata(:,:,:,1,:,:)))); 
end

% Calculating total number of readouts in the data
opt.readout=numel(find(sum(reshape(samp,size(samp,1),[]),2)));


% weights (if defined)
if isfield(data,'weightsB') && ~isempty(data.weightsB)
    weights = cat(5, data.weightsB, data.weightsX, data.weightsY, data.weightsZ);
else
    weights = logical(squeeze(abs(kdata(:,:,:,1,:,:)))); % ask this next time why weights are same as samp
end


%% Pre-processing
%---------------------------------------------------------------------------------------------------
% coil compression from actual number to opt.coil
disp('coil combining.....')
kdata = coilCombine(kdata,opt.coil,'3dt'); 

% compute time-averaged image (for walsh)
disp('computing time averaged image.....')
% first sum over frames and second over vencs
avg_k = sum(sum(kdata,5),6); 
% compute time-averaged pattern
% first sum over frames second over venc
avg_pattern = sum(sum(samp,4),5); avg_pattern_0 = logical(avg_pattern); 
% replacing 0s with inf to avoid div by zero error
avg_pattern(avg_pattern==0) = inf;
avg_k = bsxfun(@rdivide,avg_k,avg_pattern);
% 3D IFFT
avg_image = ifft3_shift(avg_k);

%% Sensitivity maps estimation

%********************************************
% disp('computing time averaged image.....')
% If you don't have sensitivity maps
% estimate walsh sensitivity maps
disp('estimating sensitivity maps.....')
p.fil = 3; % filter size
[maps,x0] = WalshCoilCombine3D(avg_image,p); %sensitivity maps & initial image x0
x0 = repmat(x0,[1,1,1,size(samp,4)]); %repeating static initial image over frames
%********************************************
%% Reconstruction via selected recon method
%---------------------------------------------------------------------------------------------------
% compute acceleration rate including partial fourier
acceleration_rate = numel(weights) / sum(weights(:).^2);
disp(['Acceleration rate: ', num2str(acceleration_rate)]);
%{ 
Note: provided sample k-space data has already been normalized by std. dev. 
of measurement noise during pre-processing of data, i.e. sig^2=1
%}
% data scaling
scale = 0.1 * max(abs(kdata(:)));

opt_scale = opt; 

for f=1:numel(recon)
t = tic; % Starting reconstruction timer
fprintf("\nRecon: "+recon(f)+"\n");
opt_scale.recon = recon(f);    % reconstruction method   
% xhat is the reconstructed image
xhat = zeros(size(squeeze(kdata(:,:,:,1,:,:))));
sname = [saveName,'_',datestr8601,'_',...
    convertStringsToChars(opt_scale.recon)];
% This loop runs over 4D Flow velocity encodings and just once for 3D cine
    for k = 1 :size(kdata, 6) % over encodings
        disp("****************Encoding: "+k+"********************");
        if opt.use_gpu
            gpuDevice(1);
        end
        [xhat(:,:,:,:,k),hist] = pMRIL14D(kdata(:,:,:,:,:,k)/scale, ...
        samp(:,:,:,:,k), weights(:,:,:,:,k), opt_scale, maps,x0/scale);
    end
fprintf('Elapsed Time = %0.2f minutes\n',toc(t)/60);

%% Magnitude Image Recovery
% Saving magnitude image in outputs struct
% and estimating magnitude image using Sum-of-squares for 4D flow
if(size(kdata, 6)>1) % for 4D flow
    
    outputs.xHat=zeros(size(xhat(:,:,:,:,1)));
for i=1:size(kdata, 6)
    outputs.xHat   = outputs.xHat+xhat(:,:,:,:,i).^2;
end
outputs.xHat = sqrt(outputs.xHat);
% Extracting just one magnitude image for 3D cine
else % for 3D cine
    outputs.xHat=xhat(:,:,:,:,1);
end

    
%% Background Phase Correction
% Bypass this section if don't intend to perform background phase
% correction on 4D flow
disp('Background phase correction...');
if(size(kdata, 6)>1) % for 4D flow only
  % Measuring correction map for x y z encodings
  % using backgroundphasecorrection(base encoding, velocity encoding)
    cmapx=backgroundCorrection3D(xhat(:,:,:,:,1),xhat(:,:,:,:,2)); %x encod
    cmapy=backgroundCorrection3D(xhat(:,:,:,:,1),xhat(:,:,:,:,3)); %y encod
    cmapz=backgroundCorrection3D(xhat(:,:,:,:,1),xhat(:,:,:,:,4)); %z encod
    % saving phase correction maps in outputs struct
    outputs.cmapx=cmapx; 
    outputs.cmapy=cmapy;
    outputs.cmapz=cmapz;
    
% Estimating flow data wrt base encoding and then applying 
% background phase correction using respective correction maps
% and saving phase images in outputs struct
    thetaXo =xhat(:,:,:,:,2).*conj(xhat(:,:,:,:,1));
    outputs.thetaX = angle(bsxfun(@times,thetaXo, exp(-1i*cmapx)));
    thetaYo =xhat(:,:,:,:,3).*conj(xhat(:,:,:,:,1));
    outputs.thetaY = angle(bsxfun(@times,thetaYo, exp(-1i*cmapy)));
    thetaZo =xhat(:,:,:,:,4).*conj(xhat(:,:,:,:,1));
    outputs.thetaZ = angle(bsxfun(@times,thetaZo, exp(-1i*cmapz)));
    
end
%%   Save reconstructed images to drive
%---------------------------------------------------------------------------------------------------

% saving the images in the sfolder path mentioned above
% the images will be saved under the name mentioned in sname
if ~exist([sfolder,'/',sname],'dir')
    mkdir([sfolder,'/',sname])
end
    
% save results
save(fullfile(sfolder,sname,[sname '.mat']), 'outputs');

% Saved image properties
gamma = 0.9; % gamma correction
fps = 10; % frames per second
clip = 2.5; % image clipping factor

%View 1 Saggital
create_GIF(abs(outputs.xHat.^gamma), fullfile(sfolder,sname, strcat(sname, '_xHat_Sag')), fps, clip);
if(size(kdata, 6)>1) %if 4D flow 
create_GIF_phase(norm_velocity(outputs.thetaX), fullfile(sfolder,sname, strcat(sname, '_thetaX_Sag')), fps, 1);
create_GIF_phase(norm_velocity(outputs.thetaY), fullfile(sfolder,sname, strcat(sname, '_thetaY_Sag')), fps, 1);
create_GIF_phase(norm_velocity(outputs.thetaZ), fullfile(sfolder,sname, strcat(sname, '_thetaZ_Sag')), fps, 1);
end

% View 2 Coronal
create_GIF(abs(permute(outputs.xHat, [1,3,2,4]).^gamma), fullfile(sfolder,sname, strcat(sname, '_xHat_Cor')), fps, clip);
if(size(kdata, 6)>1) %if 4D flow 
create_GIF_phase(norm_velocity(permute(outputs.thetaX, [1,3,2,4])), fullfile(sfolder,sname, strcat(sname, '_thetaX_Cor')), fps, 1);
create_GIF_phase(norm_velocity(permute(outputs.thetaY, [1,3,2,4])), fullfile(sfolder,sname, strcat(sname, '_thetaY_Cor')), fps, 1);
create_GIF_phase(norm_velocity(permute(outputs.thetaZ, [1,3,2,4])), fullfile(sfolder,sname, strcat(sname, '_thetaZ_Cor')), fps, 1);
end
mag_ax=permute(outputs.xHat, [2,3,1,4]);
if(size(kdata, 6)>1) %if 4D flow 
phasex_ax=permute(outputs.thetaX, [2,3,1,4]);
thetaXo_ax=permute(thetaXo, [2,3,1,4]);
cmapx_ax=permute(cmapx, [2,3,1,4]);

phasey_ax=permute(outputs.thetaY, [2,3,1,4]);
thetaYo_ax=permute(thetaYo, [2,3,1,4]);
cmapy_ax=permute(cmapy, [2,3,1,4]);

phasez_ax=permute(outputs.thetaZ, [2,3,1,4]);
thetaZo_ax=permute(thetaZo, [2,3,1,4]);
cmapz_ax=permute(cmapz, [2,3,1,4]);

end

% View 3 Tranverse
create_GIF(abs(mag_ax.^gamma), fullfile(sfolder,sname, strcat(sname, '_xHat_Ax')), fps, clip);

if(size(kdata, 6)>1)
create_GIF_phase(norm_velocity(phasex_ax), fullfile(sfolder,sname, strcat(sname, '_thetaX_Ax')), fps, 1);
create_GIF_phase(norm_velocity(phasey_ax), fullfile(sfolder,sname, strcat(sname, '_thetaY_Ax')), fps, 1);
create_GIF_phase(norm_velocity(phasez_ax), fullfile(sfolder,sname, strcat(sname, '_thetaZ_Ax')), fps, 1);
end
% View 4 
create_GIF(abs(squeeze(mag_ax(:,:,30,:)).^gamma), fullfile(sfolder,sname, strcat(sname, '_xhat_ax_V4')), fps, clip);
if(size(kdata, 6)>1)
create_GIF_phase(norm_velocity(squeeze(phasex_ax(:,:,30,:))), fullfile(sfolder,sname, strcat(sname, '_thetaX_Ax_V4')), fps, 1);
create_GIF_phase(norm_velocity(squeeze(phasey_ax(:,:,30,:))), fullfile(sfolder,sname, strcat(sname, '_thetaY_Ax_V4')), fps, 1);
create_GIF_phase(norm_velocity(squeeze(phasez_ax(:,:,30,:))), fullfile(sfolder,sname, strcat(sname, '_thetaZ_Ax_V4')), fps, 1);
end

% Saving background phase correction maps
if(size(kdata, 6)>1)
map=figure;
subplot(331); imagesc(mean(angle(thetaXo_ax(:,:,30,:)),4), [-0.1, 0.1]); axis('off','image'); colormap('jet'); title('uncorrected phase x');colorbar;
subplot(332); imagesc(mean(phasex_ax(:,:,30,:),4), [-0.1, 0.1]); axis('off','image'); colormap('jet'); title('corrected phase x');colorbar;
subplot(333); imagesc(cmapx_ax(:,:,36), [-0.1, 0.1]); axis('off','image'); colormap('jet'); colorbar;
title("X Encoding cMap Axial");

subplot(334); imagesc(mean(angle(thetaYo_ax(:,:,30,:)),4), [-0.1, 0.1]); axis('off','image'); colormap('jet'); title('uncorrected phase y');colorbar;
subplot(335); imagesc(mean(phasey_ax(:,:,30,:),4), [-0.1, 0.1]); axis('off','image'); colormap('jet'); title('corrected phase y');colorbar;
subplot(336); imagesc(cmapy_ax(:,:,36), [-0.1, 0.1]); axis('off','image'); colormap('jet'); colorbar;
title("Y Encoding cMap Axial");

subplot(337); imagesc(mean(angle(thetaZo_ax(:,:,30,:)),4), [-0.1, 0.1]); axis('off','image'); colormap('jet'); title('uncorrected phase z');colorbar;
subplot(338); imagesc(mean(phasez_ax(:,:,30,:),4), [-0.1, 0.1]); axis('off','image'); colormap('jet'); title('corrected phase z');colorbar;
subplot(339); imagesc(cmapz_ax(:,:,36), [-0.1, 0.1]); axis('off','image'); colormap('jet'); colorbar;
title("Z Encoding cMap Axial");

exportgraphics(map,fullfile(sfolder,sname, 'maps.jpg'));
end


fid = fopen(fullfile(sfolder,sname,'params.txt'),'wt');
fprintf(fid, 'Binned Filename:%s\nRecon:%s\n# of Iter:%s\nGamma:%s\nClip:%s',fname,opt_scale.recon,num2str(opt_scale.nit),num2str(gamma),num2str(clip));
fclose(fid);

writestruct(opt,fullfile(sfolder,sname,'params.xml'));


end



