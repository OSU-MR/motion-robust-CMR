%% ESPIRiT Maps Demo
% This is a demo on how to generate ESPIRiT maps. It is based on the paper
% Uecker et. al, MRM 2013 DOI 10.1002/mrm.24751. ESPIRiT is a method that
% finds the subspace of multi-coil data from a calibration region in
% k-space using a series of eigen-value decompositions in k-space and image
% space. 
function [cI, maps] = espirit_sens3D(DATA, samp, ncalib, param)
% eigThresh_1 = param.eSRT(1);
eigThresh_2 = param.eSRT;
ksize = param.fil;
ksize = min(ksize,ncalib);
map_num = param.eSmaps;

[sx,sy,sz,Nc] = size(DATA);

% ncalib = param.ACSsz;
% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.0001; % Larger value gives smoother sensitivity maps; Default: 0.02

% threshold of eigen vector decomposition in image space. 
% eigThresh_2 = 0.95; % Larger value gives more localized maps; Defalut: 0.95

% crop a calibration area
calib = crop(DATA,[ncalib,Nc]);
samp = crop(samp,ncalib(1:ndims(samp)));
%%
% Display coil images: 
im = DATA;
for n = 1:3
im = ifftc(im,n);
end

% figure;
% for j = 1 : size(im,3)
% imshow3(abs(squeeze(im(:,:,j,:,:))),[],[2,Nc]);
% title(num2str(j));
% % title('magnitude of physical coil images');
% colormap((gray(256))); %colorbar;
% pause();
% end
% 
% figure, imshow3(angle(im),[],[1,Nc]); 
% title('phase of physical coil images');
% colormap('default'); colorbar;

%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 


% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

[k,S,dim_A] = dat2Kernel3D(calib,samp,ksize);
% % estimate the noise variance
zS = max(0,(dim_A(2)-numel(S)));
S_tmp = padarray(S,[zS 0],0,'post');

seed = 100;
nStd0 = 1;
rng(100*seed); n_real = randn(dim_A);
rng(200*seed); n_im = randn(dim_A);
noise  = nStd0*complex(n_real,n_im)/sqrt(2); % noise
[~,S_noise,~] = svd(noise,'econ');
S_noise = diag(S_noise);
sigma = ((mean(S(ceil(3*end/4):end))/mean(S_noise(ceil(3*end/4):end))));
clear noise;
% 
% % soft threshold lambda
% lambda = S(end-1); %inital guess of lambda
% IS_REAL = 0;
% MSE = sure_svt(lambda, sigma, S, dim_A, IS_REAL);
% % MSEcc = zeros(numel(S),1);
% for jdx = 1:(numel(S)-1) % figure out the optimal lambda
%     testlambda = S(end-jdx);
%     testMSE = sure_svt(testlambda, (sigma), S, dim_A, IS_REAL);
% %     MSEcc(end-jdx) = testMSE; 
%     if testMSE < MSE
%         lambda = testlambda;
%         MSE = testMSE;
%     end   
% end
% % weight = max(0,S-lambda)./(abs(S-lambda)+eps);
% weight = max(0,S-lambda)./(S+eps);
% idx = max(find(weight > 0));


idx = max(find(S >= S(1)*eigThresh_1));

% idx = min(idx, 300); % to avoid memory overload

% if 'idx' is larger than 40% of the calibration matrix rows, cutoff (eigThres_1) is
% not resonable.
if ( double(idx)/numel(S(:)) >= 0.4 )
    idx = ceil(numel(S(:))*0.35);
    disp('p.eSRT(1) is not reasonable!');
end

%% 
% This shows that the calibration matrix has a null space as shown in the
% paper. 
disp([num2str(100*double(idx)/numel(S(:))) '% of singular vectors are included.']);
if ( double(idx)/numel(S(:)) >= 0.7 || double(idx)/numel(S(:)) <= 0.2)
     warning('ESPiRiT parameter estimated by SURE may not be resonable!');
     kdisp = reshape(k,[ksize(1)*ksize(2)*ksize(3)*Nc,ksize(1)*ksize(2)*ksize(3)*Nc]);
     figure, subplot(211), plot([1:ksize(1)*ksize(2)*ksize(3)*Nc],S_tmp,'LineWidth',2);
     hold on, 
     plot([1:ksize(1)*ksize(2)*ksize(3)*Nc],S_tmp(idx),'r-','LineWidth',2);
     plot([idx,idx],[0,S_tmp(1)],'g--','LineWidth',2)
     legend('signular vector value','threshold')
     title('Singular Vectors')
     subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
     xlabel('Singular value #');
     title('Singular vectors')
end

% eigThresh_2 = min(eigThresh_2, sum(S(1:idx))/sum(S));

%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
% for jdx = 1:numel(weight) %additional weight for eigen vectors
%     k(:,:,:,:,jdx) = k(:,:,:,:,jdx).*weight(jdx);
% end
[M,W] = kernelEig3D(k(:,:,:,:,1:idx),[sx,sy,sz]);

%%
% show eigen-values and eigen-vectors. The last set of eigen-vectors
% corresponding to eigen-values 1 look like sensitivity maps

% M_center=squeeze(M(:,:,ceil(size(im,3)/2),:,:));
% W_center=squeeze(W(:,:,ceil(size(im,3)/2),:));
% 
% figure, imshow3(abs(W_center),[],[1,Nc]); 
% title('Eigen Values in Image space');
% colormap((gray(256))); colorbar;
% 
% figure, imshow3(abs(M_center(:,:,:,:)),[],[Nc,Nc]); 
% title('Magnitude of Eigen Vectors');
% colormap(gray(256)); colorbar;


%%
% project onto the eigenvectors. This shows that all the signal energy
% lives in the subspace spanned by the eigenvectors with eigenvalue 1.
% These look like sensitivity maps. 

P = sum(repmat(im,[1,1,1,1,Nc]).*conj(M),4);

% figure, imshow3(abs(squeeze(P(:,:,ceil(size(im,3)),:))),[],[1,Nc]); 
% title('Magnitude of Eigen Vectors');
% colormap(sqrt(gray(256))); colorbar;


%%
% crop sensitivity maps 
maps = M(:,:,:,:,end-map_num+1:end);
cI = permute(P(:,:,:,1,end-map_num+1:end),[1,2,3,5,4]); % Coil combined image,(kx,ky,kz,set)
% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,:,end-map_num+1:end) ;
weights = (weights - eigThresh_2)./(1-eigThresh_2).* (W(:,:,:,end-map_num+1:end) > eigThresh_2);
weights = -cos(pi*weights)/2 + 1/2;
maps = maps.*sqrt(repmat(permute(weights,[1 2 3 5 4]),[1 1 1 Nc 1]));

% maps_c = squeeze(maps(:,:,ceil(size(im,3)/2),:,1));
% 
% figure, imshow3(abs(maps_c),[],[map_num,Nc]); 
% title('Absolute sensitivity maps');
% colormap((gray(256))); colorbar;
% 
% figure, imshow3(angle (maps_c),[],[map_num,Nc]); 
% title('Phase of sensitivity maps');
% colormap((jet(256))); colorbar;



