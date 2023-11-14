function [cMap] = backgroundCorrection3D(b, x )

%% options
opt.pOrd        = 3;            % polynomial order for fitting, [1, 2, 3, or 4]
opt.lam         = 5;           % regularization strength for L1-norm regularization
opt.mTh         = 0.04;         % thresholding for magnitude mask (fraction of maximum)
opt.midFOVFrac  = 0.5;          % Center fraction of pixels fit in midPE-FOV for the 0th iteration of ARTO.
opt.Kmax        = 2;            % Maximum number of ARTO iterations
opt.tau         = 3;            % Exclusion threshold for ARTO
opt.delta       = 2;            % Minimum separation b/w Gaussian components in GMM
opt.gmmComp     = 3;            % Number of components fit by GMM (editing this may break the constraints specified in gmmArto().)



%% Pre-processing

mag_cine = abs(b);
phase_cine = angle(x.*conj(b));


phase_cine = phase_cine / pi;
% Average phase cine over time
phi = mean(phase_cine, 4);
% Temporal standard deviation
sigma = std(phase_cine,0,4);
% Magnitude thresholding
mag_avg = mean(mag_cine,4);
M = ones(size(mag_avg));
% figure; imagesc(abs(squeeze(M(:,:,end/2))));
% disp(['zero fraction: ' num2str(sum(M(:)==0)/numel(M))])
M(mag_avg <= opt.mTh*(max(max(mag_avg)))) = 0;
% figure; imagesc(abs(squeeze(M(:,:,end/2))));
% disp(['zero fraction: ' num2str(sum(M(:)==0)/numel(M))])
% tmp = M.*sigma;
% figure; imagesc(abs(squeeze(tmp(:,:,end/2)))); colorbar;
M = logical(M);


[~,cMap] = wrls_arto(phi,sigma,M,opt);

cMap = cMap * pi;


end

