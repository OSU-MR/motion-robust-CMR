function [ xhatGAMP ] = GAMP4DpMRI_func(opts)
%GAMP4DPMRI Reconstructs 4D pMRI data using GAMP and 4D wavelets
%   
% fftshifts
opts.data = fftshift(fftshift(fftshift(opts.data,1),2),3);
opts.samp = fftshift(fftshift(fftshift(opts.samp,1),2),3);

% Coil combine
opts.data = coilCombine(opts.data,12,'3dt');
% Create time averaged image
avg_image = sum(opts.data,5);
avg_pattern = sum(opts.samp,4);
avg_pattern(avg_pattern==0) = inf;
avg_image = bsxfun(@rdivide,avg_image,avg_pattern);
avg_image = ifft(ifft(ifft(avg_image,[],1),[],2),[],3);
sos = sqrt(sum(abs(avg_image).^2,4));
maps = repmat(sos,[1,1,1,size(opts.data,4)]);
sos = repmat(sos,[1,1,1,size(opts.samp,4)]);
opts.GAMPopt.xhat0 = sos(:);

%% Estimate sensitivity maps
% maps = avg_image./maps;
% maps = repmat(maps,[1,1,1,1,size(opts.data,5)]);

p.mthd   = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
p.reEst  = 0; % Res-estimating sensitivities
p.fil = 3;
[maps,~] = WalshCoilCombine3D(avg_image,p);
% for ind = 1:size(maps,4)
%     figure(1)
%     subplot(3,4,ind)
%     imagesc(fftshift(fftshift(abs(maps(:,:,1,ind,1)),1),2) )
%     
%     figure(2)
%     subplot(3,4,ind)
%     imagesc( fftshift(fftshift(angle(maps(:,:,1,ind,1)),1),2))
% end
% pause

% Downsample the opts.data
opts.data = downsample_data(opts.data,opts.samp);

% Normalize The columns of A to be unit norm
R = numel(opts.samp)/length(find(opts.samp ==1));
maps = maps*sqrt(R);
opts.data = opts.data*sqrt(R);
opts.wvar = opts.wvar*R;

% Create Opterators
pMRIOp = pMRI_Op_3D_t(maps,opts.samp,1);
EstimIn = NullEstimIn(0,1);
op = ReVEAL();
op.options.SparseTrans.wname = {'db1','db1','db1','db1'};
W = ndDWTLinTrans4D(op.options.SparseTrans,size(opts.samp),1);

% dim = round(rand(1,1)*120000);
dim = 1;
x = zeros(size(opts.samp));
x(dim) = 1;
a = pMRIOp.mult(x(:));
norm(a(:))

x = zeros(size(opts.data));
x(dim) = 1;
ah = pMRIOp.multTr(x(:));
a(dim)
ah(dim)
norm(a)
norm(ah)
% Concatonate All The Linear Transform Operator Together
Op = LinTransConcat({pMRIOp;W});

% Input and outputs esimators
MeasEstimOut = CAwgnEstimOut(opts.data(:),opts.wvar,1);
opts.lambda = setLambda(size(opts.samp),opts.lambda);
AnaEstimOut1 = CplxLaplaceEstimOut(opts.lambda);
EstimOut = EstimOutConcat({MeasEstimOut;AnaEstimOut1},[pMRIOp.M,W.M]);
% size(opts.data)
% numel(opts.samp)
% pMRIOp.M
% Reconstruct
tic;
xhatGAMP = gampEst(EstimIn,EstimOut,Op,opts.GAMPopt);
display(sprintf('Reconstruction Completed in %s',num2str(toc)))
xhatGAMP = reshape(xhatGAMP,size(opts.samp));
xhatGAMP = fftshift(fftshift(fftshift(xhatGAMP,1),2),3);

end

