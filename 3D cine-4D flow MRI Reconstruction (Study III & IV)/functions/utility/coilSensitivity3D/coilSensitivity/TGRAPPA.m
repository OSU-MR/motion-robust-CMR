function [Recon] = TGRAPPA(kdata_coil,samp,R,KernelSize)

% kdata_coil: 4D data
% samp: 4D binary sampling

% s_0 = size(kdata_coil);
% option.R = R;
% option.Interleave = 1 ;

% samp = sampTest(s_0,R);
disp('Performing GRAPPA initialization ...');
tic,
kdata_coil = kdata_coil.*samp;


Ref_Img = sum(kdata_coil, 4)./(sum(samp,4)+eps); % Fix this one if number of line/frame are different
calibSize = getCalibSize(logical(sum(abs(Ref_Img),3)));
Ref_Img = crop(Ref_Img, [calibSize,size(Ref_Img,3)]);

GRAPPA_option_0.KernelSize = KernelSize;
GRAPPA_option_0.KernelPattern =  -(KernelSize(1)/2-1)*R:R:KernelSize(1)/2*R;

GRAPPA_option_0.OutPattern = [1:R-1];
[Coef, ~,~] = GRAPPA_Kernel_2D(Ref_Img, GRAPPA_option_0);

% Imgs = zeros(size(kdata_coil,1), size(kdata_coil,2), size(kdata_coil,4));
Recon=zeros(size(kdata_coil));
% for j = 1:1
for i=1:size(kdata_coil,4)
    GRAPPA_option_0.SamplingPattern = find(squeeze(max(samp(:,:,1,i),[],1))'==1);
%     Reduced_Kspace = kdata_coil(:, :, :, i);
    Recon(:,:,:,i) = GRAPPA_Recon_2D(kdata_coil(:, :, :, i), Coef, GRAPPA_option_0);
%         Imgs(:,:,i) = SoS( Kspace, ones(size(Kspace(:,:,:,i))), 0 ) ;
%         figure(2); imagesc([Imgs(:,:,i)]), title(num2str(i)), axis image, pause(0.01)
end
Recon(kdata_coil~=0) = kdata_coil(kdata_coil~=0);
disp(['GRAPPA initialization completed in ' sprintf('%0.2f', toc) ' s']);

% end

% slicer(Imgs,1);