function [ kdata_coil, noise_ch] = coilCombine( kdata_coil, n_channels)
% Reads the Data given by fName into Matlab.  Also, Compresses the Coils
% and normalizes the Data
%
% Inputs:
%           fName: Filename of the Data .mat and .h5 supported
%
% Optional Inputs:
%           samp: User supplied sampling pattern for fully sample data only
%
% Outputs:
%           kdata_coil: compressed and normalized data
%
%           samp: Sampling Pattern if Data has alread been downsampled
%           
% star_display('Loading data and compressing coils',0)
% tic,

% No User Supplied number of channels or Sampling Pattern
if nargin <2
    samp = [];
	n_channels = 16; % Number of Compressd Channels
end

Nc = [n_channels/size(kdata_coil,3), n_channels];

%% Channel compression ====================================================
n_coils = size(kdata_coil,3);
if  n_coils > n_channels
    fprintf('Compressing %d channels to %d\n',size(kdata_coil,3),n_channels)
    N=size(kdata_coil);
    ncr = min(max(round(Nc(1)*N(3)), Nc(2)), N(3));
    
    % Performing SVD for Channel Compression (optional) 
    N=size(kdata_coil);
    kdata_coil = permute(kdata_coil,[1,2,4,3]);
    kdata_coil = reshape(kdata_coil,[N(1)*N(2)*N(4),N(3)]);

    % First "ncr" orthogonal channels
    [~,S,V] = svd(kdata_coil'*kdata_coil,0);

    % Get the noise channel
    noise_ch = kdata_coil*V(:,:);
    noise_ch = noise_ch(:,end);
    noise_ch = reshape(noise_ch,[N(1),N(2),N(4),1]);
    noise_ch = permute(noise_ch,[1,2,4,3]);
    
    kdata_coil = kdata_coil*V(:,1:ncr);
    kdata_coil = reshape(kdata_coil,[N(1),N(2),N(4),ncr]);
    kdata_coil = permute(kdata_coil, [1,2,4,3]);
%     kdata_coil = kdata_coil.*repmat(permute(samp,[1,2,4,3]),[1,1,size(kdata_coil,3),1]);
end

if n_coils <= n_channels
    noise_ch = kdata_coil(:,:,end,:);
    fprintf('No Channel Compression Used.  Using last coil for noise');
end

% Normalize
% fprintf('Normalizing The Data\n')
% scale_fctr = norm(kdata_coil(:))/sqrt(numel(squeeze(kdata_coil(:,:,:,1))));
% kdata_coil = kdata_coil/scale_fctr;
% noise_ch = noise_ch/scale_fctr;
% sprintf('\n%0.10f\n',scale_fctr)
% fprintf('Data loaded and compressed in %s s\n', num2str(toc));
% star_display('',0)
% fprintf('\n')

end