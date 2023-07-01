function [ kdata_coil, noise_ch] = coilCombine( kdata_coil, n_channels, dataType)
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

if strcmpi(dataType,'2dt')

    % Channel compression =================================================
    n_coils = size(kdata_coil,3);
    if  n_coils > n_channels
        fprintf('Compressing %d channels to %d\n',size(kdata_coil,3),n_channels)

        % Performing SVD for Channel Compression (optional) 
        N=size(kdata_coil);
        N = [N,1,1];
        kdata_coil = permute(kdata_coil,[1,2,4,5,3]);
        kdata_coil = reshape(kdata_coil,[N(1)*N(2)*N(4)*N(5),N(3)]);

        [~,S,V] = svd(kdata_coil'*kdata_coil,0);

        % Get the noise channel
        noise_ch = kdata_coil*V(:,:);
        noise_ch = noise_ch(:,end);
        noise_ch = reshape(noise_ch,[N(1),N(2),N(4),N(5),1]);
        noise_ch = permute(noise_ch,[1,2,5,3,4]);

        kdata_coil = kdata_coil*V(:,1:n_channels);
        kdata_coil = reshape(kdata_coil,[N(1),N(2),N(4),N(5),n_channels]);
        kdata_coil = permute(kdata_coil, [1,2,5,3,4]);
    end

    if n_coils <= n_channels
        noise_ch = kdata_coil(:,:,end,:);
        fprintf('No Channel Compression Used.  Using last coil for noise');
    end
elseif strcmpi(dataType,'3dt')

    % Channel compression =================================================
    n_coils = size(kdata_coil,4);
    if  n_coils >= n_channels
        fprintf('Compressing %d channels to %d\n',size(kdata_coil,4),n_channels)

        % Performing SVD for Channel Compression (optional) 
        N=size(kdata_coil);
        N = [N,1];
        kdata_coil = permute(kdata_coil,[1,2,3,5,6,4]);
        kdata_coil = reshape(kdata_coil,[N(1)*N(2)*N(3)*N(5)*N(6),N(4)]);

        [~,S,V] = svd(kdata_coil'*kdata_coil,0);

        % Get the noise channel
        noise_ch = kdata_coil*V(:,:);
        noise_ch = noise_ch(:,end);
        noise_ch = reshape(noise_ch,[N(1),N(2),N(3),N(5),N(6),1]);
        noise_ch = permute(noise_ch,[1,2,3,6,4,5]);
        noise_ch = squeeze(noise_ch);

        kdata_coil = kdata_coil*V(:,1:n_channels);
        kdata_coil = reshape(kdata_coil,[N(1),N(2),N(3),N(5),N(6),n_channels]);
        kdata_coil = permute(kdata_coil, [1,2,3,6,4,5]);
        kdata_coil = squeeze(kdata_coil);
    end

    if n_coils <= n_channels
        noise_ch = kdata_coil(:,:,:,end,:,:);
        fprintf('No Channel Compression Used.  Using last coil for noise');
    end
else
    error('Unknown data type.  Use 2Dt or 3Dt');
end

end