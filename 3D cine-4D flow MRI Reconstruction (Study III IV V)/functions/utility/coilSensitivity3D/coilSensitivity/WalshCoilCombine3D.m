function [S, cI] = WalshCoilCombine3D(I, param)
%
% ------------------------ input ------------------------
% I: coil-by-coil complex images [FE, PE, PE2, Coils]. 
% param: input parameters
%        param.fil = filter size used to smoothen E[X'X] matrix
%        param.opt = type of phase correction
       

% ------------------------ output ------------------------
% cI: coil-combined image
% S: Sensitivity maps

% ---------------------- References ----------------------   
% Polarimetric techniques for enhancing SAR imagery, SPIE, 1630: 141-173 1992
% Adaptive Reconstruction of Phased Array MR Imagery, MRM 43:682-690 2000


% Check input variables
if nargin~=1 && nargin~=2
    error('Incorrect number of input arguments');
    
elseif nargin==1 
    param.fil = 9;
    param.opt = 2;
end

if ~isfield(param, 'opt')
    param.opt = 2;
elseif param.opt<1 || param.opt>3
    param.opt = 2;
    warning('param.opt is set to 2');
end

if ~isfield(param, 'fil')
    param.fil = 9;
elseif param.fil<1
    param.fil = 9;
    warning('param.fil is set to 9');
end
    


N = size(I); % Image size
FE = N(1); % Size of frequency encoding direction
PE = N(2); % Size of phase encoding direction
PE2 = N(3);% Size of the second phase encoding direction
Nc = N(4); % Number of coils
param.fil = 2*(floor(param.fil/2))+1; % To ensure it is odd

cI = zeros(FE,PE,PE2);
S  = zeros(FE,PE,PE2,Nc);


%% Bulid "covariance" matrix covI
covI = zeros(FE,PE,PE2,Nc,Nc);
for k=1:Nc
    covI(:,:,:,k,k)=I(:,:,:,k).*conj(I(:,:,:,k));
    for j=1:k-1
		covI(:,:,:,k,j)=I(:,:,:,k).*conj(I(:,:,:,j));
		covI(:,:,:,j,k)=conj(covI(:,:,:,k,j));
    end
end

% Spatially smoothen covI to generate covIs
fil3 = ones(param.fil, param.fil, 1);
covIs = zeros(FE,PE,PE2,Nc,Nc);
for i = 1:Nc
    for j=1:Nc
        covIs(:,:,:,i,j) = convn(covI(:,:,:,i,j), fil3, 'same');
%         covIs(:,:,i,j) = medfilt2(real(covI(:,:,i,j)), [13,13]) + ...
%                       1j* medfilt2(imag(covI(:,:,i,j)), [13,13]);
%         figure; imagesc(squeeze(imag(covIs(:,:,i,j))));
%         covIs(:,:,i,j) = conv2(real(covI(:,:,i,j)), fil, 'same') + 1j*conv2(imag(covI(:,:,i,j)), fil, 'same');
    end
end

% Note: This way of smoothing covI has been suggested by Peter Kellman and
% reportedly works well.


%% Performing coil combine and senstitivity estimation per Walsh et al.
for k = 1:PE2
    for i = 1:PE
        for j = 1:FE               
            covIs_2D = reshape(covIs(j, i, k, :, :), [Nc Nc]);
            [U,~] = eig((covIs_2D + covIs_2D')/2); 
            cI(j,i,k) = reshape(I(j, i, k, :), [1 Nc]) * conj(U(:,end));
            S(j,i,k,:) = reshape(U(:,end), [1 1 1 Nc]);   
    %         S(j,i,:) = reshape(abs(U(:,end)).^1.5 .* (U(:,end)./abs(U(:,end))), [1 1 Nc]);   


    %         covIs_2D = reshape(covIs(j, i, :, :), [Nc Nc]);
    %         [U,D] = eig((covIs_2D + covIs_2D')/2);        
    %         loads{1} = U*D;
    %         loads{2} = U;
    %         [~,loads] = sign_flip(loads,(covIs_2D + covIs_2D')/2);
    %         cI(j,i) = reshape(I(j, i, :), [1 Nc]) * conj(loads{2}(:,end));
    %         S(j,i,:) = reshape(loads{2}(:,end), [1 1 Nc]);    
        end
    end
end

[cI, S] = sensCorrect(cI, S, param.opt);
% S = bsxfun(@rdivide, S, sum(abs(S), 3));



function [xc, sc] = sensCorrect(x, s, option)
% Bilinear nature of pMRI introduces ambiguity in the estimation of coil
% sensitivies. Here, we devise a correction so  that image (or sensitivies)
% observe an expected behavior.

% Input ===================================================================
% x: [FE, PE, PE2]; coil combined complex image

% s: [FE, PE, PE2, Nc]; complex sensitivities

% option: '1' to achieve real image, '2' to achieve closeness to global
% phase. See "Notes" below for further explnation

% dsp: '1' to display and '0' not to display results

% Ouput ===================================================================
% xc: output image

% sc: output sensitivity maps


% Notes ===================================================================
% When option =1; do nothing

% When option = 2, each pixel is multiplied with a unit norm complex number 
% so that all pixels become real. We also multiply sensitivity maps with
% the conjugate of that complex number

% When option = 3, image is multiplied (pixel-wise) with a binary map so 
% that phase of each pixel is "closer" to the global phase. We also 
% multiply sensitivity maps with the same binary map.


if nargin<2
    error('Not enough in put arguments');
elseif nargin==2
    option = 2;
    dsp = 0;
elseif nargin==3
    dsp = 0;
end

if option~=1 && option~=2 && option~=3
    error('Incorrect value assigned to image correction option');
elseif dsp~=0 && dsp~=1
    error('Incorrect value assigned to display option');
end


[~,~,~,Nc] = size(s); % Nc: number of coils

% Option 1
if option == 1
    % do nothing
    xc = x;
    sc = s;

% Option 2
elseif option == 2
    phs = atan2(imag(x), real(x));
    map = exp(-1j*phs);
    xc = x.*map;
    sc = s.*repmat(conj(map),[1,1,1,Nc]);
    
% Option 3        
elseif option == 3
    avgPhase = atan2(imag(sum(x(:))), real(sum(x(:)))); % Global phase
    tmp = abs(x)*exp(1j*avgPhase);

    % Angle b/w global phase and the phase at each pixel
    phs = acos((real(x).*real(tmp) + imag(x).*imag(tmp))./(abs(x).*abs(tmp)));

    % Binary map which tells which pixels to "flip" so that phase at each pixel is closer to the global phase
    map = 2*(phs>(pi/2))-1; 
    % figure; imagesc(phs2); axis('image'); colorbar;

    xc = x.*map;
    sc = s.*repmat(conj(map),[1,1,1,Nc]);
end


% if dsp == 1
%     figure; subplot(131); imagesc(real(x)); axis('image'); colorbar; title('Input, real');
%             subplot(132); imagesc(imag(x)); axis('image'); colorbar; title('Input, imag');
%             subplot(133); imagesc(atan2(imag(x),real(x))); axis('image'); colorbar; title('Input, phase');
%         
%     figure; subplot(131); imagesc(real(xc)); axis('image'); colorbar; title('Output, real');
%             subplot(132); imagesc(imag(xc)); axis('image'); colorbar; title('Output, imag');
%             subplot(133); imagesc(atan2(imag(xc),real(xc))); axis('image'); colorbar; title('Output, phase');
% elseif dsp == 0
%     % Do nothing
% end
