function [xc, sc] = sensCorrect(x, s, option)
% option = 1;

% Bilinear nature of pMRI introduces ambiguity in the estimation of coil
% sensitivies. Here, we device a correction so  that image (or sensitivies)
% observe an expected behavior.

% Input ===================================================================
% x: [FE, PE]; coil combined complex image

% s: [FE, PE, Nc]; complex sensitivities

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
    option = 1;
    dsp = 0;
elseif nargin==3
    dsp = 0;
end

if option~=1 && option~=2
    error('Incorrect value assigned to image correction option');
elseif dsp~=0 && dsp~=1
    error('Incorrect value assigned to display option');
end

[~,~,Nc] = size(s); % Nc: number of coils


if option == 0
    % do nothing
    xc = x;
    sc = s;

% assign the average phase to the image
elseif option == 1
    phs = atan2(imag(x), real(x));
    map = exp(-1j*phs);
    xc = x.*map;
    sc = s.*repmat(conj(map),[1,1,Nc]);
    
% % Option 3        
% elseif option == 3
%     avgPhase = atan2(imag(sum(x(:))), real(sum(x(:)))); % Global phase
%     tmp = abs(x)*exp(1j*avgPhase);
% 
%     % Angle b/w global phase and the phase at each pixel
%     phs = acos((real(x).*real(tmp) + imag(x).*imag(tmp))./(abs(x).*abs(tmp)));
% 
%     % Binary map which tells which pixels to "flip" so that phase at each pixel is closer to the global phase
%     map = 2*(phs>(pi/2))-1; 
%     % figure; imagesc(phs2); axis('image'); colorbar;
% 
%     xc = x.*map;
%     sc = s.*repmat(conj(map),[1,1,Nc]);
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
