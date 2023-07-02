function x =  ifft2r(X)
% 2D inverse FFT from "centered" k-space to "centered" image with normalization
    x = sqrt(numel(X)) * fftshift(fftshift(ifft2(ifftshift(ifftshift(X,1),2)),2),1);
end