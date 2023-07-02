function X =  fft2r(x)
% 2D FFT from "centered" image to "centered" k-space with normalization
    X = 1/sqrt(numel(x)) * fftshift(fftshift(fft2(ifftshift(ifftshift(x,1),2)),2),1);
end
