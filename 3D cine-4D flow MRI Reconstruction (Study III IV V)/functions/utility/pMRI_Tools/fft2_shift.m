function [ y ] = fft2_shift( x )
%Unitary 2-D FFT with fftshifts
%   Detailed explanation goes here

y = 1/sqrt( (size(x,1)*size(x,2)) )*ifftshift( fft2( fftshift(x) ));
% y = 1/sqrt( (size(x,1)*size(x,2)) )*fft2(x);
end

