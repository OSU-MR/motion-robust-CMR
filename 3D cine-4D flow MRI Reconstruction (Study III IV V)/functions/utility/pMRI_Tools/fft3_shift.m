function [ y ] = fft3_shift( x )
%FFT3_SHIFT take a unitary fft on the first three dimensions of an array with
%fftshifts

y = fftshift(fftshift(fftshift(x,1),2),3);
y = fft(fft(fft(y,[],1),[],2),[],3);
y = 1/sqrt(size(x,1)*size(x,2)*size(x,3))*ifftshift(ifftshift(ifftshift(y,1),2),3);


end

