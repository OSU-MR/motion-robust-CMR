% Forward operator
function y = funA(x,s)
    y = 1/sqrt(numel(x)) * fftshift(fft2(ifftshift(x))) .* s; % fft and then downsampling
end