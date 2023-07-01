% Adjoint operator
function y = funAt(x,s)
    y = sqrt(numel(x)) * fftshift(ifft2(ifftshift(x .* s))); % downsampling and then fft
end