function y =  funA(x, sampInd, N)
% Forward operator for single-coil MRI
    y = fft2r(x);
    y = y(sampInd);
end
