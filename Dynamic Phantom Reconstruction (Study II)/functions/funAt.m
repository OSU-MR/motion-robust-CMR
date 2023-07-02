function x =  funAt(y, sampInd, N)
% Forward operator for single-coil MRI
    x = zeros(N);
    x(sampInd) = y;
    x = ifft2r(x);
end