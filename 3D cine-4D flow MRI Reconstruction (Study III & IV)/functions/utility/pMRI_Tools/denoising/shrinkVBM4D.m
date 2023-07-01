function y = shrinkVBM4D(x,sizes,sigma)
%SHRINKWAVEWIENER3D Summary of this function goes here
%   Detailed explanation goes here
%     x = real(gather(x));
    x = gather(x);
    profile   = 'lc';
%     profile   = 'np';
    do_wiener = 1; 
    sharpen = 1;  
    deflicker = 1;  
    verbose = 0;     
    x = reshape(x,sizes);
    x = fftshift(fftshift(x,1),2);
    [ y ] = vbm4d( x, sigma, profile, do_wiener, sharpen, deflicker, verbose );
    y = ifftshift(ifftshift(y,1),2);
    y = gpuArray(y(:));

end


