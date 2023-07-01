function [ w ] = shrinkWiener3D(x,sizes,noise)
%SHRINKWIENER3D Summary of this function goes here
%   Detailed explanation goes here
    
    x = reshape(x,sizes);
    x = padarray(x,[0,0,1],'replicate','both');
    x = fftshift(fftshift(x,1),2);
    w = wiener3D( x, noise );
    w = ifftshift(ifftshift(w,1),2);
    w = w(:,:,2:end-1);
    w = w(:);
end

