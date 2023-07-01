function y = shrinkWaveWiener3D(x,sizes,W,lam)
%SHRINKWAVEWIENER3D Summary of this function goes here
%   Detailed explanation goes here
%     x = fftshift(fftshift(reshape(x,sizes),1),2);
%     x = padarray(x,[0,0,3],'replicate','both');
    w = W.mult(x);
    w = reshape(w,[sizes,8]);
%     w = padarray(w,[0,0,3],'replicate','both');
%     tmp = w(:,:,:,1);
    y = wiener3D(w, lam);
%     yreal = wiener3D(real(w),lam/(2));
%     yimag = wiener3D(imag(w),lam/(2));
%     y = yreal+1j*yimag;
%     y = y(:,:,4:end-3,:);
%     y(:,:,:,1) = tmp;
    y = W.multTr(y(:));
%     y = ifftshift(ifftshift(reshape(y,sizes),1),2);
%     y = y(:,:,4:end-3,:);
    y = y(:);
end


