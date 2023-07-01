function [kernel, S, dim_A] = dat2Kernel3D(data,samp, kSize)
% kernel = dat2Kernel(data, kSize,thresh)
%
% Function to perform k-space calibration step for ESPIRiT and create
% k-space kernels. Only works for 2D multi-coil images for now.  
% 
% Inputs: 
%       data - calibration data [kx,ky,coils]
%       kSize - size of kernel (for example kSize=[6,6])
%
% Outputs: 
%       kernel - k-space kernels matrix (not cropped), which correspond to
%                the basis vectors of overlapping blocks in k-space
%       S      - (Optional parameter) The singular vectors of the
%                 calibration matrix
%
%
% See also:
%           kernelEig
%
% (c) Michael Lustig 2013



[sx,sy,sz,nc] = size(data);
imSize = [sx,sy,sz] ;

tmp = im2row3D(data,kSize); [tsx,tsy,tsz] = size(tmp);
A = reshape(tmp,tsx,tsy*tsz);

% if there is unsampled data, discard the row
samp = repmat(samp,[1 1 1 nc]);
tmp = im2row3D(samp,kSize);
index = find(sum(reshape(tmp,tsx,tsy*tsz),2) == tsy*tsz); % row index without unsampled data
A = A(index,:); 
dim_A = size(A);
% try covariance trick at some point
% [U,S,V] = svd(single(A),'econ');

[U,S,V] = svd(A);
    
kernel = reshape(V,kSize(1),kSize(2),kSize(3),nc,size(V,2));
S = diag(S);S = S(:);
