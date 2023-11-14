function [EigenVecs, EigenVals] = kernelEig3D(kernel, imSize)
% [eigenVecs, eigenVals] = kernelEig(kernel, imSize)
%
% Function computes the ESPIRiT step II -- eigen-value decomposition of a 
% k-space kernel in image space. Kernels should be computed with dat2Kernel
% and then cropped to keep those corresponding to the data space. 
%
% INPUTS:
%           kernel - k-space kernels computed with dat2Kernel (4D)
%           imSize - The size of the image to compute maps for [sx,sy]
%
% OUTPUTS:
%           EigenVecs - Images representing the Eigenvectors. (sx,sy,Num coils,Num coils)
%           EigenVals - Images representing the EigenValues. (sx,sy,numcoils )
%                       The last are the largest (close to 1)
%           
% 
% See Also:
%               dat2Kernel
% 
%
% (c) Michael Lustig 2010

nc = size(kernel,4);
nv = size(kernel,5);
kSize = [size(kernel,1), size(kernel,2), size(kernel,3)];

% "rotate kernel to order by maximum variance"
k = permute(kernel,[1,2,3,5,4]);, k =reshape(k,prod([kSize,nv]),nc);

if size(k,1) < size(k,2)
    [u,s,v] = svd(k);
else
    
    [u,s,v] = svd(k,'econ');
end

k = k*v;
kernel = reshape(k,[kSize,nv,nc]); kernel = permute(kernel,[1,2,3,5,4]);

%%
KERNEL = zeros(imSize(1), imSize(2), imSize(3), size(kernel,4), size(kernel,5));
for n=1:size(kernel,5)
    KERNEL(:,:,:,:,n) = zpad(conj(kernel(end:-1:1,end:-1:1,end:-1:1,:,n))*sqrt(imSize(1)*imSize(2)*imSize(3)), ...
        [imSize(1), imSize(2), imSize(3),size(kernel,4)]);
    for j=1:size(KERNEL,4)
    mytemp = KERNEL(:,:,:,j,n);
       for k = 1:1:3
             mytemp = fftc(mytemp,k);
       end
     KERNEL(:,:,:,j,n) = mytemp;
    end
end
 KERNEL = KERNEL/sqrt(prod(kSize));


EigenVecs = zeros(imSize(1), imSize(2), imSize(3), nc, min(nc,nv));
EigenVals = zeros(imSize(1), imSize(2), imSize(3), min(nc,nv));

for n=1:prod(imSize)
    [x,y,z] = ind2sub([imSize(1),imSize(2), imSize(3)],n);
    mtx = squeeze(KERNEL(x,y,z,:,:));

    %[C,D] = eig(mtx*mtx');
    [C,D,V] = svd(mtx,'econ');
    
    ph = repmat(exp(-1i*angle(C(1,:))),[size(C,1),1]);
    C = v*(C.*ph);
    D = real(diag(D));
    EigenVals(x,y,z,:) = D(end:-1:1);
    EigenVecs(x,y,z,:,:) = C(:,end:-1:1);
end




