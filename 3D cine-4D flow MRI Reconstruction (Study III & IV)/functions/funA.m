function [y] =funA(x,sampInd,S, nSize)
% works for 2D,3D,single and mutilple sensitivity maps
% image to k-space
% sensitivity maps
% 
%2D to pseudo 3D
if numel(nSize) == 4
    x = permute(x,[1,2,5,3,4]);% kx,ky,kz,frame,set (kz = 1)
    S = permute(S,[1,2,5,3,4]);%kx,ky,kz,coil,set (kz= 1)
end

y = sum(bsxfun(@times, S, permute(x,[1,2,3,6,5,4])),5);%kx,ky,kz,coil,set,frame; sum set dim

for n = 1:3
  y = fftc(y,n);
end
y = y(sampInd);
