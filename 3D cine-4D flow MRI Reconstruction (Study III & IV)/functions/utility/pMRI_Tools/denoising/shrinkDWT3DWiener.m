function [ y] = shrinkDWT3DWiener( x, sizes, lam)
%SHRINKDWT3D Summary of this function goes here
%   Detailed explanation goes here
    
    x = gather(reshape(x,sizes));
    w = wavedec3(x,1,'db1', 'mode','per');

    k = 5;
    lambda(1) = 0.01*lam;
    lambda(2) = lam;
    lambda(3) = lam;
    lambda(4) = lam;
    lambda(5) = k*lam;
    lambda(6) = k*lam;
    lambda(7) = k*lam;
    lambda(8) = k*lam;
    for ind = 1:8
       w.dec{ind} = wiener3D(w.dec{ind},lambda(ind));
    end
    y = waverec3(w);
    y = gpuArray(y(:));
end
