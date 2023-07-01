function w = shrinkBM3D(z,sizes,sigma)
%     z = reshape(z,sizes);
    z = real(reshape(z,sizes));
    [~,w] = BM3D(1, gather(z), sigma);
    w = gpuArray(w(:));
end