function w = shrinkMEDIAN3D(z,sizes)
    z = gather(reshape(z,sizes));
    wreal = medfilt3(real(z));
    wimag = medfilt3(imag(z));
    w = wreal(:) + 1j*wimag(:);
end