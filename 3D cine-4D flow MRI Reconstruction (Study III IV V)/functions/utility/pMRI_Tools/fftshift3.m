function y = fftshift3(x)
    y = fftshift(fftshift(fftshift(x,1),2),3);
end