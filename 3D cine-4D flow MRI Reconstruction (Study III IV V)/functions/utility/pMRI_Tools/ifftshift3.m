function y = ifftshift3(x)
    y = ifftshift(ifftshift(ifftshift(x,1),2),3);
end