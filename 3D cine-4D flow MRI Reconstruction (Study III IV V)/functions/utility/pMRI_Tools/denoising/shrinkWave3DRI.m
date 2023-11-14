function y = shrinkWave3DRI(x,sizes,W,lam)
    k = 5;
    lambda(1) = 0.01*lam;
    lambda(2) = lam;
    lambda(3) = lam;
    lambda(4) = lam;
    lambda(5) = k*lam;
    lambda(6) = k*lam;
    lambda(7) = k*lam;
    lambda(8) = k*lam;
%     x = real(x);
    w = W.mult(x);
    w = reshape(w,[sizes,8]);
    y = zeros(size(w),class(w));
    for ind = 1:8
        y(:,:,:,ind) = shrink1(real(w(:,:,:,ind)),lambda(ind),1,1e-6) +...
                       1j*shrink1(imag(w(:,:,:,ind)),lambda(ind),1,1e-6);
    end
    y = W.multTr(y(:));
    y = y(:);
end

function w = shrink1(s, alph, p, ep)
    t = abs(s);
    w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
    t(t == 0) = 1;
    w = w./t;
end