function y = shrinkWave2D(x,sizes,W,lam)
    lambda(1) = 0.01*lam;
    lambda(2) = lam;
    lambda(3) = lam;
    lambda(4) = lam;
%     x = reshape(x,sizes);
    x = real(reshape(x,sizes));
    w = W.dec(x,1);
    y = zeros(size(w),class(w));
    for ind = 1:4
        y(:,:,ind) = shrink1(w(:,:,ind),lambda(ind),1,1e-6);
    end
    y = W.rec(y);
    y = y(:);
end

function w = shrink1(s, alph, p, ep)
    t = abs(s);
    w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
    t(t == 0) = 1;
    w = w./t;
end