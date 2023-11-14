function y = shrinkWaveSNIPE3D(x,sizes,W,lam)
    k = 5;
    lambda(1) = 0.01*lam;
    lambda(2) = lam;
    lambda(3) = lam;
    lambda(4) = lam;
    lambda(5) = k*lam;
    lambda(6) = k*lam;
    lambda(7) = k*lam;
    lambda(8) = k*lam;
    w = W.mult(x);
    w = reshape(w,[sizes,8]);
    y = zeros(size(w),class(w));
    for ind = 1:8
        y(:,:,:,ind) = snipeShrink(w(:,:,:,ind),lambda(ind));
    end
    y = W.multTr(y(:));
    y = y(:);
end

function w = snipeShrink(s,lambda)

    omega = 100;
    var = 1.4*lambda^2/omega;
    snipe_thresh = SNIPEstim(omega,'isCmplx',1);
    w = snipe_thresh.estim(s,var);

end