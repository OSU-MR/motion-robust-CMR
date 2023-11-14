function y = shrinkLowRank(x,sizes,W,lam)
%SHRINKLOWRANK Summary of this function goes here
%   Detailed explanation goes here

    x = reshape(x,[sizes(1)*sizes(2),sizes(3)]);
    [U,S,V] = svd(x);
    display('blah')
end

function w = shrink1(s, alph, p, ep)
    t = abs(s);
    w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
    t(t == 0) = 1;
    w = w./t;
end