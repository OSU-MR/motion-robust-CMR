function w = shrink1(s, alph, p, ep)
t = abs(s);
w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
t(t == 0) = 1;
w = w./t;
end