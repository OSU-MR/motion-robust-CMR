function y =  sth1(x, lam)

xabs = abs(x);
y = max(xabs - lam, 0).* x;
xabs(xabs == 0) = 1;
y = y./xabs;
end