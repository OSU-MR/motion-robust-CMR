function lip=lipconst(p)
n=p.N;
A=p.A;
At=p.At;
c = randn(n);
c = c / norm(c);
thresh = 1e-3;
err = inf;
uest = inf;
while err > thresh
    c = At(A(c));
    
    %Progress
    unew = norm(c);
    err = abs(log(uest / unew));
    
    %Save the estimate
    uest = unew;
    c = c / norm(c);
end

lip = 2.05*uest;

end