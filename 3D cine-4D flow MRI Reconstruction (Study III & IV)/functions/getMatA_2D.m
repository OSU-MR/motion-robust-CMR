function [A,Af] = getMatA_2D(M,pOrd)

% Indices to fit
ind = find(M);
[Y,X] = ind2sub([size(M,1),size(M,2)],ind);
indf = find(ones(size(M)));
[Yf,Xf] = ind2sub([size(M,1),size(M,2)],indf);

% Center  so that (x,y) = (0,0) occurs at index (Lx/2+1,Ly/2+1)
% Over mask indices
X = X - (floor(size(M,2)/2)+1);
Y = Y - (floor(size(M,1)/2)+1);
Xf = Xf - (floor(size(M,2)/2)+1);
Yf = Yf - (floor(size(M,1)/2)+1);

% Transfer matrix (Ax=y)
A = polyDef(X,Y,pOrd);
Af = polyDef(Xf,Yf,pOrd);

% Normalize column norms
colNor    = diag(sqrt(Af'*Af))';
A = bsxfun(@rdivide, A, colNor); % Normalize the norm of each column of A
Af = bsxfun(@rdivide, Af, colNor); % Normalize the norm of each column of Af

end

