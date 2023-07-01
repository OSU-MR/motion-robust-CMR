function [A,Af] = getMatA(M,pOrd)

% Indices to fit
ind = find(M);
[Y,X,Z] = ind2sub([size(M,1),size(M,2),size(M,3)],ind);
indf = find(ones(size(M)));
[Yf,Xf,Zf] = ind2sub([size(M,1),size(M,2),size(M,3)],indf);

% Center  so that (x,y) = (0,0) occurs at index (Lx/2+1,Ly/2+1)
% Over mask indices
X = X - (floor(size(M,2)/2)+1);
Y = Y - (floor(size(M,1)/2)+1);
Z = Z - (floor(size(M,3)/2)+1);
Xf = Xf - (floor(size(M,2)/2)+1);
Yf = Yf - (floor(size(M,1)/2)+1);
Zf = Zf - (floor(size(M,3)/2)+1);

% Transfer matrix (Ax=y)
A = polyDef(X,Y,Z,pOrd);
Af = polyDef(Xf,Yf,Zf,pOrd);

% Normalize column norms
colNor    = diag(sqrt(Af'*Af))';
A = bsxfun(@rdivide, A, colNor); % Normalize the norm of each column of A
Af = bsxfun(@rdivide, Af, colNor); % Normalize the norm of each column of Af

end

