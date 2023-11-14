function [phi_Cor,cMap] = wrls(phi,sigma,M,opt)

% Build A matrix
if size(M, 3) > 1
    [A,Af] = getMatA(M,opt.pOrd);
else
    [A,Af] = getMatA_2D(M,opt.pOrd);
end

% Weight data
tmp = sigma(M);

A_w = bsxfun(@times, A, 1./tmp);
phi_w = phi(M).*(1./tmp);
    
AHAm = A_w'*A_w;
AHA  = @(x) AHAm*x;
AHy = A_w'*phi_w;

cLS   = AHAm \ AHy; % Initial guess

[cL1,~,~,~] = fista(AHA,AHy,cLS,opt); % Solve for coefficients via FISTA
cMap = reshape(Af*cL1,size(M));     % fitted maps 

phi_Cor = phi - cMap;
    
end



