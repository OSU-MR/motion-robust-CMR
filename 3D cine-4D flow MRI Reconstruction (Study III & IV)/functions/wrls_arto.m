function [phi_Cor,cMap] = wrls_arto(phi,sigma,M,opt)

pOrd_Orig = opt.pOrd;

[phi_Cor, useUC] = initializeWRLS_ARTO(phi, sigma, M, opt);

if useUC
    phi_Cor = phi;
    cMap = zeros(size(phi));
    return;
end

opt.pOrd = pOrd_Orig;
for k = 1 : opt.Kmax
    
    [S,mu,gamma,prob] = arto(phi_Cor,sigma,M,opt);
    
    opt.mu0 = mu; opt.sigma0 = gamma; opt.phi0 = prob;

    [phi_Cor,cMap] = wrls(phi,sigma,logical(S.*M),opt);

end



