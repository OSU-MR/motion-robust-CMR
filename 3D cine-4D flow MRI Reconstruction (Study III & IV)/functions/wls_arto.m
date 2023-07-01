function [phi_Cor,cMap] = wls_arto(phi,sigma,M,opt)

pOrd_Orig = opt.pOrd;

% Initialize exclusion operator, S
S0 = zeros(size(M));

switch opt.acqParam.InPlanePhaseEncodingDirection
    case 'ROW' % Phase encode across rows in image
        mid = floor(size(M,2)/2)+1; % Center
        extent = floor((size(M,2)*opt.midFOVFrac)/2);
        S0(:,mid-extent:mid+extent) = 1;
        S0 = logical(S0);
    case 'COL' % Phase encode across columns in image
        mid = floor(size(M,1)/2)+1; % Center
        extent = floor((size(M,1)*opt.midFOVFrac)/2);
        S0(mid-extent:mid+extent,:) = 1;
        S0 = logical(S0);
end

opt.pOrd = 1;
[phi_Cor,cMap] = wls(phi,sigma,logical(S0.*M),opt);

opt.pOrd = pOrd_Orig;
for k = 1 : opt.Kmax
    
    [S,mu,gamma,prob] = arto(phi_Cor,sigma,M,opt);
    
    opt.mu0 = mu; opt.sigma0 = gamma; opt.phi0 = prob;

    [phi_Cor,cMap] = wls(phi,sigma,logical(S.*M),opt);

end

