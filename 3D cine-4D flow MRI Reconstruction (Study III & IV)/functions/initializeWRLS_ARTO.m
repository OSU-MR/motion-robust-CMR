function [phi_Cor, useUC] = initializeWRLS_ARTO(phi, sigma, M, opt)

opt.pOrd = 1;

useUC = 0;

opt.selectInit = 0;

if opt.selectInit
    % center of the FOV
    opt.initialMaskLoc = "mid";
    S0_mid = getInitialMask(M, opt);
    [phi_Cor_mid,~] = wrls(phi,sigma,logical(S0_mid.*M),opt);
    % one edge of FOV
    opt.initialMaskLoc = "sideLo";
    S0_Lo = getInitialMask(M, opt);
    [phi_Cor_Lo,~] = wrls(phi,sigma,logical(S0_Lo.*M),opt);
    % other edge of FOV
    opt.initialMaskLoc = "sideHi";
    S0_Hi = getInitialMask(M, opt);
    [phi_Cor_Hi,cMap_Hi] = wrls(phi,sigma,logical(S0_Hi.*M),opt);
    
    % Compare weighted residuals after the three initializations
    e_mid = phi_Cor_mid ./ sigma; 
    e_Lo = phi_Cor_Lo ./ sigma;
    e_Hi = phi_Cor_Hi ./ sigma;
    
    % residual grid
    tmp_Lo = logical(S0_Lo.*M);
    tmp_mid = logical(S0_mid.*M);
    tmp_Hi = logical(S0_Hi.*M);
    I(1,1) = norm(e_Lo(tmp_Lo),2)^2;
    I(1,2) = norm(e_Lo(tmp_mid),2)^2;
    I(1,3) = norm(e_Lo(tmp_Hi),2)^2;
    I(2,1) = norm(e_mid(tmp_Lo),2)^2;
    I(2,2) = norm(e_mid(tmp_mid),2)^2;
    I(2,3) = norm(e_mid(tmp_Hi),2)^2;
    I(3,1) = norm(e_Hi(tmp_Lo),2)^2;
    I(3,2) = norm(e_Hi(tmp_mid),2)^2;
    I(3,3) = norm(e_Hi(tmp_Hi),2)^2;
    
    [V, D] = eig(I);
    A = diag(abs(V*D));
    ind = find(A == min(A));
    switch ind
        case 1
            S = S0_Lo;
            phi_Cor = phi_Cor_Lo;
        case 2
            S = S0_mid;
            phi_Cor = phi_Cor_mid;
        case 3
            S = S0_Hi;
            phi_Cor = phi_Cor_Hi;
    end
    
    convKernal = ones(10,10)*(1/(10*10));
    convPhase = conv2(phi, convKernal, 'same');
    
    
    e_Cor = phi_Cor ./ sigma;
    e_UC = phi ./ sigma;
    e_norm = norm(e_Cor(logical(S.*M)),2)^2;
    UC_norm = norm(e_UC(logical(S.*M)),2)^2;
    
    if UC_norm <= e_norm
        useUC = 1;
        phi_Cor = zeros(size(M));
        return;
    end
    
else
    % center of the FOV
    opt.initialMaskLoc = "mid";
    S0_mid = getInitialMask(M, opt);
    [phi_Cor_mid,~] = wrls(phi,sigma,logical(S0_mid.*M),opt);
    
    S0 = S0_mid;
    phi_Cor = phi_Cor_mid;    
    
end

