%% CORe - Compressive Recovery with Outlier Rejection
function [x,yo] = core(y,p)
%===========================================================================================%
% Compressive Recovery with Outlier Rejection - 'CORe' (ADMM/Split Bregman Implementation)
% Written by:
% Syed Murtaza Arshad (arshad.32@osu.edu)
% Rizwan Ahmad, PhD (ahmad.46@osu.edu)
%===========================================================================================%
% Inputs:
% y: Measured undersampled k-space data
% p: parameters
% Output:
% x: reconstructed image
%===========================================================================================%
% Extract parameters from p structure
    mu1   = p.mu1_core;
    mu2   = p.mu2_core;
    lam1  = p.lam1_core;
    lam2 = p.lam2_core;
    oIter = p.oIter_core;
    iIter = p.iIter_core;
    gStp1 = p.gStp1_core;
    gStp2 = p.gStp2_core;
    vrb   = p.vrb;
    M     = p.M;
    N     = p.N;
    W     = p.W;
    A     = p.A;
    At    = p.At;
    s     = p.s;
    
    x  = At(y); % Estimating intial image
    xmax = max(abs(x(:)));
    ymax = max(abs(y(:)));
    yo = ymax/1e3*(randn(size(y)) + 1j*randn(size(y))).*s;  % Initializing random outliers
   % Initializing v1,2 and u1,2 auxiliary variables as 0s
    v1 = xmax/10*randn(M);
    u1 = xmax/10*randn(size(v1));
    v2 = ymax/1e3*randn(1,N(2)).*s;
    u2 = ymax/1e3*randn(size(v2)).*s;
    %Walk through iterations    
    tStart = tic;  % Start iter timer
    for i = 1:oIter 
        iStart = tic;  % Start iter timer
        for j = 1:iIter
            gradA =  At(A(x) + yo - y); % Gradient of fidelity term in objective function
            gradW = mu1 * W.rec(W.dec(x,1) - v1 + u1); % Gradient of wavelet sparisty term in objective function
            x = x - gStp1*(gradA + gradW); % Taking gradient descent step to estimate true image
        end

        for j = 1:iIter
            gradA =  A(x) + yo - y;
            gradG = mu2 * (yo - bsxfun(@times, (v2-u2), bsxfun(@rdivide, yo, sqrt(sum(abs(yo).^2))+ymax*1e-6)));
            yo = yo - gStp2*(gradA + gradG); % Upldating outliers
        end
% Updating auxiliary variables
        v1 = sth2(W.dec(x,1) + u1, lam1/mu1);
        u1 = u1 + (W.dec(x,1) - v1);

        v2 = sth1(sqrt(sum(abs(yo).^2)) + u2, lam2/mu2);
        u2 = u2 + (sqrt(sum(abs(yo).^2)) - v2);

        % Displaying iteration information
        if rem(i, vrb)==0
            objA = sum(sum(abs(A(x)-y+yo).^2));
            objW = sum(sum(sum(abs(W.dec(x,1) .* permute(lam1,[3,1,2])))));
            objG = sum(lam2*sqrt(sum(abs(yo).^2)));
            fprintf('Iter = %s \tobjTOT= %s \tobjA= %s \tobjW= %s \tobjG= %s \ttime/iter = %s\n',...
                    num2str(i), num2str(objA+objW+objG,5), num2str(objA,5), num2str(objW,5), num2str(objG,5), num2str(toc(iStart),2));
        end
        
    end
    disp("Total CORe Reconstruction Time: "+toc(tStart));
end
