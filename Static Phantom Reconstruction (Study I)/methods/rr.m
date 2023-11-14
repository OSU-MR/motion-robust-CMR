%% rr- Robust Regression Algorithm
function x = rr(y,p)
%===========================================================================================%
% Robust Regression - 'RR' (ADMM/Split Bregman Implementation)
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
    mu1   = p.mu1_rr;
    mu2   = p.mu2_rr;
    lam   = p.lam_rr;
    oIter = p.oIter_rr;
    iIter = p.iIter_rr;
    gStp  = p.gStp_rr;
    vrb   = p.vrb;
    M     = p.M;
    N     = p.N;
    W     = p.W;
    A     = p.A;
    At    = p.At;
    
    x  = At(y); % Estimating intial image
    % Initializing d and b auxiliary variables as 0s
    d1 = zeros(size(y));
    d2 = zeros(M);
    b1 = zeros(size(d1));
    b2 = zeros(size(d2));
    %Walk through iterations    
    tStart = tic;  % Start recon timer
    for i = 1:oIter
        iStart = tic;  % Start iter timer
        for j = 1:iIter
            gradA = mu1 * At(A(x) - y - d1 + b1); % Gradient of fidelity term in objective function
            gradW = mu2 * W.rec(W.dec(x,1) - d2 + b2); % Gradient of wavelet sparisty term in objective function
            x = x - gStp*(gradA + gradW); % Taking gradient descent step to estimate true image
        end
        d1 = sth1(A(x) - y + b1, 1/mu1); % Updating auxiliary variables
        d2 = sth2(W.dec(x,1) + b2, lam/mu2); % Updating auxiliary variables
        b1 = b1 + (A(x) - y - d1); % Updating auxiliary variables
        b2 = b2 + (W.dec(x,1) - d2); % Updating auxiliary variables
        % Displaying iteration information
        if rem(i, vrb)==0
            objA = sum(abs(A(x)-y),'all');
            objW = sum(abs(W.dec(x,1) .* permute(lam,[3,1,2])),'all');
            fprintf('Iter = %s \tobjTOT= %s \tobjA= %s \tobjW= %s\ttime/iter = %s\n',...
                    num2str(i), num2str(objA+objW,5), num2str(objA,5), num2str(objW,5),num2str(toc(iStart),2));
        end
        
    end
    disp("Total RR Reconstruction Time: "+toc(tStart));
end