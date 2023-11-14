%% CORe - Compressive Recovery with Outlier Rejection
function [u,v] = core(y,p)
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
% u: reconstructed image
% v: rejected outliers
%===========================================================================================%
% Extract parameters from p structure
    mu1   = p.mu1_core;
    mu2   = p.mu2_core;
    lam1  = p.lam1_core;
    lam2 = p.lam2_core;
    oIter = p.oIter;
    iIter = p.iIter;
    gStp1 = p.gStp;
    gStp2 = p.gStp;
    vrb   = p.vrb;
    M     = p.M;
    N     = p.N;
    W     = p.W;
    A     = p.A;
    At    = p.At;

    xo=p.x;
    u  = At(y); % Estimating intial image
    v= zeros(size(y)); % Initializing outliers as 0s
    % Initializing v1,2 and u1,2 auxiliary variables as 0s
    d1 =  zeros(M);
    b1 = zeros(size(d1));
    d2 = zeros(size(y));
    b2 = zeros(size(d2));
    %Walk through iterations   
    tStart = tic;   % Start recon timer
    for i = 1:oIter
      iStart = tic;  % Start iter timer
        for j = 1:iIter
            gradA =  2.*At(A(u) + v - y);  % Gradient of fidelity term in objective function wrt u
            gradW = mu1 * W.rec(W.dec(u,1) - d1 + b1); % Gradient of wavelet sparisty term in objective function wrt x
            u = u - gStp1*(gradA + gradW); % Taking gradient descent step to estimate true image
        end

        for j = 1:iIter
            gradA =  A(u) + v - y;  % Gradient of fidelity term in objective function wrt v
            v=reshape(v,N(1),[]);
            gradG = mu2 .* (reshape(v,[],1) - (d2-b2).*reshape(v./(sqrt(sum(abs(v).^2))+10e-6),[],1)); % Gradient of wavelet sparisty term in objective function wrt v
            v=v(:); % Updating outliers
            v = v - gStp2*(gradA + gradG);
        end
% Updating auxiliary variables
        d1 = sth2(W.dec(u,1) + b1, lam1/mu1);
        b1 = b1 + (W.dec(u,1) - d1);
        v=reshape(v,N(1),[]);
        b2=reshape(b2,N(1),[]);
        d2 = sth1(sqrt(sum(abs(v).^2)) + b2, lam2/mu2);
        b2 = b2 + (sqrt(sum(abs(v).^2)) - d2);
       
        b2=b2(:);
        d2=d2(:);

% Displaying iteration information
        if rem(i, vrb)==0
            objA = sum(abs(A(u)+v(:)-y).^2,'all');
            objW = sum(abs(W.dec(u,1) .*permute(lam1,[3,1,2])),'all');
            objG = sum(sum(abs(v).^2,1),2).*lam2;
            fprintf('Iter = %s \tobjTOT= %s \tobjA= %s \tobjW= %s \tobjG= %s \ttime/iter = %s\n',...
                    num2str(i), num2str(objA+objW+objG,5), num2str(objA,5), num2str(objW,5), num2str(objG,5), num2str(toc(iStart),2));
        end
        v=v(:);
        
    end
    disp("Total CORe reconstruction time="+toc(tStart));
end
