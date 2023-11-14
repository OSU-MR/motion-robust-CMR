%% rr- Robust Regression Algorithm
function u = rr(y,p)
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
% u: reconstructed image
%===========================================================================================%
% Extract parameters from p structure
mu1   = p.mu1_rr;
mu2   = p.mu2_rr;
lam   = p.lam_rr;
oIter = p.oIter;
iIter = p.iIter;
gStp  = p.gStp;
vrb   = p.vrb;
M     = p.M;
N     = p.N;
W     = p.W;
A     = p.A;
At    = p.At;

xo=p.x;



u  = At(y); % Estimating intial image
% Initializing d and b auxiliary variables as 0s
d1 = zeros(size(y)); 
d2 = zeros(M);
b1 = zeros(size(d1));
b2 = zeros(size(d2));
tStart = tic;  %Start recon timer
for i = 1:oIter
    iStart = tic;  %Start iter timer
    for j = 1:iIter
        gradA = mu1 * At(A(u) - y - d1 + b1); % Gradient of fidelity term in objective function
        gradW = mu2 * W.rec(W.dec(u,1) - d2 + b2);  % Gradient of wavelet sparisty term in objective function
        u = u - gStp*(gradA + gradW); % Taking gradient descent step to estimate true image
    end
     % Updating auxiliary variables
    d1 = sth1(A(u) - y + b1, 1/mu1);
    d2 = sth2(W.dec(u,1) + b2, lam/mu2);
    b1 = b1 + (A(u) - y - d1);
    b2 = b2 + (W.dec(u,1) - d2);
    % Displaying iteration information
    if rem(i, vrb)==0
        objA = sum(abs(A(u)-y),'all');
        objW = sum(abs(W.dec(u,1) .*permute(lam,[3,1,2])),'all');
        fprintf('Iter = %s \tobjTOT= %s \tobjA= %s \tobjW= %s\ttime/iter = %s\n',...
        num2str(i), num2str(objA+objW,5), num2str(objA,5), num2str(objW,5),num2str(toc(iStart),2));
    end
       
end
 disp("Total RR reconstruction time="+toc(tStart));
end