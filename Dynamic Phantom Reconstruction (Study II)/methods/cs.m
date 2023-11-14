%% cs- Robust Regression Algorithm
function x = cs(y,p)
%===========================================================================================%
% Compressed Sensing 'CS' (ADMM/Split Bregman Implementation)
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
mu  = p.mu_cs;
lam   = p.lam_cs;
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

x  = At(y); % Estimating intial image
% Initializing d and b auxiliary variables as 0s
d = zeros(M);
b = zeros(size(d));
%Walk through iterations
tStart = tic;  % Start recon timer
for i = 1:oIter
    iStart = tic;  % Start iter timer
    for j = 1:iIter 
        gradA =  2.*At(A(x) - y); % Gradient of fidelity term in objective function
        gradW = mu * W.rec(W.dec(x,1) - d + b);  % Gradient of wavelet sparisty term in objective function
        x = x - gStp*(gradA + gradW); % Taking gradient descent step to estimate true image
    end
   
    d = sth2(W.dec(x,1) + b, lam/mu);  % Updating auxiliary variables

    b = b + (W.dec(x,1) - d);  % Updating auxiliary variables
    % Displaying iteration information
    if rem(i, vrb)==0
        objA = sum(abs(A(x)-y).^2,'all');
        objW = sum(abs(W.dec(x,1) .*permute(lam,[3,1,2])),'all');
        fprintf('Iter = %s \tobjTOT= %s \tobjA= %s \tobjW= %s\ttime/iter = %s\n',...
        num2str(i), num2str(objA+objW,5), num2str(objA,5), num2str(objW,5),num2str(toc(iStart),2));
    end
    
end
disp("Total CS reconstruction time="+toc(tStart));
end