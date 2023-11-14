%% SO - Sparse Outliers Method
function [u,v] = so(y,p)
%===========================================================================================%
% Sparse Outliers - 'SO' (ADMM/Split Bregman Implementation)
% Written by:
% Syed Murtaza Arshad (arshad.32@osu.edu)
% Rizwan Ahmad, PhD (ahmad.46@osu.edu)
%===========================================================================================%
% Inputs:
% y: Measured undersampled k-space data
% p: parameters
% Output:
% x: reconstructed image
% yo: rejected outliers
%===========================================================================================%
% Extract parameters from p structure
mu  = p.mu_so;
lam1   = p.lam1_so;
lam2    =p.lam2_so;
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
v  =zeros(size(y));   % Initializing outliers as 0s
 % Initializing d and b auxiliary variables as 0s
d = zeros(M);
b = zeros(size(d));

tStart = tic;  % Start recon timer

%Walk through iterations    
for i = 1:oIter
    iStart = tic;  % Start iter timer
    for j = 1:iIter
        gradA =  2.*At(A(u) + v - y); % Gradient of fidelity term in objective function
        gradW = mu * W.rec(W.dec(u,1) - d + b);  % Gradient of wavelet sparisty term in objective function
        u = u - gStp*(gradA + gradW);  % Taking gradient descent step to estimate true image
    end
   
    v = sth1(y-A(u),lam2); % Upldating outliers
    % Updating auxiliary variables
    d = sth2(W.dec(u,1) + b, lam1/mu);
    b = b + (W.dec(u,1) - d);
    % Displaying iteration information
    if rem(i, vrb)==0
        objA = sum(abs(A(u)+v-y).^2,'all');
        objW = sum(abs(W.dec(u,1) .*permute(lam1,[3,1,2])),'all');
        objV = sum(abs(v),'all').*lam2;

        fprintf('Iter = %s \tobjTOT= %s \tobjA= %s \tobjW= %s \tobjV= %s\ttime/iter = %s\n',...
        num2str(i), num2str(objA+objW+objV,5), num2str(objA,5), num2str(objW,5), num2str(objV,5) ...
        ,num2str(toc(iStart),2));
    end   
end
 disp("Total SO reconstruction time="+toc(tStart));
end