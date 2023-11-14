function [u,hist] = core(x0,y,A,p)
%===========================================================================================%
% Inputs:
% x0: Initial Image
% y: Measured undersampled k-space data
% A: Sensing matrix
% p: parameters
% Outputs:
% u: reconstructed image
% hist: history of objective function
%===========================================================================================%
% Compressive Recovery with Outlier Rejection 'CORe' (ADMM/Split Bregman Implementation)
% Written by:
% Syed Murtaza Arshad (arshad.32@osu.edu)
% Rizwan Ahmad, PhD (ahmad.46@osu.edu)
%===========================================================================================%
% Extract parameters from p structure
mu1  = p.mu1_l2l1g;
mu2  = p.mu2_l2l1g;
lam1   = p.lam1_l2l1g;
lam2    =p.lam2_l2l1g;
oIter = p.oIter;
iIter = p.iIter;
gStp  = p.gStp;
vrb   = p.vrb;
readout = p.readout;

% Casting the data to single precision to reduce the data on GPU
y=single(y);
u = single(x0); % u represents the current estimated image
clear x0; % intial image not needed anymore

% Initializing outliers 'v', and auxiliary variables d and b as 0s
if p.use_gpu
v=zeros(size(y),"single","gpuArray");
d1 = zeros(A.frame_size(1),A.frame_size(2),A.frame_size(3),A.Q,16, ...
    "single","gpuArray");
else
v=zeros(size(y),"single");
d1 = zeros(A.frame_size(1),A.frame_size(2),A.frame_size(3),A.Q,16, ...
    "single"); 
end
b1=d1;
d2=v;
b2=v;


compute = A.compute; % How to compute the wavelet transform

% Wavelet transform operator
W = harr_nddwt_4D('db1',[A.frame_size(1),A.frame_size(2),A.frame_size(3), ...
    A.Q],'pres_l2_norm',1,'compute',compute,'precision','single');
 mStart = tic;  % Start the reconstruction timer 
%Walk through iterations
hist=zeros(oIter,1);
for i = 1:oIter
    tStart = tic;  % Start the iteration timer
    %Counter
    for j = 1:iIter
        % Gradient of fidelity term in objective function wrt u
        gradA = 2.*A.multTr(A.mult(u) + v - y);
        % Gradient of wavelet sparisty term in objective function
        gradW = mu1 * W.rec(W.dec(reshape(u,[A.frame_size(1), ...
            A.frame_size(2),A.frame_size(3),A.Q]),1) - d1 + b1);
         % Taking gradient descent step to estimate true image
        u = u - gStp*(gradA + gradW(:));
    end
    Au=A.mult(u);
    for j = 1:iIter
        % Gradient of fidelity term in objective function wrt v
        gradA = Au + v - y;
        % Reshaping v to enforce group (readout) sparsity
        v=reshape(v,readout,[]);
        % Gradient of wavelet sparisty term in objective function wrt v
        gradW = mu2 .* (reshape(v,[],1) - (d2-b2).* ...
            reshape(v./(sqrt(sum(abs(v).^2))+10e-6),[],1));
        % Taking gradient descent step to estimate outliers
        v = v(:) - gStp*(gradA + gradW);
    end
    clear gradA;
    clear gradW;

    % Storing 16-band wavelet transform to Wdecu 
    Wdecu = W.dec(reshape(u,[A.frame_size(1),A.frame_size(2), ...
        A.frame_size(3),A.Q]),1);

    % Updating auxiliary variables
    for ind = 1:16
        d1(:,:,:,:,ind) = shrink1(Wdecu(:,:,:,:,ind)+b1(:,:,:,:,ind), ...
            lam1(ind)/mu1, 1, 1e-6);
    end
    b1 = b1 + (Wdecu - d1);
    % Reshaping v to enforce group (readout) sparsity
    v=reshape(v,readout,[]);
    b2=reshape(b2,readout,[]);
    d2 = shrink1(sqrt(sum(abs(v).^2)) + b2, lam2/mu2, 1, 1e-6);
    b2 = b2 + (sqrt(sum(abs(v).^2)) - d2);
    b2=b2(:);
    d2=d2(:);
    v=v(:);

    % Displaying iteration information
    if rem(i, vrb)==0
         objA =  round(0.5*sum(abs(Au + v - y).^2));
         objW= lam1.*reshape(Wdecu,[],16);
         objW =round(sum(abs(objW(:))));
         objV = round(sum(lam2.*sqrt(sum(abs(reshape(v,readout,[])).^2))));
         obj=objA+objW+objV;
         hist(i)=obj;
         fprintf(['CORe: Iter = %s \tobjA= %s\tobjW= %s\tobjv= ' ...
             '%s\ttotal_obj= %s\ttime/iter = %s\n'],...
         num2str(i),num2str(objA),num2str(objW),num2str(objV), ...
         num2str(obj),num2str(toc(tStart)));
    
    clear objA;
    clear objW;
    clear objV;
    clear obj;
    end
    clear Au;
    clear Wdecu;
end
     disp("CORe Total Time (mins): "+round(toc(mStart)/60));
end
   

