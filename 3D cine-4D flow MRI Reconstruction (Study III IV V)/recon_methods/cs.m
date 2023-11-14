function [u,hist] = cs(x0,y,A,p)
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
% Compressed Sensing 'CS' (ADMM/Split Bregman Implementation)
% Written by:
% Syed Murtaza Arshad (arshad.32@osu.edu)
% Rizwan Ahmad, PhD (ahmad.46@osu.edu)
%===========================================================================================%

% Extract parameters from p structure
mu    = p.mu_l2l1;
lam   = p.lam_l2l1;
oIter = p.oIter;
iIter = p.iIter;
gStp  = p.gStp;
vrb   = p.vrb;

% Casting the data to single precision to reduce the data on GPU
y=single(y);
u = single(x0); % u represents the current estimated image
clear x0; % intial image not needed anymore

% Initializing d and b auxiliary variables as 0s
if p.use_gpu
    d = zeros(A.frame_size(1),A.frame_size(2),A.frame_size(3),A.Q, ...
        16,"single","gpuArray");
else
    d = zeros(A.frame_size(1),A.frame_size(2),A.frame_size(3),A.Q, ...
        16,"single"); 
end

b=d;

compute = A.compute;  % How to compute the wavelet transform

% Wavelet transform operator
W = harr_nddwt_4D('db1',[A.frame_size(1),A.frame_size(2),A.frame_size(3), ...
    A.Q],'pres_l2_norm',1,'compute',compute,'precision','single');

mStart = tic; % Start the reconstruction timer 
%Walk through iterations
hist=zeros(oIter,1);
for i = 1:oIter
    tStart = tic;   % Start the iteration timer 
    for j = 1:iIter
        % Gradient of fidelity term in objective function
        gradA =  2.*A.multTr(A.mult(u) - y);
        % Gradient of wavelet sparisty term in objective function
        gradW = mu * W.rec(W.dec(reshape(u,[A.frame_size(1), ...
            A.frame_size(2),A.frame_size(3),A.Q]),1) - d + b);
        % Taking gradient descent step to estimate true image
        u = u - gStp*(gradA + gradW(:));
    end
    clear gradA;
    clear gradW;
        % Storing 16-band wavelet transform to Wdecu 
    Wdecu = W.dec(reshape(u,[A.frame_size(1),A.frame_size(2), ...
        A.frame_size(3),A.Q]),1);
    % Updating auxiliary variables
    for ind = 1:16
        d(:,:,:,:,ind) = shrink1(Wdecu(:,:,:,:,ind)+b(:,:,:,:,ind), ...
            lam(ind)/mu, 1, 1e-6);
    end
    b = b + (Wdecu - d);

    u=u(:);

    % Displaying iteration information
    if rem(i, vrb)==0
         objA = round(sum(abs(A.mult(u) - y).^2));
         objW= lam.*reshape(Wdecu,[],16);
         objW =round(sum(abs(objW(:))));
         obj=objA+objW;
         hist(i)=obj;
         fprintf(['CS: Iter = %s \tobjA= %s\tobjW= %s\ttotal_obj=' ...
             ' %s\ttime/iter = %s\n'],...
         num2str(i),num2str(objA),num2str(objW), num2str(obj), ...
         num2str(toc(tStart)));
    
    clear objA;
    clear objW;
    clear obj;

    end
        clear Wdecu;
end
    disp("CS Total Time (mins): "+round(toc(mStart)/60));
end
    
     

   
    

            
       

