function [xhat,maps,res] = FISTApMRI(kdata,samp,lambda,sparseTrans,varargin)
% FISTApMRI estimates and image from undersampled k-space using the FISTA
%   and the sparsifying transform specified by the user
%
%   Inputs:
%       kdata: Undersampled kspace data
%       samp: k space sampling pattern
%       lambda: regularization paramter
%       spartTrans: sparsifying transfrom object
%               
%   Optional Inputs:
%       smaps: user defined sensitivity maps
%       x0: user defined initialization
%               

uniform_var = 0;
precision = 'single';
compute = 'mat';

% force the data to be downsampled
kdata = bsxfun(@times,kdata,permute(samp,[1,2,4,3]));

%% Estimate sensitivity maps
weights_b = repmat(sum(samp,3),[1,1,size(kdata,3)]);
weights_b(find(weights_b==0)) = Inf;

% Normalize k_data by weights and take an inverse Fourier Transform
time_av_b = ifft2_shift(sum(kdata,4)./weights_b);
[x0, maps] = WalshCoilCombine(time_av_b,3);
maps = repmat(maps,[1,1,1,size(samp,3)]);
x0 = repmat(x0,[1,1,size(samp,3)]);

if mod(length(varargin),2)==0
    for ind = 1:2:length(varargin)
        switch lower(varargin{ind})
            case('smaps')
                maps = varargin{ind+1};
                disp('using user defined sensitivity maps');
            case('x0')
                x0 = varargin{ind+1};
                disp('using user defined initialization');
        end
    end
else
    error('Optional Inputs must come in pairs')
end

%% fftshift and downsample data
kdata = fftshift(fftshift(kdata,1),2);
samp = fftshift(fftshift(samp,1),2);
x0 = fftshift(fftshift(x0,1),2);
maps = fftshift(fftshift(maps,1),2);
kdata = downsample_data(kdata,samp);

% Normalize The columns of A to be unit norm
R = numel(samp)/length(find(samp ==1));
maps = maps*sqrt(R);
kdata = kdata*sqrt(R);

use_gpu = 1;
if use_gpu
   kdata = gpuArray(single(kdata));
   samp = gpuArray(single(samp));
   x0 = gpuArray(single(x0));
   maps = gpuArray(single(maps));
   compute = 'gpu';
   precision = 'single';
end

%% Create Operators =======================================================
pMRI = pMRI_Op_2D_t(maps,samp,'uniform_var',uniform_var,'precision',precision,'compute',compute);

options = FistaOpt();
options.lip = [];
options.nit = 50;
options.lam = lambda;
warning('remove change lambda back')
% global_lambda = 0.0007;
% global_lambda = 4e-06;
% options.lam = global_lambda;
options.tol = 2e-6;
options.verbose = 1;

% Run FISTA
[xhat,~,res,options] = fista_pmri(x0(:),kdata,pMRI,sparseTrans,options);


% for ind = 1:10
%     % learn the relative lambda values
%     redu_size = sum(sparseTrans.redu_factor);
%     wave = sparseTrans.mult(xhat);
%     wave = reshape(wave,[size(samp),redu_size]);
%     for ii = 1:redu_size
%         tmp = wave(:,:,:,ii);
% %         lam(ii) = 1/(mean(abs(tmp(:)))+max(abs(tmp(:)))/10);
%         lam(ii) = 1/(mean(abs(tmp(:))));
%     end
%     
%     lambda_pca = zeros(size(wave(:,:,:,1)));
%     for indpca = 1:size(wave,3)
%         tmp = squeeze(wave(:,:,indpca,1));
%         lambda_pca(:,:,indpca) = gather(1/(sum(abs(tmp(:)))+max(abs(tmp(:)))/1000));
%     end
% %     lam = global_lambda*lam/lam(1);
%     lam = global_lambda*lam;
%     fprintf('%s  ',num2str(lam))
%     fprintf('\n')
%     lambda_learn = zeros([size(samp),redu_size]);
%     for ii = 1:redu_size
%         lambda_learn(:,:,:,ii) = gather(lam(ii));
%     end
%     lambda_learn(:,:,:,1) = lambda_pca*global_lambda;
%     lambda_learn = lambda_learn(:);
% 
%     options.nit = 10;
%     options.lam = lambda_learn;
%     
%     [xhat,~,res] = fista_pmri(xhat(:),kdata,pMRI,sparseTrans,options);
% end

xhat = gather(xhat);
xhat = reshape(xhat,size(samp));
xhat = fftshift(fftshift(xhat,1),2);

end

function [x,estHist,res,options] = fista_pmri(x0,b,A,sparseTrans,options)
%FISTA Fast Iterative Soft Thresholding
%   Code to implement FISTA as described in "A Fast Iterative
%   Shrinkage-Thresholding Algorithm for Linear Inverse Problems" by Beck
%   and Teboulle. Implements the FISTA algorithm to solve the convex
%   program argmin_x ||Ax - b||_2^2 + |Lam .* x|. The implementation
%   support complex valued data.
%
%   Code inputs:
%   -x0 is the initial guess at the solution
%   -b is the measured data, i.e. b = Ax + w.
%   -A either a matrix or a linear operator defined by the LinTrans class.
%   -options a set of options of the FistaOpt class
%
%   Code Outputs:
%   -x is the estimate of the cost function's argmin.
%   -estHist history of the algorithm progress
%       estHist.xhat each iteration for the estimate


%% Input checking
t_start = tic;
% Get options
if (nargin < 5)
    options = FistaOpt();
elseif (isempty(options))
    options = FistaOpt();
end

%Check to see if history desired
saveHist = nargout > 1;


% If A is a double matrix, replace by an operator
if isa(A, 'double')
    A = MatrixLinTrans(A);
end

%Use the power iteration to estimate the Lipschitz constant if not provided
if isempty(options.lip)
    tstart= tic;
    fprintf('Calculating Lipschitz Constant...')
    
    %Initial guess
    q = randn(length(x0),1);
    q = q / norm(q);
    
    thresh = 1e-3;
    err = inf;
    uest = inf;
    while err > thresh
        q = A.multTr(A.mult(q));
        
        %Progress
        unew = norm(q);
        err = abs(log(uest / unew));
        
        %Save the estimate
        uest = unew;
        q = q / norm(q);
    end
    
    %Use a little more than twice the eigenvalue estimate
    options.lip = 2.05*uest;
    
    fprintf('  Completed in %0.2f seconds\n',toc(tstart))
end


%Convert options.lam to a vector
if length(options.lam) == 1
%     Lam = options.lam * ones(size(x0));
    Lam = options.lam;
else
    Lam = options.lam;
end


%% Iteration
%Initialize variables
z = x0;
x = x0;
t = 1;
stop = 0;
iter = 0;

%Preallocate storage if estimHist is requested
if saveHist
    estHist.xhat = zeros(length(x0),options.nit,class(x0));
end

%Walk through iterations
while ~stop
    tStart = tic; 
    %Counter
    iter = iter + 1;
    
    %Gradient step
    alphares = z - (1 / options.lip)*(2*A.multTr(A.mult(z) - b));
    alphares = sparseTrans.mult(reshape(alphares,[A.frame_size(1),A.frame_size(2),A.Q]));
    alphares = shrink1(alphares, options.lam/options.lip, 1, 1e-6);
    
    xnew = sparseTrans.multTr(alphares);
    xnew = xnew(:);
    
    %Compute the new t and y
    tnew = (1 + sqrt(1 + 4*t^2))/2;
    z = xnew + (t - 1)/tnew*(xnew - x);

    %Check max iterations
    if iter > options.nit
        stop = 1;
    else 
        %Check error tolerance
        cnew = norm(x - xnew) / norm(xnew);
        if cnew < options.tol && options.tol > 0
            stop = 1;
        end
    end
    
    %Record history if desired
    if saveHist
        estHist.xhat(:,iter) = xnew;
    end
    
    %Update on progress
    if options.verbose
        % Calculate the objective Function every 5 Iterations
        if mod(iter,5) ==0 || iter == 1
            obj = A.mult(xnew) - b;
            obj = obj'*obj;
            tmp = sparseTrans.mult(reshape(xnew,[A.frame_size(1),A.frame_size(2),A.Q]));
            obj = abs(obj) + Lam*sum(abs(tmp(:)));
%             fprintf('Iter = %s \tobjective= %s\ttime/iter = %s\n',num2str(iter), num2str(obj), ...
%                num2str(toc(tStart)));
        end
    end
      
    %Update x and t
    x = xnew;
    t = tnew;
    
end
toc(t_start);

res = (A.mult(z) - b);
%Trim the outputs if early termination occurred
if saveHist && (iter < options.nit)
    estHist.xhat = estHist.xhat(:,1:iter);
end

end


function w = shrink1(s, alph, p, ep)
t = abs(s);
w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
t(t == 0) = 1;
w = w./t;
end