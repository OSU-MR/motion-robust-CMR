function [x,estHist,res] = fista_pmri(x0,b,A,options)
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
if (nargin < 4)
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

W = nd_dwt_3D(options.wave,[A.frame_size(1),A.frame_size(2),A.Q],'pres_l2_norm',1,'compute','gpu');
k = 5;
lambda_band = zeros(8,1);
lambda_band(1) = 0.01*Lam;
lambda_band(2) = 1*Lam;
lambda_band(3) = 1*Lam;
lambda_band(4) = 1*Lam;
lambda_band(5) = k*Lam;
lambda_band(6) = k*Lam;
lambda_band(7) = k*Lam;
lambda_band(8) = k*Lam;

% lambda_band(1) = 0.01*Lam;
% lambda_band(2) = Lam;
% lambda_band(3) = Lam;
% lambda_band(4) = Lam;
% lambda_band(5) = Lam;
% lambda_band(6) = Lam;
% lambda_band(7) = Lam;
% lambda_band(8) = Lam;


%Walk through iterations
while ~stop
    tStart = tic; 
    %Counter
    iter = iter + 1;
    
    %Gradient step
    alphares = z - (1 / options.lip)*(2*A.multTr(A.mult(z) - b));
    alphares = W.dec(reshape(alphares,[A.frame_size(1),A.frame_size(2),A.Q]),1);
    
    for ind = 1:8
        alphares(:,:,:,ind) = shrink1(alphares(:,:,:,ind), lambda_band(ind)/options.lip, 1, 1e-6);
    end
    
    xnew = W.rec(alphares);
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
            tmp = W.dec(reshape(xnew,[A.frame_size(1),A.frame_size(2),A.Q]),1);
            obj = abs(obj) + Lam*sum(abs(tmp(:)));
            fprintf('Iter = %s \tobjective= %s\ttime/iter = %s\n',num2str(iter), num2str(obj), ...
               num2str(toc(tStart)));
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