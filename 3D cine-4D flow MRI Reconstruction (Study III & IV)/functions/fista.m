function [x,obj_fun,times,mse] = fista(AHA,AHy,x0,options,x_true)

% Fast Iterative Soft Thresholding (FISTA) 
%   Code to implement FISTA as described in "A Fast Iterative
%   Shrinkage-Thresholding Algorithm for Linear Inverse Problems" by Beck
%   and Teboulle. Implements the FISTA algorithm to solve the convex
%   program argmin_x ||Ax - y||_2^2 + |lam .* x|. The implementation
%   support complex valued data. The code will run a test case if no inputs
%   are provided.
%--------------------------------------------------------------------------
%
%   Code inputs:
%   -x0 is the initial guess at the solution
%   -AHA is a function handle that computes A^H*A*x
%   -AHy is the precomputed value of A^H*y
%   -x_true is the true value of the image, used to compute MSE.
%   -options is the structure of options
%       -options.L = Lipschitz constant of the operator A. This should be
%          twice the largest eigenvalue of A'A. If not provided, the code
%          will estimate it using the power iteration.
%       -options.max_iter: max allowed iterations, default 1000
%       -options.thresh: stopping criteria, ratio of cost function change
%          to cost function value. Default is 1e-6
%       -options.lam is the vector (or scalar) regularization parameter.
%       -options.p l_p norm p value
%       -option.lp_tol epsilon factor added to the shrink operator.
%          Default is 1e-6
%       -options.v, both mse and objective function values are calculated
%          after every options.v iterations. Default is 10.
%
%   Code Outputs:
%   -x_hat is the estimate of the cost function's argmin.
%   -mse is the mean-square-error, calculated when x_ture is provided
%   -obj_fun the value of the objective funtion for each iteration. Only
%            computed if option.obj=1
%   -times is the accumulative time for FISTA iterations
% -------------------------------------------------------------------------
%
% 1D FISTA code written by Jason Parker (jason.parker@wpafb.af.mil)
% Modified by Christian Austin for p<1(austinc@ece.osu.edu)
% -------------------------------------------------------------------------
%
% References:
% [1]  Amir Beck and Marc Teboulle, A Fast Iterative Shrinkage-Thresholding 
%      Algorithm for Linear Inverse Problems, SIAM J. IMAGING SCIENCES, 
%      Vol. 2, No. 1, pp. 183–202


%% Starting time
T=tic; % Total time of execution

%% Check for any missing or excessive input parameters
switch nargin
    case 0 | 1 
        error('Need at least three input arguments');
    case 2
        x0=zeros(size(AHy));
        options=[];
        x_true=[];
    case 3
        options=[];
        x_true=[];
    case 4
        x_true=[];
    case 5
        % All parameters provided, do nothing
    otherwise
%         error('Maximum number of inputs exceeded');
end

%% Verify that options are defined
if ~isfield(options,'max_iter'),    options.max_iter = 5000;    end %Check iterations
if ~isfield(options,'p'),           options.p = 1;              end %Check p value
if ~isfield(options,'thresh'),      options.thresh = 1e-8;      end %Check thresh
if ~isfield(options,'lp_tol'),      options.lp_tol = 1e-8;      end %Check lp_tol
if ~isfield(options,'lam'),         options.lam = 1e0;          end %Check lam
if (~isfield(options,'obj')),       options.obj = 0;            end %Check obj
if (~isfield(options,'v')),         options.v = 50;             end %Check v
if (~isfield(options,'lam')),       options.lam = 50;           end %Check regularization strength

%Use the power iteration to estimate the Lipschitz constant if not provided
if ~isfield(options,'L')
    %Initial guess
    q = randn(size(x0));
    q = q / norm(q);
    th = 1e-4; % threshold
    err = inf;
    uest = inf;
    while err > th
        q = AHA(q); 
        %Progress
        unew = norm(q);
        err = abs(log(uest / unew));
        %Save the estimate
        uest = unew;
        q = q / norm(q);
    end
    %Use a little more than twice the eigenvalue estimate
    options.L = 2.05*uest;
end
% Extract all parameters
lam      = options.lam;
p        = options.p;
L        = options.L;
max_iter = options.max_iter;
thresh   = options.thresh;
lp_tol   = options.lp_tol; 
obj      = options.obj;
v        = options.v;

%Verify that lam is in fact a vector
if numel(options.lam) == 1
    lam = lam * ones(size(x0));
end


%% Iteration

%Initialize variables
y = x0;
x = x0;
t = 1;
stop = 0;
iter = 0;
cchange = inf;
x_old=0;
y2=real(y'*y);
mse=zeros(max_iter,1);
obj_fun=zeros(max_iter,1);
times=zeros(max_iter,1);


%Walk through iterations
while ~stop
    ti=tic; % Time for each outer loop iteration

    %Compute the new x
    x_new = pfun(y, AHy, AHA, lam, p, L, lp_tol);
    
    %Compute the new t
    t_new = (1 + sqrt(1 + 4*t^2))/2;
    
    %Compute the new y
    y = x_new + (t - 1)/t_new*(x_new - x);

    %Counter
    iter = iter + 1;
    if rem(iter,v)==0 % Compute MSE and cost fuction for every 25th iter
%         disp(['Iteration number: ' num2str(iter)]);
        if obj==1
            obj_fun(iter/v)= real(x_new'*AHA(x_new)-...
                       2*x_new'*AHy +y2) + sum(abs(lam(:).*x_new(:)));
%               obj_fun(iter/v) = norm((AA*x_new - y),2);
        end
        if ~isempty(x_true) 
            mse(iter/v)= real((x_new-x_true)'*((x_new-x_true)))/numel(x_new); 
        end
    end
     
    %Check max iterations
    if iter >= max_iter
        stop = 1;
    else %no need to do the threshold check if we have exceeded max_iter
         %cnew = norm(A(x_new) - y)^2 + sum(lam.*((abs(x_new).^2+lp_tol).^(p/2)));
%         cnew = x_new;
        if norm(x_new)>0
            cchange = norm(x_new - x_old) / norm(x_new);
        else
            cchange = 0;
        end
        x_old=x_new;
%         if cchange < thresh
%             stop = 1;
%         end
    end
    change_delta(iter) = cchange;
    %Update x and t
    x = x_new;
    t = t_new;
    if iter==1, times(iter)=toc(ti);
    else times(iter)=times(iter-1)+toc(ti);
    end
end
mse=mse(1:floor(iter/v));
obj_fun=obj_fun(1:floor(iter/v));
times=times(v:v:iter);
run_time=toc(T);

return


%% Helper functions

function g = pfun(y,AHy,AHA,lam,p,L,lp_tol)
%% Inner function
%Compute the gradient term
alpha = y - (1 / L)*(2*(AHA(y) - AHy));
%Soft threshold
g = shrink1(alpha,lam/L,p,lp_tol);
return

function w = shrink1(s,alph, p, ep)
t = abs(s);
w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
t(t == 0) = 1;
w = w./t;
return;
