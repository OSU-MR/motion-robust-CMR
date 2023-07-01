
function [x,z, history] = admmpMRI(x0, b, A, lambda, rho, alpha)
% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

t_start = tic;

QUIET    = 0;
MAX_ITER = 50;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

% x = gpuArray(single(zeros(A.N,1)));
x = gpuArray(single(x0));
z = gpuArray(single(zeros(A.N*8,1)));
u = gpuArray(single(zeros(A.N*8,1)));
AHb = gpuArray(single(A.multTr(b)));  % precompute A^H*b

% x = x0;
% z = ((zeros(A.N*8,1)));
% u = ((zeros(A.N*8,1)));

alpha = 1;
rho = 0.1;
if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end
W = nd_dwt_3D('db1',[A.frame_size(1),A.frame_size(2),A.Q],'pres_l2_norm',1,'compute','gpu');
W_mult = @(x) reshape(W.dec(reshape(x,[A.frame_size(1),A.frame_size(2),A.Q]),1),[numel(x)*8,1]);
W_multTr = @(x) reshape(W.rec(reshape(x,[A.frame_size(1),A.frame_size(2),A.Q,8])),[numel(x)/8,1]);

% % Find Lip constant
% tstart= tic;
% fprintf('Calculating Lipschitz Constant...')
% 
% %Initial guess
% q = randn(length(x0),1);
% q = q / norm(q);
% 
% thresh = 1e-3;
% err = inf;
% uest = inf;
% while err > thresh
%     q = A.multTr(A.mult(q));
% 
%     %Progress
%     unew = norm(q);
%     err = abs(log(uest / unew));
% 
%     %Save the estimate
%     uest = unew;
%     q = q / norm(q);
% end

% %Use a little more than twice the eigenvalue estimate
% lip = 2.05*uest;

for k = 1:MAX_ITER
    
    % x-update
%     obj = @(x) 1/2*sum(abs(A.mult(x)-b).^2) + rho/2*sum(abs(W_mult(x)-z+u).^2);
    obj = @(x_d) 1/2*norm(A.mult(x_d)-b,2)^2 + rho/2*norm(W_mult(x_d)-z+u,2)^2;
%     dx_obj_old = @(x_d) A.multTr( A.mult(x_d) - b ) + rho*W_multTr( W_mult(x_d)-z+u );
    dx_obj = @(x_d) A.multTr( A.mult(x_d)) - AHb  + rho*(x_d+W_multTr(u-z));

%========================================================================== 
%     % Gradient Test
%     for ind = 1:10
%         delta = 0.01;
%         test = randn(size(x))+1j*randn(size(x));
%         dx_approx = (obj(test+delta*test)- obj(test-delta*test))/(2*delta);
%         dx_func = test'*dx_obj(test);
%         dx_approx
%         dx_func
%     end

%========================================================================== 
%     % Gradient Decent Method
%     display('Solving Data Fidelity Step')
%     for ind = 1:100
%         step = 10;
%         step_ind = 1;
%         obj_pass = 0;
%         obj_start = obj(x);
%         gradient = dx_obj(x);
%         while ~obj_pass && (step_ind < 31)
%             x_new = x - step*gradient;
%             obj_end = obj(x_new);
%             
%             if obj_start >= obj_end
%                 obj_pass = 1;
%                 x = x_new;
%             else 
%                 step = 0.7*step;
% %                 display('objective increased, decreasing step size');
%             end
%             step_ind = step_ind+1;
%         end
%         fprintf('it = %d\tstep size = %s, objective value %s\n',ind,num2str(step),num2str(obj(x)))
% %         grad_obj(ind) = obj(x);
%         if (obj_start-obj_end)/obj_start< 0.00005
%             break
%         end
%     end
    
%==========================================================================
% %     % NL conjugate Gradient
    display('Solving Data Fidelity Step')
    
    % Take one steepest descent step
    step = 1;
    step_ind = 1;
    obj_pass = 0;
    obj_start = obj(x);
    gradient = dx_obj(x);
%     x = x0;
    while ~obj_pass && (step_ind < 11)
        x_new = x - step*gradient;
        obj_end = obj(x_new);

        if obj_start >= obj_end
            obj_pass = 1;
            x = x_new;
        else 
            step = 0.5*step;
        end
        step_ind = step_ind+1;
    end
    fprintf('it = %d\tstep size = %s, objective value %s\n',1,num2str(step),num2str(obj(x)))
    nlcg_obj(1) = obj(x);
    % NL Conjugate
    if k ==1
        num_ind = 50;
    elseif k < 10
        num_ind = 25;
    else
        num_ind = 10;
    end

    s = gradient;
    stepPrev= 0.1;
    for ind = 2:num_ind
        step = 10*stepPrev;
        step_ind = 1;
        obj_pass = 0;
        obj_start = obj(x);
        gradientPrev  = gradient;
        gradient = dx_obj(x);
        b_conj  = (gradient(:)'*gradient(:))/(gradientPrev(:)'*gradientPrev(:)+eps); % Fletcher–Reeves
        s = gradient + b_conj*s; % calculate conjugate direction - ST
        while ~obj_pass && (step_ind < 11)

            x_new = x - step*s;
            obj_end = obj(x_new);

            if obj_start >= obj_end
                obj_pass = 1;
                x = x_new;
            else 
                step = 0.5*step;
            end
            step_ind = step_ind+1;
        end
%         nlcg_obj(ind) = obj(x);
        fprintf('it = %d\tstep size = %s, objective value %s\n',ind,num2str(step),num2str(obj_end))
        stepPrev = step;
        if k>0
            if (obj_start-obj_end)/obj_start<0.001
                break
            end
        end
    end

    
%==========================================================================
%     % Conjugate Gradient
%     display('Solving Data Fidelity Step')
% 
%     if k ==1
%         num_ind = 50;
%     elseif k < 10
%         num_ind = 25;
%     else
%         num_ind = 10;
%     end
%     for ind = 1:num_ind
%         obj_start = obj(x);
%         gradient = dx_obj(x);
%         r = A.multTr(A.mult(gradient)) + rho*gradient;
%         step  = abs((gradient(:)'*gradient(:))/(gradient(:)'*r(:)+eps)); % Fletcher–Reeves
% 
%         x = x + step*gradient;
%         obj_end = obj(x);
%         fprintf('it = %d\tstep size = %s, objective value %s\n',ind,num2str(step),num2str(obj(x)))
%         if k>0
%             if (obj_start-obj_end)/obj_start<0.001
%                 break
%             end
%         end
%     end

%========================================================================== 
%     % FISTA Method
%     display('Solving Data Fidelity Step')
%     t = 1;
%     zed = x;
%     obj_cur = obj(x);
%     if k ==1;
%         num_ind = 100;
%     elseif k < 10
%         num_ind = 30;
%     else
%         num_ind = 10;
%     end
%     for ind = 1:num_ind
%         gradient = dx_obj(x);
%         xnew = zed-1/lip*gradient;
%         tnew = (1 + sqrt(1 + 4*t^2))/2;
%         zed = xnew + (t - 1)/tnew*(xnew - x);
%         
%         %Update x and t
%         x = xnew;
%         t = tnew;
%         
%         obj_old = obj_cur;
%         obj_cur = obj(x);
%         fprintf('it = %d\t objective value %s\n',ind,num2str(obj_cur))
%         if k>15
%             if (obj_old-obj_cur)/obj_old<0.0001
%                 break
%             end
%         end
%     end

%==========================================================================    
%     P inverse technique
%     x = (A.multTr(A.mult(eye(A.N,A.N))) + rho*eye(A.N,A.N))^(-1)*(W_multTr(z)-W_multTr(u) + A.multTr(b));
%     x = (A.A'*A.A + rho*eye(A.N,A.N))^(-1)*(rho*W_multTr(z-u) + A.A'*b);
%     fprintf('objective value %s\n',num2str(obj(x)))

%========================================================================== 
%     % Fminunc technique
%     options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate'...
%                 ,'steepdesc','display','iter','MaxIter',10);
%     x = fminunc(obj,x,options);
    
    % z-update with relaxation
    zold = z;
    x_hat = alpha*W_mult(x) + (1 - alpha)*zold;
    z = shrink1(x_hat + u, lambda/rho,1,1e-6);
    
    % u-update
    u = u + (W_mult(x) - z);

    if mod(k,5)==0 || k ==1
        % diagnostics, reporting, termination checks
        history.objval(k)  = 1/2*sum(abs(A.mult(x)-b).^2) + sum(abs(lambda.*W_mult(x)));

        history.r_norm(k)  = norm(W_mult(x) - z);
        history.s_norm(k)  = norm(-rho*(z - zold));

    %     history.eps_pri(k) = sqrt(A.N)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    %     history.eps_dual(k)= sqrt(A.N)*ABSTOL + RELTOL*norm(rho*u);


        if ~QUIET
            fprintf('ADMM It = %3d\tprimal res%10.4f\tdual res%10.4f\tobj value %10.0f\n', k, ...
                history.r_norm(k),history.s_norm(k), history.objval(k));
        end
        
        % Update the rho parameter
        mu = 10;
        tau = 2;
        if history.r_norm(k) > mu*history.s_norm(k)
            rho = tau*rho;
            u = u/tau;
        elseif history.s_norm(k) > mu*history.r_norm(k)
            rho = rho/tau;
            u = u*tau;
        end

%         % Check for convergence
%         if k>1
%             if abs(history.objval(k-1)-history.objval(k))/history.objval(k-1) < 0.00001
%                 fprintf('ADMM It = %3d\tprimal res%10.4f\tdual res%10.4f\tobj value %10.0f\n', k, ...
%                     history.r_norm(k),history.s_norm(k), history.objval(k));
%                  break;
%             end
%         end
    end


end

if ~QUIET
    toc(t_start);
end
z = W_multTr(z);
end

function p = objective(A, b, lambda, x, z)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1) );
end

function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function [L U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end

    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

function [f,g] = obj(x, A, W, b, z, u)
    p = 1;
    f = sum(abs(A.mult(x)-b).^2) + p/2*sum(abs(W_mult(x)-z+u));
    g = A.multTr( A.mult(x) - b ) + p*W_multTr( W_mult(x)-z+ u );
end

function w = shrink1(s, alph, p, ep)
    t = abs(s);
    w = max([t - alph.*(t.^2 + ep).^(p/2 - 0.5)], 0).*s;
    t(t == 0) = 1;
    w = w./t;
end

% function [f,g] = func_obj(x_d,A,b,rho,W_mult,W_multTr,z,u)
%     f = 1/2*norm(A.mult(x_d)-b,2)^2 + rho/2*norm(W_mult(x_d)-z+u,2)^2;
%     g = gather(A.multTr( A.mult(x_d) - b ) + rho*W_multTr( W_mult(x_d)-z+u ));
% end
