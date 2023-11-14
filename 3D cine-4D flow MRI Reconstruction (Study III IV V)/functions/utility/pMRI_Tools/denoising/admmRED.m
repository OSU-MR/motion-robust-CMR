function [x,z, history] = admmRED(x0, b, A, lambda, denoiser)
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
MAX_ITER = 100;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
alpha = 1;
rho = 0.1;

x = gpuArray(single(x0));
z = gpuArray(single(x0));
u = gpuArray(single(zeros(A.N,1)));

% x = ((x0));
% z = ((x0));
% z = ((zeros(A.N,1)));
% u = ((zeros(A.N,1)));

% x = A.A'*b;
% z = A.A'*b;

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

% pinvA = (A.A'*A.A + rho*eye(A.N,A.N))^(-1);
for k = 1:MAX_ITER
    
    % x-update
    obj = @(x_d) 1/2*norm(A.mult(x_d)-b,2)^2 + rho/2*norm(x_d-z+u,2)^2;
    dx_obj = @(x_d) A.multTr( A.mult(x_d) - b)  + rho*(x_d-z+u);
    
%     %========================================================================== 
%     % Gradient Test
%     for ind = 1:10
%         delta = 0.01;
%         test = randn(size(x));
%         dx_approx = (obj(test+delta*test)- obj(test-delta*test))/(2*delta);
%         dx_func = test'*dx_obj(test);
%         dx_approx
%         dx_func
%     end
% 
%      display('Solving Data Fidelity Step')
%    

    
    %==========================================================================
    % NL conjugate Gradient
    if ~QUIET
        display('Solving Data Fidelity Step')
    end
    
    % Take one steepest descent step
    step = 1;
    step_ind = 1;
    obj_pass = 0;
    obj_start = obj(x);
    gradient = dx_obj(x);
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
    if ~QUIET
        fprintf('it = %d\tstep size = %s, objective value %s\n',1,num2str(step),num2str(obj(x)))
    end
    
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
%         b_conj  = (gradient(:)'*gradient(:))/(gradientPrev(:)'*gradientPrev(:)+eps); % Fletcher–Reeves
        b_conj =  (gradient(:)'*( gradient(:)-gradientPrev(:)))/(gradientPrev(:)'*gradientPrev(:)+eps);% Polak–Ribière
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
        
        if ~QUIET
            fprintf('it = %d\tstep size = %s, objective value %s\n',ind,num2str(step),num2str(obj_end))
        end
        stepPrev = step;
        if k>0
            if (obj_start-obj_end)/obj_start<0.001
                break
            end
        end
    end

%========================================================================== 
%     P inverse technique
%     x = (A.multTr(A.mult(eye(A.N,A.N))) + rho*eye(A.N,A.N))^(-1)*(W_multTr(z)-W_multTr(u) + A.multTr(b));
%     x = (A.A'*A.A + rho*eye(A.N,A.N))^(-1)*(rho*(z-u) + A.A'*b);
%     x = pinvA*(rho*(z-u) + A.A'*b);
%     fprintf('objective value %s\n',num2str(obj(x)))
    
    if ~QUIET
        display('Solving RED Step')
    end
    % z-update with relaxation
    zold = z;
    zStar = rho*(x+u);
    if ~QUIET
        red_cost = @(q) lambda/2*q'*(q-denoiser(q)) + rho/2*norm(q-x-u)^2;
%         fprintf('RED objective value %s\n',num2str(red_cost(z)))
    end
    
    % Solve RED Step
    if k>20
        redNit = 2;
    else
        redNit = 5;
    end
    for ind = 1:redNit
%         figure(1)
%         subplot(121)
%         tmp = reshape(z,[A.frame_size,A.Q]);
%         imagesc(abs(tmp(:,:,1)));
% %         imagesc(reshape(z,A.frame_size))
        zBar = denoiser(z);
        z = 1/(rho+lambda)*(lambda*zBar+zStar);
        if ~QUIET
%             fprintf('RED objective value %s\n',num2str(red_cost(z)))
        end
%         subplot(122)
%         tmp = reshape(z,[A.frame_size,A.Q]);
%         imagesc(abs(tmp(:,:,1)));
% %         imagesc(reshape(z,A.frame_size))
%         pause()
    end

    % u-update
    u = u + (x - z);

%     diagnostics, reporting, termination checks
    history.red(k) = lambda*x'*(x-denoiser(x));
    history.objval(k)  = 1/2*sum(abs(A.mult(x)-b).^2) + lambda*x.'*(x-denoiser(x));
%     history.objval(k)  = 1/2*sum(abs(A.mult(x)-b).^2) + lambda*x'*(x-zBar);
    history.r_norm(k)  = norm((x) - z);
    history.s_norm(k)  = norm(-rho*(z - zold));

    history.eps_pri(k) = sqrt(A.N)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(A.N)*ABSTOL + RELTOL*norm(rho*u);

    
    if ~QUIET
        fprintf('ADMM It = %3d\tprimal res%10.4f\tdual res%10.4f\tobj value %10.4f\n', k, ...
            history.r_norm(k),history.s_norm(k), history.objval(k));
    end

    % update rho to balance s_norm and r_norm
    mu = 10;
    tau = 2;
    if history.r_norm(k) > mu*history.s_norm(k)
        rho = tau*rho;
        u = u/tau;
    elseif history.s_norm(k) > mu*history.r_norm(k)
        rho = rho/tau;
        u = u*tau;
    end

%     if (history.r_norm(k) < history.eps_pri(k) && ...
%        history.s_norm(k) < history.eps_dual(k))
%          break;
%     end

end

if ~QUIET
    toc(t_start);
end

end









