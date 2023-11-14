function [x] = cgSENSEL2(x0, b, A, lam)
%

t_start = tic;

QUIET    = 0;
MAX_ITER = 100;

% x = gpuArray(single(zeros(A.N,1)));
x = gpuArray(single(x0));



% x-update
obj = @(x_d) 1/2*norm(A.mult(x_d)-b,2)^2 + lam/2*norm(x_d)^2;
dx_obj = @(x_d) A.multTr( A.mult(x_d) - b) + lam*x_d;


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
if ~QUIET
    fprintf('it = %d\tstep size = %s, objective value %s\n',1,num2str(step),num2str(obj(x)))
end

s = gradient;
stepPrev= 0.1;
for ind = 2:MAX_ITER
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

    if ~QUIET
        fprintf('it = %d\tstep size = %s, objective value %s\n',ind,num2str(step),num2str(obj_end))
    end
end

%========================================================================== 
%     P inverse technique
%     x = (A.multTr(A.mult(eye(A.N,A.N))) + rho*eye(A.N,A.N))^(-1)*(W_multTr(z)-W_multTr(u) + A.multTr(b));
%     x = (A.A'*A.A + rho*eye(A.N,A.N))^(-1)*(rho*(z-u) + A.A'*b);
%     x = pinvA*(rho*(z-u) + A.A'*b);
%     fprintf('objective value %s\n',num2str(obj(x)))
%     fprintf('objective value %s\n',num2str(obj(xt)))

if ~QUIET
    toc(t_start);
end

end









