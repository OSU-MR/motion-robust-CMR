function [lambda1,lambda2,NMSE_min] = gridSearch(fx,xTrue)
% keep_searching = 1;
% search_ind_num = 0;
% NMSE_min_prev = 100;
% range = 10;
% lambda_low = 1;
% lambda_high = 100;
% lambdaW_low = 1e-6;
% lambdaW_high = 1e-4;

num_samp = 6;
keep_searching = 1;
search_ind_num = 0;
NMSE_min_prev = 100;
range = 100;
lambda_low = 2;
lambda_high = 200;
lambdaW_low = 1e-5;
lambdaW_high = 1e-3;

while keep_searching && (search_ind_num <3)
    
    log_lambda = logspace(log10(lambda_low),log10(lambda_high),num_samp);
    log_lambdaW = logspace(log10(lambdaW_low),log10(lambdaW_high),num_samp);
    
%     log_lambda = linspace((lambda_low),(lambda_high),num_samp);
%     log_lambdaW = linspace((lambdaW_low),(lambdaW_high),num_samp);
    lambda_ind = 1;
    NMSE_admm = zeros(length(log_lambdaW),length(log_lambda));
    for lambda = log_lambda
        lambdaW_ind = 1;
        for lambdaW = log_lambdaW 

            xhat = fx(lambda,lambdaW);
            
            NMSE_admm(lambdaW_ind,lambda_ind) = 20*log10(norm( xhat(:) - xTrue(:) )/norm(xTrue(:)));
%             NMSE_admm(lambdaW_ind,lambda_ind) = xhat;
            lambdaW_ind = lambdaW_ind+1;
            
            % plot NMSE
            [x_mesh,y_mesh] = meshgrid(log_lambda,log_lambdaW);
            figure(1);
            clf
            A = axes;
            mesh(x_mesh,y_mesh,NMSE_admm)
            set(A,'XScale','log')
            set(A,'YScale','log')
            ylabel('LambdaW')
            xlabel('Lambda')
            title(['Optimizing BM3D 2D'])
            drawnow
        end
        lambda_ind = lambda_ind +1;
    end

    %%
    [x_mesh,y_mesh] = meshgrid(log_lambda,log_lambdaW);
    figure;
    clf
    A = axes;
    mesh(x_mesh,y_mesh,NMSE_admm)
    set(A,'XScale','log')
    set(A,'YScale','log')
    ylabel('LambdaW')
    xlabel('Lambda')
    title(['Optimizing BM3D 2D'])
    
    [NMSE_min,min_ind] = min(NMSE_admm(:));
    
    range = range/2;

    lambda_min = x_mesh(min_ind);
    lambdaW_min = y_mesh(min_ind);
    lambda_high = lambda_high/6 + lambda_min;
    lambdaW_high = lambdaW_high/6 + lambdaW_min;
    lambda_low = max(eps,-lambda_min/6+lambda_min);
    lambdaW_low = max(eps,-lambdaW_min/6+lambdaW_min);
    
    
    NMSE_vec(1,search_ind_num+1) = NMSE_min;
    NMSE_vec(2,search_ind_num+1) = lambda_min;
    NMSE_vec(3,search_ind_num+1) = lambdaW_min;
    
    search_ind_num = search_ind_num +1;
    
%     if search_ind_num >1
%         range = range/2;
%     end
    
%     if (NMSE_min_prev-NMSE_min)/NMSE_min_prev*100 < 1 && search_ind_num >1
%         keep_searching = 0;
    if range<=2
        keep_searching = 0;
    end
    
    NMSE_min_prev = NMSE_min;
    
end
NMSE_vec
[~,NMSE_min_col] = min(NMSE_vec(1,:))
NMSE_min= NMSE_vec(1,NMSE_min_col)
lambda1 = NMSE_vec(2,NMSE_min_col)
lambda2 = NMSE_vec(3,NMSE_min_col)
