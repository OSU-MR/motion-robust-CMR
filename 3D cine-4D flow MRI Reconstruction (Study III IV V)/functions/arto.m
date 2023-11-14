function [S,mu,gamma,prob] = arto(phi_Cor,sigma,M,opt)
% Compute weighted residuals
Wr = phi_Cor ./ sigma; % Weighted residuals
epsilon = Wr(M); 
epsilon(isnan(epsilon))=[];

epsilon = reshape(epsilon,[length(epsilon),1]);

% Gaussian mixture model

% Compute gaussian mixture model and extract central peak statistics
[muC, gammaC, mu, gamma, prob] = gmmArto(epsilon,opt);

% Upper threshold.
threshUpper = muC+opt.tau*gammaC;

% Lower threshold.
threshLower = muC-opt.tau*gammaC;

% Exclusion

S = ones(size(M));

% Loop over  1st spatial dimension (Y).
for i = 1 : size(Wr,1)
    
    % Loop over 2nd spatial dimension (X).
    for j = 1 : size(Wr,2)
        
        for k = 1 : size(Wr,3)
        
            % Exclude if outside threshold boundaries
            if Wr(i,j,k) <= threshLower || Wr(i,j,k) >= threshUpper
            
                S(i,j,k) = 0;
            
            elseif isnan(Wr(i,j,k))
            
                S(i,j,k) = 0;
            
            end
        end
    end
end

S = logical(S);

end

