function [ muC, gammaC, mu, gamma, prob ] = gmmArto( epsilon, opt )
%% gmmArto(): Gaussion Mixture Model for ARTO
%
%   The goal of this function is to find a 3-component gaussian
%   distribution using expectation maximization which fits the weighted residuals of the phase image after
%   background phase correction. The central gaussian is assumed to
%   represent the residuals of the underlying background phase only. Any side
%   peaks are instead derived from outliers such as extensive venous flow
%   or spatial wrapping. Statistics computed for the central peak are then
%   used to make a binary decision whether to keep or remove any given pixel.
%
% INPUTS:   epsilon [m x 1]: Weighted residual error vector
%           gmmComp [Scaler]: Number of peaks in mixture model (3 default)
% OUTPUTS:  mu [Scaler]: Mean of central peak
%           gamma [Scaler]: Standard deviation of central peak
% 
% Adapted from example code by Chris McCormick
% http://mccormickml.com/2014/08/04/gaussian-mixture-models-tutorial-and-matlab-code/
% =========================================================================
%% Default Options
if ~isfield(opt,'gmmComp'), opt.gmmComp     = 3;                                                    end % Number of Gaussian components
if ~isfield(opt,'delta'),   opt.delta       = 2;                                                    end % Minimum # of std's side distributions must be away from central
if ~isfield(opt,'p'),       opt.p           = 0.5;                                                  end % Minimum "strength" of central distribution
if ~isfield(opt,'Kmax'),    opt.Kmax        = 1000;                                                 end
if ~isfield(opt,'tol'),     opt.tol         = var(epsilon)*1e-4;                                    end
if ~isfield(opt,'mu0'),     opt.mu0         = [-opt.delta*std(epsilon),0,+opt.delta*std(epsilon)];  end %
if ~isfield(opt,'gamma0'),  opt.gamma0      = ones(1, opt.gmmComp) * sqrt(var(epsilon))/2;          end
if ~isfield(opt,'prob0'),   opt.prob0       = [(1-opt.p)/2,opt.p,(1-opt.p)/2];                      end

%% INITIALIZATION
epsilon = sort(epsilon);

% Define contraints on the EM algorithm

Kmax = 1000;

tol = opt.tol;

% Set 'gmmComp' to the number of clusters to find
gmmComp = opt.gmmComp;

% Distance (in the units of std. of the middle peak) from the middle peak
delta = opt.delta;

% minimum strength of the central peak, p_peak1 + p_peak2 + ... p_peakN = 1
p = opt.p;

% Choose initial values for the parameters

% Set 'm' to the number of data points
m = size(epsilon, 1);

% Use the overall variance of the dataset as the initial variance for each cluster.
gamma = opt.gamma0;

% Initialize means similar to how we expect the solution to look like
mu = opt.mu0;
if min(mu) < min(epsilon)
    mu(mu==min(mu)) = min(epsilon);
end
if max(mu) > max(epsilon)
    mu(mu==max(mu)) = max(epsilon);
end

prob = opt.prob0;

%% EXPECTATION MAXIMIZATION

% Matrix to hold the probability that each data point belongs to each
% cluster. One row per data point, one column per cluster
W = zeros(m, gmmComp);

% Loop until convergence.
for k = 1 : Kmax
    %% Expectation
    % Calculate the probability for each data point for each distribution.

    % Matrix to hold the pdf value for every data point for every cluster.
    % One row per data point, one column per cluster.
    pdf = zeros(m, gmmComp);
    
    % For each cluster...
    for j = 1 : gmmComp
        % Check if gamma is equal to zero since that will cause a division
        % by zero when evaluating the gaussian
        if gamma(j) <= std(epsilon)/50 
            % Set to lower limit
            gamma(j) = std(epsilon)/50;
        end
        % Evaluate the Gaussian for all data points for cluster 'j'.
        pdf(:, j) = gaussian1D(epsilon, mu(j), gamma(j));
    end       
    
    
    % Multiply each pdf value by the prior probability for each cluster.
    pdf_w = bsxfun(@times, pdf, prob);
    
    % Divide the weighted probabilities by the sum of weighted probabilities for each cluster.
    sum_pdf_w = sum(pdf_w, 2);
    
    % Make sure sum_pdf_w > 0 to prevent divide by zero error
    if min(sum_pdf_w(:)) == 0
        sum_pdf_wTMP = sum_pdf_w(:);
        secondMin = min(sum_pdf_wTMP(sum_pdf_wTMP > 0));
        sum_pdf_w(sum_pdf_w==0) = secondMin;
    end
    
    W = bsxfun(@rdivide, pdf_w, sum_pdf_w);

    %% Maximization
    % Calculate the probability for each data point for each distribution
    
    % Store the previous means so we can check for covergence.
    prevMu = mu;

    % For each of the clusters...
    for j = 1: gmmComp
        
        % Calculate the prior probability for cluster 'j'
        prob(j) = mean(W(:, j));
        
        % Calculate the new mean for cluster 'j' by taking the weighted
        % average of *all* data points.
        mu(j) = weightedAverage(W(:, j), epsilon);
        
        % Calculate the variance for cluster 'j' by taking the weighted
        % average of the squared differences from the mean for all data
        % points.
        variance = weightedAverage(W(:, j), (epsilon - mu(j)).^2);
        
        % Calculate gamma by taking the square root of the variance.
        gamma(j) = sqrt(variance);
        
    end
    %% CONSTRAINTS
    % Adjust mean, gamma, and prob if necessary given the defined contraints
    
    % Location of the left-most peak
    [~,mnInd] = min(mu);
    
    % Location of the right-most peak
    [~,mxInd] = max(mu); 
    
    % Location of the middle peak
    mdInd = find(ismember([1:gmmComp],[mnInd, mxInd])==0);
    
    % Push leftmost peak to the left if needed
    mu(mnInd) = min(mu(mnInd), -gamma(mdInd)*delta);
    
    % Push rightmost peak to the right if needed
    mu(mxInd) = max(mu(mxInd), +gamma(mdInd)*delta);
    
    % Make sure the central peak is the tallest
    prob(mdInd) = max(prob(mdInd), p);
    
    % Check for convergence.
    
    if (abs(mu - prevMu) < tol)
        break
    end
    
end

% Outputs

muC = mu(mdInd);
gammaC = gamma(mdInd);
