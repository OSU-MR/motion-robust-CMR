function [ lambda, theta1, theta2, sigma1, sigma2,gm] = GM_init( x, L, varargin )
%GM_INIT Initializes a GM model for use in GAMP using the EM algorithm
%
% Inputs:
%           x: A vector to fit the GMM to
%
%           L: The Number of mixture components
%
% Outputs:  
%           lambda: The estimated sparsity weighting
%
%           theta: The estimated means
%
%           sigma: The estimated variances
%
%           omega: The estimated weights
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Created:      4/26/2014
% Last update:  4/26/2014
%**************************************************************************

options = statset('Display','off');

% Plot the distributions if desired
if isempty(varargin)
    plot_out = 0;
else
    plot_out = varargin{1};
end

inits.mu = [0;0];
inits.Sigma = zeros(1,1,2);
inits.Sigma(1,1,1) = 0.001;
inits.Sigma(1,1,2) = var(x);
inits.PComponents = [0.5,0.5];
% Estimate the GM using EM
gm = fitgmdist(x,L,'Options',options);

% Find the index that corresponds to delta
[~,min_ind] = min(gm.Sigma);
[~,max_ind] = max(gm.Sigma);

% Set the values of sigma, theta, and omega
sigma1 = squeeze(gm.Sigma(min_ind)).';
theta1 = gm.mu(min_ind).';

sigma2 = squeeze(gm.Sigma(max_ind)).';
theta2 = gm.mu(max_ind).';

% Set lambda as the weight for the gaussian with smalles variance
lambda = 1-gm.PComponents(min_ind);

if plot_out
%     h = ezplot(@(x)pdf(gm,x),[min(x) max(x)]);
    t = linspace(min(x), max(x),300).';
    plot_hist_pdf(x,30)
    hold on
    h = plot(t,pdf(gm,t),'linewidth',2);
    set(h,'Color','red')
    ylim([0,0.3*max(pdf(gm,t))])
end

