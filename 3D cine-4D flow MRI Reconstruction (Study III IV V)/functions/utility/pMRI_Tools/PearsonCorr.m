function [ pearson_r ] = PearsonCorr( s1, s2, varargin)
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178.osu.edu
% Last update:  6/17/2014
%**************************************************************************
%
%       UNTITLED Plots the Pearson Correlation between two dataset s1 and 
%   s2.
% 
% Required Inputs:
%   s1:      Data Set 1 as a row vector
% 
%   s2:      Data Set 2 as a row vector
%
% Optional Inputs:
%   
% 
% Outputs:
%   pearson_r:  The Pearson Correlation coefficient
%**************************************************************************

% Init options 
opts.t1 = '';
opts.t2 = '';

% Check for input options
if length(varargin) ==1
    opts = varargin{1};
end

% Calculate Pearson r
pearson_r = corrcoef([s1.',s2.']);
pearson_r = pearson_r(2);

% Set plotting limits
% limits = [min([min(s1),min(s2)])-0.1*min([min(s1),min(s2)]) 1.1*max([max(s1),max(s2)])];
limits = [0, 1.1*max([max(s1),max(s2)])];
% Plot pearsone
% figure
hold on
scatter(s1,s2)
p1 = plot(linspace(limits(1),limits(2),300),linspace(limits(1),limits(2),300),'linewidth',2);
set(p1,'Color',[.3,0.3,0.3])
grid
xlabel(['',opts.t1],'fontsize',12)
ylabel(['',opts.t2],'fontsize',12)
xlim(limits)
ylim(limits)
end

