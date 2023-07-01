function stats = BlandAltman( s1, s2, new_fig, line_lim)
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178.osu.edu
% Last update:  6/17/2014
%**************************************************************************
%function stats = BlandAltman( s1, s2, new_fig,)
%
%       Creates a Bland-Altman plot for  two dataset s1 and s2.
% 
% Required Inputs:
%   s1:      Data Set 1 as a row vector
% 
%   s2:      Data Set 2 as a row vector
%
%   new_fig: Set to one to use a new figure.  0 uses the current figure.
%
% Optional Inputs:
%   
% 
% Outputs:
%   stats:  output object.
%
%           stats.mean_diff: Mean difference on Bland-Altman plot
%           stats.upper_lim: Upper limit of agreement 95%
%           stats.lower_lim: Upper limit of agreement 95%
%
%**************************************************************************

if nargin < 3
    new_fig = 1;
    line_lim = [];
elseif nargin < 4
    line_lim = [];
end 

% Calculate the means and differences 
means = (s1+s2)/2;
diffs = s1-s2;

% Calculate the Mean and std deviation of the differences
mean_diff(1:length(diffs)) = mean(diffs);
std_diff = std(diffs);

% Calculate the Reference Interval mean +- 1.96 std
upper_lim(1:length(diffs)) = mean_diff + 1.96*std_diff;
lower_lim(1:length(diffs)) = mean_diff - 1.96*std_diff;

% Define output object
stats.mean_diff = mean_diff(1);
stats.std_diff = std_diff(1);
stats.upper_lim = upper_lim(1);
stats.lower_lim = lower_lim(1);

% Plot Bland-Altman
if new_fig
    figure
end
hold on
p1 = scatter(means,diffs);

% If now line limits are given, make a guess
if isempty(line_lim)
    p2 = plot(linspace(min(means),max(means),length(means)),mean_diff,'linewidth',2);
    p3 = plot(linspace(min(means),max(means),length(means)),upper_lim,'--','linewidth',2);
    p4 = plot(linspace(min(means),max(means),length(means)),lower_lim,'--','linewidth',2);
else
    p2 = plot(linspace(line_lim(1),line_lime(2),length(means)),mean_diff,'linewidth',2);
    p3 = plot(linspace(line_lim(1),line_lime(2),length(means)),upper_lim,'--','linewidth',2);
    p4 = plot(linspace(line_lim(1),line_lime(2),length(means)),lower_lim,'--','linewidth',2);
end
set(p2,'Color',[0.5,0.5,0.5])
set(p3,'Color','red')
set(p4,'Color','red')

% set(p2,'Color',[.3,0.3,0.3])
% set(p3,'Color',[.3,0.3,0.3])
% set(p4,'Color',[.3,0.3,0.3])
grid
% xlabel('Mean','fontsize',12)
% ylabel('Difference','fontsize',12)
lim = 1.5*[lower_lim(1)-mean_diff(1),upper_lim(1)-mean_diff(1)] + mean_diff(1);
ylim(lim)

end

