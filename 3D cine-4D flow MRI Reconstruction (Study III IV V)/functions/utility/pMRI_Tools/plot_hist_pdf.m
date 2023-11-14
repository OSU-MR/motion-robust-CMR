function plot_hist_pdf( x, bin_length,varargin)
%PLOT_HIST_PDF Plots the histrogram of the inpute x on the same scale as
%the pdf
%
% Inputs:
%           x: vector signal to take the histogram of
%       
%           bin_length: The number of bins to use for the histogram
%
%           varargin: Optional figure number as input
if isempty(varargin)
    figure
else
    figure(varargin{1})
end
[x_hist,bins] = hist(x,bin_length);
delx = abs(bins(1)-bins(2));
x_hist = x_hist/(delx*length(x));

bar(bins,x_hist)
end

