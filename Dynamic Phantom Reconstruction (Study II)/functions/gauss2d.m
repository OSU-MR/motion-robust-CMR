function[val] = gauss2d(n,ctr,sig,row,nor)

% Generates 2D normal distribution

% n:    Size of the 2D matrix
% ctr:  Mean (central location) of distribution w.r.t floor(n/2)+1
% sig:  Std in x and y directions
% row:  Correlation coefficient
% nor:  0=don't normalize; 1=normalize area; 2 normalize peak value;
% 3=remove k from "k*exp^(-x^2)"
% ------ Rizwan Ahmad (ahmad.46@osu.edu) ---------


if nargin ==4
    nor = 1;
elseif nargin ==3
    row = 0;
    nor = 1;
elseif nargin <3
    error('At least three input arguments required');
end
[x,y] = ndgrid(1:n(1), 1:n(2));
cntr=floor(n/2)+1+ctr;
val = exp(-(1/(2*(1-row^2)))*...
      ((((x-cntr(1)).^2)/sig(1)^2 + ((y-cntr(2)).^2)/sig(2)^2) - 2*row*(x-cntr(1)).*(y-cntr(2))/(sig(1)*sig(2))));

if nor ==0
    val=1/(2*pi*sig(1)*sig(2)*sqrt(1-row^2))*val; % If FOV is large enough; area under the curve should be one
elseif nor ==1
    val=val/sum(val(:)); % For truncated FOV, returns the area under the curve to one
elseif nor ==2
    val=val/max(val(:)); % Forces the maximum value to one
elseif nor ==3
    % Do nothing
else
    error('Input argument "nar" can only assume integers between 0 and 3');
end