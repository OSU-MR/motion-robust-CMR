function [ sos ] = sos_combine( x )
%SOS_COMBINE combines multicoil data using Sum of Squares
%   Detailed explanation goes here

if length(size(x))==3
    sos = squeeze(sqrt(sum(abs(x).^2,3)));
else
    sos = squeeze(sqrt(sum(abs(x).^2,4)));
end

end

