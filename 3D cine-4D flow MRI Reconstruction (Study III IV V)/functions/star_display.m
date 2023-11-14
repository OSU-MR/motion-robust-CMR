function [ blah ] = star_display( string_in,varargin)
% Display a string inside stars to the matlab terminal
% varargin{1}   toggle wrapping the string in a set of stars
if isempty(varargin);
    wrap = 1;
    n = 80;     % star length
elseif length(varargin) ==1
    wrap = varargin{1};
    n = 80;     % star length
elseif  length(varargin) ==2
    wrap = varargin{1};
    n = varargin{2};    % star length
end

start_string = [];
mid_string = [];
string_length = length(string_in);
string_length = string_length/2;

for ind = 1:n;
    start_string = [start_string,'*'];
end

for ind = 1:n;
   if ind <= n/2-string_length
       mid_string = [mid_string,'*'];
   elseif ind == n/2+1
       mid_string = [mid_string,string_in];
   elseif ind > n/2+string_length
       mid_string = [mid_string,'*'];
   end
   
end

if wrap
    fprintf(start_string)
    fprintf(mid_string)
    fprintf([start_string,'\n'])
else
    fprintf([mid_string,'\n'])
end

end

