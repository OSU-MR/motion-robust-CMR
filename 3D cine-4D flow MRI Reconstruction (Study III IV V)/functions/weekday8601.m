function [day,str] = weekday8601(dtn,fmt)
% WEEKDAY8601 returns the Day of Week, as an integer and string.
%
% (c) 2016-2020 Stephen Cobeldick
%
% WEEKDAY8601 is a drop-in replacement for MATLAB's standard WEEKDAY
% function, but has the week starting on Monday, as per ISO 8601.
% This means an output of 1 corresponds to Monday, and 7 to Sunday.
%
%%% Syntax:
% weekday8601()
% weekday8601(dateNumbers)
% [dayNumber,dayName] = weekday8601(...)
% [dayNumber,dayName] = weekday8601(dateNumbers,nameFormat)
%
%% Input and Output Arguments
%
%%% Inputs (*=default):
% dtn = Numeric array of serial date numbers. No input uses the current date.
% fmt = CharacterVector, 'short'/'long' -> 3-character / full weekday name.
%     = LogicalScalar,    *false/true   -> 3-character / full weekday name.
%
%%% Outputs:
% day = NumericArray of the ISO 8601 day of the week (Monday=1, Sunday=7).
% str = CharacterArray of the day of the week (rows linear indexed from d).
%
% See also CALENDAR8601 DATENUM8601 DATESTR8601 WEEKDAY DATENUM DATESTR DATEVEC DATETIME DATEROUND

%% Input Wrangling
%
if nargin==0
	dtn = now;
else
	assert(isnumeric(dtn)&&isreal(dtn),...
		'SC:weekday8601:dtn:NotRealNumeric',...
		'First input <dtn> must be a real numeric array.')
	assert(all(isfinite(dtn(:))),...
		'SC:weekday8601:dtn:NotFiniteValue',...
		'First input <dtn> must contain finite values only.')
end
%
if nargin>1&&ischar(fmt)&&isrow(fmt)
	[isf,idf] = ismember(lower(fmt),{'short','long'});
	assert(isf,...
		'SC:weekday8601:fmt:NotShortLongOrLogical',...
		'Second input <fmt> must be either ''short'' or ''long''.')
else
	idf = 1+(nargin>1&&fmt);
end
%
%% Calculate Weekday
%
day = mod(fix(dtn)-3,7)+1;
nmc = {'Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'};
%
if nargout>1
	if isempty(day) % no names:
		str = '';
	elseif isscalar(day) % one name:
		str = nmc{day}(1:min(idf*idf*3,end));
	else % multiple names:
		str = char(nmc{day});
		tmp = [3,size(str,2)];
		str = str(:,1:tmp(idf));
	end
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%weekday8601