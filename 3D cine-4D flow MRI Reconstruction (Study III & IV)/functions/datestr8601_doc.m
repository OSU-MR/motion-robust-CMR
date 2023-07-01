%% DATESTR8601 Examples
% The function <https://www.mathworks.com/matlabcentral/fileexchange/34095
% |DATESTR8601|> generates ISO 8601 formatted date strings (char vectors),
% specifying the date string format using a very compact token syntax.
%
% *In general |DATESTR8601| token syntax is not compatible with |DATESTR|
% or |DATETIME|*. Note |DATESTR8601| does not parse or handle timezones.
%
% This document shows examples of |DATESTR8601| usage. Most examples
% use the date+time given by the date vector |[1999,1,3,15,6,48.0568]|.
%% Basic Usage
% Calling |DATESTR8601| without any input arguments uses the current
% time and the default token |'ymdHMS'| (i.e. basic calendar date):
datestr8601()
%% Input 1: Date Number, Date String, Date Vector, or Datetime
% The first input can be a serial date number, a date string, a date
% vector, or a datetime scalar. An empty numeric uses the current time.
datestr8601(now()-3)
datestr8601([1999,1,3])
datestr8601(datetime())
datestr8601([])
%% Inputs 2+: Timestamp Format Tokens
% Any further input arguments are tokens that define the ouput date format:
% the tokens for |DATESTR8601| are specifically designed to allow compact
% ISO 8601 timestamp definition (and are *NOT* the same as those of
% |DATESTR| or |DATETIME|). All of the date tokens follow one simple rule:
%
% * lowercase = calendar year
% * UPPERCASE = week-numbering year
%
% By default the token defines ISO8601 basic notation (i.e. no delimiter
% characters between units), extended notation (with delimiters) can be
% selected by prefixing an asterisk |'*'| to the token.
%
% *See the Mfile help for the complete timestamp token specification.*
datestr8601([],'YWDHMS') % week-numbering date, basic.
datestr8601([],'*ynHMS') % ordinal date, extended.
datestr8601([],'ymdHMS') % calendar date, basic.
datestr8601([],'ymd@H7') % reduced precision, separator, decimal fraction.
%% Inputs 2+: Single Date/Time Tokens
% The input tokens can also specify single date/time values, a few examples
% of which are shown below. All of the date tokens follow one simple rule:
%
% * lowercase = calendar year
% * UPPERCASE = week-numbering year
%
% *See the Mfile help for the complete list of single-value tokens.*
datestr8601([],'yyyy') % calendar year
datestr8601([],'YYYY') % week-numbering year
datestr8601([],'n*')   % day of the calendar year
datestr8601([],'N*')   % day of the week-numbering year
datestr8601([],'DDDD') % day of the week
datestr8601([],'XX')   % weeks remaining in week-numbering year
%% Multiple Inputs -> Separate Outputs
% If multiple tokens are provided then each token defines one output
% string: this is much faster than calling |DATESTR8601| multiple times.
[y,m,d] = datestr8601([],'yyyy','mmmm','d*')
%% Multiple Inputs -> Concatenated Output
% If more tokens are provided than outputs, then the trailing "outputs" are
% concatentated together into one string (with a space between each string):
[y,md] = datestr8601([],'yyyy','mmmm','d*')
%% Multiple Dates
% By design |DATESTR8601| accepts only *one* input date. To parse multiple
% dates simply call |DATESTR8601| using the inbuilt
% <https://www.mathworks.com/help/matlab/ref/arrayfun.html ARRAYFUN> or
% <https://www.mathworks.com/help/matlab/ref/cellfun.html CELLFUN> as appropriate:
arrayfun(@datestr8601,fix(now+(0:2)),'Uni',0)
cellfun(@(v)datestr8601(v,'*YWD'),{[2019,12,30],[2021,1,3]},'Uni',0) % first & last days of week-numbering year.
%% Bonus: ISO 8601 Calendar and Weekday Functions
% Included functions that support other ISO 8601 date functionalities,
% using ISO 8601 weekday order (with Monday as the first day of the week):
[day,name] = weekday8601(now,'long')
calendar8601()
%% Bonus: DATENUM8601 Function
% The reverse conversion (from ISO 8601 timestamps to serial date
% numbers) is easy to achieve by simply downloading my function
% <https://www.mathworks.com/matlabcentral/fileexchange/39389 DATENUM8601>:
V(1) = datenum8601('1999-01-03T15:06:48.0568'); % calendar
V(2) = datenum8601(  '1999-003T15:06:48.0568'); % ordinal
V(3) = datenum8601('1998-W53-7T15:06:48.0568'); % week
datevec(V(:))