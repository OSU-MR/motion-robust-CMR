function varargout = datestr8601(dtx,varargin)
% Convert a date vector/number/datetime to ISO 8601 formatted date strings (timestamp).
%
% (c) 2011-2020 Stephen Cobeldick
%
%%% Syntax:
%  str = datestr8601
%  str = datestr8601(dtx)
%  str = datestr8601(dtx,tkn)
%  str = datestr8601(dtx,tkn1,tkn2,...,tknN)
%  [str1,str2,...,strM] = datestr8601(dtx,tkn1,tkn2,...,tknN)
%  [...] = datestr8601([],...)
%
% Convert a date vector, serial date number, or datetime to date string/s,
% with date formats controlled by (optional) input token/s. Each token
% defines a date string with either an ISO 8601 timestamp or a single date/
% time value. If there are more input tokens than outputs, the remaining
% strings are concatenated together with a space character between each string.
%
% The ISO 8601 timestamp style options are:
% * Date in calendar, ordinal, or week-numbering notation.
% * Basic or extended format.
% * Choice of date-time separator character (one of 'T @_').
% * Full or lower precision (trailing units omitted)
% * Decimal fraction of the trailing unit.
% These style options are explained in the tables below.
%
% For date values the token's case determines the output date string's
% year-type: lowercase = calendar year, UPPERCASE = week-numbering year.
%
% Note 0: Timezones are not handled or included in the output.
% Note 1: Decimal fractions are truncated to the requested precision.
% Note 2: Some date strings use the ISO 8601 week-numbering year, where the first
%  week of the year includes the first Thursday of the year: please double check!
% Note 3: Out-of-range values are permitted in the input date vector.
%
%% Examples
%
%%% Using the date+time given by date vector [1999,1,3,15,6,48.0568].
%
% >> datestr8601()
% ans = '19990103T150648'
%
% >> datestr8601([],'yn_HM')
% ans = '1999003_1506'
%
% >> datestr8601(clock,'*ymdHMS')
% ans = '1999-01-03T15:06:48'
%
% >> [D1,D3] = datestr8601(now-2,'D','DDD')
% D1 = '5'
% D3 = 'Fri'
%
% >> datestr8601(datetime,'DDDD','d*','mmmm','yyyy')
% ans = 'Sunday 3rd January 1999'
%
% >> [da,YWD,mmyy] = datestr8601([],'d*','*YWD','mmmm','yyyy');
% >> sprintf('The %s of %s has the ISO week-date "%s".',da,mmyy,YWD)
% ans = 'The 3rd of January 1999 has the ISO week-date "1998-W53-7".'
%
%% ISO 8601 Timestamps
%
% A token consists of one letter for each of the consecutive date/time
% units in the timestamp, thus it defines the date notation (calendar,
% ordinal, or week-date) and selects either basic or extended format:
%
% Output   | Basic Format             | Extended Format (token prefix '*')
% Date     | Input  | Output Timestamp| Input   | Output Timestamp
% Notation:| <tkn>: | <str> Example:  | <tkn>:  | <str> Example:
% =========|========|=================|=========|==========================
% Calendar |'ymdHMS'|'19990103T150648'|'*ymdHMS'|'1999-01-03T15:06:48'
% ---------|--------|-----------------|---------|--------------------------
% Ordinal  |'ynHMS' |'1999003T150648' |'*ynHMS' |'1999-003T15:06:48'
% ---------|--------|-----------------|---------|--------------------------
% Week     |'YWDHMS'|'1998W537T150648'|'*YWDHMS'|'1998-W53-7T15:06:48'
% ---------|--------|-----------------|---------|--------------------------
%
% Options for reduced precision timestamps, non-standard date-time separator
% character, and the addition of a decimal fraction of the trailing unit:
%
% Omit leading and/or trailing units (reduced precision), e.g.:
% =========|========|=================|=========|==========================
%          |'DHMS'  |'7T150648'       |'*DHMS'  |'7T15:06:48'
% ---------|--------|-----------------|---------|--------------------------
%          |'mdH'   |'0103T15'        |'*mdH'   |'01-03T15'
% ---------|--------|-----------------|---------|--------------------------
% Select the date-time separator character (one of 'T',' ','@','_'), e.g.:
% =========|========|=================|=========|==========================
%          |'n_HMS' |'003_150648'     |'*n_HMS' |'003_15:06:48'
% ---------|--------|-----------------|---------|--------------------------
%          |'YWD@H' |'1998W537@15'    |'*YWD@H' |'1998-W53-7@15'
% ---------|--------|-----------------|---------|--------------------------
% Decimal fraction of the trailing date/time value by specifying digits, e.g.:
% =========|========|=================|=========|==========================
%          |'HMS4'  |'150648.0568'    |'*HMS4'  |'15:06:48.0568'
% ---------|--------|-----------------|---------|--------------------------
%          |'YW7'   |'1998W53.9471033'|'*YW7'   |'1998-W53.9471033'
% ---------|--------|-----------------|---------|--------------------------
%          |'y10'   |'1999.0072047202'|'*y10'   |'1999.0072047202'
% ---------|--------|-----------------|---------|--------------------------
%
% Note 4: The code matches single-value tokens before ISO 8601 timestamp tokens.
% Note 5: This function does not check for ISO 8601 compliance: user beware!
% Note 6: Date-time separator character must be one of 'T',' ','@','_'.
% Note 7: Date notations cannot be combined: note upper/lower case characters!
%
%% Single-Value Tokens
%
% 'q' = each quarter is three months long: Jan-Mar, Apr-Jun, Jul-Sep, Oct-Dec.
% 'Q' = each quarter is 13 weeks long (the last one may be 14).
%
% 'n'+'r' = 365 (or 366 if a leap year).
% 'N'+'R' = 7*52 or 7*53 (year dependent).
% 'W'+'X' = 52 or 53 (year dependent).
%
% Input | Output                                      | Output
% <tkn>:| <str> Date/Time Representation:             | <str> Example:
% ======|=============================================|===============
%%% Calendar Year                  (lowercase tokens) |
% 'yyyy'| year, four digit                            |'1999'
% 'n'   | day of the year, variable digits            |'3'
% 'n*'  | day of the year, with ordinal suffix        |'3rd'
% 'nnn' | day of the year, three digit, zero padded   |'003'
% 'r'   | days remaining in year, variable digits     |'362'
% 'r*'  | days remaining in year, with ordinal suffix |'362nd'
% 'rrr' | days remaining in year, three digit, padded |'362'
% 'q'   | year quarter, 3-month                       |'1'
% 'q*'  | year quarter, 3-month, with ordinal suffix  |'1st'
% 'qq'  | year quarter, 3-month, abbreviation         |'Q1'
% 'qqqq'| year quarter, 3-month, full ordinal word    |'First'
% 'm'   | month of the year, variable digits          |'1'
% 'm*'  | month of the year, with ordinal suffix      |'1st'
% 'mm'  | month of the year, two digit, zero padded   |'01'
% 'mmm' | month name, three letter abbreviation       |'Jan'
% 'mmmm'| month name, in full                         |'January'
% 'd'   | day of the month, variable digits           |'3'
% 'd*'  | day of the month, with ordinal suffix       |'3rd'
% 'dd'  | day of the month, two digit, zero padded    |'03'
% ------|---------------------------------------------|---------------
%%% Week-Numbering Year            (uppercase tokens) |
% 'YYYY'| year, four digit                            |'1998'
% 'N'   | day of the year, variable digits            |'371'
% 'N*'  | day of the year, with ordinal suffix        |'371st'
% 'NNN' | day of the year, three digit, zero padded   |'371'
% 'R'   | days remaining in year, variable digits     |'0'
% 'R*'  | days remaining in year, with ordinal suffix |'0th'
% 'RRR' | days remaining in year, three digit, padded |'000'
% 'Q'   | year quarter, 13-week                       |'4'
% 'Q*'  | year quarter, 13-week, with ordinal suffix  |'4th'
% 'QQ'  | year quarter, 13-week, abbreviation         |'Q4'
% 'QQQQ'| year quarter, 13-week, full ordinal word    |'Fourth'
% 'W'   | week of the year, variable digits           |'53'
% 'W*'  | week of the year, with ordinal suffix       |'53rd'
% 'WW'  | week of the year, two digit, zero padded    |'53'
% 'X'   | weeks remaining, variable digits            |'0'
% 'X*'  | weeks remaining, with ordinal suffix        |'0th'
% 'XX'  | weeks remaining, two digit, zero padded     |'00'
% ------|---------------------------------------------|---------------
%%% Weekday                        (uppercase tokens) |
% 'D'   | weekday number (Monday=1)                   |'7'
% 'D*'  | weekday number (Monday=1), ordinal suffix   |'7th'
% 'DD'  | weekday name, two letter abbreviation       |'Su'
% 'DDD' | weekday name, three letter abbreviation     |'Sun'
% 'DDDD'| weekday name, in full                       |'Sunday'
% ------|---------------------------------------------|---------------
%%% Time of Day                    (uppercase tokens) |
% 'H'   | hour of the day, variable digits            |'15'
% 'HH'  | hour of the day, two digit, zero padded     |'15'
% 'M'   | minute of the hour, variable digits         |'6'
% 'MM'  | minute of the hour, two digit, zero padded  |'06'
% 'S'   | second of the minute, variable digits       |'48'
% 'SS'  | second of the minute, two digit, zero padded|'48'
% 'F'   | deci-second of the second, zero padded      |'0'
% 'FF'  | centisecond of the second, zero padded      |'05'
% 'FFF' | millisecond of the second, zero padded      |'056'
% ------|---------------------------------------------|---------------
% 'MAMP'| Midnight/Ante-Meridian/Midday/Post-Meridian |'Post-Meridian'
% ------|---------------------------------------------|---------------
%
%% Input & Output Arguments
%
% Inputs:
%  dtx = datetime scalar.
%      = Date vector, [year,month,day] or [year,month,day,hour,minute,second].
%      = Serial date number, where 1 == start of 1st January of the year 0.
%      = []*, uses the current time (default).
%  tkn = String token (1xN char), chosen from the above tables (default is 'ymdHMS').
%
% Outputs:
%  str = Date string (1xN char), with date/time specified by <tkn>.
%
% See also DATENUM8601 DATEROUND CLOCK NOW DATESTR DATENUM DATEVEC DATETIME NATSORTFILES

%% Input Wrangling
%
rwf = @(szv) numel(szv)==2 && szv(1)==1;
%
% Calculate date-vector:
if nargin==0 % default:
	dtv = clock;
elseif isnumeric(dtx)
	if isempty(dtx)&&numel(size(dtx))==2 % default:
		dtv = clock;
	elseif isscalar(dtx) % datenum:
		dtv = datevecmx(dtx);
	else % datevec:
		dtn = datenummx(dtx);
		dtv = datevecmx(dtn);
	end
else % datestr & datetime:
	dtv = datevec(dtx);
end
% Check if the date vector constitutes exactly one row:
assert(rwf(size(dtv)),...
	'SC:datestr8601:dtx:NotValidDate',...
	'First input <dtx> must be a date vector, serial date number, or datetime.')
% Generate serial date number:
dtn = datenummx(dtv);
% Weekday index (Mon==1):
dtd = 1+mod(floor(dtn(1))-3,7);
% Adjust date to suit week-numbering:
dtn(2,1) = dtn(1)+4-dtd;
dtv(2,:) = datevecmx(floor(dtn(2)));
dtv(2,4:6) = dtv(1,4:6);
% Separate fraction of seconds from seconds:
dtv(:,7) = rem(dtv(1,6),1);
dtv(:,6) = floor(dtv(1,6));
% Beginning of [this,next] year:
dtb(1,:) = datenummx(dtv(1)+[0;1],1,1);
dtb(2,:) = datenummx(dtv(2)+[0;1],1,1);
dto = 3-mod(dtb(2,:),7);
dtb(2,:) = dtb(2,:)+dto;
%
varargin(1+numel(varargin):1) = {'ymdHMS'};
varargout = varargin; % preallocate output.
%
%% Create String for Each Input Token
%
raw = [...
	'00000000001111111111222222222233333333334444444444555555555566666';...
	'01234567890123456789012345678901234567890123456789012345678901234';...
	'tsnrtttttttttttttttttsnrtttttttsnrtttttttsnrtttttttsnrtttttttsnrt';...
	'htddhhhhhhhhhhhhhhhhhtddhhhhhhhtddhhhhhhhtddhhhhhhhtddhhhhhhhtddh'].';
fun = @(x) raw(1+rem(x,10)+10*any(rem(x,100)==11:13),3:4);
%
for m = 1:numel(varargin)
	% Ordinal suffix of the input token:
	sfx = sprintf('%d%s',m+1,fun(m+1));
	% Check if the input token constitutes exactly one row:
	tkn = varargin{m};
	assert(ischar(tkn)&&rwf(size(tkn)),...
		'SC:datestr8601:tkn:NotCharRowVector',...
		'%s input must be a token (1xN character).',sfx)
	tkl = numel(tkn);
	iso = strcmp(tkn(tkl),'*');
	isu = strcmp(upper(tkn),tkn);
	nrv = abs(1+floor(dtn(1))-dtb(1+isu,:));
	tkl = tkl-iso;
	switch tkn
		case {'DDDD','DDD','DD','D','D*'}
			% weekday:
			varargout{m} = ds8601Day(tkl,iso,dtd);
		case {'FFF','FF','F'}
			% deci/centi/milliseconds:
			tmp = sprintf('%#.4f',dtv(1,7));
			varargout{m} = tmp(3:2+tkl);
		case {'mm','m','m*','dd','d','d*','HH','H','MM','M','SS','S'}
			% month of the year, day of the month, hours, minutes, seconds:
			val = dtv(1,strfind('ymdHMS',tkn(1)));
			varargout{m} = raw(1+val,1+min(2-tkl,val<10):2+2*iso); %-> W.
		case {'mmmm','mmm'}
			% month of the year:
			varargout{m} = ds8601Mon(tkl,dtv(1,2));
		case {'nnn','n','n*','rrr','r','r*','NNN','N','N*','RRR','R','R*'}
			% day of the year, days remaining in the year:
			isr = strncmpi('R',tkn,1);
			varargout{m} = ds8601DoY(tkl,iso,fun,nrv(1+isr));
		case {'WW','W','W*','XX','X*','X'}
			% week of the year, weeks remaining in the year:
			isx = strncmpi('X',tkn,1);
			val = ceil((nrv(1+isx)+isx)/7)-isx;
			varargout{m} = raw(1+val,1+min(2-tkl,val<10):2+2*iso); %-> m/d/H/M/S.
		case {'qqqq','qq','q','q*','QQQQ','QQ','Q','Q*'}
			% year quarter:
			val = min(4,ceil([dtv(1,2)/3,nrv(1)/91]));
			varargout{m} = ds8601Qtr(tkl,iso,val(1+isu));
		case {'yyyy','y','YYYY','Y'}
			% year:
			varargout{m} = sprintf('%04d',dtv(1+isu));
		case 'MAMP'
			% midnight/am/noon/pm:
			val = 2+2*(dtv(1,4)>=12)-(all(dtv(1,5:7)==0)&&any(dtv(1,4)==[0,12]));
			amp = {'Midnight','Ante-Meridian','Midday','Post-Meridian'};
			varargout{m} = amp{val};
		otherwise % ISO 8601 timestamp:
			varargout{m} = ds8601ISO(tkn,sfx,dtv,dtn,dtb,dtd,nrv,raw);
	end
end
%
%% Concatenate Trailing Outputs
%
if nargin-1 > nargout
	idx = max(1,nargout);
	tmp = varargout(idx:end);
	tmp(2,:) = {' '};
	varargout{idx} = [tmp{1:end-1}];
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%datestr8601
function out = ds8601Day(tkl,iso,val)
% weekday
%
if tkl==1
	tmp = {'1st','2nd','3rd','4th','5th','6th','7th'};
	out = tmp{val}(1:1+2*iso);
else
	tmp = {'Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'};
	out = tmp{val}(1:max(tkl,end*(tkl>3))); %-> m.
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ds8601Day
function out = ds8601Mon(tkl,val)
% month
%
tmp = {'January','February','March','April','May','June','July','August','September','October','November','December'};
out = tmp{val}(1:max(tkl,end*(tkl>3))); %-> D.
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ds8601Mon
function out = ds8601Qtr(tkl,iso,val)
% year quarter
%
if tkl<3
	tmp = {'Q1st','Q2nd','Q3rd','Q4th'};
	out = tmp{val}(1+abs(tkl-2):max(2,4*iso));
else
	tmp = {'First','Second','Third','Fourth'};
	out = tmp{val};
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ds8601Qtr
function out = ds8601DoY(tkl,iso,fun,val)
% day of the year, days remaining in the year
%
if iso
	out = sprintf('%0*d%s',tkl,val,fun(val));
else
	out = sprintf('%0*d',tkl,val);
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ds8601DoY
function out = ds8601ISO(tkn,sfx,dtv,dtn,dtb,dtd,nrv,raw)
% ISO 8601 timestamp
%
% Identify format, date, separator, time, and digit characters:
isx = strncmp(tkn,'*',1);
fmt = '^(y?(m?d?|n?)|Y?W?D?)([T @_]?)(H?M?S?)(\d*)$';
tkr = regexp(tkn(1+isx:end),fmt,'tokens','once');
assert(~isempty(tkr),...
	'SC:datestr8601:tkn:NotValidToken',...
	'%s input (token) is not recognized: ''%s''. Check the help!',sfx,tkn)
assert(numel(tkr)==4,...
	'SC:datestr8601:tkn:BuggyOctaveRegexp',...
	'It appears that you are using Octave with buggy REGEXP.')
assert(isempty(tkr{2})||~isempty(tkr{3}),...
	'SC:datestr8601:tkn:TrailingSeparator',...
	'%s input (token) contains a trailing date-time separator: ''%s''',sfx,tkn)
% Identify timestamp:
stf = strfind({'ymdHMS','ynHMS','YWDHMS'},[tkr{1},tkr{3}]);
ist = find(~[isempty(stf{1}),isempty(stf{2}),isempty(stf{3})]);
assert(~isempty(ist),...
	'SC:datestr8601:tkn:UnitsOutOfSequence',...
	'%s input (token) is out of sequence or is missing units: ''%s''',sfx,tkn)
isw = ist(1)==3;
%
idb = stf{ist(1)};
cnt = numel(tkr{1})+numel(tkr{3});
idc = {1:6,[1,3:6],1:6};
idu = idc{ist(1)}(idb:idb+cnt-1);
% For calculating decimal fraction of date/time values:
idz = idu(end);
idw = false;
dtw = dtv(1+isw,:);
dtz = [1,1,1,0,0,0,0];
%
sep = [tkr{2},'T'];
sep = sep(1);
if isx % Extended-format
	dtc = {'','-','-',sep,':',':';'','','','','',''};
else % Basic-format
	dtc = {'', '', '',sep, '', '';'','','','','',''};
end
%
% hours, minutes, seconds:
for m = 4:max(idu)
	dtc{2,m} = raw(1+dtw(m),1:2);
end
%
switch ist(1)
	case 1 % Calendar
		% month, day of the month:
		for m = max(2,min(idu)):3
			dtc{2,m} = raw(1+dtw(m),1:2);
		end
	case 2 % Ordinal
		% day of the year:
		dtc{2,3} = sprintf('%03d',nrv(1));
	case 3 % Week-numbering
		% Decimal fraction of weeks, not days:
		if idz==2
			idw = true;
			dtz(3) = dtw(3)-dtd+1;
		end
		% weekday:
		if any(idu==3)
			dtc{2,3} = raw(1+dtd,2);
		end
		% week of the year:
		if any(idu==2)
			dtc{2,2} = ['W',raw(1+ceil(nrv(1)/7),1:2)];
		end
end
%
if idu(1)==1
	% year:
	dtc{2,1} = sprintf('%04d',dtw(1));
end
%
% Concatenate separator and value strings:
tmp = [idu*2-1;idu*2];
out = [dtc{tmp(2:end)}];
%
% Decimal fraction of trailing unit (decimal places):
dgt = sscanf(tkr{4},'%d');
if dgt>0
	if idz==1 % year
		dcp = 12;
		frc = (dtn(1)-dtb(isw+1,1)) ./ diff(dtb(isw+1,:));
	elseif idz==3 % day
		dcp = 10;
		frc = mod(dtn(1),1);
	elseif idz==6 % second
		dcp = 4;
		frc = dtw(7);
	elseif any(dtw(idz+1:7)>dtz(idz+1:7)) % month/week/hour/minute
		dcp = 15;
		% Floor all trailing units:
		dtw(idz+1:7) = dtz(idz+1:7);
		dtf = datenummx(dtw(1:6));
		% Increment the chosen unit:
		dtw(idz+idw) = dtw(idz+idw)+1+6*idw;
		% Decimal fraction of the chosen unit:
		frc = (dtn(isw+1)-dtf) ./ (datenummx(dtw(1:6))-dtf);
	else
		dcp = dgt;
		frc = 0;
	end
	str = sprintf('%#.*f',dcp,frc);
	str(3+dcp:2+dgt) = '0';
	out = [out,str(2:2+dgt)];
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ds8601ISO