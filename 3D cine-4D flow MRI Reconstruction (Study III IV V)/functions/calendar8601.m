function [out,y,m] = calendar8601(y,varargin)
% Displays or returns a 6*7 matrix of the calendar of the current or specified month.
%
% (c) 2016-2020 Stephen Cobeldick
%
% CALENDAR8601 is a drop-in replacement for MATLAB's standard CALENDAR
% function, but has the week starting on Monday, as per ISO 8601. This
% means the first column of the returned matrix corresponds to Monday.
%
% If no inputs: uses the current month.
% If no output: displays the calendar and marks the current date with '-'.
%
%%% Syntax:
%  calendar8601()
%  calendar8601(datetime)
%  calendar8601(year,month)
%  calendar8601(dateNumber)
%  calendar8601(dateString)
%  calendar8601(dateString,<extra arguments passed to DATEVEC>)
%  calendar8601(dateVector)
%  [out,year,month] = calendar8601(...)
%
%% Input and Output Arguments
%
%%% Inputs:
%  y = Datetime, a scalar datetime variable.
%    = NumericScalar, a valid Serial Date Number.
%    = CharacterRow,  a valid Date String.
%    = NumericVector, a valid Date Vector.
%    = NumericScalar, the year containing the desired month.
%  m = NumericScalar, the month of the year (1=January, 12=December).
%
%%% Outputs:
%  out = NumericMatrix of dates, each row corresponds to one week. Zero padded.
%  y   = NumericScalar, the year.
%  m   = NumericScalar, the month.
%
% See also WEEKDAY8601 DATENUM8601 DATESTR8601 CALENDAR DATENUM DATESTR DATEVEC DATETIME DATEROUND

%% Input Wrangling
%
dtv = clock();
erf = @(a,dn,up)isnumeric(a)&&isscalar(a)&&isreal(a)&&fix(a)==a&&(dn<a)&&(a<up);
er2 = 'Second input <%s> can be a real scalar numeric (the month).';
er1 = ['First input <%s> can be a real scalar numeric (the year),\n',...
	'a serial date number, a date vector, a date string, or a scalar datetime.'];
%
if nargin==0 % no inputs:
	y = dtv(1);
	m = dtv(2);
elseif nargin==1 && isnumeric(y)
	if numel(y)>1 % datevec:
		y = datenummx(y);
		assert(isscalar(y),'SC:calendar8601:y:NotDateVector',er1,'y')
	end % datenum:
	[y,m,~] = datevecmx(y);
elseif nargin==2 && isnumeric(y) % year & month:
	m = varargin{1};
	assert(erf(y,-Inf,Inf),'SC:calendar8601:y:NotScalarNumeric',er1,'y')
	assert(erf(m,   0, 13),'SC:calendar8601:m:NotScalarNumeric',er2,'m')
else % datestr | datetime:
	[y,m,~] = datevec(y,varargin{:});
	assert(isscalar(y),'SC:calendar8601:y:NotDatestrDatetime',er1,'y')
	assert(isscalar(m),'SC:calendar8601:y:NotDatestrDatetime',er2,'m')
end
%
% Generate Calendar
%
lpy = (mod(y,4)==0&&mod(y,100))||~mod(y,400);
dpm = [31,28+lpy,31,30,31,30,31,31,30,31,30,31];
dom = 1:dpm(m);
%
mat = zeros(7,6);
mat(dom+mod(datenummx(y,m,1)+4,7)) = dom;
%
mnc = {'January','February','March','April','May','June','July','August','September','October','November','December'};
%
if nargout
	out = mat.';
else
	ism = y==dtv(1) && m==dtv(2);
	isd = ism & mat==dtv(3);
	mat(isd) = -mat(isd);
	arr = permute(cat(3,mat,32+13*isd),[3,1,2]);
	fprintf('%*s %04d\n',18+ceil(numel(mnc{m})/2),mnc{m},y)
	fprintf('   Mon   Tue   Wed   Thu   Fri   Sat   Sun\n')
	fprintf(' %5d%c%5d%c%5d%c%5d%c%5d%c%5d%c%5d%c\n',arr)
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calendar8601