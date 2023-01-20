%--------------------------------------------------------------------------
%
%  Inputs:
%    year       year
%    mon        month
%    day        day
%    hour       universal time hour
%    minute     universal time minutes
%    sec        universal time seconds
%
%  Output:
%    Mjd        Modified julian date (UTC)
%
% Reference:
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
% New York; 3rd edition(2007).
%
%--------------------------------------------------------------------------
function Mjd = Mjday(year, mon, day, hour, minute, sec)

if (nargin < 4)
    hour = 0;
    minute = 0;
    sec = 0;
end

jd = 367 * year...
    - floor( (7 * (year + floor( (mon + 9) / 12) ) ) * 0.25 )...
    + floor( 275 * mon / 9 )...
    + day + 1721013.5...
    + ( (sec/60 + minute ) / 60 + hour ) / 24;

Mjd = jd - 2400000.5;