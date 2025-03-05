function [time] = weektow2time(week, sow, sys)

% SYNTAX:
%   [time] = weektow2time(week, sow, sys);
%
% INPUT:
%   week = GNSS week
%   sow  = GNSS seconds-of-week
%   sys  = GNSS system identifier
%
% OUTPUT:
%   time = GPS time (continuous since 6-1-1980)
%
% DESCRIPTION:
%   Conversion from GNSS time in week, seconds-of-week format to GPS time
%   in continuous format (similar to datenum).

%% ------------------------------------------------------------------------


time = week*7*86400 + sow;

BDS_mask = (sys == 'C');
if (any(BDS_mask))
    GPS_BDS_week = 1356; %GPS week number on 1st January 2006 (start of BeiDou time)
    time(BDS_mask) = (GPS_BDS_week + week(BDS_mask))*7*86400 + sow(BDS_mask);
end
