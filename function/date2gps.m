function [gps_week, gps_sow, gps_dow] = date2gps(date)

% SYNTAX:
%   [gps_week, gps_sow, gps_dow] = date2gps(date);
%
% INPUT:
%   date = date [year, month, day, hour, min, sec]
%
% OUTPUT:
%   gps_week = GPS week
%   gps_sow  = GPS seconds of week
%   gps_dow  = GPS day of week
%
% DESCRIPTION:
%   Conversion from calendar date to GPS time.
%% ------------------------------------------------------------------------


gps_start_datenum = 723186; %This is datenum([1980,1,6,0,0,0])

%number of days since the beginning of GPS time
deltat   = (datenum([date(:,1), date(:,2), date(:,3)]) - gps_start_datenum);

gps_week = floor(deltat/7);            %GPS week
gps_dow  = floor(deltat - gps_week*7); %GPS day of week
gps_sow  = (deltat - gps_week*7)*86400; 
gps_sow = gps_sow + date(:,4)*3600 + date(:,5)*60 + date(:,6); %GPS seconds of week

% %alternative way, using the Julian day
% jd = date2jd(date);
% [gps_week, gps_dow, gps_sow] = jd2gps(jd);
