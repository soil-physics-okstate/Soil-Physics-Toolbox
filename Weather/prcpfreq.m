function [statDays, vecdiff] = prcpfreq(date,prcp)
%%PRCPFREQ calculates the number of days between two consecutive
%%precipitation events.
%
% Inputs
%       date = Date in Matlab's numeric format
%
%       prcp = Precipitation in any unit. Values lower than zero are
%              converted into NaN. Cannot contain strings, such as 'Trace'.

prcp(prcp < 0) = NaN;
prcpevents = prcp > 0;
dateevents = date(prcpevents);
vecdiffdays = diff(dateevents);
vecdiffstorm = vecdiffdays(vecdiffdays>1);
statDays = [nanmedian(vecdiffstorm) nanmedian(vecdiffdays) nanmax(vecdiffdays)];



