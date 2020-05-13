function [YYYY100,occurrence] = hundredcount(date,TMAX)
%%HUNDREDCOUNT calculates monthly and yearly total number of days with TMAX
%%higher than 100 F (37.8 Celsius).

idx100 = double(TMAX >= 37.8);
[vectorDate] = datevec(date);
[YYYY100] = grpstats(idx100,vectorDate(:,1),'nansum');
[MaxVal,idxVal] = max(YYYY100);
occurrence = vectorDate(idxVal,1);
YYYY100 = [nanmedian(YYYY100) MaxVal];


