function [S]=summarizeDrainage(HC,mesoNum)
%This function reads in the HC data provided by the soilDrainage function
%and creates annual drainage summaries.

%MUST HAVE COMPLETE YEARS OF DATA
%Will work with any # of years, but must have full data for each year

%determine size of the array
[row,col]=size(HC);

%determine which year it is
startyear=mesoNum(1,1);
endyear = mesoNum(row,1);

S = nan(endyear-startyear+1,3); %set up summary array
S(:,1) = startyear:endyear; %years 
for i = startyear:endyear
    j = mesoNum(:,1)==i ;%logical index for year i
    S(i-startyear+1,2) = nansum(HC(j)); %nansum of HC for that year
    S(i-startyear+1,3) = sum(isnan(HC(j))); %count of nans **verify
end
% 
% 
% %determine how many years/iterations of for loop collecting data
% k=mesoNum(row,1)-(year-1);
% 
% %get a year's worth of data
% x=1;
% %S(k,2)=nan;
% for n=1:k
%     %is the current year a leap year?
%     leapyear=(rem(year,4)==0);
%     
%     if(leapyear)
%         S(n,1)=year;
%         year=year+1;
%         S(n,2)=nansum(HC(x:(x+365)));
%         x=x+366;
%     else
%         S(n,1)=year;
%         year=year+1;
%         S(n,2)=nansum(HC(x:(x+364)));
%         x=x+365;
%     end
% end
