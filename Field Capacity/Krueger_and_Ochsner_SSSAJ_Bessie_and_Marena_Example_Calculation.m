%Authors' names omitted
%12/5/2023

%This file loads, processes, and graphs soil moisture and precipitation 
%data from the Bessie and Marena Oklahoma Mesonet sites.  The included 
%function (peakdet_SSSAJ) identifies valid peaks in the soil moisture 
%timeseries, which are also plotted.  Data are for the 2020-2022 period.

%load soil moisture and precipitation data
Bessie = readtable('Krueger and Ochsner SSSAJ Bessie and Marena Example Data.xlsx','Sheet','Bessie');
Marena = readtable('Krueger and Ochsner SSSAJ Bessie and Marena Example Data.xlsx','Sheet','Marena');

%estimate field capacity
delta = 0.001;

[VWC_peaks_Bessie] = peakdet_SSSAJ(Bessie.VWC05,delta,Bessie.Matlab_Date,Bessie.Rain_mm);
[VWC_peaks_Marena] = peakdet_SSSAJ(Marena.VWC05,delta,Bessie.Matlab_Date,Marena.Rain_mm);

%%
%Plot soil moisture timeseries and identify soil moisture peaks

figure(1)
hold on
box on
set(gcf,'units','centimeters','Position', [10 5 18 16]); %left bottom width height

%Marena
subplot(2,1,1)
set(gca,'units','normalized','Position', [0.085 0.55 0.85 0.42]); %left bottom width height
hold on
box on

    %plot daily soil moisture data
    plot(Marena.Matlab_Date, Marena.VWC05,'k')
    %plot valid peaks
    plot(VWC_peaks_Marena.Matlab_Date, VWC_peaks_Marena.VWC_at_Peak,'marker','o','markeredgecolor','k','markerfacecolor','b','linewidth',0.5,'linestyle','none','MarkerSize',7);

    %add label
    ylabel('\bf VWC (cm^3 cm^-^3)')

    %set x axis limts
    xlim([datenum(2020,1,1) datenum(2023,1,1)])
    %set y axis limts
    ylim([0 0.48])

    %select tick mark locations
    dates_tick = datevec(Marena.Matlab_Date);
    dates_tick = unique(dates_tick(:,1:2),'rows');
    xticks(datenum(dates_tick(:,1),dates_tick(:,2),1))
    xticklabels({'2020','','','','','','','','','','','','2021','','','','','','','','','','','','2022','','','','','','','','','','','','2023'});
    xtickangle(0)
    
    %add site name
    text(datenum(2020,1,15),0.45,'Marena','fontsize',10,'FontWeight','bold')

    legend({'Daily avgerage Volumetric Soil Water Content (\theta)','Valid peak \theta'},'orientation','horizontal','Location','northeast')
    legend boxoff

%Bessie
subplot(2,1,2)
set(gca,'units','normalized','Position', [0.085 0.08 0.85 0.42]); %left bottom width height
hold on
box on
    
    %plot daily soil moisture data
    plot(Bessie.Matlab_Date, Bessie.VWC05,'k')
    %plot valid peaks
    plot(VWC_peaks_Bessie.Matlab_Date, VWC_peaks_Bessie.VWC_at_Peak,'marker','o','markeredgecolor','k','markerfacecolor','b','linewidth',0.5,'linestyle','none','MarkerSize',7);

    %add label
    ylabel('\bf VWC (cm^3 cm^-^3)')

    %set x axis limts
    xlim([datenum(2020,1,1) datenum(2023,1,1)])
    %set y axis limts
    ylim([0 0.48])
    
    %select tick mark locations
    dates_tick = datevec(Bessie.Matlab_Date);
    dates_tick = unique(dates_tick(:,1:2),'rows');
    xticks(datenum(dates_tick(:,1),dates_tick(:,2),1))
    xticklabels({'2020','','','','','','','','','','','','2021','','','','','','','','','','','','2022','','','','','','','','','','','','2023'});
    xtickangle(0)
    
    %add site name
    text(datenum(2020,1,15),0.45,'Bessie','fontsize',10,'FontWeight','bold')

%% Peak Detection Function

function [VWC_at_FC] = peakdet_SSSAJ(vwc, delta, dates, rain)
%PEAKDET Detect peaks and valleys in soil moisture timeseries data
%        This function was adapted from Eli Billauer, version 3.4.05:
%        https://billauer.co.il/blog/2009/01/peakdet-matlab-octave/
%
%        [VWC_at_FC] = PEAKDET(vwc, delta, dates, rain) finds the local
%        maxima ("peaks") in a vector of soil moisture timeseries data.
%      
%        vwc:    daily volumetric water content data
%        delta:  differential required for a peak or valley to be detected.
%        dates:  dates corresponding to each vwc measurement (in Matlab date format).
%        rain:   daily total precipitation data
% 
%   Notes:
%       This algorithm functions by looking back in time until the most 
%       recent maximum. Once VWC falls below the user identified delta, it 
%       identifies the most recent maximum as a local maximum
%       Peaks are considered calid only if they meet the following criteria:
%       1. They are within 95% of th absolute maximum VWC for the dataset.  
%          This prevents identifying days were rainfall was high but the 
%          potential for complete soil water recharge was low.
%       2. They occur during the months of December-February.  This removes
%          datapoints during times when ET is likely to be high
%       3. No rain must occur for ten days after the maxumum.  This
%          allows time for drainage without recharge. 

%% Detect Peaks
maxtab = [];
mintab = [];

%vwc = vwc(:); % Just in case this wasn't a proper vector

if nargin < 3
  dates = (1:length(vwc))';
else
  dates = dates(:);
  if length(vwc)~= length(dates)
    error('Input vectors vwc and dates must have same length');
  end
end

if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta < 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;

for i=1:length(vwc)
  vwc_today = vwc(i);     
  
  if vwc_today > mx       %if VWC today is greater than -infinity on first itteration, thereafter if VWC today is greater than yesterday.
      mx = vwc_today;     %then identify VWC today as most recent maximum
      mxpos = dates(i);
  end
  if vwc_today < mn       %if VWC today is less than infinity
      mn = vwc_today;     %then identify VWC today as most recent minimum
      mnpos = dates(i); 
  end

  if lookformax
    if vwc_today < mx-delta                 %if VWC today is less than VWC of most recent maximum - delta
      maxtab = [maxtab ; mxpos mx];         %then identify most recent maximum as local maximum
      mn = vwc_today; mnpos = dates(i);     %then identify VWC today as most recent minimum
      lookformax = 0;                       %then start to look for local minimum
    end
  else
    if vwc_today > mn+delta                 %if VWC today is more than VWC of most recent minimum + delta
      mintab = [mintab ; mnpos mn];         %then identify most recent minimum as local minimum
      mx = vwc_today; mxpos = dates(i);     %then identify VWC today as most recent maximum
      lookformax = 1;                       %then start to look for local maximum
    end
  end
end

%Save max and min values before and after applying criteria
    maxtab_pre_criteria = maxtab;
    mintab_pre_criteria = mintab;
    maxtab_post_criteria = maxtab;
    mintab_post_criteria = mintab;

%% Determine if Peaks are Valid for Field Capacity Estimation
%Select only those maximum values that meet criteria

    %criterion_max: peak within 95% of maximum for the dataset
    criterion_max = 0.95 * max((vwc));
    if ~isempty(maxtab_post_criteria)
        keep = (maxtab_post_criteria(:,2) >= criterion_max);
        maxtab_post_criteria = maxtab_post_criteria(keep,:);
    else
        maxtab_post_criteria = [NaN NaN];  %assign NaN if no peaks meet criteria (e.g., no soil moisture data for that site and depth)
    end
    
    %criterion_dorm: peak occurs during the dormant season
    criterion_dorm = [1 2 12]';     %months of January, February, or December
    [~, temp_month, ~, ~, ~, ~] = datevec(maxtab_post_criteria(:,1));
    keep = ismember(temp_month, criterion_dorm);
    maxtab_post_criteria = maxtab_post_criteria(keep,:);
    clear keep temp_month

    %criterion_rain:
    %compile rain data
        for i = 1:size(dates,1)          
            keep = ismember(dates(i), maxtab_post_criteria);
            if sum(keep) == 1
                temp_rain(i,1) = sum(rain(i+1:i+10,1)); %cumulative rain for the 10 days after maximum
            end    
        end
    %apply criterion
    keep = ismember(dates, maxtab_post_criteria);
    if sum(keep) > 0;
        temp_rain = temp_rain(keep,:);
        criterion_rain = temp_rain(:,1) < 10;
        if sum(criterion_rain) > 0;% if at least one peak passes this criterion test
            maxtab_post_criteria = maxtab_post_criteria(criterion_rain,:);
        else
            maxtab_post_criteria = NaN;
        end
    else
        maxtab_post_criteria = NaN;
    end
   
    %Retrieve VWC one day after the peak (assigned as field capacity)
    if ~isnan(maxtab_post_criteria)
        for i = 1:size(maxtab_post_criteria,1);   
            keep = dates == maxtab_post_criteria(i,1);
            %get VWC one day after peak
            temp_date = dates(keep) + 1; 
            keep = ismember(dates,temp_date);     
            if sum(keep) == 1
                VWC_at_FC(i,3) = vwc(keep,1); 
            else
                VWC_at_FC(i,3) = NaN;
            end
        end
    end
    %Assign NaN if no peaks were detected
    if exist('VWC_at_FC','var') == 0
        VWC_at_FC = [NaN];
    end

    %add dates
    VWC_at_FC(:,1) = maxtab_post_criteria(:,1);
    %add VWC at peak
    VWC_at_FC(:,2) = maxtab_post_criteria(:,2);   

    VWC_at_FC = array2table(VWC_at_FC,'VariableNames',{'Matlab_Date','VWC_at_Peak','VWC_at_FC'});

end


    

