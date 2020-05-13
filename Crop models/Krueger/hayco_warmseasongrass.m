function [Kcb, ET0, P, days, date_end_mid] = grassco_warmseason(start_grass, end_grass, growing_grass,dormant_grass,soil,theta_ini,weather,refET,pctcov,method,temp_Kcbmid)
% grassco_warmseason constructs dual crop coefficient curves for grasslands. These crop coefficient 
% curves are for the purpose of scaling reference ET to produce actual ET estimates. 
%
% v.1 ESK Jan. 21, 200
% Code modified for grasslands 25 Oct. 2016 by BMW
% Code modified (from grassco) for warm season grasses 13 Sept. 2018 ESK
%
% growing_grass = beginning of growing season date (Matlab date number format)
% dormant_grass = end of growing season (Matlab date number format)
% soil = [sand clay OM] percent sand, clay, and organic matter in soil
% theta_ini = initial volumetric water content of the surface soil
% weather    = [date doy Tmax Tmin P RHmax RHmin Rs wind]
%        date  = date in Matlab date number format
%        doy   = day of the year
%        Tmax  = maximum daily air temp (Celcius)
%        Tmin  = minimum daily air temp (Celcius)
%        P     = daily precip (mm)
%        RHmax = daily maximum relative humidity (%)
%        RHmin = daily minimum relative humidity (%)
%        Rs    = incoming solar radiation (MJ/m2/day)
%        wind  = wind speed (m/s)
%        **this array must cover the entire extent of the growing season
% refET = daily reference evapotranspiration (mm)
% pctcov = percent residue cover during initial period
%
% Kcini is calculated following:
% Allen, R.G., W.O. Pruitt, D. Raes, M. Smith, and L.S. Pereira. 2005. 
% Estimating evapotranspiration from bare soil and the crop coefficient for
% the initial period using common soils information. J. Irrig. Drain. Eng. 131:14-23.
%
% Kcmid and Kcend based on:
% Allen, R.G., L.S. Pereira, D. Raes, and M. Smith. 1998. Crop evapotranspiration:
% Guidelines for computing crop water requirements, FAO Irrigation and Drainage Paper No. 56. Rome.
%
% Kc for grass during the growing season is calculated
% following Kcini procedure above


%% Calculate the time points defining the crop coeffecient curve

h = 0.75;        % mean of maximum crop heights for forages excluding turf grass (m), Table 12, FAO-56
    
days = weather(:,1); %days in the crop coefficient period
start = start_grass;  %weather(1,1);
finish = end_grass; %weather(end,1);
index = (weather(:,1)>= start) & (weather(:,1)<=dormant_grass);
P = weather(index,5);
ET0 = refET(index);

% Values from Table 11, FAO-56 are noted (Data from various forages)
%calibrated FAO-56 parameters
if method == 1
    Lini_grass = 10; %was 30 %FAO-56, for grass pasture 10
    Ldev_grass = 45; %was 20%FAO-56, for grass pasture 20
    Lmid_grass = dormant_grass-growing_grass-10-45-90; %was 90
    Llat_grass = 90; %was 65 %FAO-56, for bermudagras for hay 35  
%parameters calibrated for FAW
elseif method == 2
    Lini_grass = 10; %was 30 %FAO-56, for grass pasture 10
    Ldev_grass = 45; %was 20%FAO-56, for grass pasture 20
    Lmid_grass = dormant_grass-growing_grass-10-45-90; %was 90
    Llat_grass = 90; %was 65 %FAO-56, for bermudagras for hay 35 
%default (uncalibrated) FAO-56 parameters
elseif method == 3 
    Lini_grass = 10; %FAO-56, for grass pasture 10
    Ldev_grass = 20; %FAO-56, for grass pasture 20
    Lmid_grass = 75; %FAO-56, for bermudagrass for hay 75
    Llat_grass = 35; %FAO-56, for bermudagras for hay 35    
end

%% Kcb values
% See table 17, FAO-56 for typical values for forages

%During spring and fall dormancy
Kcng_grass_dorm = 0.2;  %was 0.3
%Kcng_grass_dorm = Kcng_grass_dorm - Kcng_grass_dorm*pctcov/100*(.05/.1);
% Kc during the initial growth stage
Kcini_grass = 0.25;      %was 0.4
%Kcini_grass = Kcini_grass - Kcini_grass*pctcov/100*(.05/.1);
% Kc during the mid growth stage
Kcbmid_grass = temp_Kcbmid;%0.66;    %was 0.7  
% Kc at the end of the late season growth stage
Kcbend_grass = 0.2;     %was 0.2

%%
% Adjust Kcbmid and Kcbend based on weather data



%% Calculate the time points defining the crop coeffecient curve
times(1) = start;                                                               % beginning of the drainage year
times(2) =  growing_grass;                                                      % first grass growing season date
times(3) =  growing_grass + Lini_grass;                                         % date of end of initial period
times(4) =  growing_grass + Lini_grass + Ldev_grass;                            % date of end of development period
times(5) =  growing_grass + Lini_grass + Ldev_grass + Lmid_grass;               % date of end of Mid period
times(6) =  growing_grass + Lini_grass + Ldev_grass + Lmid_grass + Llat_grass;  % date of end of late period
times(7) = finish;

%Record end of mid period
date_end_mid = times(5);

%% Adjust Kcmid based on wind speed and minimum RH
index = (weather(:,1) >= times(4)) & (weather(:,1) < times(5));
u2 = nanmean(weather(index,9));  %Used nanmean to fill missing data. Would be better to 'prefill' missing data using other methods.
RHmin = nanmean(weather(index,7)); %Used nanmean to fill missing data. Would be better to 'prefill' missing data using other methods.
Kcbmid_grass = Kcbmid_grass + (0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3; %Eq. 62, FAO-56


%% Adjust Kcend based on wind speed and minimum RH during "late-season"
if times(5) < dormant_grass; % i.e., if the end of the mid season growing period occurs before harvest
    index = (weather(:,1) >= times(5)) & (weather(:,1) <= times(6));
    u2 = nanmean(weather(index,9));
    RHmin = nanmean(weather(index,7));
    Kcbend_grass = Kcbend_grass + (0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3; %Eq. 65, FAO-56
else times(5) = dormant_grass;
end;

%% Construct Kc output vector
index1 = days >= times(1) & days < times(2); %Kcng: non growing season 
index2 = days >= times(2) & days < times(3); %Kcini: initial part of growing season
index3 = days >= times(3) & days < times(4); %development period
index4 = days >= times(4) & days < times(5); %Kcmid:  middle of growing season
index5 = days >= times(5) & days < times(6); %Kcend: end of growing season or harvest (whichever comes first)
index6 = days >= times(6) & days <= times(7);%Kcng: non growing season

Kcb = ones(size(days,1),1); 
Kcb(index1) = Kcng_grass_dorm;
Kcb(index2) = Kcini_grass;
Kcb(index3) = interp1([times(3) times(4)],[Kcini_grass(end) Kcbmid_grass],days(index3));
Kcb(index4) = Kcbmid_grass;
Kcb(index5) = interp1([times(5) times(6)],[Kcbmid_grass Kcbend_grass],days(index5));
Kcb(index6) = Kcng_grass_dorm;
