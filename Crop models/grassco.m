function [Kcb, ET0, P, days] = grassco(growing_grass,dormant_grass,soil,theta_ini,weather,refET,pctcov)
% CORNCO constructs dual crop coefficient curves for grasslands. These crop coefficient 
% curves are for the purpose of scaling reference ET to produce actual ET estimates. 
%
% v.1 ESK Jan. 21, 2008
% Code modified for grasslands 25 Oct. 2016 by BMW
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

h = 0.5;        % mean of maximum crop heights for forages excluding turf grass (m), Table 12, FAO-56
Zr = 0.8;       % maximum effective rooting depth (m), based on Mesonet sensors depths
ptab = 0.57;    % mean of fraction of total available soil water that is readily available for forages except turf grass, Table 22, FAO-56

days = (weather(1,1):weather(end,1))'; %days in the crop coefficient period
start = weather(1,1);
index = (weather(:,1)>= start) & (weather(:,1)<=dormant_grass);
P = weather(index,5);
ET0 = refET(index);

% Values from Table 11, FAO-56 (Averaged for forages and rounded to
% nearest integer
Lini_grass = 10;
Ldev_grass = 20;
Lmid_grass = 24;
Llat_grass = 15;

% Kcb values from Table 17, FAO-56 (Averaged for forages and rounded to
% nearest integer
Kcbmid_grass = 0.91;     % 0.95 for Kc only (Table 12)
Kcbend_grass = 0.83;     % 0.88 for Kc only (Table 12)

% Calculate the time points defining the crop coeffecient curve
times(1) = start;                                                               % beginning of the drainage year
times(2) =  growing_grass;                                                      % first grass growing season date
times(3) =  growing_grass + Lini_grass;                                         % date of end of initial period
times(4) =  growing_grass + Lini_grass + Ldev_grass;                            % date of end of development period
times(5) =  growing_grass + Lini_grass + Ldev_grass + Lmid_grass;               % date of end of Mid period
times(6) =  growing_grass + Lini_grass + Ldev_grass + Lmid_grass + Llat_grass;  % date of end of late period

%% Kc during non-growing season:  Kcng

livecov = 'none';
[Kcng_grass ET0out Pout totaldays] = kcini(times(1),times(2)-1,soil,theta_ini,weather,refET,pctcov,livecov);

%% Kc during the initial growth period

%Calculate Kcini
livecov = 'none';
[Kcini_grass ET0out Pout totaldays] = kcini(times(2),times(3)-1,soil,theta_ini,weather,refET,pctcov,livecov);

%% Kc during mid season growth:  Kcmid

%Adjust Kcmid based on wind speed and minimum RH
index = (weather(:,1) >= times(4)) & (weather(:,1) < times(5));
u2 = mean(weather(index,9));
RHmin = mean(weather(index,7));
Kcbmid_grass = Kcbmid_grass + (0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3; %Eq. 62, FAO-56

%% Kc during late season:  Kcend

%Adjust Kcend based on wind speed and minimum RH during "late-season"
if times(5) < dormant_grass; % i.e., if the end of the mid season growing period occurs before harvest
    index = (weather(:,1) >= times(5)) & (weather(:,1) <= times(6));
    u2 = mean(weather(index,9));
    RHmin = mean(weather(index,7));
    Kcbend_grass = Kcbend_grass + (0.04*(u2-2)-0.004*(RHmin-45))*(h/3)^0.3; %Eq. 62, FAO-56
else times(5) = dormant_grass;
end;

%% Construct Kc output vector
index1 = days >= times(1) & days < times(2); %Kcng: non growing season 
index2 = days >= times(2) & days < times(3); %Kcini: initial part of growing season
index3 = days >= times(3) & days < times(4); %development period
index4 = days >= times(4) & days < times(5); %Kcmid:  middle of growing season
index5 = days >= times(5) & days <= dormant_grass; %Kcend: end of growing season or harvest (whichever comes first)

Kcb = ones(size(days,1),1);
Kcb(index1) = Kcng_grass;
Kcb(index2) = Kcini_grass;
Kcb(index3) = interp1([times(3) times(4)],[Kcini_grass(end) Kcbmid_grass],days(index3));
Kcb(index4) = Kcbmid_grass;
Kcb(index5) = interp1([times(5) times(6)],[Kcbmid_grass Kcbend_grass],days(index5));

