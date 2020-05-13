function [Kc, ET0, P, days] = cornco(planting_corn,harvest_corn,soil,theta_ini,weather,refET,pctcov)
%CORNCO constructs crop coefficient curves for corn (maize). These crop coefficient 
%curves are for the purpose of scaling reference ET to produce actual ET estimates. 
%
%v.1 ESK Jan. 21, 2008
%
%planting = planting date (Matlab date number format)
%harvest = harvest date (Matlab date number format)
%soil = [sand clay OM] percent sand, clay, and organic matter in soil
%theta_ini = initial volumetric water content of the surface soil
%weather    = [date doy Tmax Tmin P RHmax RHmin Rs wind]
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
%refET = daily reference evapotranspiration (mm)
%pctcov = percent residue cover during initial period
%
%Kcini is calculated following:
%Allen, R.G., W.O. Pruitt, D. Raes, M. Smith, and L.S. Pereira. 2005. 
%Estimating evapotranspiration from bare soil and the crop coefficient for
%the initial period using common soils information. J. Irrig. Drain. Eng. 131:14-23.
%
%Kcmid and Kcend based on:
%Allen, R.G., L.S. Pereira, D. Raes, and M. Smith. 1998. Crop evapotranspiration:
%Guidelines for computing crop water requirements, FAO Irrigation and Drainage Paper No. 56. Rome.
%
%Kc for corn (with no rye) after corn harvest until planting is calculated
%following Kcini procedure above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the time points defining the crop coeffecient curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 2; %maximum crop height (m),Table 12, FAO-56
Zr = 1; %maximum effective rooting depth (m), Table 22, FAO-56
ptab = 0.55; %fraction of total available soil water that is readily available, Table 22, FAO-56

days = (weather(1,1):weather(end,1))'; %days in the crop coefficient period
start = weather(1,1);
index = (weather(:,1)>= start) & (weather(:,1)<=harvest_corn);
P = weather(index,5);
ET0 = refET(index);

%Values from Tables 11 and 12, FAO-56
Lini_corn = 30;
Ldev_corn = 40;
Lmid_corn = 50;
Llat_corn = 50;
Kcmid_corn = 1.2;
Kcend_corn = 0.6; %for harvest of high moisture grain

%Calculate the time points defining the crop coeffecient curve
times(1) = start; % beginning of the drainage year
times(2) =  planting_corn; %corn planting date
times(3) =  planting_corn + Lini_corn; % date of end of initial period
times(4) =  planting_corn + Lini_corn + Ldev_corn; % date of end of development period
times(5) =  planting_corn + Lini_corn + Ldev_corn + Lmid_corn; % date of end of Mid period
times(6) =  planting_corn + Lini_corn + Ldev_corn + Lmid_corn + Llat_corn; % date of end of late period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during non-growing season:  Kcng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
livecov = 'none';
[Kcng_corn ET0out Pout totaldays] = kcini(times(1),times(2)-1,soil,theta_ini,weather,refET,pctcov,livecov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during the initial growth period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Kcini
livecov = 'none';
[Kcini_corn ET0out Pout totaldays] = kcini(times(2),times(3)-1,soil,theta_ini,weather,refET,pctcov,livecov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during mid season growth:  Kcmid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust Kcmid based on wind speed and minimum RH
index = (weather(:,1)>=times(4)) & (weather(:,1)<times(5));
u2 = mean(weather(index,9));
RHmin = mean(weather(index,7));
Kcmid_corn = Kcmid_corn + ( 0.04*(u2-2)-0.004*(RHmin-45) )*(h/3)^0.3; %Eq. 62, FAO-56

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during late season:  Kcend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust Kcend based on wind speed and minimum RH during
%"late-season"
if times(5) < harvest_corn; %ie, if the end of the mid season growing period occurs before harvest
    index = (weather(:,1)>=times(5)) & (weather(:,1)<=times(6));
    u2 = mean(weather(index,9));
    RHmin = mean(weather(index,7));
    Kcend_corn = Kcend_corn + ( 0.04*(u2-2)-0.004*(RHmin-45) )*(h/3)^0.3; %Eq. 62, FAO-56
else times(5) = harvest_corn;
end;

%Construct Kc output vector
index1 = days >= times(1) & days < times(2); %Kcng: non growing season 
index2 = days >= times(2) & days < times(3); %Kcini: initial part of growing season
index3 = days >= times(3) & days < times(4); %development period
index4 = days >= times(4) & days < times(5); %Kcmid:  middle of growing season
index5 = days >= times(5) & days <= harvest_corn; %Kcend: end of growing season or harvest (whichever comes first)

Kc = ones(size(days,1),1);
Kc(index1) = Kcng_corn;
Kc(index2) = Kcini_corn;
Kc(index3) = interp1([times(3) times(4)],[Kcini_corn(end) Kcmid_corn],days(index3));
Kc(index4) = Kcmid_corn;
Kc(index5) = interp1([times(5) times(6)],[Kcmid_corn Kcend_corn],days(index5));

