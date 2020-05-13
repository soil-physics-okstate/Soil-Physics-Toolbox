function [Kc, ET0, P, days] = forestco(growing_forest,dormant_forest,soil,theta_ini,weather,refET,pctcov)
%CORNCO constructs crop coefficient curves for corn (maize). These crop coefficient 
%curves are for the purpose of scaling reference ET to produce actual ET estimates. 
%
%v.1 ESK Jan. 21, 2008
%
%growing = first growing date (Matlab date number format)
%dormant = first dormant date (Matlab date number format)
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
h = 10; %maximum crop height (m),Table 12, FAO-56
%Zr = 1.0; %maximum effective rooting depth (m)
%ptab = 0.6; %fraction of total available soil water that is readily available, best estimate from Table 22, FAO-56

days = (growing_forest:growing_forest+365.25)'; %days in the crop coefficient period
start = weather(1,1);
index = (weather(:,1)>= start) & (weather(:,1)<=dormant_forest);
P = weather(index,5);

% Account for interception in forest
for i = 1:length(P);
    if P(i) > 2.0
        P(i) = P(i) - 2.0;
    else
        P(i) = 0;
    end
end

ET0 = refET(index);

%Best estimates, based on LAI curves from Johnson & Risser (1974) and Norby et al. (2003)
Lini_forest = 5;
Ldev_forest = 45;
Lmid_forest = 70;
Llat_forest = 77;
Ltrans_forest = 3;
Kcmid_forest = 1.2;       % Best guess
Kcend_forest = 0.5;     % Best guess

%Calculate the time points defining the crop coeffecient curve
times(1) = start;                                       % beginning of the drainage year
times(2) =  growing_forest;                             % beginning of growing season
times(3) =  growing_forest + Lini_forest;               % date of end of initial period
times(4) =  growing_forest + Lini_forest + Ldev_forest;                                     % date of end of development period
times(5) =  growing_forest + Lini_forest + Ldev_forest + Lmid_forest;                       % date of end of Mid period
times(6) =  growing_forest + Lini_forest + Ldev_forest + Lmid_forest + Llat_forest;         % date of end of late period
times(7) =  growing_forest + Lini_forest + Ldev_forest + Lmid_forest + Llat_forest + Ltrans_forest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during non-growing season:  Kcng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
livecov = 'none';
[Kcng_forest ET0out Pout totaldays] = kcini(times(1),times(2),soil,theta_ini,weather,refET,pctcov,livecov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during the initial growth period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Kcini
livecov = 'none';
[Kcini_forest ET0out Pout totaldays] = kcini(times(2),times(3)-1,soil,theta_ini,weather,refET,pctcov,livecov);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during mid season growth:  Kcmid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust Kcmid based on wind speed and minimum RH
index = (weather(:,1)>=times(4)) & (weather(:,1)<times(5));
u2 = mean(weather(index,9));
RHmin = mean(weather(index,7));
Kcmid_forest = Kcmid_forest + ( 0.04*(u2-2)-0.004*(RHmin-45) )*(h/3)^0.3; %Eq. 62, FAO-56

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during late season:  Kcend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust Kcend based on wind speed and minimum RH during
%"late-season"
if times(5) < dormant_forest; %ie, if the end of the mid season growing period occurs before harvest
    index = (weather(:,1)>=times(5)) & (weather(:,1)<=times(6));
    u2 = mean(weather(index,9));
    RHmin = mean(weather(index,7));
    Kcend_forest = Kcend_forest + ( 0.04*(u2-2)-0.004*(RHmin-45) )*(h/3)^0.3; %Eq. 62, FAO-56
else times(5) = dormant_forest;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kc during dormant season:  Kcdorm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
livecov = 'none';
[Kcdorm_forest ET0out Pout totaldays] = kcini(times(6),days(end),soil,theta_ini,weather,refET,pctcov,livecov);

%Construct Kc output vector
index1 = days >= times(1) & days < times(2); %Kcng: non growing season 
index2 = days >= times(2) & days < times(3); %Kcini: initial part of growing season
index3 = days >= times(3) & days < times(4); %development period
index4 = days >= times(4) & days < times(5); %Kcmid:  middle of growing season
index5 = days >= times(5) & days < times(6); %Kcend: end of growing season or harvest (whichever comes first)
index6 = days >= times(6) & days <= days(end);

Kc = ones(size(days,1),1);
Kc(index1) = Kcng_forest;
Kc(index2) = Kcini_forest;
Kc(index3) = interp1([times(3) times(4)],[Kcini_forest(end) Kcmid_forest],days(index3));
Kc(index4) = Kcmid_forest;
Kc(index5) = interp1([times(5) times(6)],[Kcmid_forest Kcend_forest],days(index5));
Kc(index6) = Kcdorm_forest;

