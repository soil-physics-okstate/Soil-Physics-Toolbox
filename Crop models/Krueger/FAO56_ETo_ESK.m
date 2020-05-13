function [date,ETo] = FAO56_ETo(weather, windheight, Lat , elevation)
%FAO56_ETO calculates daily reference evapotranspiration (ETo).  This 
%reference evapotranspiration can then be multiplied by a crop
%coefficient (Kc) for various crops to determine daily crop
%evapotranspiration (ETc)
%
%v.1 ESK Jan. 9, 2008
%
%Data required for this caluclation are as follows:
%   weather    = [date doy Tmax Tmin P RHmax RHmin Rs wind]
%        date  = date in Matlab date number format
%        doy   = day of the year
%        Tmax  = maximum daily air temp (Celcius)
%        Tmin  = minimum daily air temp (Celcius)
%        P     = daily precip (mm)
%        RHmax = daily maximum relative humidity (%)
%        RHmin = daily minimum relative humidity (%)
%        Rs    = incoming solar radiation (MJ/m2/day)
%        wind  = wind speed (m/s)
%   windheight = height at which wind measurements were taken (m)
%   latitude   = latitude of experimental site (decimal degrees)
%   elevation  = elevation above sea level of experimental site (m)
%
%ETo calculation based on:
%Allen, R.G., L.S. Pereira, D. Raes, and M. Smith. 1998. Crop evapotranspiration:
%Guidelines for computing crop water requirements, FAO Irrigation and Drainage Paper No. 56.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Date Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert Date to Day Month Year Format for output
date = datevec(weather(:,1));
date(:,4:6)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determination of ETo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETo = 0.408*delta(Rn-G) + gamma(900/T+273)u2(es-ea) / delta+gamma(1+0.34u2)     %Eq. 6, FAO-56
    %ETo   = reference evapotranspiration (mm/day)
    %Rn    = net radiation at the crop surface (MJ/m2/day)
    %G     = soil heat flux density (MJ/m2/day)
    %T     = mean daily air temperature at 2m height
    %u2    = wind speed at 2 m height (m/s)
    %es    = saturation vapor pressure (kPa)
    %ea    = actual vapor pressure (kPa)
    %es-ea = saturation vapor pressure deficit (kPa)
    %delta = slope vapor pressure curve (kPa/degree C)
    %gamma = psychrometric constant (kPa/degree C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%PSYCHROMETRIC CONSTANT: gamma
%Calculation of gamma: psychrometric constant
    %gamma = Cp*P / epsilon*lamda = 0.000665(P)                            %Eq. 8, FAO-56
        %gamma   = psychrometric constant (kPa/degree C)
        %lamda   = latent heat of vaporization, 2.45 (MJ/kg)
        %Cp      = specific heat at constant pressure (MJ/kg/degree C)
        %        = 0.001013 for average atmospheric conditions
        %epsilon = ratio of molecular weight of water vapour/dry air = 0.622
        %P       = atmospheric pressure (kPa)
        %        = 101.3(293-0.0065z/293)^5.26
        %        z = elevation above sea level (m)
        %          = 347.5 m:  Obtained from National Climate Data
        %          Center
z = elevation;
P = 101.3*((293-.0065*z)/293)^5.26;
Cp = 0.001013;
epsilon =  0.622;
lamda = 2.45;
gamma = (Cp*P)/(epsilon*lamda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WIND SPEED
%Wind speed must be corrected to 2m above the soil surface
    %u2 = uz(4.87/ln(67.8z-5.42))                                          %Eq. 47, FAO-56
        %u2 = wind speed at 2 m above ground surface (m/s)
        %uz = measured wind speed at z m above ground surface (m/s)
        %zm = height of measurement above ground surface (m)
zm = windheight;
u2 = weather(:,9);
u2 = u2*(4.87/log((67.8*zm)-5.42));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AIR TEMPERATURE
Tmax = weather(:,3);
Tmin = weather(:,4);
Tmean = (Tmax + Tmin)/2;                                                   %Eq. 9, FAO-56
T = Tmean;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AIR HUMIDITY 
%Mean saturation vapor pressure: es
    %es = (eTmax + eTmin) / 2                                              %Eq. 12, FAO-56
          %es = mean saturation vapor pressure (kPa)
          %eTmax = saturation vapor pressure at temp Tmax (kPa)
          %eTmin = saturation vapor pressure at temp Tmin (kPa)       
eTmax = 0.6108*exp(17.27*Tmax./(Tmax+237.3));                               %Eq. 11, %FAO-56
eTmin = 0.6108*exp(17.27*Tmin./(Tmin+237.3));
es = (eTmax + eTmin) / 2;

%Slope of saturation of vapor pressure curve (delta)     
%delta = 4098*(0.6108*exp(17.27*Tmean./(Tmean+237.3)))./(Tmean+237.3).^2;   %Eq. 13, %FAO-56 
         %delta = slope of saturation vapor pressure curve at air temp T (kPa/degree C)
         %Tmean = average daily air temperture
delta = 4098*(0.6108*exp(17.27*Tmean./(Tmean+237.3)))./(Tmean+237.3).^2;   
        
%Actual vapor pressure derived from relative humidity data: ea
    %ea = (eTmin(RHmax/100) + eTmax(RHmin/100)) / 2                        %Eq. 17, FAO-56
        %ea = actual vapor pressure (kPa)
        %eTmax = saturation vapor pressure at temp Tmax (kPa)
        %eTmin = saturation vapor pressure at temp Tmin (kPa)
        %RHmax = maximum relative humidity (%)
        %RHmin = minimum relative humidity (%)
RHmax = weather(:,6);
RHmin = weather(:,7);    
ea = (eTmin.*(RHmax./100) + eTmax.*(RHmin./100)) / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RADIATION
%Extraterrestrial radiation for daily periods: Ra
    %Ra = 24(60)/pi * Gsc*dr * [omega*sin(phi)sin(d) + cos(phi)cos(d)sin(omega)]  %Eq. 21, FAO-56
        %Ra     = extraterrestrial radiation (MJ / m2 /day)
        %Gsc    = solar constant
        %       = 0.0820 MJ / m2 /min
        %dr     = inverse relative distance Earth-Sun (defined below)
        %phi    = latitude (rad) (defined below)
        %d      = solar decimation (rad) (defined below) (lower case delta in FAO-56)
        %omega  = sunset hour angle (radians) (defined below)  
        %dr     = 1 + 0.033 * cos(2pi*J/365)                               %Eq. 23, FAO-56
                 %J = number of the day of the year     
        %phi    = pi/180 * decimal degrees latitude                        %Eq. 22, FAO-56      
        %d      = 0.409sin((2pi * J/365) - 1.39)
        %omega  = pi/2 - (arccos(-tan(phi)tan(d))

%dr
J = weather(:,2);
dr= 1 + 0.033 * cos(2*pi*J/365);
%phi   
%decimal degrees
phi = pi/180 * Lat;
%d  
d = 0.409*sin((2*pi * J/365) - 1.39);
%omega 
omega = (acos(-tan(phi)*tan(d)));
%Ra 
Gsc = 0.0820;
Ra = 24*(60)/pi * Gsc.*dr .* (omega*sin(phi).*sin(d) + cos(phi).*cos(d).*sin(omega));

%Clear Sky Radiation: Rso
    %Rso = calculated clear sky radiation (MJ/m2/day)
    %    = (0.75 + (2*10^-5)z)Ra                                           %Eq. 37, FAO-56
Rso =  (0.75 + (2*10^-5)*z).*Ra ;

%Solar Radiation: Rs
    %This is a measured value (MJ/m2/day)
Rs = weather(:,8);
             
%Net Radiation: Rn
    %Rn = Rns - Rnl                                                        %Eq. 40, FAO-56
        %Rn  = net radiation (MJ/m2/day)
        %Rns = (1-alpha)Rs                                                 %Eq. 38, FAO-56
             %alpha  = albedo or canopy reflection coefficient
             %       = 0.23 for hypothetical grass reference crop
        %Rnl  = (sigma((TmaxK)^4 + ((TminK)^4)/2)*(0.34 - (squareroot(ea))(1.35(Rs/Rso) - 0.35) %Eq. 39, FAO-56
             %Rnl    = net outgoing longwave radiation (MJ/m2/day)
             %sigma  = Stefan-Boltzmann constant (4.903*10^-9 MJ/K4/m2/day)
             %TmaxK  = max absolute temperature during the 24hr period
             %       K  = C + 273.16
             %TminK  = min absolute temperature during the 24hr period
             %ea     = actual vapor pressure (kPa)
             %Rs/Rso = relative shortwave radiation (limited to <= 1.0)
alpha = 0.23;
Rns = (1-alpha)*Rs; 
sigma  = 4.903*10^-9;
TmaxK = Tmax + 273.16;
TminK = Tmin + 273.16;
Rnl =  sigma*(TmaxK.^4 + TminK.^4)/2.*(0.34 - 0.14*sqrt(ea)).*(1.35*(Rs./Rso) - 0.35); 
Rn = Rns - Rnl;                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOIL HEAT FLUX DENSITY: G
    %G = soil heat flux denstiy (MJ/m2/day)
    %     = 0 for daily time steps                                          %Eq. 42, FAO-56
G = 0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETo = 0.408*delta(Rn-G) + gamma(900/T+273)u2(es-ea) / delta+gamma(1+0.34u2)  %Eq. 6, FAO-56
    %ETo   = reference evapotranspiration (mm/day)
    %Rn    = net radiation at the crop surface (MJ/m2/day)
    %G     = soil heat flux density (MJ/m2/day)
    %T     = mean daily air temperature at 2m height
    %u2    = wind speed at 2 m height (m/s)
    %es    = saturation vapor pressure (kPa)
    %ea    = actual vapor pressure (kPa)
    %es-ea = saturation vapor pressure deficit (kPa)
    %delta = slope vapor pressure curve (kPa/degree C)
    %gamma = psychrometric constant (kPa/degree C)   
ETo = (0.408*delta.*(Rn-G) + gamma*(900./(T+273)).*u2.*(es-ea)) ./ (delta+gamma*(1+0.34*u2));   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
