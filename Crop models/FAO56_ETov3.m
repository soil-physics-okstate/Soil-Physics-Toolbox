function [ETo,varargout] = FAO56_ETov3(weather, windheight, latitude , elevation,varargin)
%% Function Description
% FAO56_ETO  calculates daily reference evapotranspiration (ETo).
%
% [date,ETo] = FAO56_ETO (weather,windheight,latitude,elevation) calculates ETo using the 9 by n weather matrix (detailed below), 
% height at which wind speed was measured, latitude in decimal degrees, 
% and the elevation of the location.
%
% ETo can then be multiplied by a crop
% coefficient (Kc) for various crops to determine daily crop
% evapotranspiration (ETc).
%% Inputs description
%       weather    = [date doy Tmax Tmin P RHmax RHmin Rs wind (vpd)].
%            date  = date in Matlab date number format.
%            doy   = day of the year.
%            Tmax  = maximum daily air temp [°C].
%            Tmin  = minimum daily air temp [°C].
%            P     = daily precip [mm].
%            RHmax = daily maximum relative humidity [%].
%            RHmin = daily minimum relative humidity [%].
%            Rs    = incoming solar radiation [MJ/m2/day].
%            wind  = wind speed [m/s].
%            vpd   = daily average vapor pressure deficit [kPa] **if available**
%       windheight = height at which wind measurements were taken [m].
%       latitude   = latitude of experimental site [decimal degrees].
%       elevation  = elevation above sea level of experimental site [m].
%         varargin = whether to perform aridity adjustment [Y/N].
%
%% Outputs
%       date = Matlab date format.
%       ETo  = Resulting reference evapotrasnpiration column vector [mm/day]
%
%% Updates
%
% v.1 ESK Jan. 9, 2008.
%
% v.2 TEO Sept. 18, 2008 added option to input vpd if available.
%
% v.3 TEO July 11, 2013 adding procedures to estimate missing weather
%
% v.3 AP August 5, 2013 fixed indexing problem in actual vapor pressure (ea), and solar radiation (Rs).
%
% v.3 AP September 7, 2014 added vpd as function output.
%
% v.3 AP February 11, 2016 check for missing values. Removed date from
% being the primary output of the function. ETo and vpd are the new output
% variables.
%
% v.3 BMW & JD October 13, 2019 added aridity adjustment option
%
% Last revised on: 13-OCT-2019.
%% References
% Allen, R.G., L.S. Pereira, D. Raes, and M. Smith. 1998. Crop evapotranspiration:
% Guidelines for computing crop water requirements, FAO Irrigation and Drainage Paper No. 56.
% Allen, R.G. 1996. Assessing Integrity of Weather Data for Reference Evapotranspiration

%% Contact Information
% For more information visit our website:
% <http://soilphysics.okstate.edu/>
%
% Or contact:
% <tyson.ochsner@okstate.edu>
%
% Soil Physics. Plant and Soil Sciences. Oklahoma State University.
%

%% Set variables
matlab_date = weather(:,1);
J = weather(:,2);             % Julian day (Day of the year)
Tmax = weather(:,3);          % Maximum air temperature
Tmin = weather(:,4);          % Minimum air temperature
RHmax = weather(:,6);         % Maximum air relative humidity
RHmin = weather(:,7);         % Minimum air relative humidity
Rs = weather(:,8);            % Measured solar radiation
u2 = weather(:,9);            % Wind speed
if size(weather,2) == 10
    vpd = weather(:,10);      % Vapor pressure deficit
end

%% Check for erros in input variables

if any(matlab_date<0) || any(isnan(matlab_date))
    error('One or more dates are negative or NaN');
else
    date = datevec(matlab_date); % Date in Matlab format
    date = date(:,1:3);
end

if any(J<0) || any(isnan(J))
    error('One or more DOY are negative or NaN');
end

if any(Tmax>56.7) || any(Tmax<-89.2)
    Tmax(Tmax>56.7) = NaN; % 56.7 Celsius is the maximum AIR temperature ever recorded at ground level on earth.
    Tmax(Tmax<-89.2) = NaN;
    warning('Maximum air temperature values outside the range -89.2 to 56.7 Celsius were converted to NaN.')
end

if any(Tmin>56.7) || any(Tmin<-89.2)
    Tmin(Tmin>56.7) = NaN; % -89.2 Celsius is the lowest AIR temperature record at ground level on earth.
    Tmin(Tmin<-89.2) = NaN;
    warning('Minimum air temperature values outside the range -89.2 to 56.7 Celsius were converted to NaN.')
end

if any(RHmax>100) || any(RHmin<0)
    RHmax(RHmax>100) = NaN;
    RHmin(RHmin<0) = NaN;
    warning('Minimum or maximum air relative humidity values outside the range 0 to 100% were converted to NaN.')
end

if any(Rs<0)
    Rs(Rs<0) = NaN;
    warning('Negative solar radiation values were converted to NaN.')
end

if any(u2<0)
    u2(u2<0) = NaN;
    warning('Negative wind speed values were converted to NaN.')
end

%% General equation for ETo
% This equation corrsponds to Eq.6, FAO-56.
%
% $$ETo = 0.408\Delta(Rn-G)+\gamma(900/T+273)u2(es-ea)/\Delta+\gamma(1+0.34u2)$$  
%
% ETo   = reference evapotranspiration (mm/day)
%
% Rn    = net radiation at the crop surface (MJ/m2/day).
%
% G     = soil heat flux density (MJ/m2/day)
%
% T     = mean daily air temperature at 2 m height
%
% u2    = wind speed at 2 m height (m/s)
%
% es    = saturation vapor pressure (kPa)
%
% ea    = actual vapor pressure (kPa)
%
% es-ea = saturation vapor pressure deficit (kPa)
%
% $\Delta$ = slope vapor pressure curve (kPa/degree C)
%
% $\gamma$ = psychrometric constant (kPa/degree C)
%
%% Psychrometric constant
% This equation corresponds to Eq. 8, FAO-56. 
%
% $$\gamma = Cp*P/\epsilon*\lambda$$
%
% $\gamma$ = psychrometric constant (kPa/°C)
%
% $\lambda$ = latent heat of vaporization, 2.45 (MJ/kg)
%
% Cp = specific heat at constant pressure (MJ/kg/°C)
%
% $\epsilon$ = ratio of molecular weight of water vapour/dry air = 0.622
%
% P = atmospheric pressure (kPa)
%
% z = elevation above sea level (m)
%
%                 
z = elevation;
P = 101.3*((293-.0065*z)/293)^5.26;
Cp = 0.001013; % Approx. 0.001013 for average atmospheric conditions
epsilon =  0.622;
lamda = 2.45;
gamma = (Cp*P)/(epsilon*lamda); % Approx. 0.000665 


%% Wind speed
% Wind speed must be corrected to 2 meter above the soil surface
%
% $$u2 = uz(4.87/\ln(67.8z-5.42))$$
%
% u2 = wind speed at 2 m above ground surface (m/s)
%
% uz = measured wind speed at z m above ground surface (m/s)
%
% zm = height of measurement above ground surface (m)

zm = windheight;
u2 = u2*(4.87/log((67.8*zm)-5.42));  % Eq. 47, FAO-56
u2(isnan(u2)) = 2; %replace missing windspeed data with 2 m/s, p. 63 Allen et al. (1998) 

%% Aridity adjustment (RG Allen, 1996)

% Estimate dewpoint temp to know if ETo aridity adjustment is necessary
if nargin == 5
    if strcmp(varargin{1},'Y')
        ES = 6.1365*exp((17.502*Tmin) ./ (240.97 + Tmin));   %use daily minimum temp
        E = (RHmax/100.0).*ES; % use daily max rel humidity, which should correspond to RH at minimum temp
        Tdew = 240.97*log(E/6.1365) ./ (17.502 - log(E/6.1365));
   
        Kar = 0.5; % empirical aridity proportionment coefficient
        deltaT_clim = 2; % for semiarid, arid climates
        deltaT = Tmin - Tdew - deltaT_clim;
        deltaT(deltaT < 0) = 0;
        Tmin = Tmin - Kar*deltaT; % replace measured Tmin data (Tmin_o)
        Tmax = Tmax - Kar*deltaT; % replace measured Tmax data (Tmax_o)
    end
end

%% Air temperature
Tmean = (Tmax + Tmin)/2;   % Eq. 9, FAO-56
T = Tmean;

%% Air humidity
% Mean saturation vapor pressure as defined in Eq. 12, FAO-56
%
% $$es = (eTmax+eTmin)/2$$                                              
%
% es = mean saturation vapor pressure (kPa)
%
% eTmax = saturation vapor pressure at temp Tmax (kPa)
%
% eTmin = saturation vapor pressure at temp Tmin (kPa)
eTmax = 0.6108*exp(17.27*Tmax./(Tmax+237.3)); % Eq. 11, %FAO-56
eTmin = 0.6108*exp(17.27*Tmin./(Tmin+237.3));
es = (eTmax + eTmin) / 2;

%% Vapor pressure
% Slope of saturation of vapor pressure curve ($\Delta$) as defined in Eq. 13, FAO-56
%
% $$\Delta = 4098(0.6108\exp(17.27Tmean/(Tmean+237.3)))/(Tmean+237.3)^2$$   
%
% $\Delta$ = slope of saturation vapor pressure curve at air temp T (kPa/°C)
%
% Tmean = average daily air temperture
%
%
% Actual vapor pressure derived from relative humidity data as defined in Eq. 17, FAO-56
%
% $$ea=(eTmin(RHmax/100)+eTmax(RHmin/100))/2$$ 
%
% ea = actual vapor pressure (kPa)
%
% eTmax = saturation vapor pressure at temp Tmax (kPa)
%
% eTmin = saturation vapor pressure at temp Tmin (kPa)
%
% RHmax = maximum relative humidity (%)
%
% RHmin = minimum relative humidity (%)

delta = 4098*(0.6108*exp(17.27*Tmean./(Tmean+237.3)))./(Tmean+237.3).^2;   
 
if exist('vpd','var')
    ea = es-vpd;
else
    ea = (eTmin.*(RHmax./100) + eTmax.*(RHmin./100)) / 2;
end

% Replace missing vapor pressure data with estimate based on minimum
% temperature Eq. 48, FAO-56
ea(isnan(ea)) = 0.611*exp(17.27*Tmin(isnan(ea))./(Tmin(isnan(ea))+237.3));

%% Solar radiation
% * Extraterrestrial radiation for daily periods as defined in Eq. 21, FAO-56
%
% $$Ra = 24(60)/\pi \hspace{2mm}Gsc \hspace{2mm} dr(\omega\sin(\phi)\sin(\delta)+\cos(\phi)\cos(\delta)\sin(\omega))$$
%
% Ra = extraterrestrial radiation (MJ / m2 /day)
%
% Gsc = solar constant (MJ/m2/min)
%   
% dr = 1 + 0.033$\cos$(2$\pi$J/365)  Inverse relative distance Earth-Sun                            
%
% J = number of the day of the year     
%
% $\phi$ = $\pi$/180decimal degrees  (latitude in radians)     
%
% $\delta = 0.409\sin((2\pi J/365)-1.39)\hspace{5mm}$ Solar decimation (rad)
%
% $\omega = \pi/2-(\arccos(-\tan(\phi)\tan(\delta)) \hspace{5mm}$ sunset hour angle (radians) 

dr= 1 + 0.033 * cos(2*pi*J/365);  % Eq. 23, FAO-56 
phi = pi/180 * latitude; % Eq. 22, FAO-56   
d = 0.409*sin((2*pi * J/365) - 1.39);
omega = (acos(-tan(phi)*tan(d)));
Gsc = 0.0820; % Approx. 0.0820 
Ra = 24*(60)/pi * Gsc.*dr .* (omega*sin(phi).*sin(d) + cos(phi).*cos(d).*sin(omega));
%%
% * Clear Sky Radiation: Rso (MJ/m2/day)                                        
Rso =  (0.75 + (2*10^-5)*z).*Ra ; % Eq. 37, FAO-56
%%
% * Measured solar Radiation: Rs (MJ/m2/day)   
RsIdx = isnan(Rs);
Rs(RsIdx) = min(0.16*Ra(RsIdx).*(Tmax(RsIdx)-Tmin(RsIdx)).^0.5,Rso(RsIdx)); %replace  Rs data using Eq. 50, FAO-56

%%
% * Net Radiation: Rn (MJ/m2/day)
%
% Rn = Rns - Rnl       
%
% Rns = (1-alpha)Rs                                                 
%
% alpha = albedo or canopy reflection coefficient
%
% Rnl = (sigma((TmaxK)^4 + ((TminK)^4)/2)*(0.34 - (squareroot(ea))(1.35(Rs/Rso) - 0.35) 
%
% Rnl = net outgoing longwave radiation (MJ/m2/day)
%
% $\sigma$ = Stefan-Boltzmann constant (4.903*10^-9 MJ/K4/m2/day)
%
% TmaxK = max absolute temperature during the 24hr period (K = C+273.16)
%
% TminK = min absolute temperature during the 24hr period
%
% ea = actual vapor pressure (kPa)
%
% Rs/Rso = relative shortwave radiation (limited to <= 1.0)
alpha = 0.23; % 0.23 for hypothetical grass reference crop
Rns = (1-alpha)*Rs; % Eq. 38, FAO-56
sigma  = 4.903*10^-9;
TmaxK = Tmax + 273.16;
TminK = Tmin + 273.16;
Rnl =  sigma*(TmaxK.^4 + TminK.^4)/2.*(0.34 - 0.14*sqrt(ea)).*(1.35*(Rs./Rso) - 0.35); % Eq. 39, FAO-56
Rn = Rns - Rnl; % Eq. 40, FAO-56                  

%% Soil heat flux density  
%
% G = 0 for daily time steps  (MJ/m2/day).
%
G = 0; % Eq. 42, FAO-56
%%
% *ETo calculation* 
ETo = (0.408*delta.*(Rn-G) + gamma*(900./(T+273)).*u2.*(es-ea)) ./ (delta+gamma*(1+0.34*u2));   

if nargout == 2
    vpd = ea - es;
    varargout{2} = vpd; % Output Vapor pressure deficit (VPD)
elseif nargout == 3
    varargout{3} = date;
end

ETo(ETo<0) = NaN; % Restricts the output of negative values.
