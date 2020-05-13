function [ETc, Dr, Ke, Ks, Kcmax] = dualcropco_v2(weather,refET,Kcb,Kcmin,h,Zr,ptab,soil,rescov,thetai)
%DUALCROPCO calculates coefficients needed for the FAO-56 dual crop
%coefficient method
%
%weather    = [date doy Tmax Tmin P RHmax RHmin Rs wind]
%        1. date  = date in Matlab date number format
%        2. doy   = day of the year
%        3. Tmax  = maximum daily air temp (Celcius)
%        4. Tmin  = minimum daily air temp (Celcius)
%        5. P     = daily precip (mm)
%        6. RHmax = daily maximum relative humidity (%)
%        7. RHmin = daily minimum relative humidity (%)
%        8. Rs    = incoming solar radiation (MJ/m2/day)
%        9. wind  = wind speed (m/s)
%refET = daily reference evapotranspiration (mm)
%Kcb = basal crop coefficient curve for period of simulation
%Kcmin = single minimum Kcb value for period of simulation
%       Kcmin = 0.15 during growing season
%       Kcmin = 0 during dormant season
%h = plant height curve for period of simulation (m)
%Zr = root zone depth curve for period of simulation (m)
%ptab = depletion fraction, Table 22, FAO-56, fraction of available water
%   which can be depleted from the root zone before moisture stress occurs
%soil = [sand clay OM] percent sand, clay, and organic matter in soil
%rescov = fraction of surface covered with non-living residue or dormant
%           plants, rescov = 0 during growing season
%thetai = initial water content values for surface layer 

%Soil properties
Sa = soil(1); %percent sand in surface soil
Cl = soil(2); %percent clay in surface soil
OM = soil(3); %percent OM
[theta_fc, theta_wp] = rawls(Sa,Cl,OM); %field capacity and wilting point from Rawls ptf

%Available and evaporable water variables
Ze = 0.1 + 0.05*(1-soil(1)/100); %depth (m) of soil layer subject to evaporation
%Allen et al. (2005) recommend 0.1 coarse to 0.15 fine.  
%The equation above makes Ze scale linearly with sand content.
TEW = 1000*(theta_fc(1)-0.5*theta_wp(1))*Ze; %Eq. 73, FAO-56 total evaporable water
TEW = TEW - TEW*rescov*(.05/.1); %reduce TEW by 5% for every 10% of residue cover (Ch. 11, FAO-56)
if Sa >= 80 %Eq. 5 from Allen et al. (2005) based on Ritchie et al. (1989), readily evaporable water
    REW = 20-0.15*Sa;
elseif Cl >= 50
    REW = 11-0.06*Cl;
else
    REW = 8 + 0.08*Cl;
end
REW = min(REW,TEW); %REW is limited to be less than or equal to TEW
fw = 1; %fraction of surface wetted by precipitation

%Weather variables
days = weather(:,1);
P = weather(:,5);
RHmin = weather(:,7);
u2 = weather(:,9);
Tmax = weather(:,3);
ET0 = refET;

%Calculations for first day 
for i = 1
    Kcmax(i) = max((1.2 + [0.04*(u2(i)-2)-0.004*(RHmin(i)-45)]*(h(i)/3)^0.3),(Kcb(i) + 0.05)); %Eq. 72, FAO-56
    fc(i) = ((Kcb(i)-Kcmin(i))/(Kcmax(i)-Kcmin(i)))^(1+0.5*h(i)); %approximate crop ground cover, Eq. 76, FAO-56
    few(i) = min(1-fc(i),fw); %fraction of surface exposed and wetted
    De0 = (theta_fc(1)-thetai(1))*Ze*1000; %initial depletion of surface layer
    if De0 > REW
        Kr(i) = (TEW(i) - De0)/(TEW(i) - REW(i)); %Eq. 74, FAO-56, evaporation reduction coefficient
    else
        Kr(i) = 1;
    end
    Ke(i) = min(few(i)*Kcmax(i),Kr(i)*(Kcmax(i)-Kcb(i))); %Eq. 71, FAO-56, soil evaporation coefficient
    E(i) = Ke(i)*ET0(i); %soil evaporation (mm)
    DPe(i) = max(0,P(i)-De0); %Eq. 79, FAO-56, drainage of water from the bottom of the evaporating layer
    De(i) = De0 - P(i) + E(i)/few(i) + DPe(i); %Eq. 77, FAO-56, cumulative
    TAW(i) = 1000*(theta_fc - theta_wp)*Zr(i); %Eq. 82, FAO-56, total available water (mm)
    p = ptab + 0.04*(5-Kcb(i)*ET0(i)); %Fig. 41, FAO-56
    p = min([p 0.8]); %Fig. 41, FAO-56
    p = max([0.1 p]); %Fig. 41, FAO-56
    RAW(i) = p*TAW(i); %Eq. 83, FAO-56, readily available water (mm)
    Dr0 = (theta_fc-thetai)*Zr(i)*1000; %initial depletion of root zone
    Dr0 = min(TAW(i),Dr0); %Set Max Dr0: Dr0 cannot exceed TAW:  physically this cannot happen, but may arise when initial soil water content is low
    if Dr0 > RAW(i)
        Ks(i,1) = (TAW(i)-Dr0)/[TAW(i)-RAW(i)]; %Eq. 84, FAO-56, transpiration reduction factor
    else
        Ks(i,1) = 1;
    end
    ETc(i) = (Ks(i)*Kcb(i)+Ke(i))*ET0(i); %Eq. 80, FAO-56, total crop water use for the day (mm)
    DP(i) = max(0,P(i) - ETc(i)-Dr0); %Eq. 88, FAO-56, drainage of water from the bottom of the root zone
    Dr(i) = Dr0 - P(i) + ETc(i) + DP(i); %Eq. 85, FAO-56, root zone depletion at the end of the day
    Dr(i) = min(TAW(i),Dr(i)); %Set Max Dr(i): Dr(i) cannot exceed TAW:  physically this cannot happen, but may arise when initial soil water content is low
end

%Calculations for subsequent days
Fcum = 0; %initialize frozen water inventory
F = 0;
for i = 2:length(days)
    Kcmax(i) = max((1.2 + [0.04*(u2(i)-2)-0.004*(RHmin(i)-45)]*(h(i)/3)^0.3),(Kcb(i) + 0.05)); %Eq. 72, FAO-56
    fc(i) = ((Kcb(i)-Kcmin(i))/(Kcmax(i)-Kcmin(i)))^(1+0.5*h(i)); %approximate crop ground cover, Eq. 76, FAO-56
    few(i) = min(1-fc(i),fw); %fraction of surface exposed and wetted
    if De(i-1) > REW
        Kr(i) = (TEW(i) - De(i-1))/(TEW(i) - REW(i)); %Eq. 74, FAO-56, evaporation reduction coefficient
    else
        Kr(i) = 1;
    end
    Ke(i,1) = min(few(i)*Kcmax(i),Kr(i)*(Kcmax(i)-Kcb(i))); %Eq. 71, FAO-56, soil evaporation coefficient
    E(i,1) = Ke(i)*ET0(i); %soil evaporation (mm)
    DPe(i) = max(0,P(i)-De(i-1)); %Eq. 79, FAO-56, drainage of water from the bottom of the evaporating layer
    if Tmax(i)<=0
        F(i) = 1; %assume soil water freezes at a rate of 1 mm/day when Tmax<=0, (TEO approximation)
        Fcum(i) = F(i)+Fcum(i-1); %cumulative frozen water depth
    elseif Fcum(i-1) > 0
        F(i) = -1; %assume soil water thaws at a rate of 1 mm/day when Tmax>0
        Fcum(i) = max(0,F(i)+Fcum(i-1)); %cumulative frozen water depth
    else
        F(i) = 0;
        Fcum(i) = 0;
    end
    De(i) = De(i-1) - P(i) + E(i)/few(i) + DPe(i) + F(i); %Eq. 77, FAO-56, cumulative
    De(i) = min(TEW(i),De(i)); %Depletion cannot exceed TEW, required because F is not constrained
    TAW(i) = 1000*(theta_fc-theta_wp)*Zr(i); %Eq. 82, FAO-56, total available water (mm)
    p = ptab + 0.04*(5-Kcb(i)*ET0(i)); %Fig. 41, FAO-56
    p = min([p 0.8]); %Fig. 41, FAO-56
    p = max([0.1 p]); %Fig. 41, FAO-56
    RAW(i) = p*TAW(i); %Eq. 83, FAO-56, readily available water (mm)
    if Dr(i-1) > RAW(i)
        Dr(i-1) = min(TAW(i),Dr(i-1)); %Set Max Dr(i-1): Dr(i-1) cannot exceed TAW:  Dr(i-1) may exceed TAW when rye is going into or coming out of dormancy
        Ks(i,1) = (TAW(i)-Dr(i-1))/[TAW(i)-RAW(i)]; %Eq. 84, FAO-56, transpiration reduction factor
    else
        Ks(i,1) = 1;
    end
    ETc(i,1) = (Ks(i)*Kcb(i)+Ke(i))*ET0(i); %Eq. 80, FAO-56, total crop water use for the day (mm)
    DP(i) = max(0,P(i) - ETc(i) - Dr(i-1)); %Eq. 88, FAO-56, drainage of water from the bottom of the root zone
    Dr(i,1) = Dr(i-1) - P(i) + ETc(i) + DP(i); %Eq. 85, FAO-56, root zone depletion at the end of the day
    Dr(i) = max([Dr(i) 0]); %Eq. 85, FAO-56
    Dr(i) = min([Dr(i) TAW(i)]); %Eq. 86, FAO-56
end



