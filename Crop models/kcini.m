function [Kc ET0out Pout totaldays] = kcini(start,enddate,soil,theta_ini,weather,refET,pctcov,livecov)
%KCINI determines the crop coefficient for the initial period for the
%purpose of scaling reference ET to produce actual ET estimates.  This
%procedure is also appropriate for the intial period of perennial crops
%or outside of the normal growing season.
%
%v.1 TEO Oct. 4, 2007
%
%start = starting date (Matlab date number format)
%enddate = ending date (Matlab date number format)
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
%pctcov = percent non-living residue cover during initial period 
%livecov = 'kura' or 'none'
%
%Kcini is calculated following:
%Allen, R.G., W.O. Pruitt, D. Raes, M. Smith, and L.S. Pereira. 2005.
%Estimating evapotranspiration from bare soil and the crop coefficient for
%the initial period using common soils information. J. Irrig. Drain. Eng. 131:14-23.

totaldays = [start:enddate]'; %days in the initial period
d = length(totaldays); %number of days in initial period
n = ceil(d/30);%floor(d/30); %number of complete 30 day segments

    %Define soil specific parameters
    Sa = soil(1); % percent sand
    Cl = soil(2); %percent clay
    OM = soil(3); %percent OM
    [theta_fc theta_wp] = rawls(Sa,Cl,OM);

    if Sa >= 80 %Eq. 5 from Allen et al. (2005) based on Ritchie et al. (1989)
        REWmax = 20-0.15*Sa;
    elseif Cl >= 50
        REWmax = 11-0.06*Cl;
    else
        REWmax = 8 + 0.08*Cl;
    end

    Ze = 0.1 + 0.05*(1-Sa/100); %FAO-56 recommends 0.1 coarse to 0.15 fine.  This function makes Ze scale linearly with sand content.
    Wini = Ze*theta_ini*1000;

for j = 1:n %calculate a Kcini value for each 30 day segment or portion thereof
    if j == 1
        if d < 30
            days = totaldays;
        else
            days = [start:(start+29)]';
        end
    elseif j<n
        days = [(start+(j-1)*30):(start+j*30-1)]';
    else
        days = [(start+(j-1)*30):enddate]';
    end

    index = (weather(:,1)>= days(1)) & (weather(:,1)<=days(end));
    P = weather(index,5);
    u2 = mean(weather(index,9));
    RHmin = mean(weather(index,7));
    ET0 = refET(index);
    
    %Calc. TEWmax

    ET0mean = mean(ET0); %calculate average ET0 during initial period
    TEWmaxEq4 = 1000*(theta_fc - 0.5*theta_wp)*Ze;
    TEWmaxEq8 = TEWmaxEq4*(ET0mean/5).^0.5;
    TEWmax = min([TEWmaxEq4 TEWmaxEq8]);

    %Calc. Kcini

    nw = 0; %number of wetting events
    previous = 0; %precip on previous day
    sumP = 0; %total precip for all effective wetting events
    event = 0; %total precip for current wetting event

    for i =1:size(P,1)
        if P(i)>0
            if previous == 0
                event = P(i);
            else
                event = event + P(i);
            end
        else
            if event > ET0mean*0.2
                nw = nw+1;
                sumP = min([event TEWmax]) + sumP;
                event = 0;
            else
                event = 0;
            end
        end
        previous = P(i);
    end

    if nw > 0
        Pmean = sumP/nw; %mean depth of water added to the evaporating layer from each wetting event
        REW = REWmax*min([1 (Pmean+Wini/nw)/TEWmax]);
        TEW = min([TEWmax (Pmean+Wini/nw)]);
    else
        Pmean = 0;
        REW = REWmax*min([1 (Pmean+Wini/1)/TEWmax]);
        TEW = min([TEWmax (Pmean+Wini/1)]);
    end

    tw = (days(end)-days(1))/(nw+0.5); %avg. time between wetting events during the initial period

    if livecov == 'kura'
        Kcbini = 0.3; %Table 17, FAO-56
        h = 0.1; %height in meters 
        Kcmax = max(1.2 + [0.04*(u2-2)-0.004*(RHmin-45)]*(h/3)^0.3,Kcbini+0.05); %Eq. 72 FAO-56
    end

    if livecov == 'none'
        Kcbini = 0;
        Kcmax = 1.15;
    end

    Eso = (Kcmax - Kcbini)*ET0mean;
    t1 = REW/Eso; %time when Stage 1 drying is completed

    if tw > t1
        Kcini = (TEW - (TEW-REW)*exp(-(tw*Eso-REW)/(TEW-REW)))/(tw*ET0mean);
    else
        Kcini = Kcmax - Kcbini;
    end

    temp = ones(size(days,1),1);
    temp(:) = Kcini;

    %Save the results
    if j == 1
        Kc = temp;
        Pout = P;
        ET0out = ET0;
    else
        Kc = [Kc;temp];
        Pout = [Pout;P];
        ET0out = [ET0out;ET0];
    end
end

%Reduce Kc 5% for every 10% of ground cover
Kc = Kc*(1-0.05*pctcov/10);

