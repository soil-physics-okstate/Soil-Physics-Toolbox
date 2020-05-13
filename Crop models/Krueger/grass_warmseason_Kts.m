function [Kts] = grass_warmseason_Kts(weather_year)

%ESK
%10/19/2018

%This script was designed to calculate the temperature crop stress 
%coefficient for warm season grasses. 

%% Develop temperature stress coefficient curve
%This subroutine calculates the temperature stress coefficient for warm
%season grasses.  The stress coefficient is calculated based on the
%logistic equation (Raes et al., 2009) for upper and lower temperature
%stresses.  The logistic equation was visually compared to results of
%studies detailing the impact of air (Kakani and Reddy, 2007) or soil
%temperature (Delucia et al., 1992) on biomass accumulation in big 
%bluestem.  

    %Delucia, E.H., S.A. Heckathorn and T.A. Day. 1992. Effects of soil 
        %temperature on growth, biomass allocation and resource acquisition 
        %of Andropogon gerardii Vitman. New Phytologist. 120:543-549.

    %Kakani, V.G. and K.R. Reddy. 2007. Temperature response of C4 species 
        %big bluestem (Andropogon gerardii) is modified by growing carbon 
        %dioxide concentration. Environmental and Experimental Botany 
        %61:281-290. 
    
    %Raes D, S.P., Hsiao TC, Fereres E. 2009. AquaCrop-the FAO crop model 
        %to simulate yield response to water: II. main algorithms and 
        %software description. Agronomy journal 101:438-447.

%Calculate Q and B parameters for the logistic equation
    %B is the steepness of the curve
    %Q is an arbitrary constant
    
    %Parameter values were determined using a combination of matlab fitting
    %and visual assessment.
    
    %{
    %Matlab fit Kts_lower
        X = [0:1:28]'; 
        g = fittype( @(Q, B, x) ((1 ./ (1 + Q*exp(B*(x))))) )
        [curve2,gof2] = fit([26.3 30.8 36.1 44]',[1 0.95 0.70 0]',g,'StartPoint', [3800, -0.5]);
    %Matlab fit Kts_upper
        X = [26:1:50]'; 
        g = fittype( @(Q, B, x) (1-(1 ./ (1 + Q*exp(B*(x))))) )
        [curve2,gof2] = fit([26.3 30.8 36.1 44]',[1 0.95 0.70 0]',g,'StartPoint', [40000000, -0.4]);
    %}
    
%Final parameter values
Q_lower = 3756;
B_lower = -0.526;

Q_upper = 40000000;
B_upper = -0.463;

x_lower = [0:0.1:27]';
x_upper = [27:0.1:50]';

%Lower and upper logistic equations
y_lower = (1 ./ (1 + Q_lower*exp(B_lower*(x_lower))));  %lower temp 3800, -0.5
y_upper = 1-(1 ./ (1 + Q_upper*exp(B_upper*(x_upper))));  %upper temp 40,000,000, -0.46


Kts_x = [x_lower; x_upper]; 
Kts_y = [y_lower; y_upper]; 
%Assign Kts values between temperatures of 25 and 29 a value of 1
    keep = Kts_x > 25 & Kts_x <= 29;
    Kts_y(keep,:) = 1;
%Normalize Kts values from below 25 degrees and above 29 degrees to values
%between 0 and 1
    keep = Kts_x <= 25;
    Kts_y(keep,:) = mat2gray(Kts_y(keep,:));
    keep = Kts_x > 29;
    Kts_y(keep,:) = mat2gray(Kts_y(keep,:));
    Kts_curve = [Kts_x Kts_y];

%{    
%Plot logistic model and compare agains literature values
figure(1)
hold on
%Compare lower and upper logistic equations to literature values
    %note that y values are the fraction of biomass produced relative to
    %biomass produced at optimal temperatures.
    
    %Kakani and Reddy (2007)
    scatter([10 17 21.5 26.3 30.8 36.1 44]',[0 0.65 0.80 1 0.95 0.70 0]','ok')
    %Delucia et al. (1992)
    scatter([5.9 10.1 15.2 20.2 24.5 29.7 33.9 35.5]',[0.06 0.11 0.28 0.56 1 0.88 0.67 0.39]','xk')

%Plot upper and lower Kts
plot(Kts_curve(1:271,1),Kts_curve(1:271,2),'b');
plot(Kts_curve(272:end,1),Kts_curve(272:end,2),'r');    
    
box on
legend('Kakani and Reddy (2007)', 'Delucia et al. (1992)','Kts lower', 'Kts upper', 'location','northwest')
legend boxoff
xlabel('\bf Temperature (degrees C)')
ylabel('\bf Kts')
xlim([0 50]) 
ylim([0 1]) 
%}
    
%% Determine Kts for each day
%Fill missing data by linear interpolation
    %if missing data is < 10% of total data points
    if sum(isnan(weather_year(:,3))) < 0.1*size(weather_year(:,3),1)
        weather_year(:,3) = naninterp(weather_year(:,3));
        weather_year(:,4) = naninterp(weather_year(:,4));

        %Calculate average temperature
        Tmean = (weather_year(:,3) + weather_year(:,4)) / 2;

        %Determine Kts by comparing temperature with Kts curve
        for i = 1:size(Tmean,1)
            [~,idx]=min(abs(Tmean(i)-Kts_curve(:,1)));
            Kts(i,1) = Kts_curve(idx,2);
        end
        
    else
        Kts = ones(size(weather_year,1),1);
    end
