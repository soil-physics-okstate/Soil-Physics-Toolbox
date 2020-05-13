function [theta_10,theta_33,theta_1500,thetaS,Ks,K_theta,theta,matric] = soilwaterptf(sand,clay,om,varargin)
%% SOILWATERPTF
% Estimates soil water hydraulic properties.
%
%% Syntax
%
% |[theta_33,theta_1500,thetaS,Ks,K_theta,theta,matric] = soilwaterptf(sand,clay,om)|
%
% |[theta_33,theta_1500,thetaS,Ks] = soilwater ptf(sand,clay,om)|
% Given the fraction of |sand|,|clay|, and organic matter (|om|) the 
% function estimates soil water retained at tensions of 33 kPa(|theta_33|), 
% 1500 kPa (|theta_1500|), saturation (|thetaS|), and saturated hydraulic 
% conductivity (|Ks|). Dimensions of |theta_33|, |theta_1500|, |om|, and
% |Ks|_ are equal to inputs dimension.
%
% |[...,K_theta,theta,matric] = soilwater ptf(sand,clay,om)| estimates the 
% soil hydraulic conductivity and water retention curve from tension of 0 kPa
% to 1500 kPa. Dimensions of |K_theta|, |theta|, and |matric| are 1500 by n,
% where n is the number of columns in |sand|, |clay|, and |om|. |matric| 
% contains the soil matric potential that corresponds with each value in
% |theta| and |K_theta|.
%
%% Inputs
%       sand = row vector of soil sand content expressed as a proportion between 0 and 1.
%       clay = row vector of soil clay content expressed as a proportion between 0 and 1.
%       om   = row vector of soil organic matter expressed as a proportion between 0 and 1.
%
%% Ouputs
%      theta_33   = soil water content at tension of 33 kPa (cm3/cm3).
%      theta_1500 = soil water content at tension of 1500 kPa (cm3/cm3).
%      thetaS     = soil water content at saturation. In this function thetaS is equal to porosity (cm3/cm3).
%      Ks         = Saturated hydraulic conductivity (mm/h).
%      K_theta    = Unsaturated hydraulic conductivity (mm/h).
%      theta      = Soil water content at a specific soil matric potentials (cm3/cm3).
%      matric     = Soil matric potential corresponding to each value in K_theta (absolute value, kPa).
%      
%% Example
%
% An example to estimate the soil water characteristics of two soil is presented below.
%
% *Soil 1*: Textural class:Silt Loam Sand:15%, Clay:18%, om:1.5% (%w)
%
% *Soil 2*: Textural class:Clay Loam Sand:29%, Clay:32%, om:1.5% (%w)
%
% [theta_10,theta_33,theta_1500,thetaS,Ks,K_theta,theta,matric] = soilwaterptf([0.15 0.29],[0.18 0.32],[0.015 0.015]);
%
% figure
% [AX,H1,H2] = plotyy(theta,matric,theta,log10(K_theta)); % AX(1)-->left axes, AX(2)--> right axes.
% set(get(AX(1),'Ylabel'),'String','Tension (kPa)','FontSize',13) 
% set(get(AX(2),'Ylabel'),'String','Log_{10} K(\theta) (mm h^{-1})','FontSize',13)
% 
% set(H1(1),'LineStyle','--','LineWidth',2,'Color',[1 0 0])
% set(H1(2),'LineStyle','--','LineWidth',2,'Color',[0 0 1])
% set(H2(1),'LineStyle','-','LineWidth',1.2,'Color',[1 0 0])
% set(H2(2),'LineStyle','-','LineWidth',1.2,'Color',[0 0 1])
% 
% title('Moisture vs Tension and Conductivity','FontSize',15)
% xlabel('Volumetric water content (cm^3 cm^{-3})','FontSize',13)
% legend('Silt Loam SWRC','Clay Loam SWRC','Silt Loam K(\theta)','Clay Loam K(\theta)',...
%        'Location','NorthWest')
% legend boxoff
%% 
% <<soilwaterptf_plot1.png>>
%
% _Figure showing the relationship between moisture-tension and
% moisture-conductivity for Silt Loam and Clay Loam soil textural classes.
% This figure was reproduced based on Figs. 2 and 3 from Saxton and Rawls,
% 2006._
%
%% Prediction equations
% The following tables were obtained from Saxton and Rawls, 2006.
%
% <soilwaterptf_Tables.png Detailed description of predictive equations (Click to open)>
%
%% Definitions
%
% *Field capacity*: "The content of water, on a mass or volume basis, remaining
% in a soil 2 or 3 days after having been wetted with water and after free 
% drainage is negligible". The water content at field capacity is often, and
% arbitrarily, estimated at a tension of 33 kPa. For coarse soils tensions
% of 10 kPa and even 1 kPa have also been reported in the literature.
%
% *Permanent wilting point*:"The largest water content of a soil at which 
% indicator plants, growing in that soil, wilt and fail to recover when 
% placed in a humid chamber. Often estimated by the water content at 
% tension of 1500 kPa soil matric potential.
%
% *Plant available water (PAW)*: The difference between any soil moisture
% content and the soil water content at -1500 kPa. This variable is assumed
% to represent the amount of water currently available for plants.
%
% *Plant available water capacity*: "The amount of water released between 
% in situ field capacity and the permanent wilting point." This variable
% is assumed to represents the maximum soil moisture available for plants 
% in a given soil at a given depth.
%
% Definitions were obtained from the Soil Science Society of America Glossary.
% For more definitions, please visit:
% <https://www.soils.org/publications/soils-glossary/>
%
%% Limitations
% A total of 1722 soil samples from the Soil Characterization database 
% (USDA/NRCS) were used by Saxton and Rawls, 2006 to developed the 
% prediction equations. Samples with bulk density <1.0 or
% >1.8 g cm^{-3} were not considered. Additionally, soil samples with
% organic matter content >8% by weigth and samples with >60% clay were
% not included in the regression analysis.
% It is also worth to notice that in some soil textures the soil water
% retention at 33 kPa of suction may not well represent the concept of
% "field capacity" condition. Soil water retention at any suction can be
% obtained from the soil water retention curve.
%
%% References
% Rawls, W.J., Brakensiek, D.L. and Saxton, K.E. 1982. Estimation of Soil-Water
%        Properties. T Asae 25: 1316-&.
%
% Saxton, K.E., W.J. Rawls. 2006. Soil water characteristic estimates by
%        texture and organic matter for hydrologic solutions. Soil Sci. Soc. Am. J.
%        70:1569-1578. doi:10.2136/sssaj2005.0117
%
% Saxton, K.E., Rawls, W.J., Romberger, J.S. and Papendick, R.I. 1986. Estimating 
%        Generalized Soil-Water Characteristics from Texture. Soil Sci Soc Am J 50: 1031-1036.
%
%% Updates
%
% v.1 Created by AP March 10, 2013.
%
% v.2 AP December 5, 2013 Added density effects, moisture tension,
% moisture-conductivity equations.
%
% Last revised on: 06-Dec-2013 11:26:04
  

%%
%

% Error checking
if length(sand) ~= length(clay) || length(sand) ~= length(om)
    error('sand, clay, and om dimensions disagree');
end
if ~isrow(sand)
    sand = sand';
end
if ~isrow(clay)
    clay = clay';
end
if ~isrow(om)
    om = om';
end

%% Moisture regressions
% * Permanent wilting point.
% Calculation of soil water retention at -1500 kPa of tension (Table 1. Eq. 1)
theta_1500t = -0.024.*sand + 0.487.*clay + 0.006.*om + 0.005.*(sand .* om) - 0.013.*(clay .* om)+ 0.068.*(sand .* clay) + 0.031; 
theta_1500 = theta_1500t + (0.14.*theta_1500t - 0.02);

% * Field capacity.
% Calculation of soil water retention at -33 kPa of tension (Table 1. Eq. 2)
theta_33t= -0.251.*sand + 0.195.*clay + 0.011.*om + 0.006.*(sand .* om) - 0.027.*(clay .* om) + 0.452.*(sand .* clay) + 0.299;
theta_33 = theta_33t + 1.283.*(theta_33t).^2 - 0.374.*(theta_33t)- 0.015;

% * Saturation-Field capacity.
% Volumetric soil water content between 0 and -33 kPa of tension (Table 1. Eq. 3). 
thetaS_33t= 0.278.*sand + 0.034.*clay + 0.022.*om - 0.018.*(sand.*om)- 0.027.*(clay.*om) - 0.584.*(sand .* clay) + 0.078;
thetaS_33= thetaS_33t +(0.636.*thetaS_33t - 0.107);

% * Tension at air entry suction.
% Estimation of the tension at air entry (Table 1. Eq. 4).
psi_et = -21.67*sand - 27.93*clay - 81.97*thetaS_33 + 71.12*sand.*thetaS_33 + 8.29*clay.*thetaS_33 + 14.05*sand.*clay + 27.16;
psi_e = max(psi_et + 0.02*psi_et.^2 - 0.113*psi_et - 0.7,0.1); % For some sandy soil generates negative values. A default value of 1 kPa was assumed as a minimum air entry suction.

% * Saturation.
% Calculation of soil water retention at -0 kPa of tension (Table 1. Eq. 5)
thetaS = theta_33 + thetaS_33 - 0.097.*sand + 0.043; % Saturation at normal density

%% Density effects
if nargin == 4
    rho_DF = varargin{1}; % Allow introducing measure bulk density.
    if ~isrow(rho_DF)
        rho_DF = rho_DF';
    end
else
    rho_N= (1-thetaS)*2.65; % Particle density is assumed to be 2.65 g cm^-3. (Table 1. Eq. 6)
    DF = 1; % Compaction degree.
    rho_DF = rho_N*DF; % (Table 1. Eq. 7).
end
thetaS_DF = 1-(rho_DF/2.65); % (Table 1. Eq. 8).
theta_33_DF = theta_33 - 0.2*(thetaS-thetaS_DF); % (Table 1. Eq. 9).
thetaS_33_DF = thetaS_DF - theta_33_DF; %#ok<NASGU> % (Table 1. Eq. 10).
% thetaS_33_DF is used with equations that account for gravel, not present
% in this version. For more information see Table 1 Eqs. 19-22 in Saxton and Rawls, 2006.

%% Moisture tension
B = (log(1500) - log(33)) ./ (log(theta_33) - log(theta_1500)); % (Table 1. Eq. 15).
A = exp(log(33) + B.*log(theta_33)); % (Table 1. Eq. 14).
values = 1500;

% Create matrices for easier manipulation of data.
theta = linspacemat(thetaS_DF,theta_1500,values); % |linspacemat| is a custom 
% function that generates linear spacing in handling vectors or matrices.
matric = nan(size(theta)); % Pre-allocate matric potential array.

% First part of piecewise moisture tension 0-psi_e (air entry suction)
Idx_0_33 = bsxfun(@(x,xt1) x<=xt1,theta,thetaS_DF) &  bsxfun(@(x,xt2) x>=xt2,theta,theta_33);
Idx_33_1500 = bsxfun(@(x,xt1) x<xt1,theta,theta_33) &  bsxfun(@(x,xt2) x>=xt2,theta,theta_1500);

for i = 1:size(theta,2)
    matric(Idx_0_33(:,i),i) = 33 - ((theta(Idx_0_33(:,i),i) - theta_33(i)) ./ (thetaS_DF(i)- theta_33(i))) .* (33 - psi_e(i)); % Eq. 12.
    matric(Idx_33_1500(:,i),i) = A(i)*theta(Idx_33_1500(:,i),i).^-B(i); % Eq. 11     
end
nvalues = 5; % arbitrary number of values.
matric_0_psi_e = linspacemat(zeros(1,size(matric,2)),matric(1,:),nvalues);
matric = cat(1,matric_0_psi_e, matric);
theta = cat(1,repmat(thetaS_DF,nvalues,1),theta);

for i = 1:size(theta,2)
    if ~any(isnan(matric(:,i))) && ~any(isnan(theta(:,i)))
        theta_10(i) = interp1(matric(1:5:150,i),theta(1:5:150,i),10,'spline');
    else
        theta_10(i) = NaN;
    end
end

%% Moisture conductivity
lambda = 1./B; % (Table 1. Eq. 18).
Ks = 1930*(thetaS_DF - theta_33).^(3 - lambda); % (Table 1. Eq. 16). [mm/h]
K_theta = repmat(Ks,size(theta,1),1).*(theta./repmat(thetaS_DF,size(theta,1),1)).^(3+(2./repmat(lambda,size(theta,1),1))); % (Table 1. Eq. 17).

%% Salinity effects
% **Not enabled in this version**
%PSI_o = 36*repmat(EC,values,1); % (Table 1. Eq. 13).
%PSI_o_theta = repmat(thetaS,values,1)./theta.*PSI_o; % (Table 1. Eq. 24).
%osmotic = PSI_o_theta;

%%
% _This function is part of the Soil Physics Toolbox  created by the Soil 
% Physics team at the Plant and Soil Sciences Department, Oklahoma State University._
%