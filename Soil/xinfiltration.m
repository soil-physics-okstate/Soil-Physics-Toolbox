function [S, d, alpha, n] = xinfiltration (theta_s,theta_r,theta_ini,Ksat,distance,theta_v,time,cumInfiltration)
% XINFILTRATION
% Horizontal infiltration model.
%
%% Description
% This function calculates the sorptivity, alpha and n parameters of the 
% van Genuchten (1980) soil water retention model based on a horizontal infiltration 
% experiment using the integral solution proposed by Shao and Horton, 1998. Sorptivity 
% is calculated based on cumulative infiltration using Philip's infiltration model. 
%
%% Syntax
% |[S, d, alpha, n] = xinfiltration (theta_s,theta_r,theta_ini,Ksat,distance,theta_v,time,cumInfiltration)|
%
%% Inputs
%       |theta_s| = Saturated water content (cm3 cm-3)
%       |theta_r| = Residual water content (cm3 cm-3)
%       |theta_ini| = Initial water content (cm3 cm-3)
%       |Ksat| = Saturated hydraulic conductivity. (cm h-1)
%       |distance| = Distance of the wetting front tracked through the trasnparent tube [cm]
%       |theta_v| = Water content of each section of the horizontal tube [cm3 cm-3]
%       |time| = Time while recording cumulative infiltration [s]
%       |cumInfiltration| = Cumulative infiltration as measured from the Mariotte's bottle [cm]
%
%% Outputs
%       S = Sorptivity
%       
%       d = distance equal to Boltzmann variable (mm/s^0.5) st which
%
%       theta(d) = theta_initial. In other words we are trying to find the
%       wetting front, which is likely to be right in the transition from
%       wet soil to a soil that has not been wetted yet.
%       
%       alpha = van Genuchten SWR curve scaling parameter inversely proportional
%               to the mean pore diameter.
%       
%       n = van Genuchten SWR curve shape parameter related to the pore-size distribution.
%
%% See also
% <GreenAmpt.html |GreenAmpt|>
%
%% References
% Shao, M. and Horton, R. 1998. Integral method for estimating soil Hydraulic 
% properties. Soil Sci. Soc. Am. J. 62: 585-592.
%
%% Updates
% Created by Andres Patrignani     10-Oct-2013

Ksat = Ksat/3600; % Convert Ksat from cm/h to cm/sec

% Remove NaNs
distance(isnan(distance)) = [];
theta_v (isnan(theta_v)) = [];
time(isnan(time)) = [];
cumInfiltration(isnan(cumInfiltration)) = [];

% Step 1 Calculate sorptivity
Timeroot = time.^0.5; % Calculate the square root of time.

% Fit data forcing line through y=zero;
S=Timeroot\cumInfiltration;

% A less efficient alternative is to use the regress function.
% x = [zeros(size(Timeroot)) Timeroot]; %independent variables
% [b] = regress(CumInf,x); %determine coefficients
% For details in the statistical analysis use: [b,bint,r,rint,stats] = regress(CumInf,x); 
% S = b(2); % b(1) is the intercept which was forced to be 0, since we start at time 0.

% Step 2 Calculate wetting front distance (d)
Xf = interp1(theta_v, distance,theta_ini); % Position of the wetting front (Xf). 
% Ensure that follow the same order of the input arguments. Since we want 
% to find distance based on knowledge of distance and water content we need 
% to invert this in the function, otherwise we will have NaN as a result.
d = Xf/(Timeroot(end)); % Characteristic wetting length.

% Step 3 Calculate n
n = S/(d*(theta_s-theta_ini)-S); % Shape parameter related to the pore-size distribution. [eq 26] in Shao and Horton, 1998

% Step 4 Calculate alpha
m = 1-1/n; % fitting parameter.
Norm_theta = (theta_s-theta_ini)/(theta_s-theta_r); % Normalize water content.
alpha = (2*Ksat)/(S*d)*(1/m*Norm_theta)^(1/n); % Scaling parameter inversely 
% proportional to the mean pore diameter. [eq 25] in Shao and Horton, 1998.
% [eq 19] yields the same result as [eq 25].