function [AllParameters, theta_r, theta_s, alpha, n, m, Ksat] = rosetta(texturalClass)
%% ROSETTA
% Average soil hydraulic parameters for USDA soil textural classes. 
%
%% Syntax
%
% |[AllParameters, theta_r, theta_s, alpha, n, m, Ksat] = rosetta(texturalClass)|
%
% |[AllParameters, theta_r, theta_s, alpha, n, m, Ksat] = rosetta
% (texturalClass)| retrieves residual volumetric water content (|theta_r|),
% saturation (|theta_s|), saturated hydraulic conductivity (|Ksat|), and three
% fitting parameters (|alpha|, |n|, and |m|). AllParameters is a cell array
% containing all the previous variables and parameters with headers
% (including units).
%
%% Inputs
%       texturalClass = String containing the name of one soil textural class defined by the USDA.
%
%% Soil textural classes
%
% To use this chart use at least two out of the three soil texture components.
% Notice that numbers in the edges are tilted in the same way
% you should read the lines.
%
% <<rosetta_texture_triangle.jpg>>
%
% Graph obtained from USDA Natural ResourcesConservation Service (NRCS)
% <http://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/edu/?cid=nrcs142p2_054311>
%
%% Outputs
%       theta_r = residual volumetric water content (cm3/cm3)
%       theta_s = saturate volumetric water content (cm3/cm3)
%       alpha = fitting parameter. It is the inverse of the air entry suction (1/cm)
%       n = fitting parameter related to pore size distribution (unitless) (n>1)
%       m = fitting parameter (unitless) ($m = 1-1/n$)
%       Ksat = saturated hydraulic conductivity (cm/h)
%
%% Examples
%       [AllParameters] = rosetta (texturalClass)
%
% |AllParameters =| 
%
%    'Textural Class'    'theta_r (cm^3/cm^3)'    'theta_s (cm^3/cm^3)'    'alpha (1/cm)'    'n'         'Ksat (cm/hr)'
%    'silt loam'         [             0.0650]    [             0.4390]    [      0.0051]    [1.6600]    [      0.7600]
%
%% References
% Schaap, M.G., F.J. Leij, M.Th. van Genuchten. 2001. Rosetta: a computer 
% program for estimating soil hydraulic parameters with hierarchical
% pedotransfer functions. Journal of Hydrology 251(3–4)163–176.
%
%% Updates
%
% v.1 Created by AP Monday, ?September ?23, ?2013.
%
% Last revised on 10-Dec-2013 18:06:15
%
%%
% _This function is part of the Soil Physics Toolbox  created by the Soil 
% Physics team at the Plant and Soil Sciences Department, Oklahoma State University._
%

%%
 % Hard coded Rosetta table with average soil properties and code.
headers = {'Textural Class' 'theta_r (cm^3/cm^3)' 'theta_s (cm^3/cm^3)' 'alpha (1/cm)' 'n' 'Ksat (cm/hr)'};
ParameterList = {'clay',0.098,0.459,0.015,1.25,0.615;...
              'clay loam',0.079,0.442,0.0158,1.42,0.341;...
              'loam',0.061,0.399,0.0111,1.47,0.502
              'loamy sand', 0.049,0.390, 0.0348, 1.75,4.383;...
              'sand',0.053,0.375,0.0352,3.18,26.779;...
              'sandy clay', 0.117,0.385,0.0334,1.21,0.473;...
              'sandy clay loam',0.063,0.384,0.0211,1.33,0.549;...
              'sandy loam', 0.039,0.387,0.0267,1.45,1.595;...
              'silt',0.05,0.489,0.0066,1.68,1.823;...
              'silty clay',0.111,0.481,0.0162,1.32,0.401;...
              'silty clay loam',0.090,0.482,0.0084,1.52,0.463;...
              'silt loam',0.065,0.439,0.0051,1.66,0.760};
          
txtIdx = strcmpi(ParameterList(:,1),texturalClass); % Find match for user's input.

% Output variables
Ksat = cell2mat(ParameterList(txtIdx,6)); % Saturated hydraulic conductivity [cm/h]
theta_r = cell2mat(ParameterList(txtIdx,2)); % Residual volumetric water content [cm3/cm3]
theta_s = cell2mat(ParameterList(txtIdx,3)); % Saturated volumetric water content [cm3/cm3]
alpha = cell2mat(ParameterList(txtIdx,4)); % Inverse of bubbling pressure [1/cm]
n =cell2mat(ParameterList(txtIdx,5)); % Fitting parameter
m= 1-1/n; % Parameter as a function of n
AllParameters = [headers; ParameterList(txtIdx,1:6)]; % All parameters in a matrix.