function [AllParameters, phi, theta_e, psi_f, Ko] = greenamptpar (TexturalClass)
%% greenamptpar
% Returns the Green-Ampt model parameters according to soil textural classes.
%
%% Syntax
%
% |[AllParameters, phi, theta_e, psi_f, Ko] = greenamptpar (TexturalClass)|
%
%% Input
%       TexturalClass = one of the possible soil textural classes. Input must be a string.
%
%% Outputs
%       phi = total porosity [cm3/cm3]
%       theta_e = effective porosity [cm3/cm3]
%       psi_f = wetted front capillary pressure [cm]
%       Ko= hydraulic conductivity [cm/h]
%
%% Example
% |greenamptpar('Sand')|
%
% ans = 
%
%    'Textural Class'    'phi [cm^3/cm^3]'    'theta_e [cm^3/cm^3]'    'psi_f [cm]'    'Ko [cm/h]'
%    'Sand'              [         0.4370]    [             0.4170]    [    4.9500]    [  11.7800]
%
%% See also
% <GreenAmpt.html |GreenAmpt|>
%
%% References
% Rawls, W.J., Brakensiek, D.L. and Miller, N. 1983. Green-Ampt Infiltration
% Parameters from Soils Data. J Hydraul Eng-Asce 109: 62-70.
%

%%
%
%Create cell array
Headers = {'Textural Class' 'phi [cm^3/cm^3]' 'theta_e [cm^3/cm^3]' 'psi_f [cm]' 'Ko [cm/h]'};
ParameterList = {'Sand',0.437,0.417,4.95,11.78;...
                'Loamy Sand',0.437,0.401,6.13,2.99;...
                'Sandy Loam',0.453,0.412,11.01,1.09;...
                'Loam',0.463,0.434,8.89,0.34;...
                'Silt Loam',0.501,0.486,16.68,0.65;...
                'Sandy Clay Loam',0.398,0.330,21.85,0.15;...
                'Clay Loam',0.464,0.309,20.88,0.10;...
                'Silty Clay Loam',0.471,0.432,27.30,0.10;...
                'Sandy Clay',0.430,0.321,23.90,0.06;...
                'Silty Clay',0.479,0.423,29.22,0.05;...
                'Clay',0.475,0.385,31.63,0.03};
txtIdx = strcmpi(ParameterList(:,1),TexturalClass); % Find match for user's input.

% Prepare output
phi = cell2mat(ParameterList(txtIdx,2));
theta_e = cell2mat(ParameterList(txtIdx,3));
psi_f = cell2mat(ParameterList(txtIdx,4));
Ko = cell2mat(ParameterList(txtIdx,5));
AllParameters = [Headers; ParameterList(txtIdx,1:5)]; % All parameters in a matrix.

%%
% Copyright 2013 _This function is part of the Soil Physics Toolbox  created by the Soil 
% Physics team at the Plant and Soil Sciences Department, Oklahoma State University._
%

