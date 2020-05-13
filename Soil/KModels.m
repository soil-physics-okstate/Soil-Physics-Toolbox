function [KvG,KBC, KC] = KModels (vGparameters, variable, variableData)
%% Soil water retention models
%
%   Model       |         Description
%               |
%   'VG'        |       van Genuchten, 1980.
%   'CB'        |       Campbell, 1974
%   'BC'        |       Brooks and Corey, 1964
%
%% Model parameters (par)
%    Model VG:
%               par(1) = alpha (constant)
%               par(2) = n (constant)          
%               par(3) = porosity
%               par(4) = residual water content (theta_r).
%
%    Model CB:
%               par(1) = porosity
%               par(2) = bubbling pressure (hb).               
%               par(3) = constant (b).
%
%    Model BC:
%               par(1) = pore-size distribution (lambda)
%               par(2) = bubbling pressure (hb). 
%               par(3) = porosity.
%               par(4) = residual water content (theta_r).
%               

        

end