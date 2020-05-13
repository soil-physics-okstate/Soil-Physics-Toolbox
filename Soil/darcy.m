function [Ks, Jw] = darcy (Q, A, p1, b, z1, L)
%DARCY is a function that describes vertical macroscopic water flow in saturated porous media. 
% Assumption:
%           1) The total potential difference between both points in known.
%           2) Rigid soil.
%           3) Soil is fully saturated.
%           4) No solute membranes are present.
%           5) No other variable such as roots or temperature gradients
%           drive flow. Only the difference in hydraulic head.
% Under these conditions the total water potential is only composed by the
% hydrostatic pressure (p) and the gravitational potential (z). If head units are
% used, then the combination of both is called the hydraulic head, in order
% that H=p+z.
%
%% Inputs
%       Q = Flow. Constant in non-swelling soils [cm^3/h]. Negative when
%       downward flow.
%       A = area [cm2]
%       deltaP = hydrostatic pressure difference. Hydrostatic pressure
%       arises from an unsupported column of water [cm]. 
%       L = length of the column [cm].
%       deltaP = deltaH = Hydraulic head
%
% Upward flow is positive (i.e. evaporation), and downward flow is negative
% (i.e. drainage)
%
%% Outputs
%       Q = flow rate [cm3/d].
%       Jw = Soil water flux [cm/d].
%
% Example for a vertical soil column of 150 cm in length and area of 50 cm^2 that
% has a flow of 2000 cm^3 per hour (downward). A constant water level of 15
% cm is continuously ponded on it.
% 
% [Ks, Jw] = darcy (-2000, 50, 0, 15, 0, 150)

% Andres Patrignani. 02-Sep-2013 15:56:15

%% Calculate hydraulic head.
z2=L;
p2=b; % Hydrostatic pressure on soil surface.
H1 = p1+z1; % Hydrostatic pressure (p) + gravitational potential (z). If horizontal flow, then z=0.
H2 = p2+z2;
deltaH = H2 - H1; % Hydraulic gradient.

%% Flux form of Darcy's law
Jw = (Q ./ A); % units are cm/d
disp(['Flux is ',num2str(Jw),' cm/h']);

%% Saturated hydraulic conductivity
Ks = - Jw * L ./ deltaH ;
disp(['Ksat is ',num2str(Ks),' cm/h']);

end