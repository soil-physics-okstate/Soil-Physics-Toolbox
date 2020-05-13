function [Q, Jw, Vmax, hplot] = poiseuille (R, deltaP, L, v)
%POISEUILLE is a function based on Poiseuille's law to calculate water flow.
%
% Poiseuille calculates flow per unit
% of time (Q), and flow per unit of time and area (Jw) using
% Poiseuille's law. Q can also be called the flow rate.
%
% This law describes the water flow through a system of an ideal 
% geometry such as a capillary tube or bundle. The use of a capillary bundle is useful 
% to start understanding water flow in porous media under saturated conditions. 
%
% In short, Poiseuille's law  says that for a given pressure differential
% deltaP across a cylinder of length L, the volume of water flowing per
% unit of time (Q) is proportional to the fourth power of the radius. This
% is, when you double the size of the bundle radius, water flow will not double, but will
% increase 16 times.
%
%% Inputs
%       R = outter capillary tube radius [m]
%       deltaP = hydrostatic pressure difference across the capillary [m]
%       L = capillary bundle length [m]
%       v = coefficient of viscocity of the fluid [kg/m s]
%% Outputs
%       Q = volume fluid flow per unit of time [cm3/h].
%       Jw = volume fluid flow per unit of time and unit of area. This term
%       is also known as fluid flow rate per unit area [cm/h].
%
%% Example
% [Q, Jw, Vmax, hplot] = poiseuille ( 1e-3, 0.6, 0.05, 1);
%

% Andres Patrignani. 28-Aug-2013 11:49:11
    
%% Check user input error
%if ~isvector(R)    # IGNORE THIS LINE
%    error('R must be a vectors')    # IGNORE THIS LINE
%end    # IGNORE THIS LINE

%% Calculate parabolic velocity profile.
r = linspace(0,R,100); % Create arbitrary divisions of the bundle radius. 

Vr = (deltaP).*(R.^2-r.^2) ./ (4*L*v);  % eq 3.7. Page 76.
% Eq. 3.7 calculates the velocity of the fluid at each of the points (r)
% along the bundle radius.

%% Maximum fluid velocity Vmax through the capillary bundle 
% (r=0--> fluid velocity in the middle of the bundle).
Vmax = (deltaP).*R.^2 ./ (4*L*v);   % eq 3.8 Page 76. 

%% Poiseuille's law.
Q = (pi*R.^4.*deltaP) / (8*L*v);  % eq 3.10. Page 77.

Q = Q * 998 * 9.81; % Q * density of water in [kg/m3] * acceleration of gravity [m/s2]

Q = Q * 100.^3 *60 *60; % Q * cubic centimeters in a cubic meter * 60 sec in a min * 60 min in an hour.
% Therefore, Q has units of [cm3/h].

%% Volume flow rate per unit of area.
% Flow equations are usually expressed in volume of water flow per unit of
% time and area. For this, we can divide both sides by pi and R2 (since
% the area of a circle is pi*R2), resulting Jw=Q/pi*R^2, and therefore:

Jw = (R.^2.*deltaP) ./ (8*L*v); % 3.11. Page 77.

%% Plot of parabolic velocity distribution in the capillary tube.
radialDistance = r./R;

hplot=plot(radialDistance, log10(Vr)); % Generate plot
hold all; % Retain current graph for new graphs
xlabel('Radial Distance r / R','FontSize',14); % x axis label
ylabel('Log Water Velocity V(r) (m s^-^1)','FontSize',14); % y axis label

%legend('v = 1 kg m^{-1} s^{-1}','v = 10 kg m^{-1} s^{-1}', 'Location','SouthWest')
%legend boxoff
%saveas(hplot,'chartvisc.tif')
end