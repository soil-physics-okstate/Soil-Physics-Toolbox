function [K_Marshall, K_Jury, h, R, theta_s, theta_r] = KCapillaryModel(SoilTexturalClass, CapillaryClasses, tau)
% Capillary model for Hydraulic conductivity. T.J. Marshall. 1958.
% This model does not account for tortuosity.
%
%% Input
% SoilTexturalClass = One of the possible soil textural classes in
% vGparameter function. Input must be a string.
% CapillaryClasses = number of capillary bundles associated with the radius
% at each matric potential.
% Tau = tortuosity (constant). unitless.
%
%% Output
% K = Hydraulic conductivity [cm/h]
% h = matric potential [kPa]
% R = Maximum capillary radius filled with water. Any bundle with radius > R(i)
% will drain at h(i). [micrometers].
% theta_s = saturated volumetric water content [cm3/cm3]
% theta_r = residual volumetric water content [cm3/cm3]
%
%% Example [KCap, hCap] = KCapillaryModel('sand', 50);
%

%% Obtain van Genuchten parameters for the selected soil texture.
[~, theta_r, theta_s, alpha] = vGparameter (SoilTexturalClass);

%% Create an arbitrary set of volumetric water content linearly spaced. The model does not work at saturation.
theta_v=linspace(theta_s-1e-10,theta_r+1e-10,CapillaryClasses); % theta_s and theta_r were subtracted and added a value of 0.01 cm3/cm3 respectively to avoid infinity.

%% Obtain matric potential for the theta_v values previously created.
[output] = SWRCfun (SoilTexturalClass,'volumetric',theta_v);
h=output.h*10.197; %kPa to cm
h(h <= 1/alpha)= []; % eliminate matric potential values higher (closer to zero) than bubbling pressure.
CapillaryClasses=length(h); % updated CapillaryClasses.

%% Calculate hydraulic conductivity
K_Marshall = zeros(CapillaryClasses,1); % Pre-allocate memory.
for i=1:CapillaryClasses % Loop through every capillary class.
    multiplier=1:2:length(h(i:end))*2-1; % Generate multipliers. Multipliers come from simplification of equation 3a in Marsahall, T.J., 1958.
    K_Marshall(i,1) = 2700*3600* theta_v(i).^2 * length(h(i:CapillaryClasses)).^(-2) .* sum( (1./(h(i:CapillaryClasses).^2)).*multiplier); % Capillary model 2700 converts cm2 to cm/sec and 3600 converts cm/sec to cm/h.
    multiplier=[]; % reset multiplier in each loop since it has a different dimension in each loop (in every loop is shorter)
    R(i) = (2*0.0728)/(1000*9.81*h(i)/100) * 1e6; % *10 to obtain capillary radius in micrometers.

    delta_v=theta_v(1)-theta_v(2);
    T=20;
    viscosity = 9.62e-7*exp(2046/(273.15+T));
    K_Jury(i) = (tau*0.0728.^2 * delta_v) / (2*viscosity*1000*9.81)* sum (1./(h(i:CapillaryClasses)/100).^2)*100*3600; % converted from m/s to cm/h.
end

%% Convert matric potential from cm to kPa
h = h/10.197; % Convert back to kPa.

%% Plotting K
%figure
%plot((log10(h)),K)
%xlim([0 4]) % Set x axis limits.
%ylim([0 max(max(K),10)]) % set y axis limits.
%xlabel('Log |\psi_m| (kPa)') 
%ylabel('Hydraulic Conductivity (cm h^{-1})')

end