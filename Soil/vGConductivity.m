function [KvG, h] = vGConductivity (TexturalClass)
%% Soil parameters
[~, ~, ~, alpha, n, m,Ksat] = vGparameter (TexturalClass);
h=abs(0:-1:-1500); % absolute value of matric potential in kPa.

%% Mualem-van Genuchten K(h) Model
KvG = Ksat*((1-(alpha*h).^(n-1).*(1+(alpha*h).^n).^-m).^2) ./ (1+(alpha*h).^n).^(m/2); % Mualem-van Genuchten

%% Example
% [KvG, h] = vGConductivity (30.779,'sand')
% plot(log10(h),KvG)
end