% This script shows the relationship between water content and
% matric potential by different methods.

%% Parameters (matric potential are expresses ad absolute values).
% Parameters for a s silt loam used by Or et al. (1991).
h=0:1:1500;

% Brooks and Corey, 1964.
lambda = 0.54; % parameter related to the pore size distribution
h_bubbling_suction = 1.48; % unit is [m]. Parameter related to the soil matric potential.

% van Genuchten model, 1980.
alpha=0.417; n=1.75; theta_s=0.513; theta_r=0.05; % van Genuchten parameters.
m=1-1/n; % m parameter as a function of n.

%% Soil water retention curve (SWRC). Brooks and Corey (BC)
theta_v_BC = zeros(length(lambda),length(h)); % Pre-allocate memory for speed.
for i=length(lambda)
    theta_v_BC(i,:) = min( (h_bubbling_suction(i)./h).^lambda(i) .* (theta_s(i)-theta_r(i))+theta_r(i), theta_s);
end

%% Soil water retention curve (SWRC) van Genuchten(VG).
theta_v_VG = zeros(length(lambda),length(h)); % Pre-allocate memory for speed.
for i=1:length(alpha)
    theta_v_VG(i,:) =  ( 1 ./ (1 + (alpha(i)*h).^n(i) )).^m(i)  * (theta_s(i)-theta_r(i)) + theta_r(i);
end

%% Plot
% h(thetav)
hFig = figure; % create figure 
subplot (2,1,1) % create first subplot
plot(theta_v_VG, log10(h),'-k',theta_v_BC, log10(h),'--k');
xlim([0 0.55]); % x axis limits
title ('Silt loam soil water retention curves','FontSize',14)
legend('van Genuchten','Brooks and Corey')
legend boxoff
xlabel('Volumetric water content (cm^3 cm^-^3)','FontSize',13)
ylabel('Log \psi_m (-kPa)','FontSize',13)

% thetav(h)
subplot (2,1,2)
plot(log10(h), theta_v_VG,'-k',log10(h), theta_v_BC,'--k');
legend('van Genuchten','Brooks and Corey')
legend boxoff
xlabel('Log \psi_m (-kPa)','FontSize',13)
ylabel('Volumetric water content (cm^3 cm^-^3)','FontSize',13)

set(hFig, 'Position', [400 50 400 700])
