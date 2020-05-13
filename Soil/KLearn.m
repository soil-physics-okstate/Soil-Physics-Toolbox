%% Script that shows the hydraulic conductivity for two different soils using 
% the empirical Mualem-van Genuchten model.
% SeeFigure 3.9

%% Soil parameters
alpha=[0.050 0.01]; % 1/cm.
n=[3.1 1.2]; % unitless
theta_s=[0.375 0.519]; % cm3/cm3
theta_r=[0.05 0.09]; % cm3/cm3
Ks =[30.779 0.515]; % cm/hr
m=1-1./n; 
h=abs(0:-1:-1500); % absolute value of matric potential in kPa

%% Mualem-van Genuchten K(h) Model

for i=1:length(alpha)
    K(i,:) = Ks(i)*((1-(alpha(i)*h).^(n(i)-1).*(1+(alpha(i)*h).^n(i)).^-m(i)).^2) ./ (1+(alpha(i)*h).^n(i)).^(m(i)/2); % Mualem-van Genuchten
    intArea(i,1) = trapz(K(i,:));   % approximation using trapezoidal numerical integration
end

%% Integration methods - Analytical
%syms h;
%K(i,:) = Ks(i)*((1-(alpha(i)*h).^(n(i)-1).*(1+(alpha(i)*h).^n(i)).^-m(i)).^2) ./ (1+(alpha(i)*h).^n(i)).^(m(i)/2);
%int(K, h, 0, 15000) % Integration over an arbitrary number instead of
%infinity. No matter what the boundy of the integration is, for this case
%takes lot of time and may crush your computer!

%% Plot
hFig = figure;
subplot(1,2,1) % Create first subplot in figure
plot(log10(h), log10(K)); hold on % hold on mantains the plot active for a new graph.
legend('Sand','Clay'); % 
legend boxoff % Remove box around legend
ylim([-12 2.1]) % set limits of y axis
xlabel('Log \psi_m','FontSize',13) % x axis label and font size. \psi displays the Greek letter on screen.
ylabel('Log Hydraulic Conductivity (cm h^-^1)','FontSize',13) % y axis label
title('Log scale Y axis','FontSize',14) % Chart title.

subplot(1,2,2) % generate second subplot of the figure.
[AX]=plotyy(log10(h), K(1,:),log10(h), K(2,:)); hold on % generate plot with two y axis using plotyy function.
legend('Sand','Clay');
legend boxoff
set(hFig, 'Position',[70 100 1500 600]) % set figure position in the screen and size
title('Actual scale Y axis','FontSize',14) 
disp(['Sand Area = ', num2str(sprintf ('%.2f', intArea(1))),' cm']); % Display text in command windows.
disp(['Clay Area = ',num2str(sprintf ('%.2f', intArea(2))),' cm']);
xlabel('Log \psi_m','FontSize',13) % x axis label
ylabel(AX(1),'Hydraulic Conductivity Sand (cm h^-^1)','FontSize',13) % Left y axis label.
ylabel(AX(2),'Hydraulic Conductivity Clay (cm h^-^1)','FontSize',13) % Right y axis label.

