% Script to compare the Mualem-van Genuchten model and the capillary bundle
% model to estimate hydraulic conductivity.

%% Set parameters of script
SoilTexturalClass = 'loamy Sand';
CapillaryClasses = 1000;
tortuosity = 1;

%% Run Mualem-van genuchten model
[KvG, hvG] = vGConductivity (SoilTexturalClass);

%% Run Capillary bundle model
[KCap, KCap2, hCap, R] = KCapillaryModel(SoilTexturalClass, CapillaryClasses, tortuosity);

%% Plot Hydraulic conductivity
figure, subplot(1,2,1),
plot(log10(hvG),KvG,'-k');hold all
plot(log10(hCap),KCap,'--k');hold all
plot(log10(hCap),KCap2,'.k');hold all
title(SoilTexturalClass,'FontSize',14) % Chart title.
legend('Mualem-van Genuchten Model','Capillary Bundle Model (Marshall, 1958)', 'Capillary Bundle Model (Jury and Horton)'); % 
legend boxoff % Remove box around legend
xlabel('Log |\psi_m| (kPa) ','FontSize',13) % x axis label and font size. \psi displays the Greek letter on screen.
ylabel('Hydraulic Conductivity (cm h^-^1)','FontSize',13) % y axis label
xlim([0 4])

%% Plot Capillary radii as a function of matric potential
subplot(1,2,2),
plot(R,log10(hCap),'-k');
title([SoilTexturalClass,' capillary radii'],'FontSize',14) % Chart title.
ylabel('Log |\psi_m| (kPa) ','FontSize',13) % x axis label and font size. \psi displays the Greek letter on screen.
xlabel('Capillary radius (\mum)','FontSize',13) % y axis label
ylim([0 5])
