%% Geano's Example
% data = importdata('vario test mix.txt');
% Y = data(:,1); % Northing - Latitude
% X = data(:,2)*(-1); % Easting - Longitude
% Z = data(:,3) / 100; % Measured values
% idxnan = isnan(X) | isnan(Y) | isnan(Z);
% X(idxnan) = [];
% Y(idxnan) = [];
% Z(idxnan) = [];

% Tyson's Example 
data = xlsread('marena_hexagon.xlsx');
X = data(:,2)*(-1); %lon coordinates
Y = data(:,1); %lat coordinates
Z = data(:,3)/100; %convert vol. soil moisture from % to a fraction
[x,i] = min(X); %find the minimum value for lat and lon and the corresponding index
X(i)= -97.22365; %correct coordinate for out of range coordinate
threshold = 0.08; %minimum plausible value for the water content data
idx = Z<threshold; %find anomolous low water contents
X(idx,:)= [];
Y(idx,:)= [];
Z(idx,:)= []; %delete those rows

Y = Y*110960; %convert latitude to meters (conversion valid at 36.06 degrees lat)
X = X*90095; %convert longitude to meters (conversion valid at 36.06 degrees lat)

%% Visualize data
figure
title('Spatial variability of volumetric water content'); hold on
scatter(X,Y,[],Z,'filled')
box on
colormap('parula')
colorbar
xlabel('Longitude','FontSize',14)
ylabel('Latitude','FontSize',14)

%% Check for normality
alphaValue = 0.05;
normC = zscore(Z);
[h,p] = kstest(normC,'Alpha',alphaValue); % Kolmogorov-Smirnov test (standard normal by default).
figure
subplot(1,2,1), histfit(Z)
if h==0
    decision = 'Accepted';
else
    decision ='Rejected';
end
disp(['Results of Kolmogorov-Smirnov normality test => H0:',decision,' pvalue:',num2str(p),' alpha:',num2str(alphaValue)])
subplot(1,2,2), normplot(Z); box on
set(gcf,'Position',[100 100 1200 500])

%% Calculate correlation, covariance, and semivariance.
c1a = [0:10:20]'; %lag classes for small distances (m)
c1b = [120:300:1200]'; %lag classes for large distances (m)
c1 = [c1a;c1b]; % combined lag classes

[gamma,cvar,rho,h] = semivar(X,Y,Z,c1);
figure
subplot(3,1,1),scatter(h,gamma,'ok','filled');
subplot(3,1,2),scatter(h,cvar,'og','filled');
subplot(3,1,3),scatter(h,rho,'or','filled');
