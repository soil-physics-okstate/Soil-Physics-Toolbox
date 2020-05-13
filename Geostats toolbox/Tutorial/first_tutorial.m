%% Tutorial semivariogram
%
% Author: Andres Patrignani 21-Jun-2015
%
%--------------------------- TUTORIAL STEPS -------------------------------
% Step 1 - Load data
%
% Step 2 - Plot data
%
% Step 3 - Empirical semivariogram
%   Step 3.1 - Remove spatial trend
%   Step 3.2 - Define lags (check function nbins within function
%              semivarlags for more detail).
%   Step 3.3 - Compute empirical omnidirectional semivariogram 
%
% Step 4 - Fit semivariogram model
%
% Step 5 - Interpolate using kriging
%
% Step 6 - Plot results


%% Step 1 - Load example data (Wolfcamp aquifer, TX)
% Data is also available at: http://wiki.stat.ucla.edu/socr/index.php/SOCR_061708_NC_Data_Aquifer
data = importdata('wolfcamp_aquifer.txt');
X = data.data(:,1); % Horizontal distance in km 
Y = data.data(:,2); % Vertical distance in km
Z = data.data(:,3); % Pressure heads in meters

%% Step 2 - Plot data
figure
title('Raw data','Fontsize',14)
scatter(X,Y,[],Z,'filled')
box on
Hobs = colorbar;
set(get(Hobs,'YLabel'),'string','Pressure head (m)');
xlabel('Easting','Fontsize',14)
ylabel('Northing','Fontsize',14)

%% Step 3 - Calculate empirical semivariogram

% Step 3.1 - Remove spatial trend
variates = [ones(length(Z),1) X Y];
b = regress(Z,variates);
Ztrend = variates*b;
Zres = Z-Ztrend;

% Step 3.2 - Define lags
[lagbins,N] = semivarlags(X,Y);

% Step 3.3 - Compute empirical omnidirectional semivariogram 
[Gobs,Hobs,Ch,rho,Gcloud,Hcloud] = semivar(X,Y,Zres,lagbins);
figure, scatter(Hobs,Gobs,'ob','filled'); hold on
scatter(Hcloud,Gcloud,'.r'); hold on

%% Step 4 - Fit model to empirical data
sill0 = var(Zres);
range0 = 50;
beta0 = [sill0,range0];
% Gmodel = 'gaussV';
[Gmod,Hmod,beta,Gpars,Gmodel] = semivarfit(beta0,Hobs,Gobs);
plot(Hmod,Gmod,'-k','LineWidth',2); hold on
xlabel('Lag distance (m)','Fontsize',14)
ylabel('Semivariance','Fontsize',14)
box on
legend('Empirical semivariance','Fitted semivariance model','Location','SE')
legend boxoff
ylim([0 Gpars(3)+Gpars(3)*0.2])
plot(Hmod,Gpars(1)*ones(size(Hmod)),'--r'); hold on % Nugget
plot([Gpars(2) Gpars(2)],Gpars([1,3]),'--g'); hold on % Range
plot([0 max(Hobs)],[Gpars(3) Gpars(3)],'--b'); hold on % Sill

%% Step 5 - Interpolate pressure head values at different locations
obs = [X,Y,Z];
ngrid = 50; % Number of total point in each dimension to build the meshgrid.
[Xgrid,Ygrid] = meshgrid(linspace(min(X),max(X),ngrid), linspace(min(Y),max(Y),ngrid)); 
pred = [Xgrid(:) Ygrid(:)];
maxpoint = 5;
maxdist = 200; %Equal to range
Kmodel = 'ordinary';
[Zpred,Vpred] = krig(obs,pred,beta,maxpoint,maxdist,Gmodel,Kmodel);
Zgrid = reshape(Zpred,size(Xgrid));
Zgrid = round(Zgrid);
Vgrid = reshape(Vpred,size(Xgrid));

%% Step 6 - Plot results
figure
contourf(Xgrid,Ygrid,Zgrid,8,'ShowText','on',...
                                 'LabelSpacing',288); hold on
title('Wolfcamp Aquifer Pressure head (m) Spatial Distribution');hold on
xlabel('Easting','Fontsize',14)
ylabel('Northing','Fontsize',14)
box on

figure
contourf(Xgrid,Ygrid,Vgrid,3); hold on
title('Wolfcamp Aquifer Variance');hold on
xlabel('Easting','Fontsize',14)
ylabel('Northing','Fontsize',14)
box on
scatter3(X,Y,Z,'ow')