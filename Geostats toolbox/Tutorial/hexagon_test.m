%Script to process data from soil moisture survey at MOISST
%TEO
%Last modified 11/29/2011

%Load the data
data = xlsread('marena_hexagon.xlsx');
c = data(:,1:2); %lat-lon coordinates
z = data(:,3)/100; %convert vol. soil moisture from % to a fraction

%Visualize the statistical distribution of the data
%figure
%normplot(z); %normal probability plot for the data shows suspect data at low end
threshold = 0.08; %minimum plausible value for the water content data
j = z<threshold; %find anomolous low water contents
z(j,:)= [];%delete those rows
c(j,:)= [];
%normplot(z); %data distribution after suspect points removed

%% Visualize the data in space for quality control
figure
colorplot(c,z,'hot') %one of the measured points falls outside the study area
[x,i] = min(c); %find the minimum value for lat and lon and the corresponding index
c(i(2),2)= 97.22365; %correct coordinate for out of range coordinate
c(:,1) = c(:,1)*110960; %convert latitude to meters (conversion valid at 36.06 degrees lat)
c(:,2) = c(:,2)*90095; %convert longitude to meters (conversion valid at 36.06 degrees lat)
colorplot(c,z,'hot') %now the points are in the correct locations
axis equal %make axes proportional

%% Calculate the empirical variogram
figure
c1a = [0:10:20]'; %lag classes for small distances (m)
c1b = [120:300:1200]'; %lag classes for large distances (m)
c1 = [c1a;c1b]; % combined lag classes
method = 'kron'; %method for calculating distances
[d,V,o] = vario(c,z,c1,method,1); %calculate and plot the empirical variogram

%% Fit spherical variogram model with nugget to the empirical data
model = {'nuggetV','sphericalV'}; %specify the model
paramest = {[0.0005] [0.0015 300]}; %nugget, sill, and range (m)
[param,fval,exitflag,output]=modelfit(d,V,o,model,paramest);
hold on
d2 = (1:10:1200)'; %lag classes for model plotting (m)
modelplot(d2,model,param); %add model to empirical variogram plot

%% Make spatial predictions based on measured data and the variogram model
minc = min(c); %specify the origin of the prediction grid
dc = [30 30]; %spacing between points on the prediction grid in each direction
nc = (max(c)-min(c))./dc; %number of points on the grid in each direction
nc = floor(nc); %round down number of points to nearest whole number
ck = creategrid(minc,dc,nc); %create coordinates matrix for prediction grid
nhmax = 3; %maximum number of data points considered for spatial prediction at each location
dmax = 300; %maximum distance allowed between spatial prediction location and data point
order = 0; %order of the polynomial mean along the spatial axes
tic;
[zk,vk]=kriging(ck,c,z,model,param,nhmax,dmax,order); %make spatial predictions by kriging
toc

%% Map the spatial predictions
[ck1,ck2,zk2] = col2mat(ck,zk); %reformats data for use by pcolor
figure; %create new figure window
ck1 = ck1/110960; %convert back to degrees latitude
ck2 = ck2/90095; %convert back to degrees longitude
pcolor(ck1,ck2,zk2); %create colored plot
colormap(hot); %specify color scheme
colorbar; %add legend
axis equal %make axes proportional
title('Surface Soil Moisture, MOISST, 8/16/2011')
xlabel('Latitude')
ylabel('Longitude')


%% -------- Test Custom Functions ----------------------

% Test semivar.m
X = c(:,1);
Y = c(:,2);
Z = z;
lagbins = c1;
[gamma,h,cvar,rho] = semivar(X,Y,Z,lagbins);
figure,scatter(h,gamma);hold all

% Test semivarmod.m
modelsv = 'nuggetV';
x0 = [0.0015 150];
[Gh,beta] = semivarmod(x0,h,gamma,modelsv);
plot(Gh)

%% Test krig.m
maxpoints = 3;
maxdist = 300;
method = 'ordinary';
tic;
[Zpred] = krig([c,z],ck,beta,maxpoints,maxdist,method);
toc
[CK1,CK2,ZK2] = col2mat(ck,Zpred);

figure; %create new figure window
CK1 = CK1/110960; %convert back to degrees latitude
CK2 = CK2/90095; %convert back to degrees longitude
pcolor(CK1,CK2,ZK2); %create colored plot
colormap(hot); %specify color scheme
colorbar; %add legend
axis equal %make axes proportional
title('Surface Soil Moisture, MOISST, 8/16/2011')
xlabel('Latitude')
ylabel('Longitude')