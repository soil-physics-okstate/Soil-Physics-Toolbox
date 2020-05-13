function  examplePM10CA(tMEgrid,recalculate,spaceGrid,contMethod,showDataPoints,...
  showEstPoints,showCounties,showCaBound);

% examplePM10CA             - Example IV.2.1 : Annual geometric mean of PM10 in California (Jan 1,2001)
%
% Calculates BME estimates of the annual geometric mean of PM10 at points on
% a spatial grid for selected years, and then plots the corresponding map. 
% The data set used are soft probabilistic data for Y(s,t), the annual average
% of log(PM10), available at 191 monitoring stations and for 11 years, from
% 1987 to 1997. The value Y is defined as the annual average of log(PM10),
% where PM10 is the concentration in the air of Particulate Matters of diameter
% smaller than 10 micrometer, and expressed in microgram/cubic meters.
% The BME estimation is performed on the residual Y soft data (i.e. Y minus its
% mean trend), then the mean trend is added back to the Y residual estimate.
% Then the maps are constructed using exp(Y) which is the annual geometric mean
% of PM10, expressed in microgram/cubic meters. 
%
% SYNTAX :
% 
% examplePM10CA(tMEgrid,recalculate,spaceGrid,contMethod,showDataPoints,...
%   showEstPoints,showCounties,showCaBound);
% 
% INPUT :
% 
% tMEgrid     1 by k   Optional vector of the selected time (year) for which the 
%                      variable is estimated at all the spatial grid points
%                      The default value is [1988 1991 1994 1997]
% recalculate scalar   optional integer indicating wether to recalculate estimates
%                      1 to recalculate the BME estimates. 
%                      0 to use already calculated estimates.
%                      default value is 0.
% spaceGrid   string   Optional parameter defining the grid of estimation points
%                      'coarse' for coarse spatial grid covering California  
%                      'medium' for a grid covering Covering with medium resolution
%                      'fine' for a grid covering California on a fine resolution
%                      The default value is 'medium'
%                      Note: if spaceGrid=[], then the default value is used
% contMethod  scalar   optional integer indicating the method to use to plot contours:
%                      1 to plot contour lines, 
%                      2 to plot a map in color (using pcolor)
%                      The default is 2.
% showDataPoints 
%             scalar   optional integer indicating whether to show the monitoring stations:
%                      1 to show the monitoring stations, 0 ow
%                      The default is 1.
% showEstPoints   
%             scalar   optional integer indicating whether to show the estimation points:
%                      1 to show estimation points, 0 ow
%                      The default is 0.
% showCounties
%             scalar   optional integer indicating whether to show the county boundaries:
%                      1 to show Counties boundaries, 0 ow
%                      The default is 0 if contMethod=1 and 1 if contMethod=2.
% showCaBound scalar   optional integer indicating whether to show the state boundary:
%                      1 to show California state boundaries, 0 ow
%                      The default is 1.

%
% Set preferences
%
if nargin<1, tMEgrid=[1987 1996]; end;
if nargin<2, recalculate=0; end;
if nargin<3, spaceGrid='medium'; end;
if nargin<4, contMethod=1; end;     % plot a color map
if nargin<5, showDataPoints=1; end; % show the monitoring stations
if nargin<6, showEstPoints=0; end;  % show the estimation points
if nargin<7,                        % show counties if contMethod==2
  if contMethod==1
    showCounties=0;
  else
    showCounties=1;
  end
end
if nargin<8, showCaBound=1; end;    % show boundary of california

plotWhat=1;                         % plot the plot estimated val of annual PM10 geom mean

if exist('./examplePM10CAproba.mat')~=2
  disp('Annual geometric mean of PM10 in California. The first time you run this it will take some time.');
else
  disp('Annual geometric mean of PM10 in California.');
end;  
if nargin==0
  disp('Type help examplePM10CA for explanations about the syntax of this example.');
  disp('Note that you are currently using default options.');
  disp('The default options reproduce the plots of EXAMPLE 2.1, Chap IV of the following book:');
  disp(' Christakos, G., P. Bogaert, and M.L. Serre, ''Temporal GIS: Advanced Functions for');  
  disp(' Field-Based Applications'', Springer-Verlag, New York, NY, CD ROM included, 2001.'); 
  disp('However with this program you can easily produce different maps');
  disp('For example to produce color maps for 1988 1991 1994 1997 on a coarse grid, use:');
  disp('   >>examplePM10CA([1988 1991 1994 1997],0,''coarse'',2)');
  disp('Press any key to continue, or press Ctrl C to stop here.');
  pause;
end

%
% Check input variables
%
yr=1987:1997;
for i=1:length(tMEgrid)
  if sum(tMEgrid(i)==yr)~=1,
    error('tMEgrid is a vector specifiying year(s) that must be from 1987 to 1997')
  end
end
if recalculate~=0 & recalculate~=1
  error('recalculate must be equal to 0 or 1');
end
if contMethod~=1 & contMethod~=2
  error('contMethod must be equal to 1 or 2');
end

%
%  Set the covariance model and save it
%
covmodel{1}='exponentialC/exponentialC';
covparam{1}=[0.013 6 4];
covmodel{2}='exponentialC/exponentialC';
covparam{2}=[0.005 270 135];
save examplePM10CAcovmodel covmodel covparam;

%
%  Make a 3D plot of the covariance c(r,tau)
%
rLag=(0:0.5:10)';
nr=length(rLag);
tLag=0:0.25:4;
nt=length(tLag);
rLagMat=kron(rLag,ones(1,nt));
tLagMat=kron(tLag,ones(nr,1));
c1=[rLagMat(:) zeros(nr*nt,1) tLagMat(:)];
c0=[0 0 0];
Kval=coord2K(c0,c1,covmodel,covparam);
KvalMat=reshape(Kval,nr,nt);
figure;
mesh(rLagMat,tLagMat,KvalMat);
view([1 1 1]);
grid on;
text(7.5,4.7,'Spatial lag, r (Km)');
text(14,2,'Time lag, \tau (year)');
%text(4.7,4.1,'Spatial lag, r (Km)');      
%text(12.2,0.8,'Time lag, \tau (year)');
title('c_{\ity}(r,\tau)');


%
% Create the BME maps of the annual geometric mean of PM10 
%
for it=1:length(tMEgrid)
  disp(sprintf('\rProcessing year %d',tMEgrid(it)));
  filename=['examplePM10CAestBME' num2str(tMEgrid(it)-1900)];
  if exist(['./' filename '.mat'])~=2
    recalculate=1;
  end;
  if recalculate==1,
    disp(sprintf('Recalculating the map of BME estimates using a %s grid',spaceGrid));
    pause(0.01);
    estZsBME(spaceGrid,tMEgrid(it),filename)
  elseif recalculate==0,
    disp('Reading file with already calculated BME estimates');
  elseif recalculate==0,
    error('recalculate must be equal to 0 or 1');
  end;
  disp('Plotting the map of BME estimates');
  plotZsBME(filename,plotWhat,contMethod,showDataPoints,showEstPoints,...
    showCounties,showCaBound);
end;

%
% Display some soft data
%
cs=[];
load examplePM10CAproba;
iMS=73;
[val,valname,filetitle]=readGeoEAS('examplePM10CAcMS.dat');
cMS=val;
indx= find(cs{1}(:,1)==cMS(iMS,1) & cs{1}(:,2)==cMS(iMS,2) & nl>0);
figure;
hold on;
pdfscale=0.1;
for i=1:length(indx)
  tME=cs{1}(indx(i),3);
  y=limi(indx(i),:);
  pdf=probdens(indx(i),:);
  h(1)=plot(tME+pdf*pdfscale,y,'k');
  plot(tME+[0 0],[y(1) y(end)],'k');
end
axis([1988 1998 2.6 4]);
standard=log(30);
ax=axis;
h(2)=plot([ax(1) ax(2)],[standard standard],'-.k');
legend(h,'Soft probabilistic data','Max allowed by CA law');
xlabel('Time (years)');
ylabel('\itY');
title(sprintf('Monitoring Station %d',iMS));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function estZsBME(spaceGrid,tMEgrid,filename)

% estZsBME                  - BME estimation for maps of annual PM10 geometric mean (Jan 1, 2001)
%
% Calculates BME estimates of annual PM10 geometric mean at points on a
% spatial grid for selected years, and then plots the corresponding map. 
%
% SYNTAX :
%
% estZsBME(spaceGrid,tMEgrid,filename);
%
% INPUT :
%
% spaceGrid   string   Optional parameter defining the grid of estimation points
%                      'coarse' for coarse spatial grid covering California  
%                      'medium' for a grid covering Covering with medium resolution
%                      'fine' for a grid covering California on a fine resolution
%                      The default value is 'medium'
%                      Note: if spaceGrid=[], then the default value is used
% tMEgrid     1 by k   Optional vector of the selected time (year) for which the variable
%                      is estimated at all the spatial grid points
%                      The default value is [1997]
% filename    string   Optional name of file where to save the estimated values
%                      The default value is 'examplePM10CAestBME97'
%
% NOTES :
% 
% The number of soft data points used in the neighborhood estimation is 3, 4, or 5 
% depending on whether the grid is coarse, medium or fine grid, respectively.
%
% The estimated values saved in filename include the following:
% spaceGrid: a cell array containing two cells, with the x and y coord respectively
%   of the estimation points
% tMEgrid: a vector with the times (in years) of estimation
% XkBME: a nk by length(tMEgrid) array of estimated Z values, where nk=size(spaceGrid,1)
% XkErr: a nk by length(tMEgrid) array of stan. dev. err
%

%
% Set preferences
%
%nargin=0;
if nargin<1, spaceGrid='medium'; end;
if nargin<2, tMEgrid=[1997]; end;
if nargin<3, filename='examplePM10CAestBME97'; end;

plotErr=1;
plotDisBound=1;
plotEstPoints=1;

if isempty(spaceGrid)
  spaceGrid='medium';
end;

disp('Reading probabilistic data for log yearly PM10 in California');
%
% Load the soft data for this study.
%
if exist('./examplePM10CAproba.mat')~=2
  [cs,isST,softpdftype,nl,limi,probdens,filetitle]=readProba('examplePM10CAproba.dat');
  save examplePM10CAproba cs isST softpdftype nl limi probdens filetitle; 
else
  load examplePM10CAproba;
end

%
% Load other useful data
%
[val,valname,filetitle]=readGeoEAS('examplePM10CAcMS.dat');
cMS=val;
[val,valname,filetitle]=readGeoEAS('examplePM10CAtME.dat');
tME=val';
cs=cs{1};                   % Don't use index (one variable only) 
idxNonEmptyData=nl>0;       % Index of non-empty soft data point
                            % Remove empty soft data points
[cs,c2,nl,limi,probdens,nl2,limi2,probdens2]=probasplit(cs,...
  softpdftype,nl,limi,probdens,idxNonEmptyData);

%
% Select the grid at which estimation will be done
%
disp('Setting parameters');
if ischar(spaceGrid)
  switch spaceGrid
  case 'coarse',
    dx1=100;
  case 'medium',
    dx1=50;
  case 'fine',
    dx1=25;
  otherwise
    error('bad char string for spaceGrid')
  end;
  dy1=dx1;
  xmin1=min(cMS(:,1));
  xmax1=max(cMS(:,1));
  ymin1=min(cMS(:,2));
  ymax1=max(cMS(:,2));
  xgk1=[xmin1-dx1:dx1:xmax1+1.9*dx1];
  ygk1=[ymin1-dy1:dy1:ymax1+1.9*dy1];
  [xk1 yk1]=meshgrid(xgk1,ygk1);
  idx=(xk1>=-0.6250*yk1+3475 & xk1<=-0.6*yk1+4060 );
  xk=xk1(idx);                    % Apply filter to remove estimation
  yk=yk1(idx);                    % points that are useless for CA
  xk=[xk;cMS(:,1)];
  yk=[yk;cMS(:,2)];
else
  xk=spaceGrid(:,1);
  yk=spaceGrid(:,2);
end;

%
%  Set the covariance model
%
load examplePM10CAcovmodel;

%
% set the BME paramaters
%
nhmax=0;                % max number of hard data
switch spaceGrid
case 'coarse', 
  nsmax=3;              % max number of soft data
  maxpts=50000;         % number of function eval for integration
case 'medium',
  nsmax=4;              % max number of soft data
  maxpts=100000;        % number of function eval for integration
case 'fine',
  nsmax=5;              % max number of soft data
  maxpts=1000000;       % number of function eval for integration
end;
order=NaN;              % do not estimation the mean trend
dmax=[200 1.5 5];       % dmax(1)=200 Km, spatial search radius
                        % dmax(2)=1.5 yr, temporal search radius
                        % dmax(3)=5 Km/yr, space/time metric
rEps=0.05;              % Relative numerical error allowable
nMom=2;                 % Calculate the mean and variance
options=BMEoptions;
options(1)=1;
options(3)=maxpts;
options(4)=rEps;
options(8)=nMom;

%
% Get the S/T mean trend mst and adjust the hard and soft data
%
[val,valname,filetitle]=readGeoEAS('examplePM10CAms.dat');
ms=val;
[val,valname,filetitle]=readGeoEAS('examplePM10CAmt.dat');
mt=val';
mstgrid=stmeaninterp(cMS,tME,ms,mt,cMS,tME);
[cmst,mst]=valstg2stv(mstgrid,cMS,tME);
mst=mst(idxNonEmptyData);
[limia]=probaoffset(softpdftype,nl,limi,-mst);

XkBMEa=nan*zeros([length(xk) length(tMEgrid)]);
XkErr=nan*zeros([length(xk) length(tMEgrid)]);
XkBME=nan*zeros([length(xk) length(tMEgrid)]);

%
% For each selected year, estimate Z using BMEprobaMoments
%
for i=1:length(tMEgrid),
  tk=tMEgrid(i);
  dateyear=num2str(tk);
  %disp(sprintf('Calculating map for year %s',dateyear));
  %
  %  Estimate the Z using BMEprobaMoments
  %
  ck=[xk yk tk+0*xk];
  ch=zeros(0,3);
  zh=zeros(0,1);
  [moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limia,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  XkBMEa(:,i)=moments(:,1);
  XkErr(:,i)=sqrt(moments(:,2));
  %
  % Re-ajust XkBMEa by adding the mean to the estimated values
  %
  XkBMEtk=XkBMEa(:,i);
  XkBMEtk(isnan(XkBMEtk))=0;
  msttk=stmeaninterp(cMS,tME,ms,mt,[xk yk],tk);
  XkBME(:,i)=XkBMEtk+msttk;  
  %
  %  Re-ajust XkErr by setting NaN to max of standard dev error
  %
  XkErrtk=XkErr(:,i);
  XkErrtk(isnan(XkErrtk))=max(XkErrtk);
  if ~isreal(XkErrtk), XkErrtk=real(XkErrtk); end;
  XkErr(:,i)=XkErrtk;
end;

%
% Save the output file and plot the map
%
spaceGrid=[xk yk];
save(filename,'spaceGrid','tMEgrid','XkBME','XkErr','cMS','tME');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotZsBME(filename,plotWhat,contMethod,showDataPoints,showEstPoints,showCounties,...
  showCaBound,confLevel,griddatamethod)

% plotZsBME                 - Create figures with maps of annual PM10 geometric mean (Jan 1, 2001)
%
% Uses values of BME estimates that were already calculated in order to create the figures
% with spatial maps related to the annual PM10 geometric for selected dates. The maps created
% may be either (1) the estimated value, (2) the length of a confidence interv., (3) the 
% quantile coresponding to a given 
%
% SYNTAX :
%
% plotZsBME(idxtME,printDev,contMethod,showDataPoints,showEstPoints,showCounties,...
%   showCaBound,confLevel,griddatamethod);
%
% INPUT :
%
% filename        string   optional input with the name of data file with the estimated values.
%                          The default is 'examplePM10CAestBME97'.
% plotWhat        scalar   optional integer to indicate what to plot
%                          1 to plot estimated val, 
%                          2 to plot length of confidence interval, 
%                          where the confidence probability given by confLevel
%                          3 to plot the quantile corresponding to a given confidence probability,
%                          where the confidence probability given by confLevel
%                          4 to plot length of confidence interval for Z NORMALIZED by Z, 
%                          where the confidence probability given by confLevel
%                          Note that plotWhat=4 requires Ztrans=1
%                          The default is 1.
% contMethod      scalar   optional integer indicating the method to use to plot contours:
%                          1 to plot contour lines, 
%                          2 to plot a map in color (using pcolor)
%                          The default is 2.
% showDataPoints  scalar   optional integer indicating whether to show the monitoring stations:
%                          1 to show the monitoring stations, 0 ow
%                          The default is 1.
% showEstPoints   scalar   optional integer indicating whether to show the estimation points:
%                          1 to show estimation points, 0 ow
%                          The default is 0.
% showCounties    scalar   optional integer indicating whether to show the county boundaries:
%                          1 to show Counties boundaries, 0 ow
%                          The default is 1
% showCaBound     scalar   optional integer indicating whether to show the state boundary:
%                          1 to show California state boundaries, 0 ow
%                          The default is 1.
% confLevel       scalar   optional integer indicating a confidence level:
%                          if plotWhat=2 or 4,  confidence level of confidence interv.
%                          if plotWhat=3  confidence level of the quantile to plot
%                          The default is 0.68 if plotWhat=2 or 4, and 0.80 if plotWhat=3.
% griddatamethod string    method to smooth the contour lines. See help griddata
%                          The default is 'linear'.

%
% sets the default values
%
if nargin<1 filename='examplePM10CAestBME97'; end; % name of file with estimated data
if nargin<2 plotWhat=1; end;       % plot estimated values
if nargin<3 contMethod=1; end;     % plot a color map
if nargin<4 showDataPoints=1; end; % show the monitoring stations
if nargin<5 showEstPoints=0; end;  % show the estimation points
if nargin<6, showCounties=1; end;   % show the counties
if nargin<7 showCaBound=1; end;    % show boundary of california
if nargin<8 
  if plotWhat==2 | plotWhat==4
    confLevel=0.68;     % confidence probability of confidence interval
  elseif plotWhat==3
    confLevel=0.84;     % confidence level of quantile
  end
end
if nargin<9 griddatamethod='linear'; end;       % linear interpolation of contours

if isempty(filename), 
  filename='examplePM10CAestBME97';
end;

showLegend=1;
clabelFontSize=8;
Ztrans=1;           %1 to plot variable Z, 2 to plot log(Z)

if plotWhat==4 & Ztrans==2
  error('plotWhat=4 requires Ztrans=1');
end


%
% load the data
%
%addpath ./ZsBMEdir
load(filename);
[val]=readGeoEAS('examplePM10CAstateBorder.dat');
CAborderx=val(:,1);
CAbordery=val(:,2);
if showCounties
  pVal=readAtlasGIS('examplePM10CAcountyBorders.bna');
  for i=1:length(pVal)
    CAcountyBordersx{i}=pVal{i}(:,1);
    CAcountyBordersy{i}=pVal{i}(:,2);
  end
end

xk=spaceGrid(:,1);
yk=spaceGrid(:,2);

xmin=min([cMS(:,1);min(xk)]);
xmax=max([cMS(:,1);max(xk)]);
ymin=min([cMS(:,2);min(yk)]);
ymax=max([cMS(:,2);max(yk)]);

%
% For each selected day, plot map of estimated Z
%
for i=1:length(tMEgrid),
  tk=tMEgrid(i);
  dateyear=num2str(tk);
  %disp(sprintf('Plotting map for year %s',dateyear));
  %
  % Set new figure
  %
  figure;
  hp=[];
  hold on;
  %
  % Show Monitoring Stations
  %
  if showDataPoints==1
    hp(length(hp)+1)=plot(cMS(:,1),cMS(:,2),'hk');
    set(hp(end),'MarkerFaceColor','k');
    set(hp(end),'MarkerSize',5);
    leg{length(hp)}='Monitoring Sites';  
  end;
  %
  % Show estimation points
  %
  if showEstPoints
    hp(length(hp)+1)=plot(xk,yk,'+b');
    leg{length(hp)}='Estimation points';
  end;
  %
  % Construct the grid used to interpolate the estimation values
  %
  dx1=10;
  dy1=dx1;
  xg=[xmin-dx1:dx1:xmax+1.9*dx1];
  yg=[ymin-dy1:dy1:ymax+1.9*dy1];
  [xkm ykm]=meshgrid(xg,yg);
  %
  % Plot contour lines of estimated values if requested
  %
  Zk=XkBME(:,i);
  Zksd=XkErr(:,i);
  if plotWhat==1   % plot the estimated values
    if Ztrans==1;
      Zk=exp(Zk);
    end;
    Zkm=griddata(xk,yk,Zk,xkm,ykm,griddatamethod);
    Zkm=reshape(Zkm,size(xkm));
    maxZkm=max(max(Zkm));
    %contval=quantile(Zkm,[0.1:0.1:0.9 0.95 0.99]);
    %contval=10.^(round(log10(contval)-2)).*round(contval./10.^(round(log10(contval)-2)));
    caxis([0 50]);
    if contMethod==1
      contval=[0:5:maxZkm];
      [hc hl]=contour(xkm,ykm,Zkm,[0:5:25],'-.k');
      clabel(hc,hl,'FontSize',clabelFontSize);
      [hc hl]=contour(xkm,ykm,Zkm,[30:5:maxZkm],'-k');
      clabel(hc,hl,'FontSize',clabelFontSize);
    elseif contMethod==2
      pcolor(xkm,ykm,Zkm);
      shading interp;
      colorbar;
    end
    title(['BME map of the annual geometric mean of PM10 (\mug/m^3) in California for ' dateyear]);
  elseif plotWhat==2 % plot length of estimation interval
    if exist('norminv')==2
      confCoef=gaussinv((1+confLevel)/2,[0 1]);
    else
      F=(1+confLevel)/2;
      confCoef=sqrt(2)*erfinv(2*F-1);
    end;
    confLowBound=Zk-confCoef*Zksd;
    confUpBound=Zk+confCoef*Zksd;
    if Ztrans==1,
      confLowBound=exp(confLowBound);
      confUpBound=exp(confUpBound);
    end;
    Zkm=griddata(xk,yk,confUpBound-confLowBound,xkm,ykm,griddatamethod);
    Zkm=reshape(Zkm,size(xkm));
    maxZkm=max(max(Zkm));
    caxis([0 12]);
    if contMethod==1
      [hc hl]=contour(xkm,ykm,Zkm,'-k');
      clabel(hc,hl,'FontSize',clabelFontSize);
    elseif contMethod==2
      pcolor(xkm,ykm,Zkm);
      shading interp;
      colorbar;
    end;
    title(['Length of ' num2str(confLevel*100)  ...
        '% CI for PM10 yearly geometric mean (\mug/m^3) in ' dateyear]);
  elseif plotWhat==3 % plot quantile
    if exist('norminv')==2
      confCoef=gaussinv(confLevel,[0 1]);
    else
      F=confLevel;
      confCoef=sqrt(2)*erfinv(2*F-1);
    end;
    quantileval=Zk+confCoef*Zksd;
    if Ztrans==1
      quantileval=exp(quantileval);
    end;
    Zkm=griddata(xk,yk,quantileval,xkm,ykm,griddatamethod);
    Zkm=reshape(Zkm,size(xkm));
    maxZkm=max(max(Zkm));
    caxis([0 50]);
    if contMethod==1
      [hc hl]=contour(xkm,ykm,Zkm,'-k');
      clabel(hc,hl,'FontSize',clabelFontSize);
    elseif contMethod==2
      pcolor(xkm,ykm,Zkm);
      shading interp;
      colorbar;
    end
    title([sprintf('%.0f%% quantile for PM10 yearly geometric mean on %s, max=%.1f ',...
      confLevel*100,dateyear,maxZkm) '\mug/m^3'] );
  elseif plotWhat==4 % plot length of estimation interv normalized by est val
    if exist('norminv')==2
      confCoef=gaussinv((1+confLevel)/2,[0 1]);
    else
      F=(1+confLevel)/2;
      confCoef=sqrt(2)*erfinv(2*F-1);
    end;
    confLowBound=Zk-confCoef*Zksd;
    confUpBound=Zk+confCoef*Zksd;
    if Ztrans==1,
      confLowBound=exp(confLowBound);
      confUpBound=exp(confUpBound);
      Zk=exp(Zk);
    end;
    Zkm=griddata(xk,yk,(confUpBound-confLowBound)./Zk,xkm,ykm,griddatamethod);
    Zkm=reshape(Zkm,size(xkm));
    maxZkm=max(max(Zkm));
    caxis([0 0.3]);
    if contMethod==1
      [hc hl]=contour(xkm,ykm,Zkm,[.1:.1:.3],'-k');
      clabel(hc,hl,'FontSize',clabelFontSize);
    elseif contMethod==2
      pcolor(xkm,ykm,Zkm);
      shading interp;
      colorbar;
    end;
    title(['Normalized length of ' num2str(confLevel*100)  ...
        '% confidence interval for PM10 yearly geometric mean in ' dateyear]);
  end;  
  %
  % Show California boundary
  %
  if showCounties
    switch contMethod 
    case 1,  ltype=':k';
    case 2,  ltype='-b';
    end
    for i=1:length(CAcountyBordersx)
      htemp=plot(CAcountyBordersx{i},CAcountyBordersy{i},ltype);
    end;
    hp(length(hp)+1)=htemp;
    leg{length(hp)}='County boundaries';
  end;
  if showCaBound
    axis([1971 3157 1482 2577]);
    ax=axis;
    CApatchx=[CAborderx(1);ax(2);ax(2);ax(1);ax(1);ax(2);ax(2);CAborderx];
    CApatchy=[CAbordery(1);CAbordery(1);ax(4);ax(4);ax(3);ax(3);CAbordery(1);CAbordery];
    switch contMethod 
    case 1,
      h=fill(CApatchx,CApatchy,'w','EdgeColor','w');
      plot(CAborderx,CAbordery,'k');
    case 2,
      h=fill(CApatchx,CApatchy,'b','EdgeColor','b');
      plot(CAborderx,CAbordery,'b');
    end
  end;
  if showLegend & length(hp)>0
    legend(hp,leg);
  end
end;

xlabel('s_1 (Km)');
ylabel('s_2 (Km)');

