% BMEPROBALIBtutorial       - tutorial for the bmeprobalib directory (Jan 1,2001)

% Generates the data files needed for this tutorial 
% -------------------------------------------------
echo off;
if exist('hardS.dat')~=2 | exist('hardST.dat')~=2,
  
  %%% generate hard and soft proba data for a 2D-spatial domain
  
  rand('state',2);
  randn('state',2);
  
  [xk yk]=meshgrid([0:2:8],[0:2:10]);
  ck=[xk(:),yk(:)];
  nk=size(ck,1);
  
  nh=10;
  ns=11;
  idxperm=randperm(nk);
  ch=ck(idxperm(1:nh),:);
  cs=ck(idxperm(nh+1:nh+ns),:);
  
  cc=1;
  aa=5;
  covmodel='exponentialC';
  covparam=[cc aa];
  
  [zsim,L]=simuchol([ch;cs],covmodel,covparam,1);
  zsim=zsim+3;
  
  zh=zsim(1:nh);zh=zh(:);
  writeGeoEAS([ch zh],{'x','y','Variable 1'},'File with hard data','hardS.dat');
  
  zs=zsim(nh+1:nh+ns);zs=zs(:);
  softpdftype=2;
  NV=3*ones(ns,1);
  WidthV=[0.3*rand(ns,1) 0.7*rand(ns,1) 3.0*rand(ns,1)];
  A=0.5*sum(WidthV,2)+0.5*WidthV(:,2);
  probdensV=[zeros(ns,1) 1./A 1./A zeros(ns,1)];
  [nl,limi,probdens]=simuprobabilistic(zs,softpdftype,NV,WidthV,probdensV);
  writeProba({cs,ones(ns,1)},0,softpdftype,nl,limi,probdens,'File with proba data','probaS.dat');
  
  %%% generate hard and soft proba data for a space/time domain
  
  cc=1;
  aa=5;
  at=10;
  covmodel='exponentialC/exponentialC';
  covparam=[cc aa at];
  
  tME=[1:5];
  nME=length(tME);
  [ch,z]=valstg2stv(zeros(nh,nME),ch,tME);
  nh=size(ch,1);
  [cs,z]=valstg2stv(zeros(ns,nME),cs,tME);
  ns=size(cs,1);
  
  [zsim,L]=simuchol([ch;cs],covmodel,covparam,1);
  zsim=zsim+3;
  
  zh=zsim(1:nh);zh=zh(:);
  writeGeoEAS([ch zh],{'x','y','t','Variable 1'},'File with hard data','hardST.dat');
  
  zs=zsim(nh+1:nh+ns);zs=zs(:);
  softpdftype=2;
  NV=3*ones(ns,1);
  WidthV=[0.3*rand(ns,1) 0.7*rand(ns,1) 3.0*rand(ns,1)];
  A=0.5*sum(WidthV,2)+0.5*WidthV(:,2);
  probdensV=[zeros(ns,1) 1./A 1./A zeros(ns,1)];
  [nl,limi,probdens]=simuprobabilistic(zs,softpdftype,NV,WidthV,probdensV);
  writeProba({cs,ones(ns,1)},1,softpdftype,nl,limi,probdens,'File with proba data','probaST.dat');
  
end;

% Beginning of tutorial 
% -------------------------------------------------

%%% Clear memory content and set echo to on.

clear;
clc;
echo on;

%%% BME MAPPING USING PROBABILISTIC DATA ON A 2-DIMENTIONAL SPATIAL DOMAIN
%%% 

%%% Load the hard data and the soft data from files 

[val,valname,filetitle]=readGeoEAS('hardS.dat');
ch=val(:,1:2);
zh=val(:,3);

[cs,isST,softpdftype,nl,limi,probdens,filetitle]=readProba('probaS.dat');
cs=cs{1};

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the soft
%%% pdf for some selected soft probabilistic data points.

figure;hold on;
subplot(4,1,1);
probaplot(softpdftype,nl,limi,probdens,'-',1);
title('Soft pdf at 4 selected probabilistic data points');
ylabel('f_S(z_s)');
subplot(4,1,2);
probaplot(softpdftype,nl,limi,probdens,'-',2);
ylabel('f_S(z_s)');
subplot(4,1,3);
probaplot(softpdftype,nl,limi,probdens,'-',3);
ylabel('f_S(z_s)');
subplot(4,1,4);
probaplot(softpdftype,nl,limi,probdens,'-',4);
ylabel('f_S(z_s)');
xlabel('z_s');

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the values for
%%% the hard data using colorplot.m, with the 'hot' color
%%% map. Hard data are represented by squares, whereas soft
%%% data are represented by circles having a size which
%%% is proportional to the expected value of the soft data.

figure;hold on;
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',12,[0 0 0]};
colorplot(ch,zh,'hot',Property,Value);

Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value ={'o','k','none'};
[ms,vs]=proba2stat(softpdftype,nl,limi,probdens);
markerplot(cs,ms,[5 20],Property,Value);

title('Hard data (squares) and probabilistic soft data (circles)');
axis([0 9.4 0 10]);
xlabel('x');
ylabel('y');

%%% Type any key for continuing...
pause;
clc;

%%% Compute the BME mean estimate of the natural 
%%% variable at the node of a two dimensional grid 
%%% with regular spacing. To do the computation
%%% first set the coordinates of the estimation grid,
%%% then set the parameters for the BME estimation,
%%% and finally calculate the BME mean estimate and the 
%%% estimation error variance at all estimation points

[xk yk]=meshgrid([0:1:8],[0:1:10]);
ck=[xk(:),yk(:)];


cc=1;                       % covariance sill (variance)
aa=5;                       % covariance range
covmodel='exponentialC';    % covariance model
covparam=[cc aa];           % parameter for the covariance model
nhmax=10;                   % max number of hard data in estimation neighborhood
nsmax=2;                    % max number of soft data in estimation neighborhood
dmax=100;                   % max spatial search radius for estimation neighborhood
order=0;                    % The mean is supposed to be constant over est. neighb.
options=BMEoptions;         % Use default options


[moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
zk=moments(:,1);
vk=moments(:,2);

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the 
%%% the BME mean estimate of the natural variable
%%% using pcolor.m with a 'hot' color map.
%%% This map represents the estimate which minimizes
%%% the estimation error variance of the BME
%%% posterior pdf

figure;
[xg yg]=meshgrid([0:.5:8],[0:.5:10]);
[zg]=griddata(ck(:,1),ck(:,2),zk,xg,yg);
pcolor(xg,yg,zg);
colormap(hot);
colorbar;
shading interp;
title('Map of the BME mean estimate of the natural variable');
xlabel('x');
ylabel('y');

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the 
%%% the BME estimation error variance (i.e. 
%%% the variance of the BME posterior pdf).
%%% This map indicates that uncertainty with 
%%% the spatial estimate is highest away from
%%% the data points, it is lower at soft data
%%% points, and it is zero at the hard data points.

figure;
[xg yg]=meshgrid([0:0.5:8],[0:0.5:10]);
[vg]=griddata(ck(:,1),ck(:,2),vk,xg,yg,'v4');
pcolor(xg,yg,vg);
colormap(hot);
colorbar;
shading interp;
title('Map of the estimation error variance');
xlabel('x');
ylabel('y');

%%% Type any key to continuing...
pause;
clc;

%%% Compute the BME posterior pdf at a given
%%% estimation point, and open a new graphic
%%% window to display this pdf

ck=[5 5];                   % set the coordinate of the estimation point
[z,pdf,info]=BMEprobaPdf([],ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
figure;hold on;
plot(z,pdf);
title('BME posterior pdf of variable z_k at coordinate (x,y)=(5,5)');
xlabel('z_k');
ylabel('f_K(z_k)');

%%% Type any key to continuing...
pause;
clc;

%%% Calculate the BME confidence interval (CI) corresponding
%%% to the 68%, 90% and 95% confidence probability,
%%% and plot these interval on the graph of the BME posterior
%%% pdf. Note that the BMEprobaCI function also returns 
%%% the BME posterior pdf

options(20)=0.68;
options(21)=0.90;
options(22)=0.99;
[zlCI,zuCI,pdfCI,PCI,z,pdf]=BMEprobaCI(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);

plot([zlCI(1) zuCI(1)],[pdfCI(1) pdfCI(1)],'r');
text(0.5*(zlCI(1)+zuCI(1)),pdfCI(1)+0.01,[num2str(PCI(1)) '% CI'],'HorizontalAlignment','center','Color','r');

plot([zlCI(2) zuCI(2)],[pdfCI(2) pdfCI(2)],'r');
text(0.5*(zlCI(2)+zuCI(2)),pdfCI(2)+0.01,[num2str(PCI(2)) '% CI'],'HorizontalAlignment','center','Color','r');

plot([zlCI(3) zuCI(3)],[pdfCI(3) pdfCI(3)],'r');
text(0.5*(zlCI(3)+zuCI(3)),pdfCI(3)+0.01,[num2str(PCI(3)) '% CI'],'HorizontalAlignment','center','Color','r');

%%% Type any key to continuing...
pause;
clc;

%%% End of first part of the tutorial using a 2-D spatial domain 
%%%      Now continuing with tutorial using a space/time  domain (2D in space)...

%%% Type any key to continuing...
pause;
clc;

%%% BME MAPPING USING PROBABILISTIC DATA ON A SPACE/TIME DOMAIN, 2D in space
%%% 

%%% Load the hard data and the soft data from files 

[val,valname,filetitle]=readGeoEAS('hardST.dat');
ch=val(:,1:3);
zh=val(:,4);

[cs,isST,softpdftype,nl,limi,probdens,filetitle]=readProba('probaST.dat');
cs=cs{1};

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the values for
%%% the hard data and soft data for Day 2

th=ch(:,3);               
ts=cs(:,3);

figure;hold on;
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',12,[0 0 0]};
colorplot(ch((th==2),1:2),zh(th==2),'hot',Property,Value);

Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value ={'o','k','none'};
[ms,vs]=proba2stat(softpdftype,nl,limi,probdens);
markerplot(cs((ts==2),1:2),ms(ts==2),[5 20],Property,Value);

title('Hard (squares) and soft (circles) data for Day 2');
axis([0 9.4 0 10]);
xlabel('x');
ylabel('y');

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the values for
%%% the hard data and soft data for Day 3

th=ch(:,3);               
ts=cs(:,3);

figure;hold on;
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',12,[0 0 0]};
colorplot(ch((th==3),1:2),zh(th==3),'hot',Property,Value);

Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value ={'o','k','none'};
[ms,vs]=proba2stat(softpdftype,nl,limi,probdens);
markerplot(cs((ts==3),1:2),ms(ts==3),[5 20],Property,Value);

title('Hard (squares) and soft (circles) data for Day 3');
axis([0 9.4 0 10]);
xlabel('x');
ylabel('y');

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the values for
%%% the hard data and soft data for Day 4

th=ch(:,3);               
ts=cs(:,3);

figure;hold on;
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',12,[0 0 0]};
colorplot(ch((th==4),1:2),zh(th==4),'hot',Property,Value);

Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value ={'o','k','none'};
[ms,vs]=proba2stat(softpdftype,nl,limi,probdens);
markerplot(cs((ts==4),1:2),ms(ts==4),[5 20],Property,Value);

title('Hard (squares) and soft (circles) data for Day 4');
axis([0 9.4 0 10]);
xlabel('x');
ylabel('y');

%%% Type any key to continuing...
pause;
clc;

%%% Compute the BME mean estimate of the natural 
%%% variable at the node of a two dimensional grid 
%%% with regular spacing for day 3.
%%% To do the computation,
%%% first set the coordinates of the estimation grid,
%%% then set the parameters for the BME estimation,
%%% and finally calculate the BME mean estimate and the 
%%% estimation error variance at all estimation points

[xk yk]=meshgrid([0:1:8],[0:1:10]);
nk=prod(size(xk));
ck=[xk(:) yk(:) 3*ones(nk,1)];

cc=1;                       % covariance sill (variance)
aa=5;                       % covariance spatial range
at=10;                      % covariance temporal range
covmodel='exponentialC/exponentialC';    % covariance model
covparam=[cc aa at];        % parameter for the covariance model
nhmax=10;                   % max number of hard data in estimation neighborhood
nsmax=2;                    % max number of soft data in estimation neighborhood
dmax=[100 10 aa/at];        % dmax(1)=max spatial search radius for estimation neighborhood
                            % dmax(2)=max temporal search radius for estimation neighborhood
                            % dmax(3)=space/time metric, usually equal to aa/at
order=0;                    % The mean is supposed to be constant over the est. neighb.
options=BMEoptions;         % Use default options
options(3)=100000;
options(4)=0.002;

[moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
zk=moments(:,1);
vk=moments(:,2);

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the 
%%% the BME mean estimate of the natural variable for Day 3,
%%% using pcolor.m with a 'hot' color map.
%%% This map represents the estimate which minimizes
%%% the estimation error variance of the BME
%%% posterior pdf

figure;
[xg yg]=meshgrid([0:.5:8],[0:.5:10]);
[zg]=griddata(ck(:,1),ck(:,2),zk,xg,yg);
pcolor(xg,yg,zg);
colormap(hot);
colorbar;
shading interp;
title('Map of the BME mean estimate of the natural variable for Day 3');
xlabel('x');
ylabel('y');

%%% Type any key to continuing...
pause;
clc;

%%% Open a new graphic window and display the 
%%% the BME estimation error variance (i.e. 
%%% the variance of the BME posterior pdf), for Day 3.

figure;
[xg yg]=meshgrid([0:0.5:8],[0:0.5:10]);
[vg]=griddata(ck(:,1),ck(:,2),vk,xg,yg,'v4');
pcolor(xg,yg,vg);
colormap(hot);
colorbar;
shading interp;
title('Map of the estimation error variance for Day 3');
xlabel('x');
ylabel('y');

%%% end of tutorial
echo off;

