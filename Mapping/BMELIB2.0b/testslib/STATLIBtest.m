function STATLIBtest

% STATLIBtest             - Test for the statlib directory
%
% SYNTAX : 
%
% STATLIBtest

rand('state',0);randn('state',0);

%%% generate data for 3 fields
nh=100;
ch=10*rand(nh,2);
cc=1;
aa=5;
covmodel='exponentialC';
covparam=[cc aa];
[Z,L]=simuchol(ch,covmodel,covparam,3);
Z1=Z(:,1);
Z2=0.3*Z(:,1)+0.7*Z(:,2);
Z3=0.3*Z(:,1)+0.7*Z(:,3);
Z=[Z1 Z2 Z3];


%%% test triangularinv
disp(' ');
disp('Testing the triangularinv function');
[z]=triangularinv([.1 .3],[1 4 2]);

%%% test histscaled and histline
disp(' ');
disp('Testing the histscaled and histline function');
figure;
subplot(2,1,1);
histscaled(Z1,15,[-2 2],ones(size(Z1)));
title('histogram scaled');
subplot(2,1,2);
[hdle,zline,fline,n,x]=histline(Z1,15,[-2 2],ones(size(Z1)));
title('histogram line');
disp('Check the histscaled and histline figure to make sure it is ok');

%%% test kerneldensity, cdfest and pdf2cdf
disp(' ');
disp('Testing the kerneldensity, cdfest and pdf2cdf function');
zfileK=[-2:.02:2]';
pdfzfileK=kerneldensity(Z1,zfileK,.01);
figure;
subplot(2,1,1);
plot(zfileK,pdfzfileK,'r');
title('pdf');

subplot(2,1,2);
[zfileE,cdfzfileE]=cdfest(Z1);
[cdfzfileK]=pdf2cdf(zfileK,pdfzfileK);
plot(zfileE,cdfzfileE);
hold on;
plot(zfileK,cdfzfileK,'r');
title('cdf');
disp('Check the kerneldensity, cdfest and pdf2cdf figure to make sure it is ok');

%%% test vario, modelfit and modelplot
disp(' ');
disp('Testing the vario, modelfit and modelplot function');
cl=(0:1:5)';
figure;
optionsNS=[0,45,-45];
[dZ1,vZ1,oZ1]=vario(ch,Z1,cl,'kron',optionsNS);
plot(dZ1,vZ1,'*');
hold on;
modelZ1={'nuggetV','exponentialV'};
param0={[.1],[.9 5]};
modelplot(dZ1,modelZ1,param0);
paramZ1=modelfit(dZ1,vZ1,oZ1,modelZ1,param0);
modelplot(dZ1,modelZ1,paramZ1,'Color',[1 0 0]);
disp('Check the vario, modelfit and modelplot figure to make sure it is ok');

%%% test histscatterplot
disp(' ');
disp('Testing the histscatterplot function');
figure;
histscatterplot([Z]);
disp('Check the histscatterplot figure to make sure it is ok');

%%% test vario and coregfit
disp(' ');
disp('Testing the vario and coregfit function');
options=1;
figure;
[dT,VT,oT]=vario(ch,[Z],cl,'kron',options);
modelT={'nuggetV','exponentialV'};
paramT0={[.1],[.9 5]};
options=1;
[paramT]=coregfit(dT,VT,oT,modelT,paramT0,options);
disp('Check the vario and coregfit figure to make sure it is ok');

%%% test stmean and stcov functions
stcovtest;

%%% test complete
disp(' ');
disp('test complete');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stcovtest;

% stcovtest                 - Test for the stcov function (Jan 1, 2001) 
%
% Test the stcov function which calculates estimates of the s/t
% covariance of a Random Field Z using measurements of Z at 
% Measuring Sites with coordiantes cMS, and at Measuring Events 
% with time tME.
% This function calls first teststmean to generate values of Z,
% and calculate the s/t mean trend
%
% SYNTAX :
%
% stcovtest;

[cMS,tME,covmodel,covparam,Z,ms,mt,mss,mts]=stmeantest;

disp(' ');
disp('Testing the space/time covariance estimation function stcov');

Zmst=stmeaninterp(cMS,tME,mss,mts,cMS,tME);

Za=Z-Zmst;

rLag=   [0 2 4 6 8]';
rLagTol=[0 2 2 2 2]';
tLag=   [0 2 4 6 8 10];
tLagTol=[0 2 2 2 2 2];  
[C np]=stcov(Za,cMS,tME,Za,cMS,tME,rLag,rLagTol,tLag,tLagTol);

figure;
subplot(2,1,1);
h1=plot(rLag,C(:,1),'o');
hold on;
r=[0:0.5:8]';
Crtrue=coord2K([0 0 0],[r 0*r 0*r],covmodel,covparam);
h2=plot(r,Crtrue);
legend([h1 h2],'Calculated Cov','True Cov');
xlabel('r');
ylabel('C(r,t=0)');
title('Covariance C(r,t=0)');
ax=axis;
plot([ax(1) ax(2)],[0 0],'-k');

subplot(2,1,2);
h1=plot(tLag,C(1,:),'o');
hold on;
t=[0:0.5:10];
Cttrue=coord2K([0 0 0],[0*t' 0*t' t'],covmodel,covparam);
h2=plot(t,Cttrue);
legend([h1 h2],'Calculated Cov','True Cov');
xlabel('t');
ylabel('C(r=0,t)');
title('Covariance C(r=0,t)');
plot([ax(1) ax(2)],[0 0],'-k');

disp('Check the Covariance figure to verify that they seems correct');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cMS,tME,covmodel,covparam,Z,ms,mt,mss,mts]=stmeantest;

% stmeantest                - Test for the stmean function (Jan 1, 2001) 
%
% SYNTAX :
%
% [cMS,tME,covmodel,covparam,Z,ms,mt,mss,mts]=stmeantest;
%
% OUTPUT :
% 
% cMS       nMS by 2   matrix of spatial coordinates for the Measuring Sites
% tME       1 by nME   vector of times for the Measuring Event
% covmodel  string     name of a S/T covariance model 
% covparam  vector     parameters for covmodel
% Z         nMS by nME matrix of simulated values for Z at the nMS Measuring 
%                      Sites and nME measuring events. May contain some NaN
% ms        nMS by 1   vector of the average of Z at each Measuring Sites
% mss       nMS by 1   vector of smoothed ms
% mt        1 by nME   vector of the average of Z at each Measuring Event
% mts       nMS by 1   vector of smoothed mt

disp(' ');
disp('Testing the mean trend function stmean');

rand('state',1);
randn('state',1);

x0=0;
Lx=10;
smx=0.1;
y0=0;
Ly=15;
smy=0.05;
nMS=15;
t0=0;
t1=20;
smt=0.2;
nME=21;
covmodel='exponentialC/exponentialC';
covparam=[1 3 4];

cMS=[x0+Lx*rand(nMS,1),y0+Ly*rand(nMS,1)];
tME=t0:(t1-t0)/(nME-1):t1;
p=[kron(ones(nME,1),cMS) kron(tME',ones(nMS,1))];
[Z,L]=simuchol(p,covmodel,covparam);

Z=Z+smx*p(:,1)+smy*p(:,2)+smt*p(:,3);

Z=reshape(Z,nMS,nME);

randidx=randperm(nMS*nME);
nNaN=20;
Z(randidx(1:nNaN))=NaN;

stmeanparam=[15 4 20 3 0];
[ms,mss,mt,mts]=stmean(Z,cMS,(1:nMS)',tME,stmeanparam);


figure;
nsubplot=2;
for iMS=1:nsubplot
  subplot(nsubplot,1,iMS);
  hold on;
  ZiMS=Z(iMS,:);
  h1=plot(tME(~isnan(ZiMS)),ZiMS(~isnan(ZiMS)),'o');
  Zcalc=stmeaninterp(cMS,tME,mss,mts,cMS(iMS,:),tME);
  h2=plot(tME,Zcalc);
  Ztrue=smx*cMS(iMS,1)+smy*cMS(iMS,2)+smt*tME;
  h3=plot(tME,Ztrue,'r-');
  legend([h1 h2 h3],'data','calc mean trend','true mean trend',4);
  xlabel('t');ylabel('Z');title(sprintf('Z at Measuring Site %d',iMS));
end;


disp('Testing of stmean function complete. Check the figure of Z at different MS');

