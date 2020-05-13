function exampleBMEvsSKME(nSim,useStoredresults)

% exampleBMEvsSKME          - Comparison between BMEproba and SKME (Jan 1,2001)
% 
% In this example with probabilistic soft data the BME method 
% is compared with the Simple Kriging with Measurement Error method. The 
% steps in this example are the following:
% 1) Set the mapping situation by chosing a covariance model, and the location
%    of 2 hard data points, 2 soft data points, and 1 estimation points.
% 2) Simulate nSim realizations of the hard data, soft data, and the value ZkSim
%    of the variable at the estimatoin point
% 3) Using only the hard and soft data, estimate the varaible at the estimation
%    point for each realiation using BME and different implementation of SKME
% 4) Calculate the estimation errors as the difference for each realization 
%    between estimated and simulated (true) values.
% 5) Calculate the Mean Square Error MSE and plot the histogram of estimation 
%    errors. The better method is that with smaller MSE and greater probability
%    of producing an error of 0.
%
% SYNTAX :
% 
% exampleBMEvsSKME(nSim,useStoredresults);
%
% INPUT :
%
% nSim               scalar     optional integer indicating the number of
%                               simulations to use.
%                               default value is 100. Use higher number for better result
% useStoredresults   scalar     optional integer indicating wether to use already
%                               calculated results that were automatically stored
%                               1 tu use stored results, 0 otherwise.
%                               default value is 0.
%
% NOTE :
%
% You must first use useStoredresults=0 to calculate and store the results, 
% and then only you can use useStoredresults=1 to use already calculated results.
% This saves time the second time around so you can see the results without having 
% to recalculate them.
% If the stored results cannot be found, they are automatically (re)calculated 


if nargin<1, nSim=100; end;
if nargin<2, useStoredresults=0; end;

if useStoredresults==1 & exist('exampleBMEvsSKME.mat')~=2 
  disp('Cannot find stored results, therefore calculating them');
  useStoredresults=0;
end;

if useStoredresults==0
  
  disp(sprintf('Calculating results using %d simulations',nSim));
  
  %
  %initialization
  %
  rand('state',2);
  randn('state',13);randn(15,1);
  
  %
  % Set the coordinates of the estimated, hard and soft data points
  %
  nk=1;nh=2;ns=2;
  ck=[.4 .5];
  ch=[.9 .1;.1 .9];
  cs=[.31 .44;.54 .41];
  figure;hold on;
  hk=plot(ck(:,1),ck(:,2),'x');
  hh=plot(ch(:,1),ch(:,2),'v');
  hs=plot(cs(:,1),cs(:,2),'o');
  axis([0 1 0 1]);drawnow;
  legend([hk hh hs],'Estimation pts','Hard data points','Soft data points');
  xlabel('x');
  ylabel('y');
  title('Map of estimation and data points');
  drawnow;
  
  %
  % Sets soft probabilistic data template
  %
  incWidth=.03;
  sig=0.1;
  pickRatio=30;
  x=0:incWidth:2;
  y=exp(-(x-0.2).^2/sig^2);
  y=y+max(y)/pickRatio;
  y=y.*x.^(1/3).*(2-x).^(1/3);
  normConst=diff(x)*(y(1:length(y)-1)+0.5*diff(y))';
  y=y/normConst;
  figure;plot(x,y);title('template of soft probabilistic data');drawnow;
  softpdftype=4;
  NV=(length(x)-1)*ones(ns,1);
  WidthV=kron([incWidth],ones(ns,1));
  ProbV=kron(y,ones(ns,1));
  
  %
  % Set the covariance model and BME parameters
  %
  cNugget=0.0;
  cc=1;aa=sqrt(3);
  covmodel={'nuggetC','gaussianC'};
  covparam={cNugget,[cc aa]};
  nhmax=10; 
  nsmax=4;
  dmax=100;
  order=NaN;
  maxpts=100000;
  rEps=0.02;
  options(1)=0;
  options(3)=maxpts;
  options(4)=rEps;
  options(8)=2;
  
  %
  %  Generating simulations
  %
  [zall,L]=simuchol([ck;ch;cs],covmodel,covparam,nSim);
  zkSimAll=zall(1:nk,:);
  zhAll=zall(nk+1:nk+nh,:);
  zsAll=zall(nk+nh+1:nk+nh+ns,:);
  [Nl,LimiAll,ProbAll]=simuprobabilistic(zsAll,softpdftype,NV,WidthV,ProbV);
  
  %
  %doing the BME estimation with BMEprobaMode
  %
  disp(sprintf('Doing BME estimation with BMEprobaMode'));
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    nl=Nl;
    limi=LimiAll(:,:,iSim);
    probdens=ProbAll(:,:,iSim);
    [zk,info]=BMEprobaMode(ck,ch,cs,zh,softpdftype,...
      nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkBMEmode(1:nk,iSim)=zk;
  end;  
  tsBME=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsBME));
  
  %
  %doing the BME estimation with BMEprobaMoments
  %
  disp(sprintf('Doing BME estimation with BMEprobaMoments'));
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    nl=Nl;
    limi=LimiAll(:,:,iSim);
    probdens=ProbAll(:,:,iSim);
    [moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,...
      nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkBMEmean(1:nk,iSim)=moments(:,1);
    zkBMEvar(1:nk,iSim)=moments(:,2);
  end;  
  tsBME=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsBME));
  
  %
  %doing the krigingME estimation with proba2stat
  %
  disp('Doing krigingME estimation with proba2stat');
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    nl=Nl;
    limi=LimiAll(:,:,iSim);
    probdens=ProbAll(:,:,iSim);
    [zs,vs]=proba2stat(softpdftype,nl,limi,probdens);
    [zk,vk]=krigingME(ck,ch,cs,zh,zs,vs,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkSKMEpdf(1:nk,iSim)=zk;
    zkSKMEpdfVar(1:nk,iSim)=vk;
  end
  tsSKMEpdf=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsSKMEpdf));
  
  %
  %doing the krigingME estimation with proba2interval
  %
  disp('Doing krigingME estimation with proba2interval');
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    nl=Nl;
    limi=LimiAll(:,:,iSim);
    probdens=ProbAll(:,:,iSim);
    [a,b]=proba2interval(softpdftype,nl,limi);
    zs=0.5*(a+b);
    vs=(1/12)*(a-b).^2;
    [zk,vk]=krigingME(ck,ch,cs,zh,zs,vs,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkSKME(1:nk,iSim)=zk;
    zkSKMEvar(1:nk,iSim)=vk;
  end
  tsSKME=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsSKME));
  
  %
  %doing the SK estimation
  %
  disp('Doing SK estimation');
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    nl=Nl;
    limi=LimiAll(:,:,iSim);
    probdens=ProbAll(:,:,iSim);
    [a,b]=proba2interval(softpdftype,nl,limi);
    zs=0.5*(a+b);
    [zk,vk]=kriging(ck,[ch;cs],[zh;zs],covmodel,covparam,nhmax,dmax,order,options);
    zkSK(1:nk,iSim)=zk;
    zkSKvar(1:nk,1)=vk;
  end
  tsSK=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsSK));
  
  
  %
  %doing the SKh estimation
  %
  disp('Doing SKh estimation');
  t0=clock;
  if nh>0,
    for iSim=1:nSim,  
      zh=zhAll(:,iSim);
      nl=Nl;
      [zk,vk]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order,options);
      zkSKh(1:nk,iSim)=zk;
      zkSKhVar(1:nk,1)=vk;
    end
  else
    zkSKh=NaN*ones(nk,nSim);
    zkSKhVar=NaN*ones(nk,nSim);
  end;
  tsSKh=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsSKh));
  
  %
  %save the data
  %
  
  save exampleBMEvsSKME zkSimAll zkBMEmode zkBMEmean zkSKMEpdf zkSKME zkSK zkSKh...
    zkBMEvar zkSKMEpdfVar zkSKMEvar zkSKvar zkSKhVar nSim; 
  
elseif useStoredresults==1
  
  load exampleBMEvsSKME;
  disp(sprintf('Using stored results that were calculated using %d simulations',nSim));
  
else
  
  error('useStoredresults must be equal to either 0 or 1');
  
end;
  
%
%Calculating the estimation errors
%
zkBMEmodeErr=zkBMEvar;        
BMEmodeErr=zkBMEmode-zkSimAll;
BMEmeanErr=zkBMEmean-zkSimAll;
SKMEpdfErr=zkSKMEpdf-zkSimAll;
SKMEerr=zkSKME-zkSimAll;
SKerr=zkSK-zkSimAll;
SKhErr=zkSKh-zkSimAll;
disp(' ');
disp('                --------------------------------------------------------');
disp('                  Square root of Mean Square Error (E) for each method:');
disp('                --------------------------------------------------------');
disp('                BMEmode   BMEmean   SKMEpdf   SKME      SK        SKh'); 
disp(sprintf('Actual E     %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  ',...
  sqrt(mean(BMEmodeErr.^2)),sqrt(mean(BMEmeanErr.^2)),sqrt(mean(SKMEpdfErr.^2)),...
  sqrt(mean(SKMEerr.^2)),   sqrt(mean(SKerr.^2)),     sqrt(mean(SKhErr.^2)) ) );
disp(sprintf('Predicted E            %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  ',...
  sqrt(mean(zkBMEvar)),sqrt(zkSKMEpdfVar(1)),sqrt(zkSKMEvar(1)),sqrt(zkSKvar),...
  sqrt(zkSKhVar)));
%disp([mean(zkBMEmodeErr) mean(zkBMEvar) zkSKMEpdfVar(1) zkSKMEvar(1) zkSKvar zkSKhVar]);

%
% plotting the histogram of the Kriging methods
%
figure;
if nSim<20, Nb=11;
elseif nSim<200, Nb=20;
else Nb=50;
end;  
  
zmin=-1.5;
zmax=1.5;
h3=histline(SKMEpdfErr(1,:),Nb,[zmin,zmax]);hold on;
h4=histline(SKMEerr(1,:),Nb,[zmin,zmax]);
h5=histline(SKerr(1,:),Nb,[zmin,zmax]);
h6=histline(SKhErr(1,:),Nb,[zmin,zmax]);
set(h3,'Color','g');
set(h4,'Color','y');
set(h5,'Color','m');
set(h6,'Color','r');
set(h3,'LineStyle','-');
set(h4,'LineStyle','--');
set(h5,'LineStyle','-.');
set(h6,'LineStyle',':');
SKMEpdfLegend =sprintf('SKMEpdf, E=%.3f',sqrt(mean(SKMEpdfErr(1,:).^2)));
SKMElegend =sprintf('SKME, E=%.3f',sqrt(mean(SKMEerr(1,:).^2)));
SKlegend=sprintf('SK , E=%.3f',sqrt(mean(SKerr(1,:).^2)));
SKhLegend=sprintf('SKh , E=%.3f',sqrt(mean(SKhErr(1,:).^2)));
legend([h3 h4 h5 h6],SKMEpdfLegend,SKMElegend,SKlegend,SKhLegend,2);
ax=axis;
axis([zmin zmax ax(3) ax(4)]);
xlabel('Estimation Error');
ylabel('Frequency');

%
% plotting the histogram of the Probabilistic methods
%
figure;
Nb=11;
zmin=-1;
zmax=1;
h1=histline(BMEmodeErr(1,:),Nb,[zmin,zmax]);hold on;
h2=histline(BMEmeanErr(1,:),Nb+1,[zmin,zmax]);hold on;
h3=histline(SKMEpdfErr(1,:),Nb+2,[zmin,zmax]);
h6=histline(SKhErr(1,:),Nb+3,[zmin,zmax]);
set(h1,'Color','b');
set(h2,'Color','c');
set(h3,'Color','g');
set(h6,'Color','r');
set(h1,'LineStyle','-');
set(h2,'LineStyle','--');
set(h3,'LineStyle','-.');
set(h6,'LineStyle',':');
BMEpdf1Legend=sprintf('BMEpdfMode  , E=%.3f',sqrt(mean(BMEmodeErr(1,:).^2)));
BMEpdf2Legend=sprintf('BMEpdfMean  , E=%.3f',sqrt(mean(BMEmeanErr(1,:).^2)));
SKMEpdfLegend =sprintf('SKMEpdf, E=%.3f',sqrt(mean(SKMEpdfErr(1,:).^2)));
SKhLegend=sprintf('SKh , E=%.3f',sqrt(mean(SKhErr(1,:).^2)));
legend([h1 h2 h3 h6],BMEpdf1Legend,BMEpdf2Legend,SKMEpdfLegend,SKhLegend);
ax=axis;
axis([zmin zmax ax(3) ax(4)]);
xlabel('Estimation Error');
ylabel('Frequency');



