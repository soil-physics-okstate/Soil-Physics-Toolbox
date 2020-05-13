function exampleBMEvsIndK(nSim,useStoredresults,covtype)

% exampleBMEvsIndK          - Comparison between BMEinterval and Indicator Kriging (Jan 1,2001)
% 
% In this example with interval soft data the BME method 
% is compared with the Indicator Kriging method. The 
% steps in this example are the following:
%
% SYNTAX :
% 
% exampleBMEvsIndK(nSim,useStoredresults,covtype);
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
% covtype            string     optional input indicating the covariance model model
%                               1 for exponenetial covariance model
%                               2 for gaussian covariance model
%                               default value is 1.
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
if nargin<3, covtype=1; end;

if covtype==1, covmodel='exponentialC';
elseif covtype==2, covmodel='gaussianC';
else, error('covtype must equal to 1 or 2');
end;

if useStoredresults==1 & exist(['exampleBMEvsIndK' covmodel '.mat'])~=2 
  disp('Cannot find stored results, therefore calculating them');
  useStoredresults=0;
end;

if useStoredresults==0
  
  disp(sprintf('Calculating results using %d simulations for covmodel %s',nSim,covmodel));
  
  %
  %initialization
  %
  rand('state',16);
  randn('state',3);

  %
  % Set the coordinates of the estimated, hard and soft data points
  %
  nk=1;nh=2;ns=3;
  ck=[0.4,0.5];
  ch=[0.6,0.5
    0.43384,0.46045];
  cs=[0.92399,0.14629
    0.3023,0.97937
    0.34654,0.60146];
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
  % Sets soft interval generation parameters
  %
  Vmethod=3;
  Pkzk=[0.1 0.3 0.5 0.7 0.9];
  zkGrid=gaussinv(Pkzk,[0 1]);
  SoftInterval=zkGrid;
  Vvalue=0;

  %
  % Set the covariance model and BME parameters
  %
  cNugget=0.0;
  cc=1;
  aa=9;
  covparam=[cc aa];
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
  [aAll,bAll]=simuinterval(zsAll,zkGrid,Vmethod,Vvalue);
  
  %
  %doing the BME estimation with BMEintervalMode
  %
  disp(sprintf('Doing BME estimation'));
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    a=aAll(:,iSim);
    b=bAll(:,iSim);
    [zk,info]=BMEintervalMode(ck,ch,cs,zh,a,b,...
      covmodel,covparam,nhmax,nsmax,dmax,order,options);
    zkBME(1,iSim)=zk;
  end;  
  tsBME=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsBME));
  
  %
  %doing the SK estimation
  %
  disp('Doing Simple Kriging (SK) estimation');
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    a=aAll(:,iSim);
    b=bAll(:,iSim);
    zs=0.5*(a+b);
    vs=(1/12)*(a-b).^2;
    [zk,vk]=kriging(ck,[ch;cs],[zh;0.5*(a+b)],covmodel,covparam,nhmax,dmax,order,options);
    zkSK(1:nk,iSim)=zk;
  end
  tsSKME=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsSKME));
    
  %
  %doing the Indicator Kriging estimation
  %
  disp('Doing Indicator Kriging (IndK) estimation');
  t0=clock;
  for iSim=1:nSim,  
    zh=zhAll(:,iSim);
    zs=zsAll(:,iSim);
    for i=1:length(Pkzk)
      [Fzk(i,iSim),vk]=krigingIndGauss(ck,[ch;cs],[zh;zs],covmodel,covparam,Pkzk(i),nhmax,dmax,options);
    end
    zkIndK(1,iSim)=cdf2median(zkGrid(:),Fzk(:,iSim));
  end
  tsSKME=etime(clock,t0);
  disp(sprintf(' Computation time: %7.4gs',tsSKME));
  
  %
  %save the data
  %
  
  save(['exampleBMEvsIndK' covmodel],'zkSimAll','zkIndK','zkBME','zkSK',...
    'nk','covmodel','nSim');  
  
elseif useStoredresults==1
  
  load(['exampleBMEvsIndK' covmodel]);
  disp(sprintf('Using stored results that were calculated using %d simulations for covmodel %s',nSim,covmodel));
    
else
  
  error('useStoredresults must be equal to either 0 or 1');
  
end;
    
%
%Calculating the estimation errors
%
BMEerr=zkBME-zkSimAll;
IndKerr=zkIndK-zkSimAll;
SKerr=zkSK-zkSimAll;
disp(' ');
disp('Square Root of Mean Square Error (E) for each method:');
disp('-----------------------------------------------------');
disp('        BME        IndK       SK ');
disp(sprintf('     %8.3f   %8.3f   %8.3f',...
  sqrt(mean(BMEerr.^2)),sqrt(mean(IndKerr.^2)),sqrt(mean(SKerr.^2))));

%-------------------------------------------------------------
% Plotting the figure
%
figure;
h=[];  
switch covmodel   
case 'exponentialC',
  Nb=15;
  zmin=-1;
  zmax=1;
  h(length(h)+1)=histline(BMEerr(1,:),Nb,[zmin,zmax]);hold on;
case 'gaussianC',
  Nb=15;
  zmin=-0.5;
  zmax=0.5;
  h(length(h)+1)=histline(BMEerr(1,:),Nb);hold on;
otherwise
  error('covname should be equal to gaussianC or exponentialC');
end;
leg{length(h)}=sprintf('BME  , E=%.3f',sqrt(mean(BMEerr.^2)));

h(length(h)+1)=histline(IndKerr(1,:),Nb+1,[zmin zmax]);
set(h(end),'LineStyle','--','Color','r');
leg{length(h)}=sprintf('IndK, E=%.3f',sqrt(mean(IndKerr.^2)));

legend(h,leg);
ax=axis;
axis([zmin zmax ax(3) ax(4)]);
xlabel('Estimation Error');
ylabel('Frequency');


%-------------------------------------------------------------
function median=cdf2median(z,cdf)

% cdf2median            - finds the median for an irregular cdf
%
%
if cdf(1)>=0.5, 
  median=z(1);
elseif cdf(end)<=0.5, 
  median=z(end);
else
  i1=max((1:length(z))'.*(cdf(:)<0.5));
  median=interp1([cdf(i1) cdf(i1+1)],[z(i1) z(i1+1)],0.5);
end

