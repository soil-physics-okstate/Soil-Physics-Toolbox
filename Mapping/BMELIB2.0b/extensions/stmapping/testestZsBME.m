% testestZsBME  % tests estZsBME

% generate some data

usemeantrend=1; % 0 for no mean trend, 1 to use mean trend
harddatatype=1; % format of hard and soft data: 0 for [], 1 for cell, 2 for .mat, 3 for .dat
softdatatype=1; % format of hard and soft data: 0 for [], 1 for cell, 2 for .mat, 3 for .dat
dimspace=2;     % dimension of the spatial data (1, 2 or 3)
isST=1;         % 1 for s/t data, 0 for purely spatial
stformat=1;     % format of the space time data: 1 for stg, 2 for stg with some missing points, 3 for purely stv
estgridtype=2;  % 0 for [], 1 for pk, 2 for cell array

if dimspace~=2 & usemeantrend==1,
  disp('Since it is not a 2D case, usemeantrend is reset to 0');
  usemeantrend=0;
end
meantrend=[];

rand('state',0);

% harddata={ph,zh}; or filename.mat or filename.dat
nh=20;
ph=rand(nh,dimspace);
if dimspace~=2 
  zh=rand(nh,1);
  if isST
    th=rand(nh,1);
    ph=[ph th];
  end
elseif isST
  switch stformat
  case {1,2}
    nMS=nh;
    cMS=ph;
    nME=3;
    tME=sort(rand(1,nME));
    Zh=rand(nMS,nME);
    if usemeantrend,
      cMSm=cMS;
      tMEm=tME;
      ms=cMS(:,1);
      mt=2*tME;
      meantrend={cMSm,tMEm,ms,mt};
      [mst]=stmeaninterp(cMSm,tMEm,ms,mt,cMS,tME);
      Zh=Zh+mst;
    end
    [ph,zh]=valstg2stv(Zh,ph,tME);
    if stformat==2
      fidx=randperm(nh*3);
      fidx=fidx(1:nh);
      ph=ph(fidx,:);
      zh=zh(fidx,:);
    end
  case 3,
    th=rand(nh,1);
    ph=[ph,th];
    zh=rand(nh,1);
  otherwise
    error('stformat must be equal to 1, 2 or 3');
  end
else
  zh=rand(nh,1);
end
nh=size(ph,1);
harddata={ph,zh};


% softdata={ps,softpdftype,nl,limi,probdens}; or filename.mat or filename.dat
ns=5;
ps=rand(ns,dimspace);
if dimspace~=2 
  if isST
    ts=rand(ns,1);
    ps=[ps ts];
  end
elseif isST
  switch stformat
  case {1,2}
    [ps,dummy]=valstg2stv(ones(ns,nME),ps,tME);
    if stformat==2
      fidx=randperm(ns*nME);
      fidx=fidx(1:ns);
      ps=ps(fidx,:);
    end
  case 3,
    ts=rand(ns,1);
    ps=[ps,ts];
  end
end
ns=size(ps,1);
softpdftype=1;
NV=3*ones(ns,1);
WidthV=kron([.1 .1 .2],ones(ns,1));
probdensV=kron([1 2 1]/.5,ones(ns,1));
zs=rand(ns,1);
[nl,limi,probdens]=simuprobabilistic(zs,softpdftype,NV,WidthV,probdensV);
softdata={ps,softpdftype,nl,limi,probdens};


% estgrid=sk or {nxpix,nypix,incldatapts,inclvoronoi,ax}
nk=15;
if estgridtype==0
  estgrid=[];
elseif estgridtype==1
  sk=rand(nk,dimspace);
  estgrid=sk;
else
  estgrid={5,5,1,1,[]};
end

% esttime=tk
if ~isST
  tk=NaN;
elseif stformat==3 | dimspace~=2 
  tk=0.5;
else
  tk=tME(2);
end
esttime=tk;

if ~isST
  %covmodel='exponentialC';
  %covparam=[1 1.2];
  covmodel={'exponentialC','gaussianC'};
  covparam={[1 1.2],[0.5 5]};
else
  %covmodel='exponentialC/exponentialC';
  %covparam=[1 1.2 3];
  covmodel={'exponentialC/exponentialC','gaussianC/exponentialC'};
  covparam={[1 1.2 3],[0.5 5 4]};
end
covariance={covmodel,covparam};

switch harddatatype
case 0
  harddata=[];
case 1
  % fine
case 2
  save harddatafile ph zh;
  harddata='harddatafile.mat';
case 3
  for i=1:dimspace+isST,
    valname{i}=sprintf('s%d',i);
  end
  if isST
    valname{dimspace+isST}='time';
  end
  valname{dimspace+isST+1}='zh';
  writeGeoEAS([ph zh],valname,'hard data','harddatafile.dat');
  harddata='harddatafile.dat';
end

switch softdatatype
case 0
  softdata=[];
case 1
  % fine
case 2
  save softdatafile ps softpdftype nl limi probdens;
  softdata='softdatafile.mat';
case 3
  writeProba({ps,ones(ns,1)},isST,softpdftype,nl,limi,probdens,'soft data','softdatafile.dat');
  softdata='softdatafile.dat';
end

[estvalues]=estZsBME(meantrend,covariance,harddata,softdata,estgrid,esttime);

if isempty(estvalues)
  return;
end

sk=estvalues{1};
XkBME=estvalues{2};
XkBMEv=estvalues{3};
tk=estvalues{4};
varianceX=estvalues{5};
varname='Arsenic';
varunit='\mug/m^3';
xlab='(Km)';
ylab='(Km)';

estvalues={sk,XkBME,XkBMEv,tk,varianceX,varname,varunit,xlab,ylab};

if dimspace~=2
  return
end

transfun=[];

sh=ph(:,1:2);
ss=ps(:,1:2);

figure;
cmap=hot;
close;
cmap=cmap(end:-1:1,:);
celloptions={[],[],cmap,[]};
vecoptions=[NaN,0.5,80,40,1];

% harddata={sh,zh,shlegend,shplotfun,shProperties,shValues,shzrange,shsizelim}; or [] or sh or filename
shlegend='hard data';
shplotfun=3;
shProperties={'Marker','MarkerSize'};
shValues={'o',12};
shzrange=[];
shsizelim=[5 20];
if isST
  th=ph(:,end);
  if stformat==1 | stformat==2
    idx= (th==tME(2));
  else
    idx= (0.4<th)  & (th<0.6);
    shlegend='hard data for 0.4<t<0.6';
  end
  sh=sh(idx,:);
  zh=zh(idx);
end
if isempty(harddata)
  sh=[];
else
  harddata={sh,zh,shlegend,shplotfun,shProperties,shValues,shzrange,shsizelim};
end

% softdata={ss,zs,sslegend,ssplotfun,ssProperties,ssValues,sszrange,sssizelim}; or [] or ss or filename
sslegend='soft data';
ssplotfun=3;
ssProperties={'Marker','MarkerSize'};
ssValues={'v',12};
sszrange=[];
sssizelim=[];
if isST
  ts=ps(:,end);
  if stformat==1 | stformat==2
    idx = (ts==tME(2));
  else
    idx = (0.4<ts) & (ts<0.6);
    sslegend='soft data for 0.4<t<0.6';
  end
  ss=ss(idx,:);
  zs=zs(idx);
end
if isempty(softdata)
  ss=[];
else
  softdata={ss,zs,sslegend,ssplotfun,ssProperties,ssValues,sszrange,sssizelim};
end
% mask={maskcontour,maskfillcolor,masklinetype};
maskcontour=0.1*[1 9.99;2 2;2 9.99;9.9 9.99;9.9 .1;8 .1;.1 .1;-1 2;.1 9.99];
maskfillcolor='b';
masklinetype='m';
mask={maskcontour,maskfillcolor,masklinetype};

% areas={areascontour,areaslinetype,areaslegend};
areascontour={0.1*[1 1;4 4;4 1];0.1*[9 9;6 6;6 8]};
areaslinetype='-m';
areaslegend='areas';
areas={areascontour,areaslinetype,areaslegend};

plotZsBME(estvalues,1,2,transfun,celloptions,vecoptions,harddata,softdata,mask,areas);
plotZsBME(estvalues,4,2,transfun,celloptions,vecoptions,sh,ss,mask,areas);
