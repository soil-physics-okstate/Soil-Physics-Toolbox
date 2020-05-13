% testplotZsBME

% generate some data

plotestvalues=1;  %0 (empty) or 1 or 2 (file)
plottrans=2;      %0 , 1, 2 or 3
plotdata=1;       %0 (empty), 1, 2 (just sh,ss), or 3 (file)
plotmask=1;       %0 (empty), 1, 2 (just maskcontour), or 3 (file)
plotareas=1;      %0 (empty), 1, or 2 (file)

rand('state',0);

% harddata={sh,zh,shlegend,shplotfun,shProperties,shValues,shzrange,shsizelim}; or [] or sh or filename
nh=10;
sh=rand(nh,2);
zh=rand(nh,1)+0.01;
if plottrans>0
  zh=log(zh);
end
shlegend='hard data';
shplotfun=2;
shProperties={'Marker'};
shValues={'o'};
shzrange=[];
shsizelim=[5 20];
harddata={sh,zh,shlegend,shplotfun,shProperties,shValues,shzrange,shsizelim};

% softdata={ss,zs,sslegend,ssplotfun,ssProperties,ssValues,sszrange,sssizelim}; or [] or ss or filename
ns=6;
ss=rand(ns,2);
zs=rand(ns,1)+0.01;
if plottrans>0
  zs=log(zs);
end
sslegend='soft data';
ssplotfun=3;
ssProperties={'Marker','MarkerSize'};
ssValues={'v',15};
sszrange=[];
sssizelim=[];
softdata={ss,zs,sslegend,ssplotfun,ssProperties,ssValues,sszrange,sssizelim};

% estvalues={sk,XkBME,XkBMEv,tk,varianceX,varname,varunit,xlab,ylab}; or [] or filename
covmodel='exponentialC';
covparam=[1 1];
nhmax=10;
order=NaN;
nsmax=3;
dmax=2; 
options=BMEoptions;
[xk yk]=meshgrid([0:0.2:1],[0:0.2:1]);
sk=[[xk(:) yk(:)];sh;ss];
softpdftype=1;
[moments,info]=BMEprobaMoments(sk,[sh;ss],[],[zh;zs],softpdftype,[],[],[],covmodel,covparam,nhmax,nsmax,dmax,order,options);

XkBME=moments(:,1);
XkBMEv=moments(:,2);

tk=2;
varname='Arsenic';
varunit='\mug/m^3';
varianceX=1;
xlab='(Km)';
ylab='(Km)';

estvalues={sk,XkBME,XkBMEv,tk,varianceX,varname,varunit,xlab,ylab};

% transfun={invtransfun,transname,plottrans}; or []
% celloptions={ax,cax,cmap,contval,griddatamethod}; or []
% vecoptions=[conflevel,nonattainval,nxpix,nypix,showestpoints]; or []

if plottrans==0
  transfun=[];
else
  invtransfun='exp';
  transname='log';
  transfun={invtransfun,transname,plottrans};
end
figure;
cmap=hot;
close;
cmap=cmap(end:-1:1,:);
celloptions={[],[],cmap,[]};
vecoptions=[NaN,0.5,80,40,1];

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

switch plotestvalues
case 0,
  estvalues=[];
case 1,
  %estvalues=estvalues;
case 2,
  save estvalues sk XkBME XkBMEv tk varianceX varname varunit xlab ylab;
  estvalues='estvalues';
end

switch plotdata
case 0
  harddata=[];
  softdata=[];
case 1
  %harddata=harddata;
  %softdata=softdata;
case 2
  harddata=sh;
  softdata=ss;
case 3
  save harddata sh zh shlegend shplotfun shProperties shValues shzrange shsizelim
  harddata='harddata';
  save softdata ss zs sslegend ssplotfun ssProperties ssValues sszrange sssizelim
  softdata='softdata';
end

switch plotmask
case 0,
  mask=[];
case 1,
  %mask=mask;
case 2,
  mask=maskcontour;
case 3,
  save maskfile maskcontour maskfillcolor masklinetype
  mask='maskfile';
end

switch plotareas
case 0,
  areas=[];
case 1,
  %areas=areas;
case 2,
  save areasfile areascontour areaslinetype areaslegend
  areas='areasfile';
end

% plotZsBME(estvalues,mapwhat,mapfun,transfun,celloptions,vecoptions,harddata,softdata)

% test
disp('test');
plotZsBME(estvalues)
plotZsBME(estvalues,1,2,transfun,celloptions,vecoptions,harddata,softdata,mask,areas)
plotZsBME(estvalues,2,2,transfun,celloptions,vecoptions,sh,ss,mask,areas)
plotZsBME(estvalues,3,2,transfun,celloptions,vecoptions,sh,ss,mask,areas)
plotZsBME(estvalues,4,2,transfun,celloptions,vecoptions,sh,ss,mask,areas)
plotZsBME(estvalues,5,2,transfun,celloptions,vecoptions,sh,ss,mask,areas)
