% checksimul      - testing for simulation
%                   (please do not consider this)


clear;
clc;

load SoilsData.txt;
c=SoilsData(:,1:2);
zh=SoilsData(:,3);

cl=(0:5:50)';
method='kron';
[d,P,o]=probatablecalc(c,zh,cl,method);
dfit=(0:1:45)';
kstd=1*max(d)/length(d);
options=[0 1 1];
[Pfit]=probatablefit(dfit,d,P,o,kstd,options);

minc=[0 0];
dc=[5 5];
nc=[21 21];
[ck]=creategrid(minc,dc,nc);

nmax=3;
dmax=20;
options=[1 1 1e-3];

%%% Unconditional simulation on grid

[zhsimul]=simucat(ck,dfit,Pfit,nmax,dmax,options);

[ck1,ck2,zhsimul]=col2mat(ck,zhsimul);
figure;
pcolor(ck1,ck2,zhsimul);
colormap hot
shading flat
axis square;
xlabel('Easting');
ylabel('Northing');

%%% Conditional simulation on grid

[zhsimul]=simucatcond(ck,c,zh,dfit,Pfit,nmax,dmax,options);

[ck1,ck2,zhsimul]=col2mat(ck,zhsimul);
figure;
pcolor(ck1,ck2,zhsimul);
colormap hot
shading flat
hold on;
colorplot(c,zh);
axis square;
xlabel('Easting');
ylabel('Northing');

%%% Unconditional simulation at data points

ns=size(c,1);
zs=zeros(ns,4);
for i=1:ns,
  zs(i,zh(i))=1;    
end;

[zhsimul]=simucatcondPdf(c,c,zs,dfit,Pfit,nmax,dmax,options);