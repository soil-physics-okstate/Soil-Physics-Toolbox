% BMECATLIBtest             - Test for the bmecatlib directory
%                             (December 1, 2003)
%
% Performs a series of test for the programs in the bmecalib
% directory. These test are checkin the internal consistency
% of the results provided by the various programs in specific
% situations

%%% Load data and built all needed parameters
%%% (same as for BMECATLIBtutorial.m)

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
[Pfit]=probatablefit(dfit,d,P,o,kstd,[0 1 1]);

%%% Select a subset of data around a random estimation location
%%% for faster computations

ck=rand(1,2)*100;
[c,zh,dh,nh]=neighbours(ck,c,zh,10,100);

%%% Test 1 - BMEcatHard BMEcatSubset and BMEcatPdf should
%%%          give same results if each subset corresponds
%%%          to a single category and each Pdf has probability
%%%          one for the corresponding category and probability
%%%          zero for the others.

disp('-------- Test 1 ...');

nmax=3;
dmax=20;
options=[0 1e-2];

zs=zeros(nh,4);
for i=1:nh;
  zs(i,zh(i))=1;
end;
    
[pkHard]=BMEcatHard(ck,c,zh,dfit,Pfit,nmax,dmax,options);
[pkSubset]=BMEcatSubset(ck,c,zs,dfit,Pfit,nmax,dmax,options);
[pkPdf]=BMEcatPdf(ck,c,zs,dfit,Pfit,nmax,dmax,options);

result=[pkHard;pkSubset;pkPdf]

if sum(abs(diff(result)))==0,
  disp('results are OK');
else
  disp('something went wrong');
end;
disp(' ');

%%% Test 2 - BMEcatSubset should give similar results if locations
%%%          with zero information content (all categories are possible)
%%%          are considered or not in the initial parameters.

disp('-------- Test 2 ...');

zs(1:5,:)=1;
[pkSubset1]=BMEcatSubset(ck,c,zs,dfit,Pfit,nmax,dmax,options);
[pkSubset2]=BMEcatSubset(ck,c(6:10,:),zs(6:10,:),dfit,Pfit,nmax,dmax,options);

result=[pkSubset1;pkSubset2]

if sum(abs(diff(result)))==0,
  disp('results are OK');
else
  disp('something went wrong');
end;
disp(' ');

%%% Test 3 - BMEcatSubset and BMEcatPdf should give similar results
%%%          when probabilities for all possible categories are equal.

disp('-------- Test 3 ...');

r=floor(rand(nh,1)*4)+1;
zs=zeros(nh,4);
ps=zeros(nh,4);
for i=1:nh,
  zs(i,1:r(i))=1;
  ps(i,:)=zs(i,:)/r(i);
end; 

[pkSubset]=BMEcatSubset(ck,c,zs,dfit,Pfit,nmax,dmax,options);
[pkPdf]=BMEcatPdf(ck,c,ps,dfit,Pfit,nmax,dmax,options);

result=[pkSubset;pkPdf]

if sum(abs(diff(result)))<100*eps,
  disp('results are OK');
else
  disp('something went wrong');
end;
disp(' ');

%%% Test 4 - BMEcatPdf should provide intermediate results beween
%%%          BMEcatSubset where all categories are possible and
%%%          BMEcatSubset where only one category is possible,
%%%          with these two situations as limit cases

disp('-------- Test 4 ...');

zs=zeros(nh,4);
for i=1:nh;
  zs(i,zh(i))=1;
end;

zs(1,:)=[1 1 1 0];
[pkSubset1]=BMEcatSubset(ck,c,zs,dfit,Pfit,nmax,dmax,options);
zs(1,:)=0;
zs(1,zh(1))=1;
[pkSubset2]=BMEcatSubset(ck,c,zs,dfit,Pfit,nmax,dmax,options);

ps1=[1/3 1/3 1/3 0];
ps2=zeros(1,4);
ps2(zh(1))=1;

compt=0;
for i=0:0.05:1;
  compt=compt+1;
  disp([num2str(compt),'/',num2str(21)])
  pstemp=zs;
  pstemp(1,:)=(1-i)*ps1+i*ps2;  
  [pkPdf]=BMEcatPdf(ck,c,pstemp,dfit,Pfit,nmax,dmax,options);
  PKPdf(compt,:)=pkPdf;
end;

figure;
plot(0:0.05:1,PKPdf);
hold on;
plot(0,pkSubset1,'o');
plot(1,pkSubset2,'o');
drawnow;

%%% Test 6 - simucatcond and simucatcondPdf should provide
%%%          same results when only hard information is used.

disp('-------- Test 6 ...');

ps=zeros(nh,4);
for i=1:nh,
  ps(i,zh(i))=1;    
end;

options=[1 1 1e-3];
time=clock;
seed=sum(100*time(5:6));
Zsim1=[];Zsim2=[];

rand('state',seed);
cksim=rand(5,2)*40-20;
cksim(:,1)=cksim(:,1)+ck(1,1);
cksim(:,2)=cksim(:,2)+ck(1,2);
[zhsim1]=simucatcond(cksim,c,zh,dfit,Pfit,nmax,dmax,options);

rand('state',seed);
cksim=rand(5,2)*40-20;
cksim(:,1)=cksim(:,1)+ck(1,1);
cksim(:,2)=cksim(:,2)+ck(1,2);
[zhsim2]=simucatcondPdf(cksim,c,ps,dfit,Pfit,nmax,dmax,options);
  
result=[zhsim1';zhsim2']

if sum(abs(diff(result)))<100*eps,
  disp('results are OK');
else
  disp('something went wrong');
end;
disp(' ');

disp('======== BMECATLIBtest is over');