function BMEINTLIBtest(testtype,idxtest);

% BMEINTLIBtest             - Test for the bmeintlib directory (Jan 1, 2001) 
%
% Performs a serie of tests on BMEinterval functions (that is 
% BMEintervalMode, BMEintervalPdf, BMEintervalTMode, and BMEintervalTPdf)
%
% SYNTAX :
%
% BMEINTLIBtest(testtype,idxtest);
%
% INPUT :
%
% testtype   scalar   optional parameter indicating which function to test
%                     testtype=1 to test BMEintervalPdf and BMEintervalMode
%                     testtype=2 same as testtype=1, but WITH soft data at estimation pt
%                     testtype=3 to test BMEintervalPdf and BMEintervalMode
%                     testtype=4 same as testtype=3, but WITH soft data at estimation pt
%                     The default value is testtype=1.
% idxtest    vector   vector containing the testnumber to run, 
%                     The default value is idxtest=1:26.

if nargin<2, 
  disp('Note: Syntax is BMEINTLIBtest(testtype,idxtest);');
end;  

if nargin<1 
  testtype=1; 
  disp('Using default value testtype=1;');
end;
if nargin<2 
  idxtest=1:26; 
  disp('Using default value idxtest=1:26;');
end;

switch testtype
case 1,  
  disp('Testing functions BMEintervalPdf and BMEintervalMode');
case 2,  
  disp('Testing functions BMEintervalPdf and BMEintervalMode WITH soft data at estimation pt');
case 3,  
  disp('Testing functions BMEintervalTPdf and BMEintervalTMode');
case 4,  
  disp('Testing functions BMEintervalTPdf and BMEintervalTMode WITH soft data at estimation pt');
end;


ifail=[];
 for i=1:length(idxtest)
  itest=idxtest(i);
  figure; 
  hold on; 
  [momC,AfromPdf]=testOne(itest,testtype);
  disp(sprintf('  According to BMEprobaPdf:      mean=%10f,  var=%10f,  skewness=%10f',momC));
  disp(sprintf('                                 numerical integration of pdf=%10f',AfromPdf));
  if itest==10 | itest==20,
    disp('Take a minute to look at the figure, then press a key to close all figures and continue this test');
    pause;
    close('all');
  end;
end;


disp('test complete');

%----------------------------------------------------
function [momentsFromPdf,AfromPdf]=testOne(testnumber,testtype)
% testmoments            - Perform one test of the BMEproba functions
%
% SYNTAX :
%  [momentsFromPdf,AfromPdf]=testOne(testnumber,testtype);
% 
% INPUT (and their default values) :
% testnumber=1    a number between 1 and 26 specifying which test to run
% testtype=2      selects which function to test
%          testtype=1 to test BMEprobaCI, BMEprobaMoments, BMEprobaMode
%          testtype=2 same as testtype=1, but WITH soft data at estimation pt
%          testtype=3 to test BMEprobaTCI and BMEprobaTMode
%          testtype=4 same as testtype=3, but WITH soft data at estimation pt
%
% OUTPUT:
% momentsFromPdf   1 by 3 vector with the mean, variance and coef of skewness
%                  calculated by integration of the BME posterior pdf obtained 
%                  from BMEprobaPdf
% AfromPdf         integration value of the pdf calculated with BMEprobaPdf, 
%                  should be equal to 1.0 within numerical error

%nargin=0;
if nargin<1, testnumber=1; end;  
if nargin<2, testtype=2; end;  
  
nk=1;
nh=1;
ns=1;
n=[nk nh ns];
nSim=1;
ck=[0.5 0.5];
ch=[0 0];
cs=[1 0.9];
zh=[1.4];
zl=[0.1];
zu=[1.1];

%parameters for BMEprobaPdf and BMEprobaMoments
cNugget=0.05;               
cc=1;aa=1.0;at=1;
maxpts=1000000;aEps=0.0;
switch testtype
case {2,4}           % For test type 2 and 4 need big accuracy
  rEps=0.00002;
  nBMEpdfpts=200;
case {1,3},
  rEps=0.02;
  nBMEpdfpts=50;
end
nhmax=4; nsmax=4;
covmodel={'nuggetC','exponentialC'};
covparam={cNugget,[cc 3*aa]};
dmax=100;
order=NaN;
options(1)=0;
options(3)=maxpts;
options(4)=rEps;

%Set the transform parameters in case they are needed
%yfile=[-6:0.1:6]';
%Fyfile=gausscdf(yfile,[0 1]);
yfile=[0.01:0.01:0.09 0.1:0.1:6 10:5:25]';
Fyfile=gausscdf(log(yfile),[0 1]);

switch(testnumber)
case 1,
  disp('test 1');  
case 2,
  disp('test 2');
  nh=3;
  ch=[ch;[0.1 0.2;0.3 0.4]];
  zh=[zh;1.2;1.3];
case 3,
  disp('test 3');
  ns=2;
  cs=[cs;[0.1 0.2]];
  zl=[zl;0.1];
  zu=[zu;0.3];
case 4,
  rEps=0.02;options(4)=rEps;
  disp('test 4');
  ns=3;
  cs=[cs;[0.1 0.2;0.3 0.5]];
  zl=[zl;0.01;0.02];
  zu=[zu;1.0;1.1];
case 5,
  disp('test 5');
  ck=ch+0.001;
case 6,
  disp('test 6');
  ck=cs+0.00001;
case 7,
  disp('test 7');
  nh=0;
  ch=zeros(0,2);
  zh=zeros(0,1);
case 8,
  disp('test 8');
  ns=0;
  cs=[];
  zl=[];
  zu=[];
case 9,
  disp('test 9');
  nhmax=0;
case 10,
  disp('test 10');
  nh=3;
  ch=[ch;[10 0.2;30 0.4]];
  zh=[zh;1.2;1.3];
  nhmax=1;
case 11,
  disp('test 11');
  nsmax=0;
case 12,
  disp('test 12');
  ns=3;
  cs=[cs;[10 0.2;30 0.5]];
  zl=[zl;0.01;0.02];
  zu=[zu;1.0;1.1];
  nsmax=1;
case 13,
  disp('test 13');
  ns=3;
  cs=[cs;[10 0.2;30 0.5]];
  zl=[zl;0.01;0.02];
  zu=[zu;1.0;1.1];
  nsmax=2;
case 14,
  disp('test 14');
  it=3;covParam(2,1)=it;
  covmodel{2}='gaussianC';
  covparam{2}(2)=sqrt(3)*aa;
case 15,
  disp('test 15');
  ck=[0.5 0.5 1];
  ch=[0 0 2];
  cs=[1 0.9 3];
  it=3;covParam(2,1)=it;
  covmodel{2}='gaussianC';
  covparam{2}(2)=sqrt(3)*aa;
case 16,
  disp('test 16');
  ck=[0.5 0.5 1];
  ch=[0 0 2];
  cs=[1 0.9 3];
  dmax=[100 100 1];
  covmodel={'nuggetCST','gaussianCST'};
  covparam={cNugget,[cc 3*aa 3*at]};
case 17,
  disp('test 17');
  nh=2;
  ns=2;
  ch=[ch;[0.1 0.2]];
  cs=[cs;[0.11 0.12]];
  ck={ck,[1]};
  ch={ch,[1;2]};
  cs={cs,[1;2]};
  zh=[zh;1.2];
  zl=[zl;0.1];
  zu=[zu;1];
  nv=max([ck{2};ch{2};cs{2}]);
  nhmax=nhmax*ones(1,nv);
  nsmax=nsmax*ones(1,nv);
  covparam={ {cNugget*ones(nv,nv),[]},
             {cc*ones(nv,nv), [3*aa]} };
case 18,
  disp('test 18');
  nh=2;
  ns=2;
  ck=[0.5 0.5 1];
  ch=[0 0 2;0.1 0.2 3];
  cs=[1 0.9 3;0.11 0.12 1.5];
  ck={ck,[1]};
  ch={ch,[1;2]};
  cs={cs,[1;2]};
  zh=[zh;1.2];
  zl=[zl;0.1];
  zu=[zu;0.3];
  nv=max([ck{2};ch{2};cs{2}]);
  nhmax=nhmax*ones(1,nv);
  nsmax=nsmax*ones(1,nv);
  dmax=[100 100 1];
  covmodel={'nuggetCST','gaussianCST'};
  covparam={ {cNugget*ones(nv,nv),[]},
             {cc*ones(nv,nv), [3*aa 3*at]} };
case 22,
  disp('test 22');
  nh=3;
  ch=[ch;[0.1 0.2;0.3 0.4]];
  zh=[zh;1.2;1.3];
  order=0;
case 23,
  disp('test 23');
  nh=3;
  ch=[ch;[0.1 0.2;0.3 0.4]];
  zh=[zh;2.2;3.3];
  nsmax=0;
  order=0;
case 24,
  disp('test 24');
  nhmax=0;
  ns=2;
  cs=[1 0.9;0.11 0.12];
  zl=[zl;0.6];
  zu=[zu;2.3];
  order=0;
case 25,
  disp('test 25');
  nh=4;
  ns=2;
  ck=[0.5 0.5 1];
  ch=[0 0 2;0.1 0.2 3;1 0 1.1;0 2 1.1];
  cs=[1 0.9 3;0.11 0.12 1.5];
  ck={ck,[1]};
  ch={ch,[1;1;2;2]};
  cs={cs,[1;2]};
  zh=[zh;1.2;1.3;1.5];
  zl=[zl;0.1];
  zu=[zu;0.3];
  nv=max([ck{2};ch{2};cs{2}]);
  nhmax=nhmax*ones(1,nv);
  nsmax=nsmax*ones(1,nv);
  dmax=[100 100 1];
  covmodel={'nuggetCST','gaussianCST'};
  covparam={ {cNugget*ones(nv,nv),[]},
    {cc*ones(nv,nv), [3*aa 3*at]} };
  order=0;
  index1=(ch{2}==1);
  index2=(ch{2}==2);
case 26,
  disp('test 26');
  nh=8;
  ns=2;
  ch=[ch;[0.1 0.2;0.3 0.4;0.5 0.6;0.7 0.8;0.9 1.2;1.3 1.4;1.1 1.2]];
  zh=[zh;1.2;1.3;1.4;2.1;2.5;2.7;2.8];
  cs=[cs;0.11 0.12];
  zl=[zl;0.1];
  zu=[zu;0.3];
case 27,
  disp('test 27');
  nhmax=0;
  nsmax=0;
otherwise
  disp('No test chosen, reverting to test 1');
end;

if iscell(ck)
  maxindex=max([ck{2};ch{2};cs{2}]);
  yfilei=yfile;
  Fyfilei=Fyfile;
  yfile=[];
  Fyfile=[];
  for i=1:maxindex,
    yfile{i}=yfilei;
    Fyfile{i}=Fyfilei;
  end;
end;

hold on;
h=[];

softpdftype=1;
nl=2*ones(ns,1);
limi=[zl zu];
probdenstemplate=ones(ns,1);
if ns>0
  [probdens,norm]=proba2probdens(softpdftype,nl,limi,probdenstemplate);
end

if testtype==2 | testtype==4
  if ns>0 & nsmax>0
    if ~iscell(ck)
      ck=cs(end,:);
    else
      ck{1}=cs{1}(end,:);
      ck{2}=cs{2}(end,:);
    end
    h=probaplot(softpdftype,nl,limi,probdens,'m-',length(nl)); %plotting the soft data at the est point
    leg{1}='Soft data at est pt';
  end;
end; 

zmin=min([zh;zl;0])-0.1;
if zmin<=0, zmin=0.001; end;
zmax=max([zh;zu;0])+0.1;
zstep=(zmax-zmin)/(nBMEpdfpts-1);
z=(zmin:zstep:zmax)';

if (testtype==3 & size(cs,1)>2) | (testtype==4 & size(cs,1)-1>2)
  ns=ns-1;
  if ~iscell(cs),
    cs=cs(1:size(cs,1)-1,:);
  else
    cs{1}=cs{1}(1:size(cs{1},1)-1,:);
    cs{2}=cs{2}(1:size(cs{2},1)-1,:);
  end
  zl=zl(1:end-1);
  zu=zu(1:end-1);
end;

switch testtype
case {1,2},
  
  %zh=log(zh);
  
  disp('Calculating the mode of the BME posterior pdf using BMEintervalMode');
  [zMode,dummy]=BMEintervalMode(ck,ch,cs,zh,zl,zu,covmodel,covparam,...
    nhmax,nsmax,dmax,order,options);
  
  z=sort([z;zMode]);
  disp('Calculating BME posterior pdf using BMEintervalPdf');
  [z,pdf,info]=BMEintervalPdf(z,ck,ch,cs,zh,zl,zu,covmodel,...
    covparam,nhmax,nsmax,dmax,order,options);
  hpdf=plot(z,pdf);
  hold on;
  
  pdfzMode=pdf(z==zMode);
  if length(pdfzMode)>1, pdfzMode=max(pdfzMode); end;
  hMode=plot(zMode,pdfzMode,'r+');
  
  [meanC varC skewC AfromPdf]=pdfstat(z,pdf);
  momentsFromPdf=[meanC varC skewC];
  
case {3,4},
      
    
  disp('Calculating the mode of the BME posterior pdf using BMEintervalTMode');
  [zMode,dummy]=BMEintervalTMode(ck,ch,cs,zh,zl,zu,covmodel,covparam,...
    nhmax,nsmax,dmax,yfile,Fyfile,options);
  
  z=sort([z;zMode]);
  disp('Calculating BME posterior pdf using BMEintervalTPdf');
  [z,pdf,info]=BMEintervalTPdf(z,ck,ch,cs,zh,zl,zu,covmodel,...
    covparam,nhmax,nsmax,dmax,yfile,Fyfile,options);
  hpdf=plot(z,pdf);
  hold on;
  
  pdfzMode=pdf(z==zMode);
  if length(pdfzMode)>1, pdfzMode=max(pdfzMode); end;
  hMode=plot(zMode,pdfzMode,'r+');
  
  [meanC varC skewC AfromPdf]=pdfstat(z,pdf);
  momentsFromPdf=[meanC varC skewC];
  
end;

if ~isempty(hpdf)
  h(length(h)+1)=hpdf;
  leg{length(h)}='BME posterior pdf';
  h(length(h)+1)=hMode;
  leg{length(h)}='Mode';
  legend(h,leg,0);
end

xlabel('z');
ylabel('pdf');
title(sprintf('test %d',testnumber));
