function BMEPROBALIBtest(testtype,idxtest);

% BMEPROBALIBtest           - Test for the bmeprobalib directory (Jan 1, 2001) 
%
% Performs a serie of tests on BMEproba functions (that is 
% BMEprobaMode, BMEprobaMoments, BMEprobaCI, BMEprobaPdf,
% BMEprobaTMode, BMEprobaTCI and BMEprobaTPdf)
%
% SYNTAX :
%
% BMEPROBALIBtest(testtype,idxtest);
%
% INPUT :
%
% testtype   scalar   optional parameter indicating which function to test
%                     testtype=1 to test BMEprobaCI, BMEprobaMoments, BMEprobaMode
%                     testtype=2 same as testtype=1, but WITH soft data at estimation pt
%                     testtype=3 to test BMEprobaTCI and BMEprobaTMode
%                     testtype=4 same as testtype=3, but WITH soft data at estimation pt
%                     The default value is testtype=1.
% idxtest    vector   vector containing the testnumber to run, 
%                     The default value is idxtest=1:26.

if nargin<2, 
  disp('Note: Syntax is BMEPROBALIBtest(testtype,idxtest);');
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
  disp('Testing functions BMEprobaCI, BMEprobaMoments, BMEprobaMode');
case 2,  
  disp('Testing functions BMEprobaCI, BMEprobaMoments, BMEprobaMode WITH soft data at estimation pt');
case 3,  
  disp('Testing functions BMEprobaTCI and BMEprobaTMode');
case 4,  
  disp('Testing functions WITH soft data at estimation pt');
end;

ifail=[];
 for i=1:length(idxtest)
  itest=idxtest(i);
  figure; 
  hold on; 
  [mom,momC,AfromPdf]=testOne(itest,testtype);
  if ~isempty(mom)
    if abs(mom(1)-momC(1))/mom(1)>0.01 | abs(mom(2)-momC(2))/mom(2)>0.01 
      ifail=[ifail itest];
    end;
    merr(i)=abs(mom(1)-momC(1))/mom(1);
    sterr(i)=abs(mom(2)-momC(2))/mom(2);
    disp(sprintf('  According to BMEprobaMoments:  mean=%10f,  var=%10f,  skewness=%10f',mom));
  end;
  disp(sprintf('  According to BMEprobaPdf:      mean=%10f,  var=%10f,  skewness=%10f',momC));
  disp(sprintf('                                 numerical integration of pdf=%10f',AfromPdf));
  if itest==10 | itest==20,
    disp('Take a minute to look at the figure, then press a key to close all figures and continue this test');
    pause;
    close('all');
  end;
end;

if ~isempty(mom)
  
  %disp('itest  Mean Err   St Dev Err');
  %for i=1:length(idxtest)
  %  disp(sprintf('%2d  %10.4f  %10.4f',idxtest(i),merr(i),sterr(i)));
  %end;
  %disp(sprintf('ave %10.4f  %10.4f',mean(merr),mean(sterr)));
  
  if isempty(ifail)
    disp('It appears that all tests were successful');
  else
    disp('Attention!!: The following tests might have failed:');
    disp(sprintf('%d ',ifail));
  end;
  
end

disp('test complete');

%----------------------------------------------------
function [moments,momentsFromPdf,AfromPdf]=testOne(testnumber,testtype)
% testmoments            - Perform one test of the BMEproba functions
%
% SYNTAX :
%  [moments,momentsFromPdf,AfromPdf]=testOne(testnumber,testtype);
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
% moments          1 by 3 vector with the mean, variance and coef of skewness
%                  calculated with BMEproMoments
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
softpdftype=1;
nl=[4];
limi=[0.1 0.3 0.7 1.1];
probdens=[1 1.5 0.5];

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
PCI=[0.25 0.5 0.75 0.9];
order=NaN;
options(1)=0;
options(3)=maxpts;
options(4)=rEps;
options(6)=nBMEpdfpts;
options(7)=0.001;
options(8)=3;
options(20:20+length(PCI)-1)=PCI;

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
  nl=[nl;2];
  limi=[limi;0.1 0.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
case 4,
  rEps=0.02;options(4)=rEps;
  disp('test 4');
  ns=3;
  cs=[cs;[0.1 0.2;0.3 0.5]];
  softpdftype=2;
  nl=[4;3;3];
  %limi=[-0.9 -0.3 0.2 1.0;-1.0 -0.1 0.9 NaN;-0.8 0.1 1.1 NaN];
  limi=[0.01 0.03 0.2 1.0;0.01 0.09  0.9 NaN;0.02 0.1  1.1 NaN];
  probdens=[0 1 1 0;0 1 0 NaN;0 1 0 NaN];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
  options(6)=25;  
case 5,
  disp('test 5');
  ck=ch+0.001;
case 6,
  disp('test 6');
  ck=cs+0.00001;
case 7,
  disp('test 7');
  nh=0;
  ch=[];
  zh=[];
case 8,
  disp('test 8');
  ns=0;
  cs=zeros(0,2);
  nl=zeros(0,1);
  limi=zeros(0,1);
  probdens=zeros(0,1);
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
  nl=[nl;2;3];
  limi=[limi;0.1 0.3 NaN NaN;0.1 0.3 0.7 NaN];
  probdens=[probdens;5 NaN NaN;4 0.5 NaN];
  nsmax=1;
case 13,
  disp('test 13');
  ns=3;
  cs=[cs;[10 0.2;30 0.5]];
  nl=[nl;2;3];
  limi=[limi;0.1 0.3 NaN NaN;0.1 0.3 0.7 NaN];
  probdens=[probdens;5 NaN NaN;4 0.5 NaN];
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
  nl=[nl;2];
  limi=[limi;0.1 0.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
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
  nl=[nl;2];
  limi=[limi;0.1 0.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
  nv=max([ck{2};ch{2};cs{2}]);
  nhmax=nhmax*ones(1,nv);
  nsmax=nsmax*ones(1,nv);
  dmax=[100 100 1];
  covmodel={'nuggetCST','gaussianCST'};
  covparam={ {cNugget*ones(nv,nv),[]},
             {cc*ones(nv,nv), [3*aa 3*at]} };
case 19,
  disp('test 19');
  softpdftype=2;
  nl=[4];
  limi=[0.1 0.3 0.7 1.1];
  probdens=[0 2 2 1];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
case 20,
  disp('test 20');
  softpdftype=3;
  nl=[4];
  limi=[2 1 5];
  probdens=[1 2 3];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
case 21,
  disp('test 21');
  softpdftype=4;
  nl=[4];
  limi=[2 1 5];
  probdens=[0 2 2 1];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
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
  nl=[nl;2];
  limi=[limi;0.6 2.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
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
  nl=[nl;2];
  limi=[limi;0.1 0.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
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
  %zh=[zh;2.2;3.3;3.4;0.1;0.5;0.7;0.8];
  cs=[cs;0.11 0.12];
  nl=[nl;2];
  limi=[limi;0.1 0.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
  order=2;
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

if ns==1
  rEps=rEps/1000;
  options(4)=rEps;
end;

hold on;
h=[];

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
 
switch testtype
case {1,2},
  
  %zh=log(zh);
  
  disp('Calculating BME posterior pdf and confidence intervals using BMEprobaCI');
  [zlCI,zuCI,fCI,PCI,zc,pdfc]=BMEprobaCI(ck,ch,cs,zh,softpdftype,...
    nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  z=zc{1};
  pdf=pdfc{1};
  hpdf=plot(z,pdf);
  hold on;
  
  %y=exp(z);
  %y=y(y<4);
  %ypdf=interp1(z,pdf,log(y))./y;
  %plot(y,ypdf,'r');

  [meanC varC skewC AfromPdf]=pdfstat(z,pdf);
  momentsFromPdf=[meanC varC skewC];
  for iCI=1:length(PCI)
    hCI=plot([zlCI(1,iCI) zuCI(1,iCI)],[fCI(1,iCI) fCI(1,iCI)],'r');
    text(0.5*(zlCI(1,iCI)+zuCI(1,iCI)),fCI(1,iCI),...
      sprintf('%.2f %% CI',PCI(iCI)),'HorizontalAlignment','center');
  end;
  
  disp('Calculating the moments of the BME posterior pdf using BMEprobaMoments');
  [moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
  disp('Calculating the cross validation moments of the BME posterior pdf using BMEprobaMomentsXvalidation');
  [momentsXval,info,MSE,MAE,ME]=BMEprobaMomentsXvalidation(1,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options);
  [momentsXval,info,MSE,MAE,ME]=BMEprobaMomentsXvalidation(2,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
  disp('Calculating the mode of the BME posterior pdf using BMEprobaMode');
  [zMode,dummy]=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  [dummy,fzMode,dummy]=BMEprobaPdf(zMode,ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  hMode=plot(zMode,fzMode,'r+');
  
case {3,4},
  
  disp('Calculating BME posterior pdf and confidence intervals using BMEprobaTCI');
  [zlCI,zuCI,fCI,PCI,zc,pdfc]=BMEprobaTCI(ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,Fyfile);
  z=zc{1};
  pdf=pdfc{1};
  hpdf=plot(z,pdf);
  hold on;
  [meanC varC skewC AfromPdf]=pdfstat(z,pdf);
  momentsFromPdf=[meanC varC skewC];
  for iCI=1:length(PCI)
    hCI=plot([zlCI(1,iCI) zuCI(1,iCI)],[fCI(1,iCI) fCI(1,iCI)],'r');
    text(0.5*(zlCI(1,iCI)+zuCI(1,iCI)),fCI(1,iCI),...
      sprintf('%.2f %% CI',PCI(iCI)),'HorizontalAlignment','center');
  end;
  
  disp('Calculating the mode of the BME posterior pdf using BMEprobaTMode');
  [zMode,dummy]=BMEprobaTMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,Fyfile);
  [dummy,fzMode,dummy]=BMEprobaTPdf(zMode,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,Fyfile);
  hMode=plot(zMode,fzMode,'r+');
  
  moments=[];
end;

if ~isempty(hpdf)
  h(length(h)+1)=hpdf;
  leg{length(h)}='BME posterior pdf';
  h(length(h)+1)=hMode;
  leg{length(h)}='Mode';
  h(length(h)+1)=hCI;
  leg{length(h)}='CI';
  legend(h,leg,0);
end

xlabel('z');
ylabel('pdf');
title(sprintf('test %d',testnumber));
