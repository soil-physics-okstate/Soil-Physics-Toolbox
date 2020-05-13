function IOLIBtest(testtype,idxtest);

% IOLIBtest                 - Test for the iolib directory (Jan 1, 2001)
%
% Performs a serie of tests of the writeBMEproba and readBMEproba 
% functions. The tests consists in generating the input BME parameters
% and data for the BMEproba functions, writing them in a file, reading
% the file, and comparing the orginal data with that obtained after
% the write/read procedure.
%
% SYNTAX :
%
% IOLIBtest(testtype,idxtest);
%
% INPUT :
%
% testtype   scalar   optional parameter indicating which function to test
%                     testtype=1 to test BMEprobaMode
%                     testtype=2 same as testtype=1, but WITH soft data at estimation pt
%                     testtype=3 to test BBMEprobaTMode
%                     testtype=4 same as testtype=3, but WITH soft data at estimation pt
%                     testtype=5 to test BMEprobaMoments
%                     testtype=6 same as testtype=5, but WITH soft data at estimation pt
%                     The default value is testtype=1.
% idxtest    vector   vector containing the testnumber to run, 
%                     The default value is idxtest=1:26.

if nargin<2, 
  disp('Note: Syntax is IOLIBtest(testtype,idxtest);');
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
  disp('Testing IOLIB functions with BMEprobaMode');
case 2,  
  disp('Testing IOLIB functions with BMEprobaMode and WITH soft data at estimation pt');
case 3,  
  disp('Testing IOLIB functions with BMEprobaTMode');
case 4,  
  disp('Testing IOLIB functions with BMEprobaTMode and WITH soft data at estimation pt');
case 5,  
  disp('Testing IOLIB functions with BMEprobaMoments');
case 6,  
  disp('Testing IOLIB functions with BMEprobaMoments and WITH soft data at estimation pt');
end;


ifail=[];
for i=1:length(idxtest)
  itest=idxtest(i);
  ifail(i)=testwritereadOne(itest,testtype);
end;

idxfail=find(ifail>0);
if isempty(idxfail)
  disp('It appears that all tests were successful');
else
  disp('Attention!!: The following tests might have failed:');
  disp(sprintf('%d ',idxfail));
end;

disp('test complete');

%-----------------------------------------------------
function [ifail]=testwritereadOne(testnumber,testtype)
% testwritereadOne  - Perform one test of writeBMEproba and readBMEproba
%
% SYNTAX :
%  [ifail]=testreadwriteOne(testnumber,testtype);
% 
% INPUT (and their default values) :
% testnumber=1    a number between 1 and 26 specifying which test to run
% testtype=1      selects which function to test betwwen write and read
%          testtype=1 to test BMEprobaMode
%          testtype=2 same as testtype=1, but WITH soft data at estimation pt
%          testtype=3 to test BBMEprobaTMode
%          testtype=4 same as testtype=3, but WITH soft data at estimation pt
%
% OUTPUT:
% ifail    reports whether the test fail or succeded         
%          0- test succeded
%          1- BME calculation was different between a original data and 
%             data obtained after write and read of original data
%          2- Data file changed between first data file written and 
%             the second data file obtained after write and read  of the
%             first data file.

%nargin=0;
if nargin<1, testnumber=25; end;  
if nargin<2, testtype=1; end;  

ifail=0;

fileNames1={'BME1.par';'BMEhard1.dat';'BMEproba1.dat';'BMEest1.dat';'BME1.out';'BME1.dbg'};

fileNames2={'BME2.par';'BMEhard2.dat';'BMEproba2.dat';'BMEest2.dat';'BME2.out';'BME2.dbg'};

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
case {2,4,6}           % For test type 2 and 4 need big accuracy
  rEps=0.00002;
  nBMEpdfpts=200;
case {1,3,5},
  rEps=0.02;
  nBMEpdfpts=50;
end
nhmax=4; nsmax=4;
covmodel={'nuggetC','exponentialC'};
covparam={cNugget,[cc 3*aa]};
%covmodel='exponentialC';
%covparam=[cc 3*aa];
dmax=100;
PCI=[0.25 0.5 0.75 0.9];
order=NaN;
options(1)=0;
options(2)=1e-04;
options(3)=maxpts;
options(4)=rEps;
options(6)=nBMEpdfpts;
options(7)=0.001;
options(8)=3;
options(14)=100;
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
  ch=zeros(0,2);
  zh=zeros(0,1);
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
  %probdens=[1 2 2];
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
  cs=[cs;0.11 0.12];
  nl=[nl;2];
  limi=[limi;0.1 0.3 NaN NaN];
  probdens=[probdens;5 NaN NaN];
  probdens=proba2probdens(softpdftype,nl,limi,probdens);
  order=2;
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

if testtype==2 | testtype==4 | testtype==6
  if ns>0 & nsmax>0
    if ~iscell(ck)
      ck=cs(end,:);
    else
      ck{1}=cs{1}(end,:);
      ck{2}=cs{2}(end,:);
    end
  end;
end; 

switch testtype,
case {1,2},
  
  BMEmethod='BMEprobaMode';

  [zMode1,dummy]=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
  writeBMEproba(BMEmethod,fileNames1,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
  [zMode2,dummy]=BMEprobaFromFile(fileNames1{1});
  
  if abs((zMode1-zMode2)/zMode1)>5e-4
    ifail=1; 
    disp(sprintf('warning: Different mode estimates %f and %f',zMode1,zMode2));
  end; 
  
  clear ck ch cs zh softpdftype nl limi probdens covmodel covparam;
  clear nhmax nsmax dmax order options;
  
  [BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,...
      covparam,nhmax,nsmax,dmax,order,options]=readBMEproba(fileNames1{1});
  
  writeBMEproba(BMEmethod,fileNames2,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
case {3,4},
  
  BMEmethod='BMEprobaTMode';
  
  [zMode1,dummy]=BMEprobaTMode(ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,...
    yfile,Fyfile);
  
  writeBMEproba(BMEmethod,fileNames1,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,...
    yfile,Fyfile);
  
  [zMode2,dummy]=BMEprobaFromFile(fileNames1{1});
  
  if abs((zMode1-zMode2)/zMode1)>5e-4
    ifail=1; 
    disp(sprintf('warning: Different mode estimates %f and %f',zMode1,zMode2));
  end; 
  
  clear ck ch cs zh softpdftype nl limi probdens covmodel covparam;
  clear nhmax nsmax dmax order options;
  
  [BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,...
      covparam,nhmax,nsmax,dmax,order,options,yfile,Fyfile]=...
    readBMEproba(fileNames1{1});
  
  writeBMEproba(BMEmethod,fileNames2,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,...
    yfile,Fyfile);
  
case {5,6},
  
  BMEmethod='BMEprobaMoments';

  [zMoments1,dummy]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);

  writeBMEproba(BMEmethod,fileNames1,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
  [zMoments2,dummy]=BMEprobaFromFile(fileNames1{1});
  
  if abs((zMoments1(1)-zMoments2(1))/zMoments1(1))>5e-4
    ifail=1; 
    disp(sprintf('warning: Different mode estimates %f and %f',zMoments1(1),zMoments2(1)));
  end; 
  
  clear ck ch cs zh softpdftype nl limi probdens covmodel covparam;
  clear nhmax nsmax dmax order options;
  
  [BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,...
      covparam,nhmax,nsmax,dmax,order,options]=readBMEproba(fileNames1{1});
  
  writeBMEproba(BMEmethod,fileNames2,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  

  
end;

%%%%% Compare the data files

fidx=[2 3 4];
for ifidx=1:length(fidx)
  i=fidx(ifidx);
  fsame=fcmp(fileNames1{i},fileNames2{i});
  if fsame==-1,
    warning(sprintf('Could not open file %s or %s',fileNames1{i},fileNames2{i}));
    ifail=2;
  elseif fsame~=1   
    warning(sprintf('files %s and %s are different',fileNames1{i},fileNames2{i}));
    ifail=2;
  end;
end;

%%%%% Compare the parameters files

fileNames2bis={};
[fileNames2bis{1:6,1}]=deal(fileNames2{1},fileNames1{2:6});
switch testtype,
case {1,2,5,6},
  writeBMEproba(BMEmethod,fileNames2bis,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
case {3,4},
  writeBMEproba(BMEmethod,fileNames2bis,ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,...
    yfile,Fyfile);
end;
fsame=fcmp(fileNames1{1},fileNames2{1});
if fsame==-1,
  warning(sprintf('Could not open file %s or %s',fileNames1{1},fileNames2{1}));
  ifail=2;
elseif fsame~=1   
  warning(sprintf('files %s and %s are different',fileNames1{1},fileNames2{1}));
  ifail=2;
end;