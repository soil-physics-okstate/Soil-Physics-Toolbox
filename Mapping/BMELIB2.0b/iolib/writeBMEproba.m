function []=writeBMEproba(BMEmethod,fileNames,ck,ch,cs,zh,softpdftype,nl,limi,...
  probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile);

% writeBMEproba             - Writes data files for the BMEproba functions (Jan 1, 2001)
%
% Writes to files all the parameters and data needed for a BME estimation
% with hard data and soft data of probabilistic type (i.e. estimation using 
% either the BMEprobaMode, BMEprobaTMode, BMEprobaMoment, or BMEprobaTMoment 
% functions).
%
% SYNTAX :
%
% writeBMEproba(BMEmethod,fileNames,ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
%   covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile);
%
% INPUT :
%
% BMEmethod      string      string indicating the BME method to use. The string must  
%                            be one of the following: 'BMEprobaMode', 'BMEprobaTMode',
%                            'BMEprobaMoment', or 'BMEprobaTMoment'. Use the help for
%                            any of these methods to learn more about it.
% fileNames     cell         6 by 1 cell array of strings with the following filenames
%                            fileNames{1}   file with BME parameters
%                            fileNames{2}   file with hard data
%                            fileNames{3}   file with probabilistic data
%                            fileNames{4}   file with coordinates of estimation points
%                            fileNames{5}   file for BME estimated values
%                            fileNames{6}   file for BME estimation messages
%                            See ioparamsyntax for explanations of the syntax of the
%                            BME parameter file
% ck,ch,etc.                 same arguments as for the BME estimation function specified by
%                            BMEmethod. See the help for that function to get a complete 
%                            description of the input arguments.
% 
% NOTE :
%
% After you execute this function, the BME parameters are stored in fileNames{1}
% and the data is stored in the files fileNames{2}, fileNames{3} and fileNames{4}.
% Once written to these files, you can read the data back by simply using
% the function readBMEproba(paramfile), where paramfile=fileNames{1}.
% See ioparamsyntax for explanations of the syntax of the BME parameter files,
% and see the help for readBMEproba for more explanations on reading BME data. 
% Using writeBMEproba and readBMEproba is particularly useful as a 
% platform to exchange BME data between different people.

%%%%%%%%% Checking data
if ~iscell(fileNames) error('fileNames must be a cell'); end;
if ~prod(double(size(fileNames)==[6 1])), error('fileNames must be a 6 by 1 cell'); end;
[ch,zh,cs,nl,limi,probdens]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
  covmodel,covparam,nhmax,nsmax,dmax);
[isST,isSTsep,modelS,modelT]=isspacetime(covmodel);
defaultOptions=BMEoptions;
if nargin<17,
  options=defaultOptions
else
  options=[options defaultOptions(length(options)+1:end)];
end;
if nargin<19          % If yfile is not specified then there will be no transformation
  yfile=[];           %    in which case set yfile=[]
  cdfyfile=[];
else
  FyTransformCheckArgs(yfile,cdfyfile,ck,ch,cs,zh,...
    softpdftype,nl,limi,probdens);
end;

%%%%%% Build the general data structure for the singular cases

noindex=(iscell(ck)==0);
nocell=(iscell(covmodel)==0);

if nocell==1 & noindex==1,          % case for 1 model and 1 variable
  nk=size(ck,1);                    % compute the number of coordinates
  nh=size(ch,1);
  ns=size(cs,1);
  ck={ck,ones(nk,1)};               % create an index with values 1
  ch={ch,ones(nh,1)};
  cs={cs,ones(ns,1)};
  nm=1;                             % set number of models equal to 1
  nv=1;                             % set number of variables equal to 1
  covmodel={covmodel};                    % create a cell array with model
  np=length(covparam);                 % create a cell array with param
  covparam={{covparam(1),covparam(2:np)}};
end;

if nocell==1 & noindex==0,          % case for 1 model and nv variables
  nk=size(ck{1},1);                 % compute the number of coordinates
  nh=size(ch{1},1);
  ns=size(cs{1},1);
  nm=1;                             % set number of models equal to 1
  nv=max([ck{2};ch{2};cs{2}]);            % compute the number of variables 
  covmodel={covmodel};                    % create a cell array with model
  covparam={covparam};
end;

if nocell==0 & noindex==1,          % case for nm models and 1 variable
  nk=size(ck,1);                    % compute the number of coordinates
  nh=size(ch,1);
  ns=size(cs,1);
  ck={ck,ones(nk,1)};               % create an index with values 1
  ch={ch,ones(nh,1)};
  cs={cs,ones(ns,1)};
  nm=length(covmodel);                 % compute the number of models 
  nv=1;                             % set the number of variables to 1 
  for i=1:nm,                       % create cell arrays with param vectors 
    np=length(covparam{i});
    covparam{i}={covparam{i}(1),covparam{i}(2:np)};
  end;
end;

if nocell==0 & noindex==0;          % case for nm models and nv variables
  nk=size(ck{1},1);                 % compute the number of coordinates
  nh=size(ch{1},1);
  ns=size(cs{1},1);
  nm=length(covmodel);                 % compute the number of models
  nv=max([ck{2};ch{2};cs{2}]);            % compute the number of variables 
end;

d=size(ck{1},2)-isST;

dotransform=~isempty(yfile);
if noindex==1
  yfile={yfile};
  cdfyfile={cdfyfile};
end;

%%%%% Write the parameter file
fidParam=fopen(fileNames{1},'w');
str30='                                ';     %string with 30 empty char used for padding
fprintf(fidParam,'         Parameters for BMEproba \r');
fprintf(fidParam,'            ---------------------\r');
fprintf(fidParam,'START OF PARAMETERS:\r');
str=BMEmethod;str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\BME method\r',str);
str=sprintf('%s',fileNames{2});str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\file for hard data\r',str);
str=sprintf('%d',nv);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\nv, number of variables primary + others\r',str);
column=[1:d zeros(1,3-d) isST*(d+1) d+isST+1:d+isST+nv];
colString=sprintf('%d ',column);
str=sprintf('%s',colString);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\column for x,y,z,t and variable(s)\r',str);
str=sprintf('%s',fileNames{3});str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\file for probabilistic soft data\r',str);
dataMin=min([zh;limi(:,1)]);
if dataMin>0 
  tmin=0;
  outOfRange=-1;
else
  tmin=-1e30;   
  while tmin>dataMin
     tmin=2*tmin;
  end;
  outOfRange=2*tmin;
end;
switch softpdftype,
  case {1,2}, b=limi((1:ns)'+(nl-1)*ns);    
  case {3,4}, b=limi(:,3);
end;
dataMax=max([zh;b]);
tmax=1e30;
while tmax<dataMax
   tmax=2*tmax;
end;
str=sprintf('%d %d',tmin,tmax);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\data trimming limits\r',str);
str=sprintf('%s',fileNames{4});str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\file with the estimation points\r',str);
str=sprintf('%s',fileNames{5});str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\output file for BME results\r',str);
str=sprintf('%s',fileNames{6});str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\output file for BME messages\r',str);
str=sprintf('%d',nm);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\number of covariance models\r',str);
for im=1:nm,
  str=covmodel{im};str=char(str,str30);str=str(1,:);
  fprintf(fidParam,'%s \\covariance model %d\r',str,im);
  if nv==1
    str=sprintf('%g',covparam{im}{1});str=char(str,str30);str=str(1,:);
    fprintf(fidParam,'%s \\variance for model %d\r',str,im);
    str=sprintf('%g ',covparam{im}{2});str=char(str,str30);str=str(1,:);
    fprintf(fidParam,'%s \\additional parameters for model %d\r',str,im);
  else
    for iv=1:nv
      str=sprintf('%g ',covparam{im}{1}(iv,:));str=char(str,str30);str=str(1,:);
      fprintf(fidParam,'%s \\row %d of the nv by nv matrix of cross-variances for model %d\r',str,iv,im);
    end
    str=sprintf('%g ',covparam{im}{2});str=char(str,str30);str=str(1,:);
    fprintf(fidParam,'%s \\additional parameters for model %d\r',str,im);
  end;
end;
str=sprintf('%d ',nhmax);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\max number of hard data\r',str);
str=sprintf('%d ',nsmax);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\max number of soft data\r',str);
str=sprintf('%d ',dmax);str=char(str,str30);str=str(1,:);
column=dmax(:)';
colString=sprintf('%d ',column);
str=sprintf('%s',colString);str=char(str,str30);str=str(1,:);
if length(dmax)==1
  fprintf(fidParam,'%s \\max search radius\r',str);
else
  fprintf(fidParam,'%s \\max search radius, max search time, S/T metric\r',str);
end;  
if isnan(order),
  order=-1;
end;
str=sprintf('%d',order);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\order of mean trend (use -1 for no mean trend)\r',str);
str=sprintf('%d %g',options(3),options(4));str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\maxpts and rEps for numerical integration\r',str);
str=sprintf('%g %d',options(2),options(14));str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\relative tolerance and maxeval for fminbnd\r',str);
str=sprintf('%g %g',options(6:7));str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\options(6) and options(7) to plot BME posterior pdf\r',str);
str=sprintf('%g %g',options(8));str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\options(8) for BMEprobaMoments, number of BME moments to calculate\r',str);
PCI=options(20:end);
PCI=PCI(PCI>0);
str=sprintf('%g ',PCI);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\Confidence Interval Probabilities, up to 10 values\r',str);
str=sprintf('%d',dotransform);str=char(str,str30);str=str(1,:);
fprintf(fidParam,'%s \\dotransform, 0 for no transformation, 1 otherwise\r',str);
for iv=1:dotransform*nv,
  str=sprintf('%d',length(yfile{iv}));str=char(str,str30);str=str(1,:);
  fprintf(fidParam,'%s \\number of yfile points for variable %d\r',str,iv);
  for iy=1:length(yfile{iv})
    str=sprintf('%g %g',yfile{iv}(iy),cdfyfile{iv}(iy));str=char(str,str30);str=str(1,:);
    fprintf(fidParam,'%s \\value %d of yfile and cdfyfile for variable %d\r',str,iy,iv);
  end;
end;
st=fclose(fidParam);
if st==-1,
  warning('could not close parameter file properly');
end;


%%%%%%% Writing hard data file
fidHard=fopen(fileNames{2},'w');
fprintf(fidHard,'BME hard data\r');
fprintf(fidHard,'%d\r',d+isST+nv);
switch d,
  case 1, fprintf(fidHard,'s1\r');
  case 2, fprintf(fidHard,'s1\rs2\r');
  case 3, fprintf(fidHard,'s1\rs2\rs3\r');
end;
if isST
  fprintf(fidHard,'time\r');
end;
fprintf(fidHard,'primary variable\r');
for i=2:nv
   fprintf(fidHard,'variable %d\r',i);
end;
for i=1:nh,
  str=sprintf('%15g ',ch{1}(i,:)); 
  for j=1:ch{2}(i)-1,
    str=[str sprintf('%12g ',outOfRange)];
  end;
  str=[str sprintf('%12g ',zh(i))];
  for j=ch{2}(i)+1:nv
     str=[str sprintf('%12g ',outOfRange)];
  end;
  fprintf(fidHard,'%s\r',str);
end;
st=fclose(fidHard);
if st==-1,
  warning('could not close hard data file properly');
end;

%%%%%%% Writing soft probabilistic data file
fidProba=fopen(fileNames{3},'w');
fprintf(fidProba,'BME Probabilistic data\r');
fprintf(fidProba,'%d    NOTE: The actual number of columns varies from line to line\r',d+isST+5);
switch d,   
  case 1, fprintf(fidProba,'s1\r');
  case 2, fprintf(fidProba,'s1\rs2\r');
  case 3, fprintf(fidProba,'s1\rs2\rs3\r');     
end;
if isST
  fprintf(fidProba,'time\r');
end;  
if nv>1, fprintf(fidProba,'code for the variable (between 1 and %d)\r',nv);
else fprintf(fidProba,'code for the variable (always equal to one)\r'); end;
fprintf(fidProba,'Type of soft pdf (always equal to %d, corresponding to ', softpdftype);
switch softpdftype
  case 1, fprintf(fidProba,'histogram)\r');
  case 2, fprintf(fidProba,'linear)\r');
  case 3, fprintf(fidProba,'histogram on regular grid)\r');
  case 4, fprintf(fidProba,'linear on regular grid)\r');
end;
fprintf(fidProba,'number of limit values, nl\r');
switch softpdftype
  case {1,2}, fprintf(fidProba,'limits of intervals (nl values)\r');
  case {3,4}, fprintf(fidProba,'limits of intervals (minval,step,maxval)\r');
end;
switch softpdftype
  case {1,3}, fprintf(fidProba,'probability density (nl-1 values)\r');
  case {2,4}, fprintf(fidProba,'probability density (nl values)\r');
end;
for i=1:ns
  str=sprintf('%15g ',cs{1}(i,:)); 
  str=[str sprintf('%3d ',cs{2}(i))];
  str=[str sprintf('%2d ',softpdftype)];
  str=[str sprintf('%4d ',nl(i))];
  switch softpdftype
    case {1,2}, str=[str sprintf('%9g ',limi(i,1:nl(i)))];
    case {3,4}, str=[str sprintf('%9g ',limi(i,:))];
  end;
  str=[str '          '];  
  switch softpdftype
    case {1,3}, str=[str sprintf('%g ',probdens(i,1:nl(i)-1))];
    case {2,4}, str=[str sprintf('%g ',probdens(i,1:nl(i)))];
  end;
  fprintf(fidProba,'%s\r',str);
end;
st=fclose(fidProba);
if st==-1,
  warning('could not close Proba data file properly');
end;

%%%%%%% Writing estimation points file
fidEst=fopen(fileNames{4},'w');
fprintf(fidEst,'BME estimation points\r');
fprintf(fidEst,'%d\r',d+isST+1);
switch d,   
  case 1, fprintf(fidEst,'s1\r');
  case 2, fprintf(fidEst,'s1\rs2\r');
  case 3, fprintf(fidEst,'s1\rs2\rs3\r');     
end;
if isST
  fprintf(fidHard,'time\r');
end;
if nv>1, fprintf(fidEst,'code for the variable (between 1 and %d)\r',nv);
else fprintf(fidEst,'code for the variable (always equal to one)\r'); end;
for i=1:nk
  str=sprintf('%15g ',ck{1}(i,:)); 
  str=[str sprintf('%3d ',ck{2}(i))];
  fprintf(fidEst,'%s\r',str);
end;  
st=fclose(fidEst);
if st==-1,
  warning('could not close estimation points data file properly');
end;

