function [BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile]=...
    readBMEproba(paramfile);

% readBMEproba              - Reads data files for the BMEproba functions (Jan 1, 2001)
%
% Reads from files all the parameters and data needed for a BME estimation
% with hard data and soft data of probabilistic type (i.e. estimation using 
% either the BMEprobaMode, BMEprobaTMode, BMEprobaMoment, or BMEprobaTMoment 
% functions).
%
% SYNTAX :
%
% [BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
%   covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile]=...
%   readBMEproba(paramfile);
%
% INPUT :
%
% BMEparamfile   string      filename of a BME parameter file specifying 
%                            the BME parameter and the BME data files.
%                            See ioparamsyntax for detailed explanations.
%
% OUTPUT :
%
% BMEmethod      string      string indicating the BME method to use. The string must be one 
%                            of the following: 'BMEprobaMode', 'BMEprobaTMode',
%                            'BMEprobaMoment', or 'BMEprobaTMoment'. Use the help for
%                            any of these methods to learn more about it.
% fileNames      cell        2 by 1 cell array of strings, with the names of the BME result 
%                            file, the BME error messages and warnings file.
% ck,ch,etc.                 same arguments as for the BME estimation function specified by
%                            BMEmethod. See the help for that function to get a complete 
%                            description of the input arguments.

%%%%%%%%% Checking data
if ~ischar(paramfile) error('paramfile must be a character string'); end;

%%%%%%%%%% read the parameter file
fidParam=fopen(paramfile,'r');
if fidParam==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',paramfile));
end;
fgetl(fidParam);
fgetl(fidParam);
fgetl(fidParam);
str=fgetl(fidParam);      % Reading the BME method
BMEmethod=sscanf(str,'%s',1);
str=fgetl(fidParam);      % Reading the filename for hard data
hardfile=sscanf(str,'%s',1);
str=fgetl(fidParam);      % Reading the number of variables, nv
nv=sscanf(str,'%d',1);
if length(nv)~=1 | nv<=0        
  fclose('all');
  error('nv must be a value greater than zero'); 
end;
str=fgetl(fidParam);      % Reading the column numbers
colnum=sscanf(str,'%d')';
if length(colnum)~=nv+4, 
   fclose('all');
   error(sprintf('line 7 of file %s should have %d values',paramfile,4+nv)); 
end;
if sum(colnum(1:3)>0)==0
   fclose('all');
   error(sprintf('line 7 of file %s, at least one of the first 3 values must be non zero',paramfile)); 
end;
if sum(colnum(4+1:4+nv)>0)~=nv
   fclose('all');
   error(sprintf('line 7 of file %s, the column number for variable(s) must be non zero',paramfile)); 
end;
d=sum(colnum(1:3)>0);
isST=(colnum(4)>0);
str=fgetl(fidParam);      % Reading the filename for soft data
softfile=sscanf(str,'%s',1);
str=fgetl(fidParam);      % Reading the data trimming limits
trimlim=sscanf(str,'%g');
if length(trimlim)~=2, 
   fclose('all');
   error(sprintf('line 9 of file %s should have 2 values',paramfile)); 
end;
estfile=fscanf(fidParam,'%s',1); % Reading the filename for the estimation points
fgetl(fidParam);
resfile=fscanf(fidParam,'%s',1); % Reading the filename for BME output results
outputfiles{1}=resfile;
fgetl(fidParam);
debugfile=fscanf(fidParam,'%s',1); % Reading the filename for BME output messages
outputfiles{2}=debugfile;
fgetl(fidParam);
str=fgetl(fidParam);      %Reading the covariance models
linenumber=13;
nm=sscanf(str,'%d');
if length(nm)~=1, 
   fclose('all');
   error(sprintf('line %d of file %s should have 2 values',linenumber,paramfile)); 
end;
for im=1:nm
  str=fgetl(fidParam);
  linenumber=linenumber+1;;
  covmodel{im}=sscanf(str,'%s',1);
  for iv=1:nv
    str=fgetl(fidParam);
    linenumber=linenumber+1;
    tempval=sscanf(str,'%g');
    if length(tempval)~=nv
      fclose('all');
      error(sprintf('line %d of file %s should have %d values',linenumber,paramfile,nv)); 
    end;
    covparam{im}{1}(iv,1:nv)=tempval;
  end;
  str=fgetl(fidParam);
  linenumber=linenumber+1;
  tempval=sscanf(str,'%g');
  if length(tempval)==0
    covparam{im}{2}=[];
  else
    covparam{im}{2}=tempval(:)';
  end;
end;
str=fgetl(fidParam);      %Reading the maximum number of hard data
linenumber=linenumber+1;
nhmax=sscanf(str,'%d');
if length(nhmax)~=nv, 
   fclose('all');
   error(sprintf('line %d of file %s should have %d values',linenumber,paramfile,nv)); 
end;
str=fgetl(fidParam);      %Reading the maximum number of soft data
linenumber=linenumber+1;
nsmax=sscanf(str,'%d');
if length(nsmax)~=nv, 
   fclose('all');
   error(sprintf('line %d of file %s should have %d values',linenumber,paramfile,nv)); 
end;
str=fgetl(fidParam);      %Reading the maximum search radius
linenumber=linenumber+1;
dmax=sscanf(str,'%g');
if ~isST & length(dmax)~=1
  error(sprintf('line %d of file %s should have 1 value (spatial case)',linenumber,paramfile)); 
elseif isST & length(dmax)~=3 
  error(sprintf('line %d of file %s should have 3 values (space/time case)',linenumber,paramfile)); 
end;
str=fgetl(fidParam);      %Reading the order of the mean trend
linenumber=linenumber+1;
tempval=sscanf(str,'%d');
if length(tempval)~=1,
  fclose('all');
  error(sprintf('line %d of file %s should have 1 value',linenumber,paramfile));
end;
if tempval==-1,
  order=NaN;
else
  order=round(tempval);
end;
str=fgetl(fidParam);      %Reading options(3) and options(4)
linenumber=linenumber+1;
tempval=sscanf(str,'%g');
if length(tempval)~=2,
  fclose('all');
  error(sprintf('line %d of file %s should have 2 values',linenumber,paramfile));
end;
options(3)=round(tempval(1));
options(4)=tempval(2);
str=fgetl(fidParam);      %Reading options(2) and options(14)
linenumber=linenumber+1;
tempval=sscanf(str,'%g');
if length(tempval)~=2, 
   fclose('all');
   error(sprintf('line %d of file %s should have 2 values',linenumber,paramfile)); 
end;
options(2)=tempval(1);
options(14)=round(tempval(2));
str=fgetl(fidParam);      %Reading options(6) and options(7)
linenumber=linenumber+1;
tempval=sscanf(str,'%g');
if length(tempval)~=2, 
   fclose('all');
   error(sprintf('line %d of file %s should have 2 values',linenumber,paramfile)); 
end;
options(6)=tempval(1);
options(7)=tempval(2);
str=fgetl(fidParam);      %Reading options(8)
linenumber=linenumber+1;
tempval=sscanf(str,'%g');
if length(tempval)~=1, 
   fclose('all');
   error(sprintf('line %d of file %s should have 2 values',linenumber,paramfile)); 
end;
options(8)=round(tempval);
str=fgetl(fidParam);      %Reading options(20) up to options(29)
linenumber=linenumber+1;
tempval=sscanf(str,'%g');
if length(tempval)>10, 
   fclose('all');
   error(sprintf('line %d of file %s should have less than 10 values',linenumber,paramfile)); 
end;
if min(tempval)<0 | max(tempval)>1 
   fclose('all');
   error(sprintf('line %d of file %s should have values between 0 and 1',linenumber,paramfile)); 
end;
options(20:20+length(tempval)-1)=tempval(:)';
str=fgetl(fidParam);        %Reading dotransform
linenumber=linenumber+1;
tempval=sscanf(str,'%d');
if length(tempval)~=1, 
   fclose('all');
   error(sprintf('line %d of file %s should have 1 number',linenumber,paramfile)); 
end;
dotransform=tempval;
if dotransform==0
  yfile=[];
  cdfyfile=[];
end
for iv=1:dotransform*nv 
  str=fgetl(fidParam);       %Reading ny=length(yfile{iv})
  linenumber=linenumber+1;
  tempval=sscanf(str,'%d');
  if length(tempval)~=1, 
    fclose('all');
    error(sprintf('line %d of file %s should have 1 number',linenumber,paramfile)); 
  end;
  ny=tempval;
  for iy=1:ny
    str=fgetl(fidParam);     %Reading iy-th value of yfile,cdfyfile
    linenumber=linenumber+1;
    if str==-1, 
      fclose('all');
      error(sprintf('line %d of file %s, problems reading that line',linenumber,paramfile)); 
    end;
    tempval=sscanf(str,'%g');
    if length(tempval)~=2, 
      fclose('all');
      error(sprintf('line %d of file %s should have 2 values',linenumber,paramfile)); 
    end;
    yfile{iv}(iy,1)=tempval(1);
    cdfyfile{iv}(iy,1)=tempval(2);
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read hard data file
fidHard=fopen(hardfile,'r');
if fidHard==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',hardfile));
end;
fgetl(fidHard);
str=fgetl(fidHard);
nc=sscanf(str,'%d',1);
if nc>max(colnum)
  error(sprintf('%s requires %d columns of hard data, but only %d are found in %s',paramfile,max(colnum),nv,hardfile));   
end;
for ic=1:nc,
   fgetl(fidHard);
end;
idxcoord(1,1:d)=colnum(find(colnum(1:3)));
if isST,
   idxcoord=[idxcoord colnum(4)];
end;
idxdata=colnum(4+1:4+nv);
ih=0;
linenumber=nc+2;
tempval=[];
while ~feof(fidHard)
  str=fgetl(fidHard);
  tempval=sscanf(str,'%g'); 
  linenumber=linenumber+1;
  if length(tempval)~=nc, 
    fclose('all');
    error(sprintf('line %d of file %s should have %d values',linenumber,hardfile,nc)); 
  end;
  coord=tempval(idxcoord);
  data=tempval(idxdata);
  for iv=1:nv
    if data(iv)>trimlim(1) & data(iv)<trimlim(2),
      ih=ih+1;
      ch{1}(ih,1:d+isST)=coord;
      ch{2}(ih,1)=iv;
      zh(ih,1)=data(iv);
    end;
  end;
end;
if ih==0,
   ch{1}=zeros(0,d+isST);
   ch{2}=zeros(0,1);
   zh=zeros(0,1);   
end;
st=fclose(fidHard);
if st==-1,
  warning(sprintf('could not close %s properly',hardfile));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read soft probabilistic data file
fidSoft=fopen(softfile,'r');
if fidSoft==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',softfile));
end;
fgetl(fidSoft);
str=fgetl(fidSoft);
nc=sscanf(str,'%d',1);
for ic=1:nc,
   fgetl(fidSoft);
end;
idxcoord=1:d+isST;
idxcode=d+isST+1;
idxspt=d+isST+2;
idxnl=d+isST+3;
linenumber=nc+2;
is=0;
tempval=[];
while ~feof(fidSoft)
  str=fgetl(fidSoft);
  tempval=sscanf(str,'%g');
  tempval=tempval(:)';
  linenumber=linenumber+1;
  is=is+1;  
  if length(tempval)<idxnl, 
    fclose('all');
    error(sprintf('line %d of file %s should have at least %d values',linenumber,softfile,idxnl)); 
  end;
  cs{1}(is,1:d+isST)=tempval(idxcoord);
  cs{2}(is,1)=tempval(idxcode);
  if is==1
     softpdftype=tempval(idxspt);
  elseif softpdftype~=tempval(idxspt),
     fclose('all');
     error(sprintf('On line %d of file %s, the softpdftype should be the same as in previous lines',linenumber,softfile)); 
  end;
  nl(is,1)=tempval(idxnl);  
  switch softpdftype
    case {1}, nlimi=nl(is);nprobdens=nl(is)-1;
    case {2}, nlimi=nl(is);nprobdens=nl(is);
    case {3}, nlimi=3;nprobdens=nl(is)-1;
    case {4}, nlimi=3;nprobdens=nl(is);
    otherwise,
      fclose('all');
      error(sprintf('line %d of file %s softpdftype must equal to 1, 2, 3 or 4',linenumber,softfile)); 
  end;
  nval=d+isST+3+nlimi+nprobdens;  
  if length(tempval)<nval, 
    fclose('all');
    error(sprintf('line %d of file %s should have at least %d values',linenumber,softfile,nval)); 
  end;
  limi(is,1:nlimi)=tempval(d+isST+3+1:d+isST+3+nlimi);
  if softpdftype==3 | softpdftype==4
    limivec=limi(is,1):limi(is,2):limi(is,3);
    if nl(is)~=length(limivec) 
      fclose('all');
      error(sprintf('on line %d of file %s nl=%d, but limi would requires that nl=%d',...
        linenumber,softfile,nl(is),length(limivec))); 
    end;
  end
  probdens(is,1:nprobdens)=tempval(d+isST+3+nlimi+1:d+isST+3+nlimi+nprobdens);
end;
if is==0,
   cs{1}=zeros(0,d+isST);
   cs{2}=zeros(0,1);
   softpdftype=1;
   nl=zeros(0,1);
   limi=zeros(0,1);
   probdens=zeros(0,1);
end;
st=fclose(fidSoft);
if st==-1,
  warning(sprintf('could not close %s properly',softfile));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read coordinates of estimation points
fidEst=fopen(estfile,'r');
if fidEst==-1,
   fclose('all');
   error(sprintf('Problem opening file %s',estfile));
end;
fgetl(fidEst);
str=fgetl(fidEst);
nc=sscanf(str,'%d',1);
for ic=1:nc,
   fgetl(fidEst);
end;
idxcoord=1:d+isST;
idxcode=d+isST+1;
linenumber=nc+2;
ik=0;
tempval=[];
while ~feof(fidEst)
  str=fgetl(fidEst);
  tempval=sscanf(str,'%g');
  tempval=tempval(:)';
  linenumber=linenumber+1;
  ik=ik+1;  
  if length(tempval)~=d+isST+1, 
    fclose('all');
    error(sprintf('line %d of file %s should have %d values',linenumber,estfile,d+isST+1)); 
  end;
  ck{1}(ik,1:d+isST)=tempval(idxcoord);
  ck{2}(ik,1)=tempval(idxcode);
end;
if ik==0,
   ck{1}=zeros(0,d+isST);
   ck{2}=zeros(0,1);
end;
st=fclose(fidEst);
if st==-1,
  warning(sprintf('could not close %s properly',estfile));
end;

[ch,zh,cs,nl,limi,probdens]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
  covmodel,covparam,nhmax,nsmax,dmax,order);

st=fclose(fidParam);
