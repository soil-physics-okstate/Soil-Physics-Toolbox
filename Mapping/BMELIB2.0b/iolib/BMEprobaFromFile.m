function [BMEest,info]=BMEprobaFromFile(BMEparamfile);

% BMEprobaFromFile          - BME estimation using input from a file (Jan 1, 2001)
%
% Computation of the BME estimate of a variable at a set of estimation
% points, using both hard data and soft data of probabilistic type
% for an arbitrary number of different variables.
% This function reads the BME parameters from a parameter file, and
% the soft and hard data from data files. It calculates the BME estimate
% using the estimation method specified in the parameter file (i.e. either
% the mode or mean estimate, using transformation if requested), and
% it stores the estimate values in an output file.
%
% SYNTAX :
%
% [BMEest,info]=BMEprobaFromFile(BMEparamfile);
%
% INPUT :
%
% BMEparamfile string        filename of a BME parameter file specifying 
%                            the BME parameter and the BME data files.
%                            See ioparamsyntax for detailed explanations.
%
% OUTPUT :
%
% BMEest       nk by q       BME estimate, either BMEmode or BMEmean  
%                            depending on the BME estimation method BMEmethod specified 
%                            in the BME parameter file. 
%                            if BMEmethod='BMEprobaMode' or 'BMEprobaTMode'
%                                then BMEest is a nk by 1 vector of BMEmode estimates 
%                                at the nk est points (see BMEprobaMode)
%                            if BMEmethod='BMEprobaMoments'
%                                then BMEest is a nk by 3 matrix of BME moments 
%                                estimates at the nk est points (see BMEprobaMoment)
%                                The first column of BMEest is the BMEmean estimate
%                                The second column of BMEest (if calculated )is the 
%                                   standard deviation of the BME posterior pdf
%                                The third column of BMEest (if calculated )is the 
%                                   coefficient of skewness of the BME posterior pdf
% info         nk by q       information status obtained by the BME estimation.
%                            for BMEmode q=1, for BMEmean q=3. The value of info are  
%                            info=NaN : no computation, no hard or soft data
%                            info=0   : computation using BME with soft data 
%                            info=1   : dubious results, integration error above tolerance
%                                       when integrating over the soft pdf. Try to increase 
%                                       the value of options(3)
%                            info=3   : computation using kriging, no soft data
%                            info=4   : computation provides hard data value at estimation point
%                            info=10  : dubious results from integration routine. Try reducing
%                                       the number of soft data points ns, or option(3) or option(4)
%
% NOTES :
%
% The outputs from this program, i.e. BMEest and info, are also stored in the
% output file in GEO format. Once you use this function, look for a new file created
% with the output filename given in the parameter file. In this output file you
% will find the BME estimate output as well as the information parameter.


%%%%%% Error messages

if ~ischar(BMEparamfile)
  error('BMEparamfile must be a char array');
end

[BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile]=...
    readBMEproba(BMEparamfile);

switch(BMEmethod)
case 'BMEprobaMode'
  [BMEest,info]=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options);
  filetitle='BMEmode estimate calculated using BMEprobaMode';
    
case 'BMEprobaMoments'
  [BMEest,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  filetitle='BMEmean estimate calculated using BMEprobaMoments';
  
case 'BMEprobaTMode'
  [BMEest,info]=BMEprobaTMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
    covmodel,covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile);
  filetitle='BMEmode estimate calculated using BMEprobaTMode';
    
otherwise
  error(sprintf('BMEmethod=%s is not a valid choice for BMEprobaFromFile',BMEmethod));
  
end;

[isST,isSTsep,modelS,modelT]=isspacetime(covmodel);
if iscell(ck)
  val=ck{1};
else
  val=ck;
end
d=size(val,2)-isST;
switch d,
  case 1, valname={'x'};
  case 2, valname={'x','y'};
  case 3, valname={'x','y','z'};
end;
if isST
  valname{d+1}='time';
end;
switch(BMEmethod)
case {'BMEprobaMode','BMEprobaTMode'}
  val(:,d+isST+1)=BMEest;
  valname{d+isST+1}='BMEmode';
  val(:,d+isST+2)=info;
  valname{d+isST+2}='infoBMEmode';
case {'BMEprobaMoments','BMEprobaTMoments'}
  val(:,d+isST+1:d+isST+3)=BMEest;
  valname{d+isST+1}='BMEmean';
  valname{d+isST+2}='BMEvar';
  valname{d+isST+3}='BMEskewness';
  val(:,d+isST+4:d+isST+6)=info;
  valname{d+isST+4}='infoBMEmean';
  valname{d+isST+5}='infoBMEvar';
  valname{d+isST+6}='infoBMEskewness';
end;

colnum=1:size(val,2);
minval=min(min(val));
if minval>-1
  outofrangeval=-1;
elseif min>-1e30
  outofrangeval=-1e30;
elseif min>-1e300
  outofrangeval=-1e300;
elseif sum(isnan(BMEest(:)))>0
  warning('could not set a out of range value, setting it to 0');
  outofrangeval=0;
end

writeGeoEAS(val,valname,filetitle,outputfiles{1},colnum,outofrangeval);