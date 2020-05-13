function [yk,info]=BMEprobaTMode(ck,ch,cs,yh,ysoftpdftype,ynl,ylimi,yprobdens,covmodel,...
  covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile);

% BMEprobaTMode             - BME mode prediction with transform and probabilistic data (Jan 1, 2001)
%
% BME computation of the estimated mode at a set of estimation 
% points, using both hard data, soft data of probabilistic type
% and a data transformation toward a multivariate Gaussian 
% distribution for an arbitrary number of different variables.
% The transformation of the original Random Field Y to the
% transformed Random Field Z, Y->Z, is such that
% the maximum entropy distribution for Z is supposed to be multivariate
% normal. The solution corresponds to the mode of
% the posterior pdf on the original scale Y.
%
% SYNTAX :
%
% [yk,info]=BMEprobaTMode(ck,ch,cs,yh,ysoftpdftype,ynl,ylimi,yprobdens,covmodel,...
%   covparam,nhmax,nsmax,dmax,order,options,yfile,cdfyfile);
%
% INPUT :
%
% ck          nk by d     matrix of coordinates for the estimation locations.
%                         A line corresponds to the vector of coordinates at
%                         an estimation location, so the number of columns
%                         corresponds to the dimension of the space. There is
%                         no restriction on the dimension of the space.
% ch          nh by d     matrix of coordinates for the hard data locations,
%                         with the same convention as for ck.
% cs          ns by d     matrix of coordinates for the soft data locations,
%                         with the same convention as for ck.
% yh          nh by 1     vector of the values for the hard data for Y
%ysoftpdftype scalar      indicates the type of soft pdf for the variable Y   
%                         representing the probabilitic soft data at the 
%                         coordinates specified in cs.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
% ynl         ns by 1     vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. (see probasyntax for more explanations)
% ylimi       ns by l     matrix of interval limits for variable Y, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         soft data. (see probasyntax for more explanations)
% yprobdens   ns by p     matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
% covmodel   string       string that contains the name of the covariance model 
%                         for the Gaussian-transformed Z data.
%                         (see the modelslib directory for a list of available models)
%                         Variogram models are not available for this function.
%                         See modelsyntax for more explanations on covariance models.
% covparam   1 by k       vector of values for the parameters of covmodel, according
%                         to the convention for the corresponding covariance model.
% nhmax      scalar       maximum number of hard data values that are considered
%                         for the estimation at the locations specified in ck.
% nsmax      scalar       maximum number of soft data values that are considered for
%                         the estimation at the locations specified in ck. As the
%                         computation time is exponentially increasing with nsmax,
%                         it is not advised to use more than few soft data locations.
%                         In any case, nsmax should be lower than 20 in order to
%                         avoid numerical computation problems.
% dmax       scalar       maximum distance between an estimation location and
%                         existing hard/soft data locations. All hard/soft data
%                         locations separated by a distance smaller than dmax from an
%                         estimation location will be included in the estimation process
%                         for that location, whereas other data locations are neglected.
% order      scalar       order of the polynomial drift along the spatial axes at the
%                         estimation locations. For the zero-mean case, NaN (Not-a-Number)
%                         is used.
%                         order=Nan for a zero mean
%                              =0 for a constant mean
%                              =1 for a constant+linear mean along each axis
%                              =2 for a constant+linear+quadratic mean along each axis, etc.
% options    1 by 1 or 14 vector of optional parameters that can be used if default
%                         values are not satisfactory (otherwise this vector can simply be
%                         omitted from the input list of variables), where :
%                         options(1), options(2) and options(14) are values used by the
%                         fminbnd.m MATLAB optimization routine that finds the mode of the
%                         probability distribution function (default values are the same
%                         as for fminbnd.m),
%                         options(3) specifies the maximum number of evaluation that can
%                         be done by the FORTRAN77 subroutines for the integrals (default
%                         value is 50 000 ; this value should be increased if a warning
%                         message appears on the screen during the computation),
%                         options(4) specifies the maximum admissible relative error on the
%                         estimation of these integrals (default value is 1e-4).  
% yfile      k by 1       vector of values for Y sorted in ascending order
% cdfyfile   k by 1       vector of values of the cumulative distribution for Y
%                         The values of cdfyfile must be strickly increasing, 
%                         between 0 and 1, and cdfyfile(1)<0.001 and cdfyfile>0.999 so that 
%                         yfile covers most of the range of Y values.
%
% OUTPUT :
%
% yk         nk by 1      vector of mode estimates for Y using kriging or BME with data transformation
% info       nk by 1      vector of information about the computation in the local 
%                         neighbourhood around the estimation point :
%                         info=NaN : no computation at all, no hard or soft data
%                         info=0   : computation using BME using soft data 
%                         info=1   : dubious results, integration error above tolerance
%                                    when integrating over the soft pdf. Try to increase 
%                                    the value of options(3)
%                         info=2   : dubious results, did not converge in the maximum 
%                                    number of iterations corresponding to options(14)
%                                    Try to increase options(14)
%                         info=3   : computation using kriging, no soft data
%                         info=4   : computation provides hard data value at estimation point 
%                         info=10  : dubious results from integration routine. Try reducing
%                                    the number of soft data points ns, or option(3) or option(4)
%
% NOTE :
%
% All the specific conventions for specifying nested models, multivariate 
% or space-time cases are the same as for BMEprobaMode.m.
% In the multivariate case (more than one field variables), yfile and cdfyfile
% must be cell arrays with nv vectors, where each vector correspond to one
% of the nv field variables.

%%%%%% Error messages

[ch,yh,cs,ynl,ylimi,yprobdens]=fixempty(ck,ch,yh,cs,ynl,ylimi,yprobdens);
BMEprobaCheckArgs(ck,ch,cs,yh,ysoftpdftype,ynl,ylimi,yprobdens,...
  covmodel,covparam,nhmax,nsmax,dmax);
if nargin<17          % If yfile is not specified then there will be no transformation
  yfile=[];           %    in which case set yfile=[]
  cdfyfile=[];
else
  FyTransformCheckArgs(yfile,cdfyfile,ck,ch,cs,yh,...
    ysoftpdftype,ynl,ylimi,yprobdens);
end;

%%%%%% Initialize the parameters

global INFOINTEG      % information message from the Fortran subroutine
if nargin<15,         % initialize options with default values
  options=BMEoptions;
  options1=0;
else
  options1=options(1);
  options(1)=0;
end;

if length(options)<2 | ~options(2)
  options(2)=1e-4;
end  
TolX=options(2);
if length(options)<14 | ~options(14)
  options(14)=500;
end  
MaxFunEvals=options(14);
fminbndoptions = optimset('Display','off','TolX',TolX,'MaxFunEvals',MaxFunEvals);


noindex=~iscell(ck);       % test if there is an index for the variables
if noindex==1,
  nk=size(ck,1);           % nk is the number of estimation points
  nh=size(ch,1);           % nh is the number of hard data
  ns=size(cs,1);           % ns is the number of soft data
else
  nk=size(ck{1},1);
  nh=size(ch{1},1);
  ns=size(cs{1},1);
end;

if options1==1,
  num2strnk=num2str(nk);
end;

zk=zeros(nk,1)*NaN;
info=zeros(nk,1)*NaN;

%%%%%% Transform the y values to Gaussian-distributed values

yk=zk;                        % initialize yk to zeros
dotransform=~isempty(yfile);  % do the y-z transformation if yfile is not empty

if dotransform
  if noindex==1,
    indexk=ones(nk,1);
    indexh=1;
    indexs=1;
    zh=other2gauss(yh,yfile,cdfyfile);
    [zsoftpdftype,znl,zlimi,zprobdens]=probaother2gauss(ysoftpdftype,ynl,ylimi,yprobdens,yfile,cdfyfile);
    yfile={yfile};
    cdfyfile={cdfyfile};
  else
    indexk=ck{2};
    indexh=ch{2};
    indexs=cs{2};
    zh=other2gauss({yh,indexh},yfile,cdfyfile);
    [zsoftpdftype,znl,zlimi,zprobdens]=probaother2gauss(ysoftpdftype,{ynl,indexs},ylimi,yprobdens,yfile,cdfyfile);
  end;
else
  zh=yh;
  zsoftpdftype=ysoftpdftype;
  znl=ynl;
  zlimi=ylimi;
  zprobdens=yprobdens;
end;

%%%%%% Compute the derivative of the transforms

if dotransform
  maxindex=max([indexk;indexh;indexs]);
  for i=1:maxindex,
    dgfile{i}=transformderiv(yfile{i},cdfyfile{i});
    zfile{i}=other2gauss({yfile{i},i},yfile,cdfyfile);
  end;
end;

%%%%%% Main loop starts here

for i=1:nk,
  
  if noindex==1,
    ck0=ck(i,:);
  else
    ck0={ck{1}(i,:),ck{2}(i)};
  end;

  %%%%%% Select the local neighbourhood for all variables

  [chlocal,zhlocal,dh,sumnhlocal]=neighbours(ck0,ch,zh,nhmax,dmax);
  [cslocal,znllocal,zlimilocal,zprobdenslocal,ds,sumnslocal]=...
    probaneighbours(ck0,cs,zsoftpdftype,znl,zlimi,zprobdens,nsmax,dmax);
  
  %%%%%% Test if local neighbourhood is empty
  
  iscomputed=0;
  if isempty(zhlocal) & isempty(znllocal)
    zk(i)=NaN;
    info(i)=NaN;
    iscomputed=1;
  end;

  %%%%%% Test if there is a hard data at estimation point

  [index]=findpairs(ck0,chlocal);
  if iscomputed==0 & ~isempty(index),            % if there is a hard data at estimation point,
    zk(i)=zhlocal(index(2));                     % return the hard data value as the estimate
    info(i)=4;                                   % set the value of info to 4 
    iscomputed=1;                                % specify that the pdf has been computed
  end;

  %%%%%% Do the BME estimation if iscomputed==0
  
  if iscomputed==0
                 %%%%%% Test if there is a soft data at estimation point

    isduplicate=0;
    [index]=findpairs(ck0,cslocal);                % test if there is a soft data at estimation point
    if ~isempty(index),                            % and split the soft data into two sets
      [csest,cslocal,znlest,zlimiest,zprobdensest,znllocal,zlimilocal,zprobdenslocal]=...
        probasplit(cslocal,zsoftpdftype,znllocal,zlimilocal,zprobdenslocal,index(2));
      sumnslocal=sumnslocal-1;
      isduplicate=1;                               % specify that there is a soft data at estimation point
    end;

                 %%%%%%% Calculate the covariance matrices
    Kk=coord2Kinterface(ck0,ck0,covmodel,covparam);           % compute variance for estimation point
    Kh=coord2Kinterface(chlocal,chlocal,covmodel,covparam);   % compute covariance matrix for hard data
    Ks=coord2Kinterface(cslocal,cslocal,covmodel,covparam);   % compute covariance matrix for soft data
    Kk_h=coord2Kinterface(ck0,chlocal,covmodel,covparam);     % compute cross-covariance vector estimation/hard data
    Kk_s=coord2Kinterface(ck0,cslocal,covmodel,covparam);     % compute cross-covariance vector estimation/soft data
    Ks_h=coord2Kinterface(cslocal,chlocal,covmodel,covparam); % compute cross-covariance matrix soft/hard data
    Kkh=[[Kk,Kk_h];[Kk_h',Kh]];                      % build estimation+hard data covariance matrix
    Ks_kh=[Kk_s',Ks_h];                              % build cross-covariance matrix soft/estimation+hard data
    invKkh=inv(Kkh);
    BsIFkh=Ks_kh*invKkh;
    KsIFkh=Ks-BsIFkh*Ks_kh';                         % compute the conditional covariance matrix soft given k and h
    invKh=inv(Kh);
    BkIFh=Kk_h*invKh;
    KkIFh=Kk-BkIFh*Kk_h';
    invKkIFh=inv(KkIFh);

                 %%%%%%% Calculate the mean trend and remove it from the data
    [zmslocal,zvslocal]=proba2stat(zsoftpdftype,znllocal,zlimilocal,zprobdenslocal);
    [zmkest,zmhest,zmsest]=localmeanBME(ck0,chlocal,cslocal,zhlocal,zmslocal,zvslocal,Kh,Ks_h,Ks,order);
    if sumnhlocal>0  
      xhlocal=zhlocal-zmhest;                       % substract the mean from the hard data
    else
      xhlocal=zeros(0,1);
    end;
    xsoftpdftype=zsoftpdftype;
    if sumnslocal>0                                 % substract the mean from the soft data
      xnllocal=znllocal;
      [xlimilocal]=probaoffset(zsoftpdftype,znllocal,zlimilocal,-zmsest);
      xprobdenslocal=zprobdenslocal;
    else
      xnllocal=zeros(0,1);
      xlimilocal=zeros(0,1);
      xprobdenslocal=zeros(0,1);
    end;
    if isduplicate                                  % substract mean from soft data at est pt
      xnlest=znlest;
      [xlimiest]=probaoffset(zsoftpdftype,znlest,zlimiest,-zmkest);
      xprobdensest=zprobdensest;
    end;
    
                  %%%%%%% Transform the data if needed
    if dotransform                                % if transforming the data 
      zfilei=zfile{indexk(i)};                    %   select the proper transformation file for the 
      dgfilei=dgfile{indexk(i)};                  %   estimation point
    else                                          % else set the files to []
      zfilei=[];
      dgfilei=[];
    end
    
                  %%%%%%% Find the mode of the BME posterior pdf
    if isduplicate==1,                            %   if there is a soft data at the estimation point
                                                  %      use the minimum and maximum value for the soft data 
                                                  %      at the estimation point as min/max for fminbnd
      [xkmin,xkmax]=proba2interval(xsoftpdftype,xnlest,xlimiest);
      [xk(i),FVAL,EXITFLAG,OUTPUT]=fminbnd(@fminBMEprobaDupliTMode,xkmin,xkmax,...
        fminbndoptions,xhlocal,xsoftpdftype,xnllocal,xlimilocal,xprobdenslocal,...
        BkIFh,invKkIFh,BsIFkh,KsIFkh,options,xnlest,xlimiest,xprobdensest,...
        zfilei,dgfilei,zmkest);
      options(10)=OUTPUT.funcCount;
      zk(i)=xk(i)+zmkest;                         %      and add the mean to the mode estimate
    else                                          %   if there is no soft data at the estimation point
      xkmin=0;                                    %      use the minimum and maximum of hard data, lower bound of 
      xkmax=0;                                    %      soft data and zero value as min/max for fminbnd       
      if sumnhlocal>0
        xkmin=min([xkmin;xhlocal]);
        xkmax=max([xkmax;xhlocal]);
      end
      if sumnslocal>0
        [xalocal,xblocal]=proba2interval(xsoftpdftype,xnllocal,xlimilocal);
        xkmin=min([xkmin;xalocal]);
        xkmax=max([xkmax;xblocal]);
      end
      if dotransform
        xkmin=min(zfilei)-zmkest;
        xkmax=max(zfilei)-zmkest;
      end
      [xk(i),FVAL,EXITFLAG,OUTPUT]=fminbnd(@fminBMEprobaTMode,xkmin,xkmax,...
        fminbndoptions,xhlocal,xsoftpdftype,xnllocal,xlimilocal,xprobdenslocal,...
        BkIFh,invKkIFh,BsIFkh,KsIFkh,options,zfilei,dgfilei,zmkest);
      options(10)=OUTPUT.funcCount;
      zk(i)=xk(i)+zmkest;                            %      and add the mean to the mode estimate
    end;

                   %%%%%% Select the appropriate info value
    if options(10)>=options(14),
       info(i)=2;
    else
      if INFOINTEG==0,info(i)=0;end;
      if INFOINTEG==1,info(i)=1;end
      if INFOINTEG>1,info(i)=10;end
    end;
  end;
  
  %%%%%%% Transform back the zk estimate to a yk estimate
  
  if dotransform
    yk(i)=interp1(zfile{indexk(i)},yfile{indexk(i)},zk(i),'linear');
  else
    yk(i)=zk(i);
  end;
  
  %%%%%% Display the point number if options(1)==1
  
  if options1==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;