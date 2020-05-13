function [zk,info]=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,...
  covparam,nhmax,nsmax,dmax,order,options);

% BMEprobaMode              - BME mode prediction with probabilistic data (Jan 1, 2001)
%
% BME computation of the mode of the posterior pdf at a set
% of estimation points, using both hard data and soft probabilistic data
% that can reasonably be assumed as Gaussian distributed. The
% function is intendend to be as general as possible, covering
% various situations, like multi field variables, nested models, 
% space-time estimations, non-stationarity of the mean, etc. 
% Depending on the case, specific format are needed for the input 
% variables. This function is the most general function of space/time
% estimation in BMELIB as it has the potential to process both hard 
% data and soft probability data when both are available, but it 
% reduces to the kriging estimate when only hard data is used.
%
% SYNTAX :
%
% [zk,info]=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,...
%   covparam,nhmax,nsmax,dmax,order,options);
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
% zh          nh by 1     vector of values for the hard data at the coordinates
%                         specified in ch.
% softpdftype scalar      indicates the type of soft pdf representing the  
%                         probabilitic soft data at the coordinates specified in cs.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
% nl          ns by 1     vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. (see probasyntax for more explanations)
% limi        ns by l     matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         soft data. (see probasyntax for more explanations)
% probdens    ns by p     matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
% covmodel   string       string that contains the name of the covariance model
%                         that is used for the estimation (see the modelslib directory)
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
%                         estimation of these integrals (default value is 1e-4). The values
%                         for options(5) to options(13) are not used.
%
% OUTPUT :
%
% zk         nk by 1      vector of estimated values at the estimation locations. A value
%                         coded as NaN means that no estimation has been performed at that
%                         location due to the lack of available data.
% info       nk by 1      vector for information about the computation of estimated values.
%                         Different possible values are :
%                         info=NaN if there is no computation at all (no hard and soft data
%                         are available around or at the estimation location),
%                         info=0 when computation is made using BME with at least 1 soft data,
%                         info=1 when the estimation is dubious due to an integration error
%                         above the tolerance specified in options(4),
%                         info=2 when the estimation is dubious as the result did not converge
%                         within the maximum number of iterations specified in options(14),
%                         info=3 when the computation is made using kriging with no soft data,
%                         info=4 when there is a hard data value at the estimation location,
%                         info=10 when dubious results from integration routine. Try reducing
%                               the number of soft data points ns, or option(3) or option(4)
%
% NOTE :
%
% All the specific conventions for specifying nested models, multivariate 
% or space-time cases are the same as for BMEprobaMoments.m.

%%%%%% Error messages

[ch,zh,cs,nl,limi,probdens]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
  covmodel,covparam,nhmax,nsmax,dmax);

%%%%%% Initialize the parameters

global INFOINTEG      % information message from the Fortran subroutine
if nargin<15,         % initialize options with default values
  options=BMEoptions;
  options1=1;
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

%%%%%% Main loop starts here

for i=1:nk,
  
  if noindex==1,
     ck0=ck(i,:);
  else
     ck0={ck{1}(i,:),ck{2}(i)};
  end;
   
  %%%%%% Select the local neighbourhood for all variables
    
  [chlocal,zhlocal,dh,sumnhlocal]=neighbours(ck0,ch,zh,nhmax,dmax);
  [cslocal,nllocal,limilocal,probdenslocal,ds,sumnslocal]=...
    probaneighbours(ck0,cs,softpdftype,nl,limi,probdens,nsmax,dmax);
  
  %%%%%% Test if local neighbourhood is empty
  
  iscomputed=0;
  if isempty(zhlocal) & isempty(nllocal)
    zk(i)=NaN;
    info(i)=NaN;
    iscomputed=1;
  end;

  %%%%%% Test if there is a hard data at estimation point
  
  [index]=findpairs(ck0,chlocal);                
  if iscomputed==0 & ~isempty(index),            % if there is a hard data at estimation point,
    zk(i)=zhlocal(index(2));                     % return the hard data value as the estimate
    info(i)=4;                                   % set the value of info to 4 
    iscomputed=1;                                % specify that the value has been computed
  end;
  
  %%%%%% Do the BME estimation if iscomputed==0
  
  if iscomputed==0
                 %%%%%% Test if there is a soft data at estimation point
                 
    isduplicate=0;
    [index]=findpairs(ck0,cslocal);                % test if there is a soft data at estimation point           
    if ~isempty(index),                            % and split the soft data into two sets
      [cest,cslocal,nlest,limiest,probdensest,nllocal,limilocal,probdenslocal]=...
        probasplit(cslocal,softpdftype,nllocal,limilocal,probdenslocal,index(2));
      sumnslocal=sumnslocal-1;
      isduplicate=1;                               % specify that there is a soft data at estimation point  
    end;
    
                 %%%%%% Returns the cokriging estimate when sumnslocal==0 and isduplicate==0
  
    % This next block can be turned on or off.  
    % Turning the block ON will will use universal kriging when their is no soft data
    % Turning the block OFF will treat universal kriging as a special case of BME
    useUniversalKriging=0;
    if useUniversalKriging==1
      if (sumnhlocal>0)&(sumnslocal==0)&(iscomputed==0)&(isduplicate==0),% test if there are no soft data
        K=coord2K(chlocal,chlocal,covmodel,covparam);   % built the left-hand side matrix
        k=coord2K(chlocal,ck0,covmodel,covparam);       % built the right-hand side vector
        [X,x]=krigconstr(chlocal,ck0,order);            % built the constraints matrix and vector
        nx=size(X,2);
        Kadd=[[K,X];[X',zeros(nx)]];
        kadd=[k;x];
        lam=Kadd\kadd;                                  % compute the kriging weights lam
        lam=lam(1:sumnhlocal);                          % remove the Lagrangians from the solution
        lamt=lam';
        zk(i)=lamt*zhlocal;                             % compute the kriging estimate zk(i)
        info(i)=3;
        iscomputed=1;  
      end;
    end;
    
  end;
  
                 %%%%%% Returns the BME estimate if we still have iscomputed==0,
  if iscomputed==0
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
    [mslocal,vslocal]=proba2stat(softpdftype,nllocal,limilocal,probdenslocal);
    [mkest,mhest,msest]=localmeanBME(ck0,chlocal,cslocal,zhlocal,mslocal,vslocal,Kh,Ks_h,Ks,order);
    if sumnhlocal>0  
      zhlocal=zhlocal-mhest;                          % substract the mean from the hard data
    end;  
    if sumnslocal>0                                   % substract the mean from the soft data
      [limilocal]=probaoffset(softpdftype,nllocal,limilocal,-msest);
    end;
    if isduplicate                                  % substract mean from soft data at est pt
      [limiest]=probaoffset(softpdftype,nlest,limiest,-mkest);
      [aest,best]=proba2interval(softpdftype,nlest,limiest);
    end;
    
                  %%%%%%% Find the mode of the BME posterior pdf
    if isduplicate==1,                              %   if there is a soft data at the estimation point
      zkmin=aest;                                   %      use the minimum and maximum value for the soft data 
      zkmax=best;                                   %      at the estimation point as min/max for fminbnd
      [zk(i),FVAL,EXITFLAG,OUTPUT]=fminbnd(@fminBMEprobaDupliMode,zkmin,zkmax,...
        fminbndoptions,zhlocal,softpdftype,nllocal,limilocal,probdenslocal,...
        BkIFh,invKkIFh,BsIFkh,KsIFkh,options,nlest,limiest,probdensest);
      options(10)=OUTPUT.funcCount;
      zk(i)=zk(i)+mkest;                            %      and add the mean to the mode estimate
    else                                            %   if there is no soft data at the estimation point
      [alocal,blocal]=proba2interval(softpdftype,nllocal,limilocal);
      zkmin=min([zh;alocal;0]);                     %      use the minimum and maximum of hard data, lower bound of 
      zkmax=max([zh;blocal;0]);                     %      soft data and zero value as min/max for fmin       
      [zk(i),FVAL,EXITFLAG,OUTPUT]=fminbnd(@fminBMEprobaMode,zkmin,zkmax,...
        fminbndoptions,zhlocal,softpdftype,nllocal,limilocal,probdenslocal,BkIFh,...
        invKkIFh,BsIFkh,KsIFkh,options);
      options(10)=OUTPUT.funcCount;
      zk(i)=zk(i)+mkest;                            %      and add the mean to the mode estimate
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
  
  %%%%%% Display the point number if options(1)==1
  
  if options1==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;
