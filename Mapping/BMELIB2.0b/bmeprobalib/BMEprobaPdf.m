function [z,pdf,info]=BMEprobaPdf(z,ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
  covmodel,covparam,nhmax,nsmax,dmax,order,options);

% BMEprobaPdf               - BME pdf prediction with proba data (Jan 1, 2001)
%
% BME computation of the posterior probability density function
% at an estimation point, using both hard data and soft probabilistic data
% that can reasonably be assumed as Gaussian distributed. The
% function is intendend to be as general as possible, covering
% various situations, like multi field variables, nested models, 
% space-time estimations, non-stationarity of the mean, etc. 
% Depending on the case, specific format are needed for the input 
% variables. This function has the potential to process both hard 
% data and soft probability data when both are available, but it 
% reduces to the kriging posterior pdf when only hard data is used.
%
% SYNTAX :
%
% [z,pdf,info]=BMEprobaPdf(z,ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
%   covmodel,covparam,nhmax,nsmax,dmax,order,options);
%
% INPUT :
%
% z          nz by 1      vector of the values for which the posterior pdf is 
%                         computed. z may be empty, in which case a grid of
%                         z values is automatically constructed (see options(6)
%                         and options(7) below)
% ck          1 by d      row vector of coordinates for the estimation locations.
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
%                         estimation of these integrals (default value is 1e-4).  
%                         options(6) is used  when z=[] to determine how fine of a grid of z 
%                         values toconstruct (default value is 25 ; increase this value to  
%                         construct a finer grid),
%                         options(7) is used  when z=[] to determine the extend of the grid of 
%                         z values to construct (default value is 0.001 ; decrease this value 
%                         to widen the range of the grid of z values)
%
% OUTPUT :
%
% z          nz by 1      vector of the values for which the posterior
%                         pdf is computed
% pdf        nz by 1      vector of posterior pdf values using BME with probas
% info       nz by 1      vector of information about the computation in the local 
%                         neighbourhood around the estimation point :
%                         info=NaN : no computation at all, no hard or soft data
%                         info=0   : computation using BME with at least 1 soft data
%                         info=1   : dubious results, integration error above tolerance
%                                    when integrating over the soft pdf. Try to increase 
%                                    the value of options(3)
%                         info=3   : computation using kriging, no soft data
%                         info=4   : computation provides hard data value at estimation point
%                         info=10  : dubious results from integration routine. Check the
%                                    routine used for the numerical integration.
%
% NOTE :
%
% All the specific conventions for specifying nested models, multivariate 
% or space-time cases are the same as for BMEprobaMoments.m.

%%%%%% Error messages

[ch,zh,cs,nl,limi,probdens]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
  covmodel,covparam,nhmax,nsmax,dmax);
if length(order)>1, error('order must be a scalar'); end;
noindex=~iscell(ck);  % test if there is an index for the variables
if noindex==1,
  nk=size(ck,1);      % nk is the number of estimation points
else
  nk=size(ck{1},1);
end;
if nk~=1, error('nk must be equal to 1'); end;

%%%%%% Initialize default value of some input parameters

if nargin<16,         % initialize options with default values if not provided as input
  options=BMEoptions;
end;
if ~isempty(z), z=z(:); end;

%%%%%% Select the local neighbourhood

ck0=ck;
[chlocal,zhlocal,dh,sumnhlocal]=neighbours(ck0,ch,zh,nhmax,dmax);
[cslocal,nllocal,limilocal,probdenslocal,ds,sumnslocal]=...
  probaneighbours(ck0,cs,softpdftype,nl,limi,probdens,nsmax,dmax);
[chlocal,zhlocal,cslocal]=fixempty(ck0,chlocal,zhlocal,cslocal);

%%%%%% Test if local neighbourhood is empty

if isempty(zhlocal) & isempty(nllocal)
  nz=length(z);         % nz is the number of values for pdf computation
  pdf=zeros(nz,1)*NaN;
  info=zeros(nz,1)*NaN;
  iscomputed=1;
  return;
end;

%%%%%% If z is empty, then construct of grid of z values to plot the posterior pdf

zempty=isempty(z);
if zempty
  [index]=findpairs(ck,chlocal);                         % test if there is a hard data at estimation point
  if ~isempty(index),  
    z=zhlocal(index(2));                                      % return the hard data value as the z vector
    pdf=Inf;                                                  % set the corresponding pdf value to +Inf
    info=4;                                                   % set the value of info to 4
    return;
  end;
  if nargin<16 | length(options)<7                  % if there is no hard data data at the estimation point
    zNPTS=25;                                       % then construct the grid of z values using the 
    frEps=0.001;                                    % parameters options(6) and options(7)
  else
    zNPTS=options(6);
    frEps=options(7);
  end 
  if zNPTS<=2 warning('options(6) seems to be too small, you should use at least 5'); end;
  if frEps<=0 warning('options(7) seems to be suspiciously small'); end;
  
  zMode=BMEprobaMode(ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...      % get the mode of the posterior pdf
    probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  
  if isnan(zMode)
    z=[];
    pdf=[];
    info=NaN;
    return;
  end;
  
  [index]=findpairs(ck,cslocal);                                   % set the min and max of z values
  if ~isempty(index), 
    [zMin,zMax]=proba2interval(softpdftype,nllocal(index(2)),...   
      limilocal(index(2),:)); 
  else
    [a,b]=proba2interval(softpdftype,nl,limi);
    zMin=min([zh;a]);                                         
    zMax=max([zh;b]);                                         
    if zMode-zMin<=0
      if zMax-zMode>0
        zMin=zMode-(zMax-zMode);
      else
        zMin=zMode-1; 
      end;
    end;
    if zMax-zMode<=0
      if zMode-zMin>0
        zMax=zMode+(zMode-zMin);
      else
        zMax=zMode+1; 
      end;
    end;
    [dummy,fzMode]=BMEprobaPdf(zMode,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...
      probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    [dummy,fzMin]=BMEprobaPdf(zMin,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...
      probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    while fzMin/fzMode<0.2*frEps                                            % Added Nov 14, 2000
      zMin=zMin+0.4*(zMode-zMin);                                           % Added Nov 14, 2000  
      [dummy,fzMin]=BMEprobaPdf(zMin,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...    % Added Nov 14, 2000
        probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);         % Added Nov 14, 2000
    end;                                                                    % Added Nov 14, 2000
    while fzMin/fzMode>frEps
      zMin=zMin-0.5*(zMode-zMin);  
      [dummy,fzMin]=BMEprobaPdf(zMin,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...
        probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options); 
    end;
    [dummy,fzMax]=BMEprobaPdf(zMax,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...
      probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    while fzMax/fzMode<0.2*frEps                                            % Added Nov 14, 2000
      zMax=zMax-0.4*(zMax-zMode);                                           % Added Nov 14, 2000
      [dummy,fzMax]=BMEprobaPdf(zMax,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...    % Added Nov 14, 2000
        probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);         % Added Nov 14, 2000
    end;                                                                    % Added Nov 14, 2000
    while fzMax/fzMode>frEps
      zMax=zMax+0.5*(zMax-zMode);
      [dummy,fzMax]=BMEprobaPdf(zMax,ck,chlocal,cslocal,zhlocal,softpdftype,nllocal,limilocal,...
        probdenslocal,covmodel,covparam,nhmax,nsmax,dmax,order,options);
    end;
  end;
  zStep=(zMax-zMin)/zNPTS;                                     % Use the min and max of z and the mode to 
  z=(zMin:zStep:zMax)';                                        % determine the grid of z values.
  zStep=0.2*zStep;
  zNPTS=ceil(zNPTS/2);
  z=[z; ((zMode-zNPTS*zStep):zStep:(zMode+zNPTS*zStep))'];
  z=sort(z);
end;

%%%%%% Initialize the parameters

nz=length(z);         % nz is the number of values for pdf computation

maxpts=options(3);
aEps=0;
rEps=options(4);

pdf=zeros(nz,1)*NaN;
info=zeros(nz,1)*NaN;

%%%%%% Test if there is a hard data at estimation point

iscomputed=0;  
[index]=findpairs(ck0,chlocal);          % test if there is a hard data at estimation point
if ~isempty(index),                      %   If there is hard data at the estimation point
  idxzEqHard=find(z==zhlocal(index(2))); %   set the value of the pdf to Inf when z=hard val        
  if ~isempty(idxzEqHard)
    pdf(idxzEqHard)=Inf;
  end;
  idxzNonEqHard=find(z~=zhlocal(index(2)));
  if ~isempty(idxzNonEqHard)             %   else set the pdf to zero
    pdf(idxzNonEqHard)=0;
  end;
  info(:)=4;
  iscomputed=1;
  return;  
end;

%%%%%% Test if there is a soft data at estimation point

isduplicate=0;
[index]=findpairs(ck0,cslocal);                % test if there is a soft data at estimation point
if ~isempty(index),                            % and split the soft data into two sets
  [cest,cslocal,nlest,limiest,probdensest,nllocal,limilocal,probdenslocal]=...
    probasplit(cslocal,softpdftype,nllocal,limilocal,probdenslocal,index(2));
  sumnslocal=sumnslocal-1;
  isduplicate=1;                               % specify that there is a soft data at estimation point
end;

%%%%%% Returns the kriging pdf when sumnslocal==0

% This next block can be turned on or off.  
% Turning the block ON will will use universal kriging when their is no soft data
% Turning the block OFF will treat universal kriging as a special case of BME
% Note: When order=NaN (simple kriging case) then kriging and BME give the same error variance.
%       However when order=0, 1, or 2, then there is a difference between the error
%       variance returned by BME and Universal Kriging. The difference comes from the fact 
%       that BME neglects the uncertainty associated with the estimated mean trend.

useUniversalKriging=0;
if useUniversalKriging==1
 if (sumnhlocal>0)&(sumnslocal==0)&(isduplicate==0)&(iscomputed==0), % test if there are no soft data
  K=coord2K(chlocal,chlocal,covmodel,covparam);    % built the left-hand side
  k=coord2K(chlocal,ck0,covmodel,covparam);        % built the right-hand side
  Kk=coord2K(ck0,ck0,covmodel,covparam);           % 
  [X,x]=krigconstr(chlocal,ck0,order);              % built the constraints matrix and vector
  nx=size(X,2);
  Kadd=[[K,X];[X',zeros(nx)]];
  kadd=[k;x];
  lammu=Kadd\kadd;                                  % compute the kriging weights and Lagrangian coef
  lam=lammu(1:sumnhlocal);                          % get the kriging weights from the solution
  mu=lammu(sumnhlocal+1:end);                       % get the Lagrangian coef from the solution
  lamt=lam';
  zk=lamt*zhlocal;                                  % compute the universal kriging estimate zk
  vk=Kk-lamt*k;                                 % and the universal kriging variance vk
  if ~isempty(mu)
    vk=vk-mu'*x;
  end;
  pdf=gausspdf(z,[zk vk]);                          % compute the value of the Gaussian pdf
  info=ones(nz,1)*3;                                % set the value of info to 3
  if isduplicate
    Kh=K;                                           % if duplicate soft data at est
    Ks=zeros(0,0);                                  % point, then construct all these
    Kk_h=k';                                        % cov matrix for later use in
    Kk_s=zeros(1,0);                                % convulution with soft info at
    Ks_h=zeros(0,1);                                % est point
    Kkh=[[Kk,Kk_h];[Kk_h',Kh]];
    Ks_kh=[Kk_s',Ks_h];
    invKh=inv(Kh);  
  end;
  iscomputed=1;  
 end;
end;

%%%%%% Calculate and returns the BME pdf in the general case

ismeanadjusted=0;
if (iscomputed==0),                % test if the pdf has not already been computed 

  pdf=zeros(nz,1)*NaN;                            % initialize the pdf values
  info=zeros(nz,1)*NaN;                           % initialize the info values
  Kk=coord2Kinterface(ck0,ck0,covmodel,covparam);         % compute variance for estimation point
  Kh=coord2Kinterface(chlocal,chlocal,covmodel,covparam); % compute covariance matrix for hard data
  Ks=coord2Kinterface(cslocal,cslocal,covmodel,covparam); % compute covariance matrix for soft data
  Kk_h=coord2Kinterface(ck0,chlocal,covmodel,covparam);     % compute cross-covariance vector estimation/hard data
  Kk_s=coord2Kinterface(ck0,cslocal,covmodel,covparam);     % compute cross-covariance vector estimation/soft data
  Ks_h=coord2Kinterface(cslocal,chlocal,covmodel,covparam); % compute cross-covariance matrix soft/hard data
  Kkh=[[Kk,Kk_h];[Kk_h',Kh]];                      % build estimation+hard data covariance matrix
  Ks_kh=[Kk_s',Ks_h];                                % build cross-covariance matrix soft/estimation+hard data
  
                                    %%%%%%% Remove the mean trend from the data
  ismeanadjusted=1;
  [mslocal,vslocal]=proba2stat(softpdftype,nllocal,limilocal,probdenslocal);
  [mkest,mhest,msest]=localmeanBME(ck0,chlocal,cslocal,zhlocal,mslocal,vslocal,Kh,Ks_h,Ks,order);
  z=z-mkest;                                      % substract the mean from the z variable
  if sumnhlocal>0  
    zhlocal=zhlocal-mhest;                        % substract the mean from the hard data
  end;  
  if sumnslocal>0                                 % substract the mean from the soft data
    [limilocal]=probaoffset(softpdftype,nllocal,limilocal,-msest);
  end;
  if isduplicate                                  % substract mean from soft data at est pt
    [limiest]=probaoffset(softpdftype,nlest,limiest,-mkest);
    [aest,best]=proba2interval(softpdftype,nlest,limiest);
  end;
  
                                    %%%%%%% Calculate fden, the denominator of the pdf
                                          % if no duplicate (i.e. no soft data at est point)
  if ~isduplicate                         %   fden=mvnpdf(Xh) Int[dXs fS(Xs) mvnpdf(Xs|h)]
    if sumnhlocal==0,                             % if no hard data use the soft covariance 
      ms=0*nllocal;
      [fden err info]=mvPro(softpdftype,nllocal,limilocal,probdenslocal,ms,Ks,maxpts,aEps,rEps);
    else                                          % else use the conditional soft covariance
      invKh=inv(Kh);
      cstden=1/(((2*pi)^(sumnhlocal/2))*sqrt(det(Kh)));
      fden=cstden*exp(-0.5*zhlocal'*invKh*zhlocal);
      if sumnslocal>0
        BsIFh=Ks_h*invKh;
        KsIFh=Ks-BsIFh*Ks_h';                     % compute the conditional covariance matrix soft given h
        msIFh=BsIFh*zhlocal;
        [Pden Perr Pinfo]=mvPro(softpdftype,nllocal,limilocal,probdenslocal,msIFh,KsIFh,maxpts,aEps,rEps);
        fden=fden*Pden;
      end;
    end;  
  else                                    % else in duplicate case,                         
    Kks=[[Kk,Kk_s];[Kk_s',Ks]];           %  fden=mvnpdf(Xh) Int[dXks fS(Xks) mvnpdf(Xks|h)]
    Kks_h=[Kk_h;Ks_h];
    [dummy,nlks,limiks,probdensks]=...    %combine soft data at est point and soft non-est points
      probacat(softpdftype,nlest,limiest,probdensest,softpdftype,nllocal,limilocal,probdenslocal);
    if sumnhlocal==0,                               % if no hard data use cov Kks
      mks=[0*nlks];
      [fden err info]=mvPro(softpdftype,nlks,limiks,probdensks,mks,Kks,maxpts,aEps,rEps);
    else                                            % else use the conditional cov KksIFh
      invKh=inv(Kh);
      cstden=1/(((2*pi)^(sumnhlocal/2))*sqrt(det(Kh)));
      fden=cstden*exp(-0.5*zhlocal'*invKh*zhlocal);
      BksIFh=Kks_h*invKh;
      KksIFh=Kks-BksIFh*Kks_h';
      mksIFh=BksIFh*zhlocal;
      [Pden Perr Pinfo]=mvPro(softpdftype,nlks,limiks,probdensks,...
        mksIFh,KksIFh,maxpts,aEps,rEps);
      fden=fden*Pden;
    end;
  end;
                                   %%% Calculate fnum, the numerator of the pdf
  fnum=zeros(nz,1);                %%% fnum=mvnpdf(Xkh) Int[dXs fS(Xs) mvnpdf(Xs|kh)] 
  info=zeros(nz,1); 
  if fden==0,                                     % if denominator=0, set pdf values to 0 
    pdf=zeros(nz,1); 
  else                                            % else compute the numerator values
    cstnum=1/(((2*pi)^((sumnhlocal+1)/2))*sqrt(det(Kkh)));
    invKkh=inv(Kkh);
    BsIFkh=Ks_kh*invKkh;
    KsIFkh=Ks-BsIFkh*Ks_kh';
    for i=1:nz,
      msIFkh=BsIFkh*[z(i);zhlocal];
      fnum(i)=cstnum*exp(-0.5*[z(i);zhlocal]'*invKkh*[z(i);zhlocal]);
      info(i)=0;
      if sumnslocal>0,
        [Pnum err infointeg]=mvPro(softpdftype,nllocal,limilocal,probdenslocal,msIFkh,KsIFkh,maxpts,aEps,rEps);
        fnum(i)=fnum(i)*Pnum;
        if infointeg==0,info(i)=0;end;
        if infointeg==1,info(i)=1;end;
        if infointeg>1,info(i)=10;end;
      end;
    end;
    if isduplicate                               % if duplicate, fnum=fnum*fS(Xk)
      fnum=proba2val(z,softpdftype,nlest,limiest,probdensest).*fnum;
    end;
  end;
  
  if fden>0
    pdf=(1/fden)*fnum;                  %%%% Apply the normalization constant to fnum
  end;
    
                                       %%%%%%  If duplicate add bounds to the pdf
  if zempty & isduplicate
    idxinside=(aest<=z & z<=best);
    z=[aest;z(idxinside);best];            % add the bounds of the interval to the z vector
    pdf=[0;pdf(idxinside);0];              % add 0 as value of the pdf at the bounds
    ninfo=length(info);                    % compute the length of the info vector
    info=[info(1);info(idxinside);info(ninfo)];       % add first and last info values at the bounds
  end;
  
  if ismeanadjusted  %%%%%% Add the mean at the estimation point
    z=z+mkest;
  end;
end;


