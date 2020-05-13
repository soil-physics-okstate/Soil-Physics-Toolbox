function [moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);

% BMEprobaMoments           - BME Moments prediction with probabilistic data (Jan 1, 2001)
% 
% BME computation of the statistical moments of the posterior pdf at a set
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
% [moments,info]=BMEprobaMoments(ck,ch,cs,zh,softpdftype,nl,limi,probdens,...
%   covmodel,covparam,nhmax,nsmax,dmax,order,options);
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
%                         fmin.m MATLAB optimization routine that finds the mode of the
%                         probability distribution function (default values are the same
%                         as for fmin.m),
%                         options(3) specifies the maximum number of evaluation that can
%                         be done by the FORTRAN77 subroutines for the integrals (default
%                         value is 50 000 ; this value should be increased if a warning
%                         message appears on the screen during the computation),
%                         options(4) specifies the maximum admissible relative error on the
%                         estimation of these integrals (default value is 1e-4). 
%                         options(8) number of moments to calculate (1, 2 or 3)
%                                 =1 to calculate the mean of the BME posterior pdf,
%                                 =2 for the mean and estimation error variance
%                                 =3 for the mean, estimation variance and coef of skewness.
%                         The default value for options(8) is 2
%
% OUTPUT :
%
% moments    nk by 3      matrix of the moments of the cBME estimate at each of the nk estimation points.
%                         The first column are the BME mean, 
%                         The second column has the BME estimation error variance (if options(8)>=2)
%                         The third column has the coeff of skewness (if options(8)>=3)
% info       nk by 3      vector of information about the computation in the local 
%                         neighbourhood around the estimation point :
%                         info=NaN : no computation, no hard or soft data
%                         info=0   : computation using BME with soft data 
%                         info=1   : dubious results, integration error above tolerance
%                                    when integrating over the soft pdf. Try to increase 
%                                    the value of options(3)
%                         info=3   : computation using only hard data
%                         info=4   : computation provides hard data value at estimation point
%                         info=10  : dubious results from integration routine. Try reducing
%                                    the number of soft data points ns, or option(3) or option(4)
%
% NOTE :
%
% 0- In the case that there only soft probabilistic data,
% ch and zh can be entered as the empty [ ] matrices. Similarly when
% there is only hard data, then cs, nl limi and probdens can be
% entered as empty [ ] matrices, and the results are the kriging
% estimates.
%
% 1- It is important to make sure that the number of points in each spatial
% neighbourhood which is determined by nhmax, nsmax and dmax is sufficient with
% respect to the selected order value. E.g., for order=1, there must be
% at least 1+d points in each neighbourhood, for order=2, there must be
% at least 1+2*d points, etc., where d is the dimension of the space.
%
% 2- These are the generic input and output variables. The specific conventions
% for some of these input are depending on the case that is considered, so
% that different situations can be easily covered by the same function. These
% specific conventions have been adopted in a standard way whenever it was
% possible to do so for the functions in this package. The various cases are
% the following :
%
% i) The univariate case
%    -------------------
%
% The simplest use of the BMEprobaMoments.m function is in the case where there is
% only one field variable that is used in the estimation process (univariate case).
% When the covariance model is unique (no nested models), covmodel
% is a simple string that contains the name of the covariance model,
% e.g., covmodel='exponentialC'. In that case, covparam corresponds to the vector
% of parameters for the corresponding covmodel, e.g., covparam=[1,4] for the 
% exponential model with sill (i.e. variance) and range equal to 1 and 4, respectively.
% 
% For nested models, one need to specify the names of the various models and
% their corresponding parameters. This can be done by using the MATLAB cell
% array notation. E.g., for a nested model including also a nugget effect, the
% notation becomes covmodel={'nuggetC','exponentialC'}. The corresponding covparam
% vector is now, e.g., covparam={[0.2],[1,4]}, where 0.2 is the sill of the nugget
% effect model. Using this notation, one can use as many nested models as wanted.
%
% ii) The multivariate case
%     ---------------------
%
% It is also possible to process several field variables at the same time 
% (multivariate case). It is then needed to specify additional tags in the ck, 
% ch and cs matrices. These tags are provided as vectors of values that refer to 
% the field variable, the values in these vectors ranging from 1 to nv, where nv is 
% the number of field variables. E.g., if there are 3 natural field variables, the indexk, 
% indexh and indexs vectors must be defined, having same number of rows as ck, ch and cs, 
% respectively, and they contain values equal to 1, 2 or 3. These vectors are grouped 
% with the coordinates matrices using the MATLAB cell array notation, so that 
% ck={ck,indexk}, ch={ch,indexh} and cs={cs,indexs} are now the correct input variables.
% If indexk is a single scalar value, the same field variable specified by indexk will be
% estimated at all ck coordinates. For the multivariate case, the nhmax variable
% is now a 1 by nv vector specifying the maximum number of hard data values of each
% field variable to be used for the estimation. Using the same logic, order is now a nv
% by 1 vector specifying the order of the polynomial drift along the spatial axes
% for each field variable.
%
% When there is only one covariance model to specify, covmodel is still a string 
% that contains the name of the covariance model. For covparam, it is now needed 
% to specify the covariance matrix between the field variables. E.g., for 
% an exponential model having a range of 1.4, we would have covparam={C,[1.4]}, 
% where C is the symmetric nv by nv covariance matrix between the nv field
% variables. Note that covariance models that are characterized by a
% single parameter (e.g., nuggetC.m) require that the second cell of the param 
% variable is the empty [ ] matrix. E.g., for a nugget effect covariance model, 
% param={C,[ ]}.
%
% When the covariance model is a nested model, the covmodel cell array
% is written in the same way than for nested models in the univariate case, but
% covparam is now a cell array containing as many cells as there are variables. E.g.,
% for a nested model including also a nugget effect, the notation becomes covparam=
% {{Cn,[ ]},{Ce,[1.4]}}, where Cn and Ce are the symmetric nv by nv covariance
% matrices associated with the nugget effect and the exponential models (see also
% the vario.m, covario.m and coregfit.m functions for explanations about estimating
% and modeling these covariance matrices). Nested models in the multivariate case
% are thus defined as a Linear Model of Coregionalization. Using this notation, one
% can use as many nested models as wanted.
%
% iii) The space-time case
%      -------------------
%
% For space/time data, the convention is that the last column of the ck, ch and cs
% matrices of coordinates corresponds to the time axis. As space/time data are not
% necessarily processed using the same covariance models than spatial
% or temporal data, it is necessary to precise if the data are space/time or not
% using specific conventions for the name of the variogram or covariance models.
%
% Separable covariance models, which are products of spatial and temporal covariance
% models, are coded as model='covmodelS/covmodelT', where covmodelS refers to the
% name of the spatial covariance model and covmodelT refers to the name of the
% temporal covariance model. The function will automatically detect the occurence
% of the '/' character inside the string and will process the data accordingly. Any
% combination of valid covariance models can be used. E.g., if the spatial covariance
% model is exponential and the temporal covariance model is spherical, model=
% 'exponentialC/sphericalC'. For the covparam vector, the first value in the vector
% always corresponds to the space/time variance. The remaining values are the
% remaining parameters for the spatial and the temporal covariance models, in this
% order. E.g., param=[1,1.4,0.7] if the variance of the space/time model is equal to 1,
% the range of the exponential spatial covariance model is equal to 1.4 and the range
% of the spherical temporal covariance model is equal to 0.7. For non-separable
% variogram or covariance models, the name contained in the string must end up with the
% ST characters. The occurence of these characters in the string is automatically
% detected by the function. For the param vector, the first value must always
% correspond to the space/time variance. For the multivariate case or for nested
% covariance models, the conventions are exactly the same as described for spatial data,
% whether the space/time model is separable or not.
%
% For space/time data, the dmax variable becomes a 1 by 3 vector instead of a single
% scalar value. dmax(1) is the maximum spatial distance between an estimation location
% and existing data locations, and dmax(2) is the maximum temporal lag between an
% estimation location and existing data locations. Only data locations that respect
% both conditions are considered for the estimation. dmax(3) refers to a space/time
% metric, such that the space/time distance=spatial distance+dmax(3)*temporal
% distance. The definition of this space/time distance is needed for selecting the
% nhmax closest data locations from the estimation locations.
%
% Is is possible to specify a different order for the polynomial along the spatial
% axes and the temporal axis. For the univariate case, order is a 1 by 2 vector,
% where order(1) is the order of the polynomial along the spatial axes and order(2)
% is the order of the polynomial along the temporal axis. For the multivariate case
% where there are nv different variables, order is a nv by 2 matrix, where the first
% and second columns of order contain the order for the spatial and temporal
% polynomials for each one of the nv variables, respectively. If in that case order is
% entered as a 1 by 2 vector, the same spatial order corresponding to order(1) and
% the same temporal order corresponding to order(2) will be used for all the variables.
%
% iv) The geometric anisotropy case
%     -----------------------------
%
% In some circumstances, the covariance or variogram model may have parameters that
% change according to the direction of the distance vector. E.g., the range can be
% maximum in a specific direction and minimum in the perpendicular direction. If the
% way the range parameter changes according to the angle can be described by an ellipse
% for two dimensional coordinates or by an ellipsoid for three dimensional coordinates,
% this anisotropy is called a geometric anisotropy. The BMEprobaMoments.m function has not been
% written in a way that it will specifically process an anisotropic covariance.
% The reason is that it is always possible to do a coordinate transformation  
% using the aniso2iso.m function, so that ck and ch coordinates are mapped
% into an isotropic space. The ck and ch vectors are then simply substituted in the
% BMEprobaMoments.m function with the vectors aniso2iso(ck,angle,ratio) and aniso2iso(ch,angle,
% ratio), respectively (see the aniso2iso.m function for more information about the
% definition of angle and ratio). The value for the range of the covariance
% function model specified in param is always the range of this model along the principal
% axis of the ellipse or ellipsoid. It is worth noting that for nested covariance 
% functions (i.e., when covmodel and covparam are cell arrays), the same anisotropy
% will hold for all of the specified models.

%%%%%% Error messages

[ch,zh,cs,nl,limi,probdens]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax);

%%%%%% Initialize the parameters

global INFOINTEG      % information message from the Fortran subroutine
if nargin<15,         % initialize options with default values
  options=BMEoptions;
  options1=1;
else
  options1=options(1);
  options(1)=0;
end;

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

moments=zeros(nk,3)*NaN;
info=zeros(nk,3)*NaN;

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
    moments(i,:)=NaN;
    info(i,:)=NaN;
    iscomputed=1;
  end;

  %%%%%% Test if there is a hard data at estimation point
  
  [index]=findpairs(ck0,chlocal);                
  if iscomputed==0 & ~isempty(index),            % if there is a hard data at estimation point,
    moments(i,:)=0;
    moments(i,1)=zhlocal(index(2));              % return the hard data value as the mean estimate
    info(i,:)=4;                                 % set the value of info to 4 
    iscomputed=1;                                % specify that the value has been computed
  end;
  
  %%%%%% Do the BME estimation if there is no hard data at the estimation point
  
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
    
                 %%%%%%% Calculate the covariance matrices
    Kk=coord2Kinterface(ck0,ck0,covmodel,covparam);           % compute variance for estimation point
    Kh=coord2Kinterface(chlocal,chlocal,covmodel,covparam);   % compute covariance matrix for hard data
    Ks=coord2Kinterface(cslocal,cslocal,covmodel,covparam);   % compute covariance matrix for soft data
    Kk_h=coord2Kinterface(ck0,chlocal,covmodel,covparam);     % compute cross-covariance vector estimation/hard data
    Kk_s=coord2Kinterface(ck0,cslocal,covmodel,covparam);     % compute cross-covariance vector estimation/soft data
    Ks_h=coord2Kinterface(cslocal,chlocal,covmodel,covparam); % compute cross-covariance matrix soft/hard data
    
                 %%%%%%% Calculate the mean trend and remove it from the data
    [mslocal,vslocal]=proba2stat(softpdftype,nllocal,limilocal,probdenslocal);
    [mkest,mhest,msest]=localmeanBME(ck0,chlocal,cslocal,zhlocal,mslocal,vslocal,Kh,Ks_h,Ks,order);
    if sumnhlocal>0                                 % substract the mean from the hard data
      zhlocal=zhlocal-mhest;
    end;
    if sumnslocal>0                                 % substract the mean from the soft data
      [limilocal]=probaoffset(softpdftype,nllocal,limilocal,-msest);
    end;
    if isduplicate                                  % substract mean from soft data at est pt
      [limiest]=probaoffset(softpdftype,nlest,limiest,-mkest);
    end;
      
                  %%%%%%% Calculate the statistical moments
    if ~isduplicate
      if sumnhlocal==0
        BsIFh=[];
        KsIFh=Ks;
      else
        invKh=inv(Kh);
        BsIFh=Ks_h*invKh;
        KsIFh=Ks-BsIFh*Ks_h';
      end;
      Khs=[[Kh,Ks_h'];[Ks_h,Ks]];
      Kk_hs=[Kk_h Kk_s];
      invKhs=inv(Khs);
      BkIFhs=Kk_hs*invKhs;
      KkIFhs=Kk-BkIFhs*Kk_hs';
      [BMEmean,stdDev,skewCoef,info(i,1:3)]=momentsFun(zhlocal,softpdftype,...
        nllocal,limilocal,probdenslocal,options,BsIFh,KsIFh,BkIFhs,KkIFhs);
      moments(i,1)=BMEmean; 
      moments(i,2)=stdDev^2; 
      moments(i,3)=skewCoef; 
    else
      Kks=[[Kk,Kk_s];[Kk_s',Ks]];
      if sumnhlocal==0
        BksIFh=[];
        KksIFh=Kks;
      else
        invKh=inv(Kh);
        Kks_h=[Kk_h;Ks_h];
        BksIFh=Kks_h*invKh;
        KksIFh=Kks-BksIFh*Kks_h';
      end;
      [BMEmean,stdDev,skewCoef,info(i,1:3)]=momentsDupFun(zhlocal,softpdftype,...
        nllocal,limilocal,probdenslocal,options,BksIFh,KksIFh,nlest,limiest,probdensest);
      moments(i,1)=BMEmean; 
      moments(i,2)=stdDev^2; 
      moments(i,3)=skewCoef; 
    end;
    
                  %%%%%%% Add the mean trend back to the BME mean
    moments(i,1)=moments(i,1)+mkest;
  end;
  
  %%%%%% Display the point number of options(1)==1
  if options1==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;
