function  [zlCI,zuCI,pdfCI,PCI,z,pdf]=BMEprobaCI(ck,ch,cs,zh,softpdftype,...
  nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);

% BMEprobaCI                - BME Confidence Interval and posterior pdf (Jan 1, 2001)
%
% BME computation of the confidence interval (CI) and the BME posterior pdf 
% at a set of estimation points, using both hard and soft probabilistic data 
% that can reasonably be assumed as Gaussian distributed. The
% function is intendend to be as general as possible, covering
% various situations, like multi field variables, nested models, 
% space-time estimations, non-stationarity of the mean, etc. 
% Depending on the case, specific format are needed for the input 
% variables. This function has the potential to process both hard 
% data and soft probability data when both are available, but it 
% reduces to the kriging CI when only hard data is used. In the most
% simple case the CI is just one interval with lower bound zlCI and upper
% bound zuCI, but in general the CI may be composed of the union of several 
% disconnected intervals, in which case zlCI and zuCI store the minimum of the
% lower bounds and maximum of the upper bounds of these intervals, respectively. 
%
% SYNTAX :
%
% [zlCI,zuCI,pdfCI,PCI,z,pdf]=BMEprobaCI(ck,ch,cs,zh,softpdftype,...
%   nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
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
% options    1 by 1 or 29 vector of optional parameters that can be used if default
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
%                         options(6) is used to determine how fine of a grid of z values to
%                         construct (default value is 25 ; increase this value to construct 
%                         a finer grid),
%                         options(7) is used to determine the extend of the grid of z values
%                         to construct (default value is 0.001 ; decrease this value to widen
%                         the range of the grid of z values), 
%                         options(20:29) confidence probability for which to compute the confidence 
%                         interval. The CI are calculated only for options(20:29) which 
%                         are >=0.01 and <=0.99. The default value is options(20)=0.68 and 
%                         options(21:29)=0
%
% OUTPUT :
%
% zlCI       nk by nCI    matrix storing the lower bound of the CI, where nk 
%                         is the number of est pts, and nCI the number of confidence probabilities
% zuCI       nk by nCI    matrix storing the upper bounds of the CI.
% pdfCI      nk by nCI    matrix of pdf value at the bounds of the CI
% PCI        1 by nCI     vector of the probability of the confidence intervals
% z          nk by 1 cell of vectors storing the grid of z value defining the 
%                         BME posterior pdf at each of the nk estimation point
% pdf        nk by 1 cell of vectors storing the BME pdf value for each of 
%                         z grid
%
% NOTE :
%
% All the specific conventions for specifying nested models, multivariate 
% or space-time cases are the same as for BMEprobaMoments.m.

%%%%%% Error messages

[ch,zh,cs,nl,limi,probdens]=fixempty(ck,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,...
  nhmax,nsmax,dmax);

%%%%%% Get the number of estimation points and the probabilities of CI

noindex=~iscell(ck);       % test if there is an index for the variables
if noindex==1,
  nk=size(ck,1);           % nk is the number of estimation points
else
  nk=size(ck{1},1);
end;

if length(options)<20, 
  error('options(20) must define the confidence probability'); 
end;

PCI=options(20:end);
PCI=PCI(~isnan(PCI));
PCI=PCI( 0.01<=PCI & PCI<=0.99 );
nCI=length(PCI);

%%%%%% Calculate the BME posterior pdf and the CI

for ik=1:nk
  if noindex==1,
    cki=ck(ik,:);           % ik-th estimation point
  else
    cki={ck{1}(ik,:),ck{2}(ik,:)};
  end;
  [z{ik,1},pdf{ik,1},info]=BMEprobaPdf([],cki,ch,cs,zh,softpdftype,nl,limi,...
    probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  [zlCIik,zuCIik,pdfCI(ik,1:nCI),lCI]=pdf2CI(z{ik},pdf{ik},PCI);
  if nCI==1, 
    zlCI(ik,1)=min(zlCIik);
    zuCI(ik,1)=max(zuCIik);
  else
    for iCI=1:nCI,
      zlCI(ik,iCI)=min(zlCIik{iCI});
      zuCI(ik,iCI)=max(zuCIik{iCI});
    end;
  end;
end;
