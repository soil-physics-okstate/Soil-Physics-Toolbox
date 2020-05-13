function [BMEmean,stdDev,skewCoef,info]=momentsFun(zh,softpdftype,nl,limi,probdens,...
  options,BsIFh,KsIFh,BkIFhs,KkIFhs)

% momentsFun                - Moments of the BME posterior pdf (Jan 1, 2001)
%
% Calculates the mean , the standard deviation and the coefficient 
% of skewness of the BME posterior pdf.
% The number of moments calculated is controled by nMom=options(8) (1 for 
% mean, 2 for mean and sdt dev, 3 for all three moments)
%
% SYNTAX :
%
% [BMEmean,stdDev,skewCoef,info]=momentsFun(zh,softpdftype,nl,limi,probdens,...
%   PDFparam,BsIFh,KsIFh,BkIFhs,KkIFhs) 
%
% INPUT :
%
% zh         nh by 1   vector of the values for the hard data
% softpdftype scalar   integer for the type of soft pdf, representing the 
%                      probabilitic soft data. These soft pdf types are as follow:
%                      1-Histogram, 2-Linear, 3-Grid histogram, 4-Grid Linear. 
%                      (see also probasyntax)
% nl         ns by 1   vector of number of interval limits. nl(i) is 
%                      the number of interval limits used to define the soft pdf 
%                      for soft data point i. (see also probasyntax)
% limi       ns by l   matrix of interval limits, where l is equal to
%                      either max(nl) or 3 depending of the softpdftype.
%                      limi(i,:) are the limits of intervals for the i-th 
%                      soft data. (see also probasyntax)
% probdens   ns by p   matrix of probability density values, where p is 
%                      equal to either max(nl)-1 or max(nl), depending on the 
%                      softpdftype. probdens(i,:) are the values of the probability 
%                      density corresponding to the intervals for the i-th soft data 
%                      defined in limi(i,:). (see also probasyntax)
% options    vector of optional parameters, where :
%                      options(3) is the maximum number of evaluation of the integral.
%                                 The default value is 50000
%                      options(4) is the relative error on the estimation of the integral.
%                                 The default value is equal to 1e-4 
%                      options(8) number of moments to calculate (1, 2 or 3)
%                       =1 to calculate the mean of the BME posterior pdf,
%                       =2 for the mean and standard deviation,
%                       =3 for the mean, standard deviation and coef of skewness.
% BsIFh      ns by nh  matrix, BsIFh=Ks_h*inv(Kh)
% KsIFh      ns by ns  matrix, KsIFh=BsIFh*Ks_h' is the cov mat of soft data 
%                                               conditional to hard data
% BkIFhs     1 by nh+ns matrix, BkIFhs=Kk_hs*inv(Khs)
% KkIFhs     1 by 1    matrix, KkIFhs=BkIFhs*Kk_hs' is the cov mat of est point
%                                               conditional to hard+soft data     
%
% OUTPUT :
%
% BMEmean    scalar    the mean of the BME posterior pdf 
% stdDev     scalar    the standard deviation of the BME posterior pdf
% skewCoef   scalar    the coefficient of skewness of the BME posterior pdf
%
% info       1 by 3    vector of information about the numerical computation of
%                      the 3 moments
%                      info=NaN : no computation, no hard or soft data
%                      info=0   : computation using BME with soft data 
%                      info=1   : dubious results, integration error above tolerance
%                                 when integrating over the soft pdf. Try to increase 
%                                 the value of options(3)
%                      info=3   : computation using simple kriging, no soft data
%                      info=10  : dubious results from integration routine. Try reducing
%                                 the number of soft data points ns, or option(3) or option(4)

%%%%%%% Set some parameters

nh=length(zh);
ns=length(nl);
maxpts=options(3);
aEps=0;
rEps=options(4);
nMom=options(8);

if (nMom~=1) & (nMom~=2) & (nMom~=3) 
  error('Illegal value for nMom');
end;

%
% Case with no soft data points
%
% m1=meank=BkIFhs*zh
% m2=vark=KkIFhs
% stdDev=sqrt(vark)
% m3=0;
% skewCoef=m3/sdtDev^3;

if (ns==0),
  normConstant=1;
  if nh==0
    BMEmean=0;
  else
    BMEmean=BkIFhs*zh;
  end;
  stdDev=sqrt(KkIFhs);
  skewCoef=0;
  info=[0 0 0];
  return;
end;

%
% Case with soft data points
%
% A=Int[ dXs fS(Xs) mvnpdf(Xs|h)]
% m1=meank=1/A * Int[dXs fS(Xs) Bk|hs*Xsh mvnpdf(Xs|h)]
% m2=vark=Kk|sh + 1/A * Int[dXs fS(Xs) (Bk|hs*Xhs-meank)^2 mvnpdf(Xs|h)]
% stdDev=sqrt(vark)

BsMean(1,1)=1;
if (nh==0),
  msIFh=zeros(ns,1);
  BsMean(1,2)=0;
else
  msIFh=BsIFh*zh;
  BsMean(1,2)=BkIFhs(:,1:nh)*zh;
end;

As=zeros(ns,1);
As(1:ns,2)=BkIFhs(nh+1:nh+ns)';
Pmean=[1,1];nMomMean=2;
[Val Err fInfo]=mvMomVec(softpdftype,nl,limi,probdens,msIFh,KsIFh,nMomMean,As,BsMean,Pmean,maxpts,aEps,rEps);
if (fInfo(end)~=0)
  disp(sprintf('momentsFun warning: fInfo=%d, rError for mean=%9.3g',fInfo(end),max(abs(Err./Val))));
end;

normConstant=Val(1);
BMEmean=Val(2)/normConstant;
info(1)=fInfo(end);

if (nMom==1),
  stdDev=NaN; 
  skewCoef=NaN;  
  return; 
end;

As=kron(As(:,2),ones(1,nMom-1));


Bs(1,1)=BsMean(2)-BMEmean;
P(1)=2;

if nMom==3,
  Bs(1,2)=Bs(1,1);
  P(1,2)=3;
end;

[Val Err fInfo]=mvMomVec(softpdftype,nl,limi,probdens,msIFh,KsIFh,nMom-1,As,Bs,P,maxpts,aEps,rEps);
if (fInfo(1)~=0)
  disp(sprintf('momFun warning: fInfo=%d, rError for moments 2 and 3=%9.3g',fInfo(1),max(Err./Val)));
end;
stdDev=sqrt(KkIFhs+Val(1)/normConstant);
info(2)=fInfo(1);

if nMom==2,
  skewCoef=NaN; 
  info(3)=NaN;  
else
  skewCoef=(Val(2)/normConstant)/stdDev^3;
  info(3)=fInfo(end);
end;

