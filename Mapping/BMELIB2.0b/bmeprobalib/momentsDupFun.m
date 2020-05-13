function [BMEmean,stdDev,skewCoef,info]=momentsDupFun(zh,softpdftype,nl,limi,probdens,...
  options,BksIFh,KksIFh,nlest,limiest,probdensest)

% momentsDupFun             - Moments of the BME posterior pdf with soft data at est pt (Jan 1, 2001)
%
%   Calculates the mean , the standard deviation and the coefficient 
%   of skewness of the BME posterior pdf when there is a soft data 
%   at the estimation point.
%   The number of moments calculated is controled by nMom=options(8) (1 for 
%   mean, 2 for mean and sdt dev, 3 for all three moments)
%
% SYNTAX :
%
% [BMEmean,stdDev,skewCoef,info]=momentsDupFun(zh,softpdftype,nl,limi,probdens,...
%   PDFparam,BksIFh,KksIFh,nlest,limiest,probdensest) 
%
% INPUT :
%
% zh          nh by 1     vector of the values for the hard data
% softpdftype scalar      integer for the type of soft pdf, representing the 
%                         probabilitic soft data. These soft pdf types are as follow:
%                         1-Histogram, 2-Linear, 3-Grid histogram, 4-Grid Linear. 
%                         (see also probasyntax)
% nl, limi, probdens      are matrices with ns rows representing the soft
%                         probablistic data at the ns soft data points
%                         (see probasyntax for more explanation on soft data)
% options     1 by k       vector of optional parameters, where :
%                         options(3) is the maximum number of evaluation of the integral.
%                         The default value is 50000
%                         options(4) is the relative error on the estimation of the integral.
%                         The default value is equal to 1e-4 
%                         options(8) number of moments to calculate (1, 2 or 3)
%                          =1 to calculate the mean of the BME posterior pdf,
%                          =2 for the mean and standard deviation,
%                          =3 for the mean, standard deviation and coef of skewness.
% BksIFh      ns+1 by nh  matrix, BksIFh=Kks_h*inv(Kh)
% KksIFh      ns+1 by nh  matrix, KksIFh=BksIFh*Kks_h' is the cov mat of  
%                                   est point and soft data given the hard data
% nlest, limiest and probdensest are matrices with one row representing the soft
%                          probablistic data at the estimation point
%
% OUTPUT :
%
% BMEmean     scalar      the mean of the BME posterior pdf 
% stdDev      scalar      the standard deviation of the BME posterior pdf
% skewCoef    scalar      the coefficient of skewness of the BME posterior pdf
% info        1 by 3      vector of information about the numerical computation of
%                         the 3 moments
%                         info=NaN : no computation, no hard or soft data
%                         info=0   : computation using BME with soft data 
%                         info=1   : dubious results, integration error above tolerance
%                                    when integrating over the soft pdf. Try to increase 
%                                    the value of options(3)
%                         info=3   : computation using simple kriging, no soft data
%                         info=10  : dubious results from integration routine. Try reducing
%                                    the number of soft data points ns, or option(3) or option(4)

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


%%%%%%% Concatenate soft data at the estimation point and at the soft data points

%nlks=[nlest;nl];
%[s1,s2]=size(limiest);   limiks(1,1:s2)=limiest;
%[s1,s2]=size(limi);      limiks(1+(1:s1),1:s2)=limi;
%[s1,s2]=size(probdensest); probdensks(1,1:s2)=probdensest;
%[s1,s2]=size(probdens);    probdensks(1+(1:s1),1:s2)=probdens;
[dummy,nlks,limiks,probdensks]=...
  probacat(softpdftype,nlest,limiest,probdensest,softpdftype,nl,limi,probdens);


%%%%%%%% Calculate the normalization constant and the mean of the BME posterior pdf

if nh==0
  mksIFh=zeros(ns+1,1);
else
  mksIFh=BksIFh*zh;
end;

Bs(1,1)=1;
As(1:ns+1,1)=0;
Bs(1,2)=0;
As(1,2)=1;
As(2:ns+1,2)=0;
P=[1,1];
[Val Err fInfo]=mvMomVec(softpdftype,nlks,limiks,probdensks,mksIFh,KksIFh,...
  length(P),As,Bs,P,maxpts,aEps,rEps);
if (fInfo(end)~=0)
  disp(sprintf('momentsFun warning: fInfo=%d, rError for mean=%9.3g',fInfo(end),max(Err./Val)));
end;

normConstant=Val(1);
BMEmean=Val(2)/normConstant;
info(1)=fInfo(end);

if (nMom==1),
  stdDev=NaN; 
  skewCoef=NaN;  
  return; 
end;

%%%%%%%% Calculate the moments of order 2 and 3 if requested

Bs=[];As=[];P=[];
Bs(1,1)=-BMEmean;
As(1,1)=1;
As(2:ns+1,1)=0;
P(1,1)=2;

if nMom==3,
  Bs(1,2)=-BMEmean;
  As(1,2)=1;
  As(2:ns+1,2)=0;
  P(1,2)=3;
end;

[Val Err fInfo]=mvMomVec(softpdftype,nlks,limiks,probdensks,mksIFh,KksIFh,nMom-1,As,Bs,P,maxpts,aEps,rEps);
if (fInfo(1)~=0)
  disp(sprintf('momFun warning: fInfo=%d, rError for moments 2 and 3=%9.3g',fInfo(1),max(Err./Val)));
end;
stdDev=sqrt(Val(1)/normConstant);
info(2)=fInfo(1);

if nMom==2,
  skewCoef=NaN; 
  info(3)=NaN;  
else
  skewCoef=(Val(2)/normConstant)/stdDev^3;
  info(3)=fInfo(end);
end;

