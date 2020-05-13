function [mlogfy,INFOINTEG]=fminBMEprobaTMode(xk,xh,xsoftpdftype,xnl,xlimi,...
  xprobdens,BkIFh,invKkIFh,BsIFkh,KsIFkh,options,zfile,dgfile,zmk);

% fminBMEprobaTMode         - fmin subroutine for BMEprobaTMode.m (Jan 1, 2001)
%
% Objective function used by the function BMEprobaTMode
% for the BME case in order to identify the mode of the pdf
% on the original scale.
% This function calculates the -log(fY(x)), where X is a zero mean 
% multivariate Gaussian stationary/homogeneous S/TRF, Z(p)=X(p)+mZ(p) 
% has the Expected value mZ=E[Z] added, and Y(p)=g(Z(p)) is a transform 
% of Z using the one to one function g(.).
%
% SYNTAX :
%
% [mlogfy,INFOINTEG]=fminBMEprobaTMode(xk,xh,xsoftpdftype,xnl,xlimi,...
%   xprobdens,BkIFh,invKkIFh,BsIFkh,KsIFkh,options,zfile,dgfile,zmk);
%
% INPUT :
%
% xk          scalar    value of the variable X for which the minus log of the 
%                       Y posterior pdf is computed 
% xh          nh by 1   vector of conditioning hard values for X
%                       NOTE:  z2=[zk;zh] is a nk+nh by 1 vector of conditioning values
% xsoftpdftype scalar   integer for the type of soft pdf, representing the 
%                       probabilitic soft data for X. These soft pdf types are as follow:
%                       1-Histogram, 2-Linear, 3-Grid histogram, 4-Grid Linear. 
%                       (see also softPdfType)
% xnl         ns by 1   vector of number of interval limits for X. nl(i) is 
%                       the number of interval limits used to define the soft pdf 
%                       for soft data point i. (see also softPdfType)
% xlimi       ns by l   matrix of interval limits for X, where l is equal to
%                       either max(nl) or 3 depending of the softpdftype.
%                       limi(i,:) are the limits of intervals for the i-th 
%                       soft data. (see also softPdfType)
% xprobdens   ns by p   matrix of probability density values for X, where p is 
%                       equal to either max(nl)-1 or max(nl), depending on the 
%                       softpdftype. probdens(i,:) are the values of the probability 
%                       density corresponding to the intervals for the i-th soft data 
%                       defined in limi(i,:). (see also softPdfType)
% BkIFh      nk by 1    vector equal to Kk_h*inv(Kh), where
%                       Kk_h and Kh are covariance matrices for X
% invKkIFh   nk by nk   matrix for the inverse of KkIFh, where KkIFh is the 
%                       covariance matrix for xk given xh for X
% BsIFkh     ns by 1    vector equal to Ks_kh*inv(Kkh), where
%                       Ks_kh and Kkh are covariance matrices
% KsIFkh     ns by ns   matrix of conditional covariances for xs given xkh
% options    vector of parameters as described in BMEprobaMode.m
% zfile and dgfile are two optional vectors of same size defining the 
%                      derivative dg=dy/dz of the transformation function y=g(z)
%                      the default value is [] (identity transformation)
% zmk        scalar    optional input specifing the expected value for Z, 
%                      i.e. mzk=E[Zk], such that Z=X+zmk.
%
% OUTPUT :
%
% mlogfy      scalar   a value of -log(fY(x)) up to some constants
% INFOINTEG   scalar   the information returned by the Fortran integration
%                      subroutine mvProAG2

%Note: The following line performing error checking of input arguments
%      was removed to speed up the calculation
%softpdftypeCheckArgs(xsoftpdftype,xnl,xlimi,xprobdens);  

global INFOINTEG;      % declares INFOTEG as global, so it can be used
				         % in the main program
                     
if nargin<13
  zfile=[];
  dgfile=[];
end
if nargin<14
  zmk=0;
end

if length(xnl)>0       % Check to see if there are soft data point
  maxpts=options(3);   %   if there are soft data point, calculate the integral
  aEps=0;
  rEps=options(4);
  xkh=[xk;xh];
  xmsIFkh=BsIFkh*xkh;  %        compute the conditional mean
  [P Perr INFOINTEG]=mvPro(xsoftpdftype,xnl,xlimi,xprobdens,...
    xmsIFkh,KsIFkh,maxpts,aEps,rEps);
  Pa=max([P,1e-323]);  %        make sure that log(P) does not yield -Inf
  logPa=log(Pa);
else
  INFOINTEG=0;
  logPa=0;             %   however if there are no soft data, assign logPa=0
end;
  
if length(xh)>0
  xmkIFh=BkIFh*xh;
else
  xmkIFh=0;
end;
xka=xk-xmkIFh;
mlogfx=0.5*(xka'*invKkIFh*xka)-logPa;   % compute the value of the -log posterior pdf

if (INFOINTEG~=0)
   disp(sprintf('fminBMEprobaTMode warning: fInfo=%d, P=%9.3g,Perr=%9.3g',...
      INFOINTEG,P,Perr));
end;

dotransform=~isempty(zfile);
if dotransform
  zk=xk+zmk;
  dg=interp1(zfile,dgfile,zk,'linear');
  mlogfy=mlogfx+log(dg);
else
  mlogfy=mlogfx;
end;
