function [m,v,s,A]=pdfstat(z,pdf);

% pdfstat                   - statistics of a pdf  (Jan 1,2001)
% 
% Compute the mean, the variance and the mediane of a
% distribution based on a discrete definition of a
% univariate probability distribution function.
%
% SYNTAX :
%
% [m,v,s,A]=pdfstat(z,pdf);
%
% INPUT :
%
% z      n by 1   vector of value.
% pdf    n by 1   vector of the probability distribution
%                 function values at the z values.
%
% OUTPUT :
%
% m      scalar   mean of the distribution.
% v      scalar   variance of the distribution.
% s      scalar   skewness of the distribution.
% A      scalar   normalization constant corresponding to a trapezoidal
%                 integration over z of the probability distribution
%                 function values given in pdf. If the probability
%                 distribution is sufficiently finely discretized, A
%                 should be close to the theoretical value of 1.

%%%% If the pdf has some NaN, try to remove them

isnanpdf=isnan(pdf);
if sum(~isnanpdf)<=1 & sum(isinf(pdf))==0
  m=NaN;
  v=NaN;
  s=NaN;
  A=NaN;  
  return;
elseif sum(isnanpdf)>0
  warning('Removing NaN values from the pdf when calculating the moments');
  z=z(~isnanpdf);
  pdf=pdf(~isnanpdf);
end;

%%%%% If the pdf has an Inf value, then the RV has a deterministic value = z(isinf(pdf))

idxInf=find(isinf(pdf));
if length(idxInf)>1, 
  error('pdf must have at most one infinite value'); 
end;
if length(idxInf)==1,
  zDeterministic=z(idxInf);
  m=zDeterministic;
  v=0;
  s=0;
  A=1;
  return;
end;

%%%% Get size, min and max of z

n=length(z);
zmin=min(z);
zmax=max(z);

%%%% Do calculation of moments

[m]=trapezint(z,z.*pdf,zmin,zmax);
if nargout>=2,
  [v]=trapezint(z,((z-m).^2).*pdf,zmin,zmax);
end;
if nargout>=3,
  [s]=trapezint(z,((z-m).^3).*pdf,zmin,zmax)/(sqrt(v)^3);
end;
if nargout>=4,
  [A]=trapezint(z,pdf,zmin,zmax);
end;
