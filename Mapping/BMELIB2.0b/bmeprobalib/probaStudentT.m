function [softpdftype,nl,limi,probdens]=probaStudentT(zave,zvar,nobvs)
% probaStudentT      - creates student-t soft probabilistic data  at a set of data points
%
% This program generates the probabilistic soft pdf of the expected value m (i.e. the  
% long run arithmetic average) of a set of observation values at each data point, 
% with data points in space/time vector format. At each data point the estimated 
% arithmetic average zave of nobvs observation values Xi is zave=1/nobvs*Sum(Xi), 
% and the estimated variance is var=1/(nobvs-1)*Sum((Xi-Z)^2).
% The soft probrabilistic distribution is given by the student-t pdf of 
% degree n=nobvs-1, tpdf_n, as follow
% fS(m)=1/sn*tpdf_n( (m-zave) /sn );  where sn=sqrt(var)/sqrt(nobvs)
%
% SYNTAX :
%
% [softpdftype,nl,limi,probdens]=probaStudentT(zave,zvar,nobvs);
%
% INPUT :
%
% zave     ns by 1   vector of the mean value at each of ns soft data point
% zvar     ns by 1   vector of the variance at each of ns soft data point
% nobs     ns by 1   vector of the variance the number of observations at each of ns soft data point
%
% OUTPUT :
%
% softpdftype scalar=2   indicates the type of soft pdf representing the probabilitic soft data.  
% nl          ns by 1    vector of number of limits  nl (see probasyntax)
% limi        ns by l    matrix representing the limits  values (see probasyntax)
% probdens    ns by p    matrix representing the probability densities (see probasyntax)
%
% NOTE : 
% 
% This program needs the statistics toolbox of MATLAB because it uses the 
% tpdf and tinv functions
%
% SEE ALSO probaGaussian, probaUniform, probasyntax.m and stgridsyntax.m 

if exist('tinv')~=2
  error('Function tinv.m missing. You need to have the MATLAB stats toolbox to run probaStudentT')
end


sig=sqrt(zvar);
[ns]=length(zave);

softpdftype=2;
nl=13*ones(ns,1);

limi=NaN*zeros(ns,13);
probdens=NaN*zeros(ns,13);
for is=1:ns
  xnorm=gaussinv([0.00001 0.001 0.01 0.05 0.3 .4 .5 .6 .7 0.95 0.99 0.999 0.99999],[0 1]);
  if nobvs(is)>1
    xnormT=tinv([0.001 0.01 0.02 0.05 0.3 .4 .5 .6 .7 0.95 0.98 0.99 0.999],nobvs(is)-1);
    xnorm=sqrt(xnormT.*xnorm).*sign(xnorm);
  end
  sn=sig(is)/sqrt(nobvs(is));
  x=sn*xnorm+zave(is);
  if nobvs(is)>1
    f=tpdf((x-zave(is))/sn,nobvs(is)-1)/sn;
  else
    f=gausspdf(x,[zave(is) sn^2]);   % because the student-t for degree 0 does not exist
  end
  limi(is,1:13)=x;
  probdens(is,1:13)=f;
  if mod(is,500)==0, disp(sprintf('probaStudentT is %d / %d = %.1f%% done',is,ns,100*is/ns)); end
end

norm=sum(0.5*(probdens(:,1:end-1)+probdens(:,2:end)).*diff(limi,1,2),2);
probdens=probdens./kron(norm,ones(size(xnorm)));
