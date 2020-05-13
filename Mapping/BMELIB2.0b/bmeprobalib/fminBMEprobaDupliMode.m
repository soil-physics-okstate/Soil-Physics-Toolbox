function [f,INFOINTEG]=fminBMEprobaDupliMode(zk,zh,softpdftype,nl,limi,probdens,...
  BkIFh,invKkIFh,BsIFkh,KsIFkh,options,nlest,limiest,probdensest);

% fminBMEprobaDupliMode     - fmin subroutine for BMEprobaMode.m when there is soft data at est pt (Jan 1, 2001)
%
% Objective function used by the function BMEprobaMode for the BME  
% case in order to identify the mode of the pdf when there is a duplicate
% soft data at the estimation point
%
% SYNTAX :
%
% [f,INFOINTEG]=fminBMEprobaDupliMode(zk,zh,softpdftype,nl,limi,probdens,...
%   BkIFh,invKkIFh,BsIFkh,KsIFkh,options,nlest,limiest,probdensest);
%
% INPUT :
%
% zk         scalar       value of variable at estimation point for which the posterior pdf 
%                         is computed 
% zh         nh by 1      vector of conditioning hard values
%                         NOTE:  z2=[zk;zh] is a nk+nh by 1 vector of conditioning values
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
% BkIFh      nk by 1      vector equal to Kk_h*inv(Kh), where
%                         Kk_h and Kh are covariance matrices
% invKkIFh   nk by nk     matrix for the inverse of KkIFh, where KkIFh is the 
%                         covariance matrix for zk given zh
% BsIFkh     ns by 1      vector equal to Ks_kh*inv(Kkh), where
%                         Ks_kh and Kkh are covariance matrices
% KsIFkh     ns by ns     matrix of conditional covariances for zs given zkh
% options    vector of parameters as described in BMEprobaMode.m
% nlest,limiest and probdensest define the soft data at the estimation point
%
% OUTPUT :
%
% f          scalar       a value of -log(pdf) up to some constants
% INFOINTEG  scalar       the information returned by the Fortran integration
%                         subroutine mvProAG2

%softpdftypeCheckArgs(softpdftype,nl,limi,probdens);
%warning('need to remove softPdfTypeCheckArgs test for speed');

global INFOINTEG;      % declares INFOTEG as global, so it can be used
				         % in the main program
                     
                     
if length(nl)>0        % Check to see if there are soft data point
  maxpts=options(3);   %   if there are soft data point, calculate the integral
  aEps=0;
  rEps=options(4);
  zkh=[zk;zh];
  msIFkh=BsIFkh*zkh;             % compute the conditional mean
  [P Perr INFOINTEG]=mvPro(softpdftype,nl,limi,probdens,msIFkh,KsIFkh,maxpts,aEps,rEps);
  Pa=max([P,1e-323]);            % make sure that log(P) does not yield -Inf
  logPa=log(Pa);
else
  INFOINTEG=0;
  logPa=0;             %   however if there are no soft data, assign logPa=0
end;

if length(zh)>0
  mkIFh=BkIFh*zh;
else
  mkIFh=0;
end;
yk=zk-mkIFh;
spdf=proba2val(zk,softpdftype,nlest,limiest,probdensest);
spdfa=max([spdf,1e-323]);
f=0.5*(yk'*invKkIFh*yk)-logPa-log(spdfa);   % compute the value of the -log posterior pdf

if (INFOINTEG~=0)
   disp(sprintf('fminBMEprobaDupliMode warning: fInfo=%d, P=%9.3g,Perr=%9.3g',...
      INFOINTEG,P,Perr));
end;