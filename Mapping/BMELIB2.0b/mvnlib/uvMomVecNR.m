function [Val,Err,Info]=uvMomVecNR(softpdftype,nl,limi,probdens,Mean,C,As,Bs,P,aEps,rEps);

% uvMomVecNR                - Uni-Variate integrals arising from calculation of Moments (Jan 1, 2001)
%
% Uses a numerical integration to calculate a univariate integral arising
% from the calculation of moments in statistics. The form of the integral
% to calculate is as follow
%  b
% I dXs fs(Xs) (As'*Xs+Bs).^P mvnpdf(Xs,Mean,Cov),
%  a
% where Xs is a variable, fs(.) is a univariate function,
% and mvnpdf(.) is a uni-variate normal pdf with mean vector 
% Mean, and covariance matrix Cov. The integration limits are given by 
% a and b. As and Bs are user defined variables, while P is a vector with the 
% order of the moments to calculate. Note that fs is a pdf defined 
% by the soft probabilistic data variables softpdftype, nl, limi and probdens
% (see probasyntax for more explanation on soft probabilistic data)
% This function is vectorised in that there are nMom sets of 
% As, bs and p coeffients defined in As, Bs and P. 
%
% A simple interative trapezoidal method is used (see, for example,
% "numerical recipies").
%
% SYNTAX :
%
%   [Val,Err,Info]=uvMomVecNR(softpdftype,nl,limi,probdens,Mean,C,As,Bs,P,aEps,rEps);
%
% INPUT :
%
% softpdftype scalar      indicates the type of soft pdf representing the  
%                         function fs.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
% nl          N by 1      vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for the i-th 
%                         dimention of fs. (see probasyntax for more explanations)
%                         NOTE THAT FOR unMomVecNR, N=1 (univariate case).
% limi        N by l      matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         dimension of fs. (see probasyntax for more explanations)
%                         NOTE THAT FOR unMomVecNR, N=1 (univariate case).
% probdens    N by q      matrix of probability density values fs, where q is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th dimension 
%                         defined in limi(i,:). (see probasyntax for more explanations)
%                         NOTE THAT FOR unMomVecNR, N=1 (univariate case).
% Mean        N by 1      vector of ns means    
%                         NOTE THAT FOR unMomVecNR, N=1 (univariate case).
% C           N by N      covariance matrix
%                         NOTE THAT FOR unMomVecNR, N=1 (univariate case).
% As          N by nMom   matrix of coeffients
% Bs          1 by nMom   vector of nMom scalar coefficients bs
% P           1 by nMom   vector of nMom power coefficients p
% maxpts      scalar      maximum number of function values allowed.
% aEps        scalar      absolute error tolerance.
% rEps        scalar      relative error tolerance.
%
% OUTPUT :
%
% Val         1 by nMom   estimated values for the univariate integrals
% Err         1 by nMom   estimated absolute error.
% Info        1 by nMom   vector of end of computation status
%                         info = 0 for normal exit, when Error <=  aEps or
%                                  Error <=  abs(Value)*rEps 
%                         info = 1 too many iterations necessary to obtain 
%                                  the required accuracy.
%                         info = 2 if softpdftype is not equal to 1, 2, 3, or 4
%                         info = 2 if softpdftype is not equal to 1, 2, 3, or 4

nMom=length(As);
if length(Bs)~=nMom
  error('As and Bs must have same length');
end
if length(P)~=nMom
  error('As and P must have same length');
end

for i=1:nMom
  [Val(i) Err(i) Info(i)]=uvMomNR(softpdftype,nl,limi,probdens,Mean,C,...
    As(i),Bs(i),P(i),aEps,rEps);
end