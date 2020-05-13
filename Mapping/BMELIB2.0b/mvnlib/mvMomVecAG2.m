function [Val,Err,Info]=mvMomVecAG2(softpdftype,nl,limi,probdens,Mean,C,...
  nMom,As,Bs,P,maxpts,aEps,rEps);

% mvMomVecAG2               - Multi-Variate integrals arising from calculation of Moments (Jan 1, 2001)
%
% Uses a numerical integration to calculate a multivariate integral arising
% from the calculation of moments in statistics. The form of the integral
% to calculate is as follow
%  b
% I dXs fs(Xs) (As'*Xs+Bs).^P mvnpdf(Xs,Mean,Cov),
%  a
% where Xs is a vector of N variables, fs(.) is a N-dimentional function,
% and mvnpdf(.) is a N-dimentional multi-variate normal pdf with mean vector 
% Mean, and covariance matrix Cov. The integration limits are given by 
% a and b. As and Bs are user defined variables, while P a vector with the 
% order of the different moments to calculate. Note that fs is a pdf defined 
% by the soft probabilistic data variables softpdftype, nl, limi and probdens
% (see probasyntax for more explanation on soft probabilistic data)
% This function is vectorised in that there are nMom sets of 
% As, bs and p coeffients defined in As, Bs and P. 
%
% If the number of dimensions for the integral must be greater or equal to 2,
% The AG2 (Adaptative Genz 2) algorithm used is the subregion 
% Adaptative integration program dcuhre for multiple integral. See:
% Berntsen, J., O. Espelid and A. Genz, "An Adaptive Multidimensional
% Integration Routine for a Vector of Integrals".
%
% SYNTAX :
%
% [Val,Err,Info]=mvMomVecAG2(softpdftype,nl,limi,probdens,Mean,C,...
%   nMom,As,Bs,P,maxpts,aEps,rEps);
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
% limi        N by l      matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         dimension of fs. (see probasyntax for more explanations)
% probdens    N by q      matrix of probability density values fs, where q is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th dimension 
%                         defined in limi(i,:). (see probasyntax for more explanations)
% Mean        N by 1      vector of N means    
% C           N by N      covariance matrix
% nMom        scalar      number of moments to calculate
% As          N by nMom   matrix of coeffients
% Bs          1 by nMom   vector of nMom scalar coefficients bs
% P           1 by nMom   vector of nMom power coefficients p
% maxpts      scalar      maximum number of function values allowed.
% aEps        scalar      absolute error tolerance.
% rEps        scalar      relative error tolerance.
%
% OUTPUT :
%
% Val         nMom by 1  estimated values for the multiple integrals
% Err         nMom by 1  estimated absolute error, with 99% confidence level.
% Info        scalar      Info = 0, normal completion with ERR < EPS;
%                         Info = 1, completion with ERR > EPS and MAXPTS function 
%                             values used; increase MAXPTS to decrease ERROR;
%                         Info = 2 to 11: see IFAIL in source code dcuhre.f                    

error('It seems that you did not compile the mvnlib directory (see mvnlibcompile)');
