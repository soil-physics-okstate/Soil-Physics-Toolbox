function    [Val,Err,info]=uvProNR(softpdftype,nl,limi,probdens,Mean,C,aEps,rEps);

% uvProNR                   - Uni-Variate integral arising from calculation of probability (Jan 1, 2001)
%
% Uses a numerical integration to calculate a univariate integral arising
% from the calculation of a probability in statistics. The form of the integral
% to calculate is as follow
%  b
% I dXs fs(Xs) uvnpdf(Xs,Mean,Cov),
%  a                 
% where Xs is a variable, fs(.) is a univariate function,
% and mvnpdf(.) is a uni-variate normal pdf with mean vector 
% Mean, and covariance matrix Cov. The integration limits are given by 
% a and b. As and bs are user defined variables, while p is a value with the 
% order of the moment to calculate. Note that fs is a pdf defined 
% by the soft probabilistic data variables softpdftype, nl, limi and probdens
% (see probasyntax for more explanation on soft probabilistic data)
%
% A simple interative trapezoidal method is used (see, for example,
% "numerical recipies").
%
% SYNTAX :
%
% [Val,Err,Info]=uvProNR(softpdftype,nl,limi,probdens,Mean,C,aEps,rEps);
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
%                         NOTE THAT FOR unProNR, N=1 (univariate case).
% limi        N by l      matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         dimension of fs. (see probasyntax for more explanations)
%                         NOTE THAT FOR unProNR, N=1 (univariate case).
% probdens    N by q      matrix of probability density values fs, where q is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th dimension 
%                         defined in limi(i,:). (see probasyntax for more explanations)
%                         NOTE THAT FOR unProNR, N=1 (univariate case).
% Mean        N by 1      vector of ns means    
%                         NOTE THAT FOR unProNR, N=1 (univariate case).
% C           N by N      covariance matrix
%                         NOTE THAT FOR unProNR, N=1 (univariate case).
% maxpts      scalar      maximum number of function values allowed.
% aEps        scalar      absolute error tolerance.
% rEps        scalar      relative error tolerance.
%
% OUTPUT :
%
% Value       scalar      estimated value for the integral
% Error       scalar      estimated error for the integral
% info        scalar      end of computation status
%                         info = 0 for normal exit, when Error <=  aEps or
%                                  Error <=  abs(Value)*rEps 
%                         info = 1 too many iterations necessary to obtain 
%                                  the required accuracy.
%                         info = 2 if softpdftype is not equal to 1, 2, 3, or 4

error('It seems that you did not compile the mvnlib directory (see mvnlibcompile)');
