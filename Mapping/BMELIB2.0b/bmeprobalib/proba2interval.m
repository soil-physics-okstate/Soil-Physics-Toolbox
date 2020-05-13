function [a,b]=proba2interval(softpdftype,nl,limi)

% proba2interval            - Transform probabilistic soft data to interval soft data (Jan 1, 2001)
%
% Transform probabilistic soft data to interval soft data. 
%
% SYNTAX :
%
% [a,b]=proba2interval(softpdftype,nl,limi)
%
% INPUT:
%
% softpdftype scalar      indicates the type of soft pdf representing the  
%                         probabilitic soft data.  
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
%
% OUTPUT :
%
% a          ns by 1      vector with the lower bounds of the domains where
%                         the soft probabilistic data is defined
% b          ns by 1      vector with the upper bounds of the domains where
%                         the soft probabilistic data is defined

ns=length(nl);
if ns==0
  a=NaN;
  b=NaN;
else  
  a=limi(:,1);
  switch softpdftype
  case {1,2}, b=limi((1:ns)'+(nl-1)*ns);    
  case {3,4}, b=limi(:,3);
  end;
end;

