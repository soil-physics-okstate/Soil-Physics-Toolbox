function [limioffset]=probaoffset(softpdftype,nl,limi,offset)

% probaoffset               - Offset the probabilistic soft data (Jan 1, 2001)
%
% Offsets the limi argument of some probabilistic soft data. 
%
% SYNTAX : 
%
% [limioffset]=probaoffset(softpdftype,nl,limi,offset)
%
% INPUT :
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
% offset      ns by 1     vector specifying the offset for each row of limi.
%                         Note: offset may be a scalar, in which case the
%                         same offset is used for each row
%
% OUTPUT :
%
% limioffset  ns by l     matrix of same size as limi with the new
%                         offseted limi values

ns=length(nl);
if length(offset)==1
  offset=offset*ones(ns,1);
end;

switch softpdftype
case {1,2},
  limioffset=limi+kron(offset,ones(1,size(limi,2)));
case {3,4},
  limioffset=limi;
  limioffset(:,[1 3])=limi(:,[1 3])+kron(offset,[1 1]);
otherwise
  error('wrong value for softpdftype');
end



