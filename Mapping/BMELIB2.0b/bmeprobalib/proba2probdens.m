function [probdens,norm]=proba2probdens(softpdftype,nl,limi,probdenstemplate);

% proba2probdens            - Normalizes the probability density function (Jan 1, 2001)
%
% Calculates the norm (area under the curve) of a function template, and returns
% a normalized pdf. This function uses the syntax of probabilistic data (see 
% probasyntax) to define the pdf.
%
% SYNTAX :
%
% [probdens,norm]=proba2probdens(softpdftype,nl,limi,probdenstemplate);
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
% probdenstemplate        ns by p matrix of non-normalized probability density values,  
%                         where p is equal to either max(nl)-1 or max(nl), depending on 
%                         the softpdftype. probdenstemplate(i,:) are the values of the  
%                         non-normalized probability density corresponding to the intervals  
%                         for the i-th soft data defined in limi(i,:). (see probasyntax for 
%                         more explanations)
%
% OUTPUT :
%
% probdens    ns by p     matrix of normalized probability density values,
%                         each row of probadens is equal to the corresponding row
%                         of probadenstemplate divided by its normalization constant
% norm        ns by 1     vector of the normalization constants (area under the curve)
%                         of each row of probadenstemplat, e.g. norm(is) is the
%                         normalization constant for probdenstemplate(is,:)

probdens=probdenstemplate;
ns=length(nl);
for is=1:ns
  switch softpdftype
  case 1,
    norm(is,1)=sum(probdens(is,1:nl(is)-1).*diff(limi(is,1:nl(is))));
    probdens(is,1:nl(is)-1)=probdens(is,1:nl(is)-1)/norm(is,1);
  case 2,
    norm(is,1)=sum(0.5*(probdens(is,1:nl(is)-1)+probdens(is,2:nl(is)))...
           .*diff(limi(is,1:nl(is))));
    probdens(is,1:nl(is))=probdens(is,1:nl(is))/norm(is,1);
  case 3,
    norm(is,1)=sum(probdens(is,1:nl(is)-1)*limi(is,2));
    probdens(is,1:nl(is)-1)=probdens(is,1:nl(is)-1)/norm(is,1);
  case 4,
    norm(is,1)=sum(0.5*(probdens(is,1:nl(is)-1)+probdens(is,2:nl(is)))...
      *limi(is,2));
    probdens(is,1:nl(is))=probdens(is,1:nl(is))/norm(is,1);
  end
end

