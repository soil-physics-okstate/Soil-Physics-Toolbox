function [softmean,softvar]=proba2stat(softpdftype,nl,limi,probdens)

% proba2stat                - Transform probabilistic soft data to mean/variance data (Jan 1, 2001) 
%
% Calculates the mean and variance of some probabilistic soft data
%
% SYNTAX :
%
% [softmean,softvar]=proba2stat(softpdftype,nl,limi,probdens)
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
% probdens    ns by p     matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
%
% OUTPUT :
%
%  softmean  ns by 1      vector of mean values of the soft probabilistic data
%  softvar   ns by 1      vector of variance of the soft probabilistic data

%%%%%% Check the input data

softpdftypeCheckArgs(softpdftype,nl,limi,probdens)

%%%%%% Calculate the mean and variance of the soft probabilistic data

ns=length(nl);

if ns==0
  softmean=[];
  softvar=[];
  return;
end;

for is=1:ns,
  if (softpdftype==1) | (softpdftype==2),
    L1=limi(is,1:nl(is)-1);
    L2=limi(is,2:nl(is));
  elseif (softpdftype==3) | (softpdftype==4),
    L=limi(is,1):limi(is,2):limi(is,3);
    L1=L(1:nl(is)-1);
    L2=L(2:nl(is));
  end;
  P1=probdens(is,1:nl(is)-1);
  if (softpdftype==1) | (softpdftype==3),
    XsMean(is,1)=(1/2)*P1*(L2.^2-L1.^2)';
    Xs2Mean(is,1)=(1/3)*P1*(L2.^3-L1.^3)';
  elseif (softpdftype==2) | (softpdftype==4)
    P2=probdens(is,2:nl(is));
    fsp=(P2-P1)./(L2-L1);
    fso=P1-L1.*fsp;
    XsMean(is,1)= (1/2)*(fso*(L2.^2-L1.^2)')+(1/3)*(fsp*(L2.^3-L1.^3)');
    Xs2Mean(is,1)=(1/3)*(fso*(L2.^3-L1.^3)')+(1/4)*(fsp*(L2.^4-L1.^4)');
  end;
end;
softmean=XsMean;
softvar=Xs2Mean-XsMean.^2;
