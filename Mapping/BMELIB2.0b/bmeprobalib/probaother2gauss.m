function [zsoftpdftype,znl,zlimi,zprobdens]=probaother2gauss(ysoftpdftype,ynl,ylimi,...
  yprobdens,yfile,Fyfile);

% probaother2gauss          - Transform soft data for arbitrary RV Y to soft data for Gaussian RV Z (Jan 1, 2001)
%
% Transforms the soft probabilitic data for an arbitrary random
% variable Y to the corresponding soft probabilistic data for 
% the random variable Z, where Y and Z are related by the one to
% one function Y=g(Z). 
% The one-to-one function g(.) is such that it transforms a zero 
% mean unit variance random variable Z to a random variable Y=g(Z) 
% with an arbitrary cumulative distribution function F(Y). Hence the 
% function g(.) is completely defined by the general knowledge (at the 
% prior stage) of the cdf F(Y) for the random variable Y.
% 
% SYNTAX :
%
% [zsoftpdftype,znl,zlimi,zprobdens]=probaother2gauss(ysoftpdftype,ynl,ylimi,...
%   yprobdens,yfile,Fyfile);
%
% INPUT :
%
% ysoftpdftype scalar     indicates the type of soft pdf for the variable Y   
%                         representing the probabilitic soft data at the 
%                         coordinates specified in cs.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
% ynl          ns by 1    vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. (see probasyntax for more explanations)
% ylimi        ns by l    matrix of interval limits for variable Y, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         soft data. (see probasyntax for more explanations)
% yprobdens    ns by p     matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
% y            n by 1     vector of values for the Y variable to be transformed
% yfile        k by 1     vector of values for Y sorted in ascending order
% Fyfile       k by 1     vector of values of the cumulative distribution for Y
%                         The values of Fyfile must be strickly increasing, 
%                         between 0 and 1, and Fyfile(1)<0.001 and Fyfile>0.999 so that 
%                         yfile covers most of the range of Y values.
%
% OUTPUT :
%
% zsoftpdftype,znl,zlmi and zprobdens represent the soft probabilistic data for the
%                         variable Z with zero-mean unit variance Gaussian distributed values.
%
% NOTE :
%
% 1- If the value of the cumulative distribution function
% associated with a y-value is outside the definition of
% the cumulative distribution defined by (yfile,Fyfile),
% the returned z-value will be equal to NaN. The minimum
% and maximum possible values for z are determined by the 
% gaussinv.m function.
%
% 3- It is possible to specify an additional index for ynl,
% taking integer values from 1 to nv, which specifies which
% cumulative distribution function should be used. The
% vector ynl is then a cell array with the first cell containing
% the ynl values and the second cell containing the index. This
% index can be a vector of the same size than the ynl vector
% or a single integer if the same cumulative distribution
% function must be used for all the values. If an index is
% specified, both yfile and Fyfile are cell arrays  where
% each of the nv cells corresponds to the definition of a
% specific cumulative distribution function

if ( ~iscell(ynl) & (iscell(yfile) | iscell(Fyfile) ) ) ...
    | ( iscell(ynl) & (~iscell(yfile) | ~iscell(Fyfile) ) ),
  error('ynl, yfile and Fyfile must be either all non-cells, or all cells');
end;

if ~iscell(ynl),
  [zsoftpdftype,znl,zlimi,zprobdens]=scalarproba2gauss(ysoftpdftype,ynl,ylimi,yprobdens,yfile,Fyfile);
else
  index=ynl{2};
  if length(index)==0,
    zsoftpdftype=ysoftpdftype;
    znl=ynl;
    zlimi=ylimi;
    zprobdens=yprobdens;
  elseif length(index)==1,
    [zsoftpdftype,znl,zlimi,zprobdens]=scalarproba2gauss(ysoftpdftype,ynl{1},...
      ylimi,yprobdens,yfile{index},Fyfile{index});
  else
    ns=length(ynl{1});
    nv=max(index);
    for i=1:nv,
      indexi=find(index==i);
      [zsoftpdftype,znli,zli,zpdi]=scalarproba2gauss(ysoftpdftype,ynl{1}(indexi),...
        ylimi(indexi,:),yprobdens(indexi,:),yfile{i},Fyfile{i});
      znl(indexi,1)=znli;
      zlimi(indexi,1:size(zli,2))=zli;
      zprobdens(indexi,1:size(zpdi,2))=zpdi;
    end;
  end;
end;


function [zsoftpdftype,znl,zlimi,zprobdens]=scalarproba2gauss(ysoftpdftype,ynl,ylimi,yprobdens,yfile,Fyfile);

% scalarproba2gauss               - transform arbitrary proba soft data to Gaussian (April 20, 2000)
%
% Transform a soft probabilitic data for a variable with an 
% arbitrary cumulative distribution function to soft 
% probabilistic data for a variable with zero-mean unit 
% variance Gaussian distributed values
%
% This function perform the work when ynl is not a cell

ns=length(ynl);
switch ysoftpdftype
case 1,
  zsoftpdftype=ysoftpdftype;
  znl=ynl;
  zlimi=NaN*ylimi;
  zprobdens=NaN*yprobdens;
  for is=1:ns,
    yl=ylimi(is,1:ynl(is));
    yld=diff(yl);
    ypd=yprobdens(is,1:ynl(is)-1);
    Fyl=interp1(yfile,Fyfile,yl,'linear');
    zl=gaussinv(Fyl,[0 1]);
    zld=diff(zl);
    zpd=yld.*ypd./zld;
    zlimi(is,1:znl(is))=zl;
    zprobdens(is,1:znl(is)-1)=zpd;
  end
case 2,
  dgfile=transformderiv(yfile,Fyfile);
  zsoftpdftype=ysoftpdftype;
  znl=ynl;
  zlimi=NaN*ylimi;
  zprobdens=NaN*yprobdens;
  for is=1:ns,
    yl=ylimi(is,1:ynl(is));
    ypd=yprobdens(is,1:ynl(is));
    dg=interp1(yfile,dgfile,yl,'linear');
    Fyl=interp1(yfile,Fyfile,yl,'linear');
    zl=gaussinv(Fyl,[0 1]);
    zpd=ypd.*dg;
    [zpd,norm]=proba2probdens(zsoftpdftype,znl(is),zl,zpd);
    zlimi(is,1:znl(is))=zl;
    zprobdens(is,1:znl(is))=zpd;
  end
case 3,
  zsoftpdftype=1;
  znl=ynl;
  zlimi=NaN*ones(ns,max(znl));
  zprobdens=NaN*ones(ns,max(znl)-1);
  for is=1:ns,
    yl=ylimi(is,1):ylimi(is,2):ylimi(is,3);
    if length(yl)~=ynl(is), 
      error(['mismatch on row ' num2str(is) ' of ynl and ylimi']); 
    end;
    yld=diff(yl);
    ypd=yprobdens(is,1:ynl(is)-1);
    Fyl=interp1(yfile,Fyfile,yl,'linear');
    zl=gaussinv(Fyl,[0 1]);
    zld=diff(zl);
    zpd=yld.*ypd./zld;
    zlimi(is,1:znl(is))=zl;
    zprobdens(is,1:znl(is)-1)=zpd;
  end
case 4,
  dgfile=transformderiv(yfile,Fyfile);
  zsoftpdftype=2;
  znl=ynl;
  zlimi=NaN*ones(ns,max(znl));
  zprobdens=NaN*ones(ns,max(znl));
  for is=1:ns,
    yl=ylimi(is,1):ylimi(is,2):ylimi(is,3);
    if length(yl)~=ynl(is), 
      error(['mismatch on row ' num2str(is) ' of ynl and ylimi']); 
    end;
    ypd=yprobdens(is,1:ynl(is));
    dg=interp1(yfile,dgfile,yl,'linear');
    Fyl=interp1(yfile,Fyfile,yl,'linear');
    zl=gaussinv(Fyl,[0 1]);
    zpd=ypd.*dg;
    [zpd,norm]=proba2probdens(zsoftpdftype,znl(is),zl,zpd);
    zlimi(is,1:znl(is))=zl;
    zprobdens(is,1:znl(is))=zpd;
  end
end    


