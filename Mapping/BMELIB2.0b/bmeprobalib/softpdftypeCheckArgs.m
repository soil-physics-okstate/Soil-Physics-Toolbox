function softpdftypeCheckArgs(softpdftype,nl,limi,probdens)

% softpdftypeCheckArgs      - Check soft pdf arguments (Jan 1, 2001)
%
% This function checks the arguments defining a soft pdf for probabilisitic
%  soft data. If an error in the arguments is function, an error message is
%  generated explaining what the problem is.
%
% SYNTAX :
%
% softpdftypeCheckArgs(softpdftype,nl,limi,probdens)
%
% INPUT :
%
% softpdftype scalar   integer for the type of soft pdf, representing the 
%                      probabilitic soft data. These soft pdf types are as follow:
%                      1-Histogram, 2-Linear, 3-Grid histogram, 4-Grid Linear. 
%                      (see also probasyntax)
% nl         ns by 1   vector of number of interval limits. nl(i) is 
%                      the number of interval limits used to define the soft pdf 
%                      for soft data point i. (see also probasyntax)
% limi       ns by l   matrix of interval limits, where l is equal to
%                      either max(nl) or 3 depending of the softpdftype.
%                      limi(i,:) are the limits of intervals for the i-th 
%                      soft data. (see also probasyntax)
% probdens=  ns by p   matrix of probability density values, where p is 
%                      equal to either max(nl)-1 or max(nl), depending on the 
%                      softpdftype. probdens(i,:) are the values of the probability 
%                      density corresponding to the intervals for the i-th soft data 
%                      defined in limi(i,:). (see also probasyntax)

%%%%% Check the types of the arguments
if ~isnumeric(softpdftype), error('softpdftype must be a numeric integer'); end;
if ~isnumeric(nl) error('nl must be a numeric vector'); end;  
if ~isnumeric(limi) error('limi must be a numeric matrix'); end; 
if ~isnumeric(probdens) error('probdens must be a numeric matrix'); end;

%%%%% Check size of softpdftype
if length(softpdftype)~=1, error('softpdftype must be a 1 by 1 integer'); end;
if softpdftype~=1 & softpdftype~=2 & softpdftype~=3 & softpdftype~=4 
  error('softpdftype must be equal to 1, 2, 3 or 4');
end;

%%%%% Check size of nl and the number of rows of limi and probdens
ns=size(nl,1);
if ns>0 & size(nl,2)~=1, error('nl must be ns by 1 vector'); end;
if size(limi,1)~=ns | size(probdens,1)~=ns
   error('nl, limi and probdens must have the same number of rows');
end;

%%%%% Check size of limi
if size(limi,1)~=ns error('cs and limi must have the same number of rows'); end;
if ns>0
  nlMax=max(nl);
  if ndims(limi)>2, error('length(size((limi)) cannot be more than 2'); end;
  sl=size(limi);
  switch softpdftype,
    case {1,2}
      if ~prod(double(sl==[ns nlMax])), error('limi must be a ns by max(nl) matrix'); end;
    case {3,4}
      if ~prod(double(sl==[ns 3])), error('limi must be a ns by 3 matrix'); end;
  end;
end;  
  
%%%%% Check size of probdens
if size(probdens,1)~=ns error('cs and probdens must have the same number of rows'); end;
if ns>0
  if ndims(probdens)>2, error('length(size((probdens))  cannot be more than 2'); end;
  sp=size(probdens);  
  switch softpdftype,
    case {1,3}
      if ~prod(double(sp==[ns nlMax-1])), error('probadens must be a ns by max(nl)-1 matrix'); end;
    case {2,4}
      if ~prod(double(sp==[ns nlMax])), error('probadens must be a ns by max(nl) matrix'); end;
  end;
end;

%%%%% If softpdftype>2, make sure that nl is consistent with limi(1):limi(2):limi(3)
if softpdftype>2
  for is=1:ns
    nlcalc=length(limi(is,1):limi(is,2):limi(is,3));
    if nl(is)~=nlcalc, 
      error(sprintf('row %d of limi requires that nl(%d)=%d',is,is,nlcalc));
    end
  end
end
