function FyTransformCheckArgs(yfile,Fyfile,ck,ch,cs,yh,ysoftpdftype,ynl,ylimi,yprobdens);

% FyTransformCheckArgs      - Check arguments for the BMEprobaT functions (Jan 1, 2001)
%
% Checks the inputs arguments for data tranformation Y=g(Z) in the BMEprobaT
% functions (such as BMEprobaTMode). This function checks that the cdf for 
% the Y variable is properly described using yfile and Fyfile, and that the
% cdf covers a sufficient range of Y values to cover all the hard and soft data
% 
% Returns an error message if a problem is found.
%
% SYNTAX :
%
% FyTransformCheckArgs(yfile,Fyfile,ck,ch,cs,yh,ysoftpdftype,ynl,ylimi,yprobdens);
%
% INPUT : Similar to those of the BMEprobaT functions (such as BMEprobaTMode).
%
% yfile      k by 1   vector of values for Y sorted in ascending order
% Fyfile     k by 1   vector of values of the cumulative distribution for Y
% ck         nk by p  matrix of coordinates for the estimation point,
%                     where d is the dimension of the space
% ch         nh by d  matrix of coordinates for the hard data 
% cs         ns by d  matrix of coordinates for the soft data
% yh         nh by 1  vector of the values for the hard data for Y
% ysoftpdftype scalar integer for the type of soft pdf, representing the 
%                     probabilitic soft data for Y. These soft pdf types are as follow:
%                     1-Histogram, 2-Linear, 3-Grid histogram, 4-Grid Linear. 
%                     (see also probasyntax)
% ynl        ns by 1  vector of number of interval limits. nl(i) is 
%                     the number of interval limits used to define the soft pdf 
%                     for soft data point i. (see also probasyntax)
% ylimi      ns by l  matrix of interval limits, where l is equal to
%                     either max(nl) or 3 depending of the softpdftype.
%                     limi(i,:) are the limits of intervals for the i-th 
%                     soft data. (see also probasyntax)
% yprobdens= ns by p  matrix of probability density values, where p is 
%                     equal to either max(nl)-1 or max(nl), depending on the 
%                     softpdftype. probdens(i,:) are the values of the probability 
%                     density corresponding to the intervals for the i-th soft data 
%                     defined in limi(i,:). (see also probasyntax)
%
% NOTE :
%
% 2- It is possible to specify an additional index for ck,
% ch and cs, taking integer values from 1 to nv. This index
% specifies which variable is known at each coordinate.
% In that case, ck, ch, cs are cell arrays, where the first
% cell is the matrix of coordinates and the second cell is
% the vector of index values. 
% If an index is specified, both yfile and Fyfile are cell
% arrays of size nv, where each cell corresponds to the
% definition of a specific cumulative distribution function.
% The nhmax and nsmax variables are then vectors of nv
% elements, where sum(nsmax)<=20.

%%%% Find every variables and call the checking function for each variable

noindex=~iscell(ck);
if noindex
  noindexCheckArgs(yfile,Fyfile,yh,ysoftpdftype,ynl,ylimi,yprobdens);
else
  indexk=ck{2};
  indexh=ch{2};
  indexs=cs{2};
  nv=max([indexk;indexh;indexs]);
  if ~iscell(yfile) | ~iscell(Fyfile)
    error('If ck is a cell, then yfile and Fyfile must be cells as well');
  end;
  if length(yfile)~=nv | length(Fyfile)~=nv
    error('yfile and Fyfile must be cell arrays of length nv');
  end;
  for i=1:nv,
    [csi,c2,ynli,ylimii,yprobdensi,nl2,limi2,probdens2]=...
      probasplit(cs,ysoftpdftype,ynl,ylimi,yprobdens,find(indexs==i));
    noindexCheckArgs(yfile{i},Fyfile{i},yh(indexh==i),...
      ysoftpdftype,ynli,ylimii,yprobdensi);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The checking of the arguments start here for each variables

function noindexCheckArgs(yfile,Fyfile,yh,ysoftpdftype,ynl,ylimi,yprobdens);

%%%% No checking necessary if yfile is empty 
if isempty(Fyfile)
  return;
end;

%%%% No checking necessary if there is no hard nor soft data
nh=length(yh);
ns=length(ynl);
if nh+ns==0
  return;
end

%%%% Check the size of yfile and Fyfile
if size(yfile,2)~=1, error('yfile must be a column vector'); end;
if size(Fyfile,2)~=1, error('Fyfile must be a column vector'); end;
if length(yfile)~=length(Fyfile), 
  error('yfile and Fyfile must be vectors of same length');
end;

%%%% Check that yfile and Fyfile define a proper cdf
sortyfile=sort(yfile);
if sum(sortyfile==yfile)~=length(yfile)
  error('yfile file must be a vector of sorted increasing values');
end;
diffFyfile=diff(Fyfile);
if sum(diffFyfile==0)>0
  error('Fyfile cannot have equal consecutive cdf values')
end;
if sum(diffFyfile<0)>0,
  error('The cdf values in Fyfile must be increasing');
end;
if sum(Fyfile<0 | Fyfile>1)>0,
  error('Fyfile must values between 0 and 1'); 
end;
if Fyfile(1)>0.001, error('Fyfile(1) must be smaller than 0.001'); end;
if Fyfile(end)<0.999, error('Fyfile(end) must be greater than 0.999'); end;

%%%% Make sure that the range covered by yfile includes all hard and soft data
if ns>0,
  [ysmin,ysmax]=proba2interval(ysoftpdftype,ynl,ylimi);
end;
if nh==0,
  ymin=ysmin;
  ymax=ysmax;
elseif ns==0,
  ymin=min(yh);
  ymax=max(yh);
else
  ymin=min([yh;ysmin]);
  ymax=max([yh;ysmax]);
end
if yfile(1)>ymin
  error('There is a hard or soft data that is smaller than the min yfile value'); 
end;
if yfile(end)<ymax
  error('There is a hard or soft data that is greater than the max yfile value'); 
end;
