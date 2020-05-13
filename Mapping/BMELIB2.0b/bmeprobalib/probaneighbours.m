function [csub,nlsub,limisub,probdenssub,dsub,nsub,index]=...
  probaneighbours(c0,c,softpdftype,nl,limi,probdens,nmax,dmax);

% probaneighbours           - Radial neighbourhood selection (Jan 1, 2001)
%
% Select a subset of coordinates and probabilistic soft data based on
% their distances to an origine point with coordinate c0
%
% SYNTAX :
%
% [csub,nlsub,limisub,probdenssub,dsub,nsub,index]=...
%   probaneighbours(c0,c,nl,limi,probdens,nmax,dmax);
%
% INPUT :
%
% c0          1 by d      coordinate of origina point, where d is the dimension of the space
% c           n by d      coordinates of n points with soft probabilistic data
% softpdftype scalar      indicates the type of soft pdf representing the  
%                         probabilitic soft data at the coordinates specified in c.  
%                         softpdftype may take value 1, 2, 3 or 4, as follow:
%                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
%                         and 4 for Grid Linear. (see probasyntax for more explanations)
% nl          n by 1      vector of the number of interval limits. nl(i) is the number  
%                         of interval limits used to define the soft pdf for soft data 
%                         point i. (see probasyntax for more explanations)
% limi        n by l      matrix of interval limits, where l is equal to
%                         either max(nl) or 3 depending of the softpdftype.
%                         limi(i,:) are the limits of intervals for the i-th 
%                         soft data. (see probasyntax for more explanations)
% probdens    n by p      matrix of probability density values, where p is 
%                         equal to either max(nl)-1 or max(nl), depending on the 
%                         softpdftype. probdens(i,:) are the values of the probability 
%                         density corresponding to the intervals for the i-th soft data 
%                         defined in limi(i,:). (see probasyntax for more explanations)
% nmax        scalar      integer indicating the maximum number of points with coordinates
%                         c that must be kept. 
% dmax        scalar      maximum Euclidian distance between c0 and c. 
%
% OUTPUT :
%
% csub        m by d      matrix which is a subset of lines of the c matrix (m<=n)
% nlsub       m by 1      matrix which is a subset of lines of the nl matrix
% limisub     m by l      matrix which is a subset of the limi matrix
% probdenssub m by d      matrix which is a subset of the probdens matrix
% dsub        m by d      matrix which is a vector of distances between csub and c0 
% nsub        scalar      length of the dsub vector, i.e. nsub=m
% index       m by 1      vector giving the ordering of the lines in csub with
%                         respect to the initial matrix c.
%
% NOTE :
%
% It is possible to specify an additional index for c0
% and c, taking integer values from 1 to nv. Thix index
% specifies which variable is known at each coordinate.
% In that case, c0 and c are cell arrays, where the first
% cell is the matrix of coordinates and the second cell is
% the vector of index values. The nmax variable is then a
% vector of nv elements.

[csub,probasub,dsub,nsub,index]=neighbours(c0,c,[nl limi probdens],nmax,dmax);

if nsub>0
  nlsub=probasub(:,1);
  switch softpdftype
  case {1,2}, l=max(nl); lsub=max(nlsub);    
  case {3,4}, l=3; lsub=3;  
  end;   
  switch softpdftype
  case {1,3}, psub=max(nlsub)-1;    
  case {2,4}, psub=max(nlsub);  
  end;   
  limisub=probasub(:,2:1+lsub);  
  probdenssub=probasub(:,2+l:1+l+psub);
else
  nlsub=[];
  limisub=[];  
  probdenssub=[];
end;

