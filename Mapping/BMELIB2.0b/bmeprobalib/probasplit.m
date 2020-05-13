function [c1,c2,nl1,limi1,probdens1,nl2,limi2,probdens2]=probasplit(c,softpdftype,...
  nl,limi,probdens,index);

% probasplit                - Splits probabilitic data in two sets (Jan 1, 2001)
%
% Split the matrix of coordinates and probabilistic soft data
% into two sets according to the vector index that refers
% to the line numbers.
%
% SYNTAX :
%
% [c1,c2,nl1,limi1,probdens1,nl2,limi2,probdens2]=probasplit(c,softpdftype,nl,...
%   limi,probdens,index);
%
% INPUT :
%
% c           n by d      coordinates of n points with soft probabilistic data.
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
% index       ni by 1     matrix of index for the c1 points, where ni<=n 
%
% OUPUT :
%
% c1          ni by d     matrix of the ni selected coordinates (as specified by index)
% c2      (n-ni) by d     matrix of the reminder of coordinates
% nl1,limi1,probdens1     soft probabilistic data corresponding to the points c1
% nl2,limi2,probdens2     soft probabilistic data corresponding to the points c2
%
% NOTE :
%
% It is possible to specify an additional index for c,
% taking integer values from 1 to nv. Thix index
% specifies which variable is known at each coordinate.
% In that case, c is a cell array, where the first cell
% is the matrix of coordinates and the second cell is
% the vector of index values. The output matrices c1 and
% c2 are then cell arrays too, with the index splitted
% accordingly.

[c1,c2,nl1,nl2]=split(c,nl,index);
[dummy1,dummy2,limi1,limi2]=split(c,limi,index);
if softpdftype==1 | softpdftype==2
  limi1=limi1(:,1:max(nl1));
  limi2=limi2(:,1:max(nl2));
end;                              
[dummy1,dummy2,probdens1,probdens2]=split(c,probdens,index);
switch softpdftype
case {1,3}, 
  probdens1=probdens1(:,1:max(nl1)-1); 
  probdens2=probdens2(:,1:max(nl2)-1); 
case {2,4}, 
  probdens1=probdens1(:,1:max(nl1)); 
  probdens2=probdens2(:,1:max(nl2)); 
end;