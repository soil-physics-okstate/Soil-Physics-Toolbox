function [K]=coord2Kinterface(c1,c2,model,param);

% coord2Kinterface          - Interface for coord2K (Jan 1, 2001)
%
% Interface for coord2K.  Return K=zeros(n1,n2) if n1 or n2
% equal 0, else call coord2K which computes the covariance or 
% variogram matrix based on the distances between two sets of 
% coordinates
%
% SYNTAX :
%
% [K]=coord2Kinterface(c1,c2,model,param);
%
% INPUT :
%  same as coord2K
%
% OUTPUT :
% K       n1 by n2 covariance or variogram matrix
%

%%%%%% find n1 and n2

noindex=(iscell(c1)==0);

if noindex==1,          % case for 1 variable
  n1=size(c1,1);
  n2=size(c2,1);
else                    % case for nv variables
  n1=size(c1{1},1);
  n2=size(c2{1},1);
end;

%%%%%% Calculate K

if n1==0 | n2==0,               %%%Return K=zeros(n1,n2) if n1==0 or n2==0
   K=zeros(n1,n2);
else                            %%%Else call coor2K
  K=coord2K(c1,c2,model,param);
end;



