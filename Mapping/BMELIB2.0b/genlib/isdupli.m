function [c]=isdupli(p)
%  isdupli                - True if there are duplicate coordinates             
%
%  SYNTAX :
% 
%  [c]=isdupli(p);
% 
%  INPUT :
%
%  p          n by d     matrix of the coordinates of n points in a space of dimension d
% 
%  OUTPUT :
%
%  c         logical     True if there are duplicates, False otherwise


c = size(p,1) > size(unique(p,'rows'),1);


