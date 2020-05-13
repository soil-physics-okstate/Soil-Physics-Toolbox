function [S]=sumoverone(P,dim);

% sumoverone             - sum along one dimension of a multidimensional table
%                          (December 1, 2003)
%
% Sum the values of a multidimensional table along a specified dimension.
%
% SYNTAX : 
%
% [S]=sumoverone(P,dim);
%
% INPUT :
%
% P       nc by ... by nc   ndim-dimensional table with nc elements along
%                           each one of the ndim dimensions.
% dim     scalar            dimension along which the summation must
%                           operate.
%
% OUTPUT :
%
% S       nc by ... by nc   multidimensional table with nc elements along
%                           each one of the ndim-1 dimensions

S=sum(P,dim);

 
