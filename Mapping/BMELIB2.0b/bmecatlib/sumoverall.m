function [S]=sumoverall(P);

% sumoverall             - sum along all dimensions of a multidimensional table
%                          (December 1, 2003)
%
% Sum the values of a multidimensional table along all of its dimensions.
%
% SYNTAX : 
%
% [S]=sumoverall(P);
%
% INPUT :
%
% P       nc by ... by nc   ndim-dimensional table with nc elements along
%                           each one of the ndim dimensions.
%
% OUTPUT :
%
% S       scalar          that corresponds to the sum of all the values
%                         contained in P.

P=P(:),
S=sum(P);
 
