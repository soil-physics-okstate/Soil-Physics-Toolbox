function [P]=createones(ndim,nc);

% createones             - creates a table filled with values equal to one
%                          (December 1, 2003)
%
% Create a multidimensional table with all values equal to one
%
% SYNTAX : 
%
% [P]=createones(ndim,nc);
%
% INPUT :

% ndim     scalar            dimension of the table to be created
% nc       scalar            number of elements along each dimension
%
% OUTPUT :
%
% P        nc by ... by nc   ndim-dimensional table with nc elements along
%                            each dimension, all values being equal to one.

tab=ones(1,ndim)*nc;
P=ones(tab);
