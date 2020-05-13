function [P]=createzeros(ndim,nc);

% createzeros            - creates a table filled with values equal to zero
%                          (December 1, 2003)
%
% Create a multidimensional table with all values equal to zero
%
% SYNTAX : 
%
% [P]=createzeros(ndim,nc);
%
% INPUT :
%
% ndim     scalar            dimension of the table to be created
% nc       scalar            number of elements along each dimension
%
% OUTPUT :
%
% P        nc by ... by nc   ndim-dimensional table with nc elements along
%                            each dimension, all values being equal to zero.

tab=ones(1,ndim)*nc;
P=zeros(tab);
