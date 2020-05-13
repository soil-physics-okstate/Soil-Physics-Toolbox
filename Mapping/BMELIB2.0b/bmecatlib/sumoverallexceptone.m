function [S]=sumoverallexceptone(P,dim);

% sumoverallexceptone    - sum along all dimensions of a multidimensional table except one
%                          (December 1, 2003)
%
% Sum the values of a multidimensional table along all of its dimensions
% except one of them.
%
% SYNTAX : 
%
% [S]=sumoverallexceptone(P,dim);
%
% INPUT :
%
% P         nc by ... by nc   ndim-dimensional table with nc elements along
%                             each one of the ndim dimensions.
% dim       scalar            dimension along which the summation does
%                             not operate.
%
% OUTPUT :
%
% S         nc by 1           vector obtained after summation over all dimensions
%                             of table P except over dimension dim.

ndim=length(size(P));
S=P;

for i=1:ndim,
  if dim~=i,
    S=sum(S,i);
  end;
end;
S=S(:);
