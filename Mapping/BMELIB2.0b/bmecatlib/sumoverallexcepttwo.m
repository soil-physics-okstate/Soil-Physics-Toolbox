function [S]=sumoverallexcepttwo(P,dim1,dim2);

% sumoverallexcepttwo    - sum along all dimensions of a multidimensional table except two
%                          (December 1, 2003)
%
% Sum the values of a multidimensional table along all of its dimensions
% except two of them.
%
% SYNTAX : 
%
% [S]=sumoverallexcepttwo(P,dim1,dim2);
%
% INPUT :
%
% P         nc by ... by nc   ndim-dimensional table with nc elements along
%                             each one of the ndim dimensions.
% dim1      scalar            first dimension along which the summation does
%                             not operate.
% dim2      scalar            second dimension along which the summation does
%                             not operate.
%
% OUTPUT :
%
% S         nc by nc          two-dimensional table obtained by summation over P
%                             except over dimensions dim1 and dim2.

ndim=length(size(P));
S=P;

for i=1:ndim,
  if (dim1~=i)&(dim2~=i),
    S=sum(S,i);
  end;
end;
