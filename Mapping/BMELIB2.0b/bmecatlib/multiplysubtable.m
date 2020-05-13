function [P]=multiplysubtable(P,dim,index,ptheor,pest);

% multiplysubtable       - multiplication of a subtable by a ratio of scalars
%                          (December 1, 2003)
%
% Multiply a bidimensional marginal subtable by a ratio of
% scalars. This function is used by iterativerescaling.m
%
% SYNTAX :
%
% [P]=multiplysubtable(P,dim,index,ptheor,pest);
%
% INPUT :
%
% P         nc by ... by nc   ndim-dimensional table with nc elements along
%                             each one of the ndim dimensions.
% dim       1 by 2            vector that refers to the dimensions to be
%                             considered in the nc by nc subtable.
% index     1 by 2            vector that refers to the location of the value
%                             to be considered in the nc by nc subtable
% ptheor    scalar            theoretical value for the bivariate probability
%                             in the marginal subtable at the position specified
%                             by index.
% pest      scalar            estimated value for the bivariate probability
%                             in the marginal subtable at the position specified
%                             by index. The marginal subtable is obtained by
%                             summation over all dimensions of P except for 
%                             dimensions dim(1) and dim(2).
%
% OUTPUT :
%
% P         nc by ... by nc   ndim-dimensional table with nc modified elements
%                             along each one of the ndim dimensions.

ndim=ndims(P);

strindex=[];
for i=1:ndim;
  if (dim(1)~=i)&(dim(2)~=i),
    strindex=[strindex,':,'];
  else
    if dim(1)==i,
      strindex=[strindex,num2str(index(1)),','];
    end;
    if dim(2)==i,
      strindex=[strindex,num2str(index(2)),','];
    end;
  end;
end;

strindex=strindex(1:length(strindex)-1);

if ptheor==0,
  eval(['P(',strindex,')=0;']);
else
  eval(['P(',strindex,')=P(',strindex,')*ptheor/pest;']);
end;
