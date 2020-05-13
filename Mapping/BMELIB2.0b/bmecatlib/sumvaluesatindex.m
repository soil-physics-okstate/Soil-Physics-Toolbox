function [S]=sumvaluesatindex(P,index);

% sumvaluesatindex       - sum over several values in a multidimensional table
%                          (December 1, 2003)
%
% Sum the values contained at specified positions in a multidimensional table.
%
% SYNTAX : 
%
% [S]=sumvaluesatindex(P,index);
%
% INPUT :
%
% P           nc by ... by nc   ndim-dimensional table with nc elements along
%                               each one of the ndim dimensions.
% index       ndim by n         matrix that specify the position of the values
%                               over which summation should occur. Each column
%                               corresponds to the position of an element in P.
%
% OUTPUT :
%
% S           scalar            that corresponds to the sum of the values of P
%                               located at the n positions specified by index.

ndim=ndims(P);

strindex=[];
for i=1:ndim;
  strindex=[strindex,'[',num2str(find(index(i,:)==1)),'],'];
end;

strindex=strindex(1:length(strindex)-1);
eval(['Pindex=P(',strindex,');']);
S=sum(Pindex(:));

