function [Pindex]=extractvalueatindex(P,index);

% extractvalueatindex    - extract a single value from a multidimensional table
%                          (December 1, 2003)
%
% Extract a single value from a multidimensional table at
% a specified position.
%
% SYNTAX : 
%
% [Pindex]=extractvalueatindex(P,index);
%
% INPUT :
%
% P          nc by ... by nc   ndim-dimensional table with nc elements along
%                              each of the ndim dimensions.
% index      ndim by 1         vector of value that specify the position of
%                              the value to be extracted.
%        
% OUTPUT :
%
% Pindex     scalar            value of P at the position specified in index

ndim=ndims(P);

strindex=[];
for i=1:ndim;
  strindex=[strindex,num2str(index(i)),','];
end;

strindex=strindex(1:length(strindex)-1);
eval(['Pindex=P(',strindex,');']);

