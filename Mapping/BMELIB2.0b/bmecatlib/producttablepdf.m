function [Pprod]=producttablepdf(P,pdf);

% productablepdf         - product between ndim margins and a ndim-dimensional table
%                         (December 1, 2003)
%
% Product between ndim margins and a ndim-dimensional table.
% This function is used by BMEcategPdf.m for taking into account
% a probability distribution function at an estimation location
%
% SYNTAX : 
%
% [Pprod]=producctablepdf(P,pdf);
%
% INPUT :
%
% P           nc by ... by nc   ndim-dimensional table with nc elements along
%                               each one of the ndim dimensions.
% pdf         ndim by nc        matrix of probability values. Each line correspond
%                               to a specific margin and each column corresponds to
%                               a category, so that each line must sum up to 1.
%        
% OUTPUT :
%
% Pprod       nc by ... by nc   ndim-dimensional table obtained by multiplying P by
%                               each line of pdf along its 1-dimensional margins.

ndim=ndims(P);
nc=size(pdf,2);

strindex=[];
for i=1:ndim;
  strindex=[strindex,':,'];
end;
strindex=strindex(1:length(strindex)-1);

Pprod=P;
for i=1:ndim,
  for j=1:nc,
    strindextemp=strindex;
    strindextemp(2*i-1)=num2str(j);
    eval(['Pprod(',strindextemp,')=Pprod(',strindextemp,')*pdf(i,j);']);
  end;
end;


