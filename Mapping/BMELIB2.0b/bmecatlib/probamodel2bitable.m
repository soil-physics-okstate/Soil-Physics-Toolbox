function [Pd]=probamodel2bitable(d,dmodel,Pmodel);

% probamodel2bitable     - build a bivariate probability table from a model
%                          (December 1, 2003)
%
% Build the bivariate probability table at a specified distance from
% the set of bivariate probability tables Pmodel available at distances
% dmodel. The probabilities are obtained using a linear interpolation method.
%
% SYNTAX :
%
% [Pd]=probamodel2bitable(d,dmodel,Pmodel);
%
% INPUT :
%
% d         scalar     distance for which the bivariate probability table
%                      has to be computed.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of bivariate
%                      probability values between two categories at distances
%                      specified in dmodel.
%
% OUTPUT :
%
% Pd        nc by nc   table of bivariate probabilities at distance d

nc=size(Pmodel,1);
Pd=zeros(nc,nc);

for k=1:nc,
  for l=1:nc,
    Pd(k,l)=interp1(dmodel,Pmodel{k,l},d,'linear');
  end;
end;

