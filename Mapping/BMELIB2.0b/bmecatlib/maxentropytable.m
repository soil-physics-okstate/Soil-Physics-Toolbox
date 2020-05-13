function [Pfit,niter]=maxentropytable(c,dmodel,Pmodel,tol);   

% maxentropytable        - estimate the maximum entropy multivariate probability table
%                          (December 1, 2003)
%
% Estimate the maximum entropy multivariate probability table from
% the knowledge of the bivariate probabilities between categories
% for given distances. The algorithm which is used is the iterative
% rescaling procedure.
%
% SYNTAX :
%
% [Pfit,niter]=maxentropytable(c,dmodel,Pmodel,tol);
%
% INPUT :
%
% c         n  by d    matrix of coordinates for the locations where the
%                      multivariate probabilities for the categories have
%                      to be estimated. A line corresponds to the vector
%                      of coordinates at a location, so the number of columns
%                      is equal to the dimension of the space. There is no
%                      restriction on the dimension of the space.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of
%                      bivariate probability values between two categories at
%                      distances specified in dmodel.
% tol       scalar     value of the stopping criterion for the iterative scaling,
%                      that corresponds to the maximum of the absolute differences
%                      between joint probabilities estimates for two successive
%                      iterations. Default value is equal to 1e-3.
%
% OUTPUT :
%
% Pfit   n by ... by n   n-dimensional table of estimated joint probability values
%                        with nc elements along each of the n dimensions.
% niter  scalar          number of iterations for reaching the stopping criterion.

%%%%%% Initialize the parameters

if nargin<4,
  tol=1e-3;
end;

ncat=size(Pmodel,1);
ndim=size(c,1);

%%%%%% Build the bivariate probability tables

Pbiv=cell(ndim,ndim);
for i=1:ndim,
  for j=i+1:ndim,
    d=coord2dist(c(i,:),c(j,:));
    Pbiv{i,j}=probamodel2bitable(d,dmodel,Pmodel);
  end;
end;

%%%%%% Estimate the maximum entropy table

[Pfit,niter]=iterativerescaling(Pbiv,ndim,ncat,tol);



