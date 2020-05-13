function [Pfit,niter]=maxentropytablek(c,dmodel,Pmodel,bitable,tol);   

% maxentropytablek       - estimate the maximum entropy multivariate probability table
%                          (December 1, 2003)
%
% Estimate the maximum entropy multivariate probability table from
% the knowledge of the bivariate probabilities between categories
% for given distances. The algorithm that is used is the iterative
% rescaling procedure.
% This function is exclusively used by BMEcatHardk.m
%
% SYNTAX :
%
% [Pfit,niter]=maxentropytablek(c,dmodel,Pmodel,bitable,tol);
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
% bitable   nc by nc   table of bivariate probabilities between categories
%                      of the primary and secondary variables at a same location
% tol       scalar     optional parameter that can be used if the default value
%                      is not satisfactory (otherwise it can simply be omitted
%                      from the input list of variables). tol is the value of the
%                      stopping criterion, that corresponds to the maximum of the
%                      absolute differences between estimated probabilities for
%                      two successive iterations. Default value is equal to 1e-3.
%
% OUTPUT :
%
% Pfit   n by ... by n   n-dimensional table of estimated joint probability values
%                        with nc elements along each of the n dimensions.
% niter  scalar          number of iterations for reaching the stopping criterion.

%%%%%% Initialize the parameters

if nargin<5,
  tol=1e-3;
end;

ncat=size(Pmodel,1);
ndim=size(c,1);

%%%%%% Build the bivariate probability tables

Pbiv=cell(ndim+1,ndim+1);
for i=1:ndim,
  for j=i+1:ndim,
    d=coord2dist(c(i,:),c(j,:));
    Pbiv{i,j}=probamodel2bitable(d,dmodel,Pmodel);
  end;
end;
Pbiv{1,ndim+1}=bitable;
for i=2:ndim,
  Pbiv{i,ndim+1}=zeros(ncat)*NaN;
end;

%%%%%% Estimate the maximum entropy table

[Pfit,niter]=iterativerescaling(Pbiv,ndim+1,ncat,tol);



