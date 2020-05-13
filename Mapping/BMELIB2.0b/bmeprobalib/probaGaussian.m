function [softpdftype,nl,limi,probdens]=probaGaussian(zm,zv)
% probaGaussian     - creates Gaussian soft probabilistic data at a set of data points
%
% This program generates Gaussian probabilistic soft pdf at a set of data points.
% At each data point the mean and variance are known to be zm and the variance is zv.
% The soft probrabilistic distribution is given by the Gaussian pdf with mean zm and 
% variance zv
%
% SYNTAX :
%
% [softpdftype,nl,limi,probdens]=probaGaussian(zm,zv);
%
% INPUT :
%
% zm     ns by 1       vector of the mean value at each of ns soft data point
% zv     ns by 1       vector of the variance at each of ns soft data point
%
% OUTPUT :
%
% softpdftype scalar=2   indicates the type of soft pdf representing the probabilitic soft data.  
% nl          ns by 1    vector of number of limits  nl (see probasyntax)
% limi        ns by l    matrix representing the limits  values (see probasyntax)
% probdens    ns by p    matrix representing the probability densities (see probasyntax)
%
% SEE ALSO probaStundentT, probaUniform and probasyntax.m 

zm=zm(:);
zv=zv(:);
sig=sqrt(zv);
ns=length(zm);

softpdftype=2;
nl=13*ones(ns,1);

xnorm=gaussinv([0.0001 0.001 0.01 0.05 0.3 .4 .5 .6 .7 0.95 0.99 0.999 0.9999],[0 1]);
limi=kron(sig,xnorm)+kron(zm,1+0*xnorm);
fnorm=gausspdf(xnorm,[0 1]);
probdens=kron(1./sig,fnorm);

norm=sum(0.5*(probdens(:,1:end-1)+probdens(:,2:end)).*diff(limi,1,2),2);
probdens=probdens./kron(norm,ones(size(xnorm)));
