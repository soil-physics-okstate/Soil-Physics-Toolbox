function [softpdftype,nl,limi,probdens]=probaUniform(zlow,zup)
% probaUniform      - creates Uniform soft probabilistic data  at a set of data points
%
% This program generates uniform probabilistic soft pdf at a set of data points.
% At each data point the lower and upper bound of the z values are known to be zlow and zup.
% The soft probrabilistic distribution is given by the uniform pdf as follow
% pdf(z)=1/(zup-low) for zlow<z<zup, 0 otherwise.
%
% SYNTAX :
%
% [softpdftype,nl,limi,probdens]=probaUniform(zlow,zup);
%
% INPUT :
%
% zlow     ns by 1       vector of the lower bound value at ns soft data point
% zup      ns by 1       vector of the upper bound value at ns soft data point
%
% OUTPUT :
%
% softpdftype scalar=2   indicates the type of soft pdf representing the probabilitic soft data.  
% nl          ns by 1    vector of number of limits  nl (see probasyntax)
% limi        ns by 2    matrix representing the limits  values (see probasyntax)
% probdens    ns by 2    matrix representing the probability densities (see probasyntax)
%
% SEE ALSO probaGaussian, probaStundentT and probasyntax.m 

softpdftype=2;
ns=length(zlow);
nl=2*ones(ns,1);
limi=[zlow(:) zup(:)];
probdens=kron(1./diff(limi,1,2),[1 1]);

