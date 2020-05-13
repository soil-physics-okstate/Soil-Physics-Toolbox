function [ypdf]=pdfgauss2other(y,z,zpdf,yfile,Fyfile);

% pdfgauss2other            - Transform pdf for Gaussian RV Z to pdf for arbitrary RV Y (Jan 1, 2001)
%
% Transforms the posterior pdf for a random variable Z to the 
% corresponding posterior pdf for an arbitrary random variable Y,
% where Z and Y are related by the one to one function Y=g(Z). 
% The one-to-one function g(.) is such that it transforms a zero 
% mean unit variance random variable Z to a random variable Y=g(Z) 
% with an arbitrary cumulative distribution function F(Y). Hence the 
% function g(.) is completely defined by the general knowledge (at the 
% prior stage) of the cdf F(Y) for the random variable Y.
% 
% SYNTAX :
%
% [ypdf]=pdfgauss2other(y,z,zpdf,yfile,Fyfile);
%
% INPUT :
%
% y                 ny by 1  vector with the y values for which the 
%                            posterior pdf for Y needs to be calculated
% z and zpdf        nz by 1  vectors with the z and pdf values defining the
%                            posterior pdf for Z
% yfile and Fyfile  k by 1   vectors of with y and F(y) defining the cumulative
%                            distribution F(y) for Y given by the general knowledge
%                            (at the prior stage). yfile values must be sorted by ascending
%                            order
%
% OUTPUT :
%
% ypdf              ny by 1  vectors with the pdf values defining the
%                            posterior pdf for Y corresponding to input y values
%
% NOTE :
%
% 1- If a y value is outside the definition domain
% of the cumulative distribution defined by (yfile,Fyfile),
% the returned ypdf value will be equal to NaN. 

ypdf=NaN*ones(size(y));
[dgfile]=transformderiv(yfile,Fyfile);
[zfile]=other2gauss(yfile,yfile,Fyfile);
nzfile=length(zfile);
isinterp=(min(yfile)<=y & y<=max(yfile));
ytemp=y(isinterp);
ztemp=z(isinterp);
ypdftemp=zpdf(isinterp);
dgtemp=interp1(zfile,dgfile,ztemp,'linear');
if ~isempty(ypdftemp)
  ypdftemp=ypdftemp./dgtemp;
end;
ypdf(isinterp)=ypdftemp;
ypdf(~isinterp)=NaN;
