function [q]=quantest(Z,p)

% quantest        - quantiles estimated from raw data
%
% Compute the quantiles associated with probabilities using
% a set of raw data values representing the distribution
% of interest
%
% SYNTAX :
%
% [q]=quantest(Z,p);
%
% INPUT :
%
% Z     n1 by n2  matrix of values of the raw data
% p     m by k    matrix of probability values for which the quantiles
%                 must be computed.
%                 default is [0.1:0.1:0.9]
%
% OUTPUT :
%
% q      m by k   matrix of quantile values associated with the probabilities
%                 specified in the p matrix.

if nargin<2 p=[0.1:0.1:0.9]; end;

if isempty(Z) | sum(~isnan(Z(:)))==0
  q=p*NaN; 
else
  z=Z(~isnan(Z));
  if length(unique(z(:)))<2
    q=p*NaN;
    q(0<=p & p<=1)=unique(z(:));
  else
    z=sort(z(:));
    n=length(z);
    F=[0:n-1]'/n;
    idx=[logical(1);diff(z)>0];
    z=z(idx);
    F=F(idx);
    z=[z(1);z(2:end)-diff(z)/2;z(end)];
    F=[F;1];
    q=interp1(F,z,p);
  end;
end

