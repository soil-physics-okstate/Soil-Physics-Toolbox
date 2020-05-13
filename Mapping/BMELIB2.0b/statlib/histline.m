function [hdle,zline,fline,n,x]=histline(z,nbins,bounds,w);

% histline                  - density-scaled histogram plotted using a line (Jan 1, 2001)
%
% Plot with line a density scaled histogram using regularly spaced bins
% between specified minimum and maximum values. The density
% scaled histogram is an histogram that has been rescaled so
% that the histogram has the meaning of a probability distribution
% function, i.e. the sum of the displayed rectangular surfaces is
% equal to 1. 
%
% SYNTAX :
%
% [hdle,zline,fline,n,x]=histline(z,nbins,bounds,w);
%
% INPUT :
%
% z        n by 1      column vector of values.
% nbins    scalar      number of bins to be used for drawing the histogram.
% bounds   1 by 2      optional vector of values for the lower and upper limits of
%                      the first and last bins, respectively.
%                      default is bounds=[min(z) max(z)];
% w        n by 1      optional vector of weights associated with the values
%                      in z (see decluster.m).
%
% OUTPUT :
%
% hdle     scalar      handle for the histogram line, obtained as follow: hdle=plot(zline,fline);
% zline    1 by k      vector of z values for the histogram line
% fline    1 by k      vector of f values for the histogram line
% n        nbins by 1  vector for the values of the density scaled histogram.
% x        nbins by 1  vector of the bins centers for which the values in n
%                      have been computed.
%
% NOTE :
%
% 1- Not-a-Number (NaN) values for z are automatically stripped out when
% drawing the histogram.
%
% 2- Values in z that are outside the interval specified by the bounds
% vector are automatically counted in the first or last bin if they
% are below the lower limit or above the upper limit, respectively.

z=z(:);
if nargin<3, bounds=[min(z) max(z)]; end;
if nargin<4, w=ones(size(z))/length(z); end;
zmin=bounds(1);
zmax=bounds(2);

[frequency,zcenter]=histscaled(z,nbins,bounds,w);
nbins=length(frequency);
diffzcenter=diff(zcenter);
zbegin(1)=zcenter(1)-diffzcenter(1)/2;
zbegin(2:nbins)=zcenter(2:nbins)-diffzcenter/2;
zend=zcenter(1:nbins-1)+diffzcenter/2;
zend(nbins)=zcenter(nbins)+diffzcenter(nbins-1)/2;
indexBegin=(1:2:2*nbins) + 1;
indexEnd=(2:2:2*nbins) + 1;
zline(1)=zbegin(1);Cbin(1)=0;fline(1)=0;
zline(indexBegin)=zbegin;
zline(indexEnd)=zend;
fline(indexBegin)=frequency;
fline(indexEnd)=frequency;
zline(2*nbins+2)=zend(nbins);Cbin(2*nbins+2)=0;fline(2*nbins+2)=0;
hdle=plot(zline,fline);

n=frequency;
x=zcenter;
