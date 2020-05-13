function [n,x]=histscaled(z,nbins,bounds,w);

% histscaled                - density-scaled histogram  plotted using bars (Jan 1,2001)
%
% Plot with bars a density scaled histogram using regularly spaced bins
% between specified minimum and maximum values. The density
% scaled histogram is an histogram that has been rescaled so
% that the histogram has the meaning of a probability distribution
% function, i.e. the sum of the displayed rectangular surfaces is
% equal to 1. 
%
% SYNTAX :
%
% [n,x]=histscaled(z,nbins,bounds,w);
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
%
% 3- When output variables are specified, the graphic is not displayed
% and histscaled.m simply returns the values for these variables instead.
% The optional ouput variables have the same meaning than for the hist.m
% function.

if nargin<3, bounds=[min(z(:)) max(z(:))]; end;
if nargin<4,
  w=ones(size(z))/length(z);
end;

w=w(find(~isnan(z)));
z=z(find(~isnan(z)));

class=bounds(1):(bounds(2)-bounds(1))/nbins:bounds(2);
x=(class(2:nbins+1)+class(1:nbins))'/2;
n=zeros(size(x));
for i=1:nbins,
  index=find(z>class(i) & z<=class(i+1));
  n(i)=sum(w(index));
end;

n(1)=n(1)+sum(w(z<=class(1)));
n(nbins)=n(nbins)+sum(w(z>class(nbins+1)));
n=n/(x(2)-x(1));

if nargout==0,
  test=(ishold==1);
  bar(x,n);
  if test==0,
    hold off;
  end;
end;


