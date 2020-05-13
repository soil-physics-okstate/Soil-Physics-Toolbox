function [lagbins,mergedN,edges,N] = semivarlags(X,Y)
%%SEMIVARLAGS estimates the number of bins and bin boundaries to create
%%empirical semivariograms.
%
% Inputs
%       X = Column vector with spatial coordinates in the X direction. 
%
%       Y = Column vector with spatial coordinates in the Y direction.
%
% Outputs
%       lagbins = edges of bins with more than 50 elements.
%
%       mergedN = number of elements in each merged bin
%
%       edges = edges of bins, including those with < 50 elements
%
%       N = number of elements in each bin, including those with < 50 elements
%
%
% Author: Andres patrignani 21-Jun-2015
% Modified the outputs by JD 3-July-2015
% Modified the outputs by JD 14-July-2015

% Eliminate NaN
nanidx = isnan(X) | isnan(Y); % identify nans in any of the three input vectors
X = X(~nanidx); % X vector without nans
Y = Y(~nanidx); % Y vector without nans
n = length(X);

%% Compute Euclidean distance from each point to all points in the set. 
df = kron(ones(n,1),[X,Y]) - kron([X,Y],ones(n,1)); % Compute differences in each dimension using Kronecker product.
d = sqrt(sum(df.^2,2)); % Calculate Euclidean distance using Pythagoras' Theorem.
d = reshape(d,n,n); % Create n by n square matrix
loweridx = logical(tril(ones(size(d)),-1));
hcloud = d(loweridx);
[~,N,edges] = nbins(hcloud); % Find optimal number of bins. N is the number of elements in each bin.
[sortedh] = sort(hcloud); 
idx = N>50;
csN = cumsum(N);
csNidx = csN(idx);
lagbins = sortedh([1 csNidx]); % Find edges of bins with >50 elements
mergedN = (csNidx - [0 csNidx(1:end-1)])'; % Calculate the number of elements in each merged bin.

