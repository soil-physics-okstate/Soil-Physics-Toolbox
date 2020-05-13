function [D] = distmat(X,Y)
%DISTMAT computes matrix of Euclidean distances
%
% Inputs
%       X = column vector with X map coordinates
%
%       Y = column vector with Y map coordinates
%
% Outputs
%       D = square matrix of distances
%
% Andres Patrignani 13-Jul-2015

n = length(X);
A = ones(n,1);
B = [X Y];

df = kronfast(A,B) - kronfast(B,A); % Compute differences in each dimension using Kronecker product.
D = sqrt(sum(df.^2,2)); % Compute Euclidean distance.
D = reshape(D,n,n); % Generate n by n square matrix

