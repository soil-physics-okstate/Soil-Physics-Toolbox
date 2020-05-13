function [K] = kronfast(A,B)
%KRONFAST is an alternative to Matlab's kron function to avoid running out
%of memory.
%
% Obtained from Peter J. Acklam

[ma,na] = size(A);
[mb,nb] = size(B);

K = reshape( permute( reshape( B(:)*A(:)', ...
         [mb nb ma na] ),...
         [1 3 2 4] ),...
         [ma*mb na*nb] );