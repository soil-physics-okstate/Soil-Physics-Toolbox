function [C]=nuggetCST(Ds,Dt,param);

% nuggetCST                 - Nugget non separable S/T covariance model (Jan 1, 2001)
%
% Nugget space-time covariance model (i.e. model with nugget
% effect in space and time)
%
% SYNTAX :
%
% [C]=nuggetCST(Ds,Dt,param);
%
% INPUT :
%
% Ds     n by m   distance matrix in space, having real values >=0.
% Dt     n by m   distance matrix in time with same dimensions as Ds
% param  1 by 1   vector of parameters such that :
%                 param(1) is the sill of the space-time model,
%
% OUTPUT :
%
% C      n by m   space/time covariance matrix with same size as Ds and Dt.

Cs=nuggetC(Ds,1);
Ct=nuggetC(Dt,1);
C=param(1)*Cs.*Ct;
