function    [P,Err,Info]=mvnAG1(a,b,C,maxpts,aEps,rEps);

% mvnAG1                    - Multi Variate Normal probability, using the AG1 algorithm (Jan 1, 2001)
%
% Calculates the probability P of a multivariate normal random variable with
% zero mean and covariance matrix C to have a value in a hyper rectangle
% define by the lower bound a and upper bound b.
% The numerical algorithm used is AG1 (Adaptative Genz 1) which uses a
% uses a subregion adaptative method to do the multiple integration of the 
% transformed mvn integral.
% Reference: Alan Genz, "Numerical Computation of Multivariate Normal 
% Probabilities", J. of Computational and Graphical Stat., 1(1992), pp. 141-149
%
% SYNTAX :
%
% [P,Err,Info]=mvnAG1(a,b,C,maxpts,aEps,rEps);
%
% INPUT :
%
% a       N by 1     vector of lower integration limits of length N (N<20)
% b       N by 1     vector of upper integration limits of length N.
% C       N by N     covariance matrix of size N by N
% maxpts  scalar     maximum number of function values allowed. Start with
%                    maxpts = 1000*N, and increase maxpts if Err is too large.
% aEps    scalar     absolute error tolerance.
% rEps    scalar     relative error tolerance.
%
% OUTPUT :
%
% P       scalar     estimated value for the integral
% Err     scalar     estimated absolute error, with 99% confidence level.
% Info    scalar     if INFO = 0, normal completion with ERR < EPS;
%                    if INFO = 1, completion with ERR > EPS and MAXPTS 
%                    function values used; increase MAXPTS to decrease ERROR;

error('It seems that you did not compile the mvnlib directory (see mvnlibcompile)');
