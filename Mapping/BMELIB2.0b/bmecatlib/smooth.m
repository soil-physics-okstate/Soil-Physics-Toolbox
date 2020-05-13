function [pk]=smooth(dk,d,p,o,kstd,order,weighting);

% smooth                 - smoothing using a Gaussian kernel regression method
%                          (December 1,2003)
%
% Implementation of the regression.m function in a moving
% neighbourhood context. Regression is conducted locally
% at a set of coordinates using a least squares estimation
% procedure for a linear regression model, where the
% deterministic part of the linear model is a polynomial of
% a given order. Instead of using ordinary least squares,
% the function uses a diagonal covariance matrix, where the
% variances are inversely proportional to the weights provided
% by a Gaussian kernel, as well as inversely proportional to
% the number of observations. 
%
% SYNTAX :
%
% [pk]=smooth(dk,d,p,o,kstd,order,weighting); 
%
% INPUT :
%
% dk          nk by 1   vector of distances for the smoothed estimates.
% d           n by 1    vector of distances for which the bivariate 
%                       probabilities have been estimated.
% p           n by 1    vector of probabilities at distances specified in d.
% o           n by 1    vector of number of pairs separated by distances
%                       specified in d. 
% kstd        scalar    standard deviation of the Gaussian kernel function
%                       which is used in the kernel regression smoothing
% order       scalar    order of the polynomial used in the kernel regression
%                       smoothing
% weighting   scalar    weigthing is equal to 1 or 0, depending if the user wants
%                       or does not want to use the number of observations
%                       as weights.
%
% OUTPUT :
%
% pk          nk by 1   vector of estimated probabilities at distances dk.

nk=length(dk);
pk=zeros(nk,1)*NaN;
nmax=Inf;

for i=1:nk,
  dmaxi=4*kstd;
  [dsub,posub,trash,nsub]=neighbours(dk(i),d,[p o],nmax,dmaxi);
  if nsub>0,
    dsub=dsub-dk(i);
    w=gausspdf(dsub,[0 kstd^2]);
    w=1./w;
    if weighting==1,
      w=w.*posub(:,2);
    end;
    K=diag(w);
    [best,Vbest]=regression(dsub,posub(:,1),order,K);
    pk(i)=best(1);
  else
    disp(['Warning: No estimation at distance ',num2str(dk(i)),' - increase smoothing parameter']);
  end;
end;

