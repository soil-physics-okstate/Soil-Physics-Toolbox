% BMELIB analysis of space/time variability, statitics and estimation of variograms (Jan 1, 2001)
%
%
% Analyzing the distribution of a variable: sample histogram, pdf and cdf and other stats :
%
%   decluster                 - spatial declustering weights 
%   histline                  - density-scaled histogram plotted using a line 
%   histscaled                - density-scaled histogram  plotted using bars 
%   kerneldensity             - pdf estimation using a Gaussian kernel 
%   pdf2cdf                   - compute the cdf from the pdf 
%   cdfest                    - experimental cumulative distribution 
%   cdf2pdf                   - compute the pdf from the cdf 
%   quantest                  - quantiles estimated from raw data
%   skewness                  - experimental skewness coefficient 
%   kurtosis                  - experimental kurtosis coefficient 
%   spearman                  - experimental Spearman's correlation matrix 
%   values2rank               - transform a matrix of values into a matrix of rank 
%
% Is your data Gaussian distributed? Transformation to and from the Gaussian distribution :
%
%   gaussplot                 - Gaussian probability plot
%   other2gauss               - transform from an arbitrary pdf to a Gaussian pdf 
%   gauss2other               - transform from a Gaussian pdf to an arbitrary pdf 
%
% Operations on discretized pdf :
%
%   pdf2CI                    - Compute confidence set from a pdf 
%   pdfconvol                 - convolution of a pdf with a Gaussian pdf 
%   pdfprod                   - product of two pdf's. 
%   pdfstat                   - statistics of a pdf  
%   quantile                  - quantiles computed from the pdf  
%   entropy                   - Shannon's entropy of a pdf  
%
% Estimation of sample variograms from experimental data :
%
%   pairsplot                 - display pairs of points separated by a given distance interval
%   pairsindex                - finds pairs of points separated by a given distance interval
%   vario                     - multivariate variogram estimation 
%   covario                   - multivariate covariance function estimation 
%   crosscovario              - single cross covariance function estimation 
%   crossvario                - single cross variogram estimation 
%   crosscovarioST            - space/time cross covariance estimation for s/t vector data
%   crossvarioST              - space/time cross variogram estimation for s/t vector data
%   stcov                     - space/time cross covariance estimation for s/t grid data 
%   stmean                    - space/time mean estimation for s/t grid data
%   stmeaninterp              - interpolation of calculated s/t grid mean estimates
%   stmeaninterpstv           - Interpolate s/t mean onto s/t vector data
%
% Fitting a variogram model to an estimated sample variogram :
%
%   coregfit                  - fitting of the coregionalization model 
%   modelfit                  - single variogram/covariance least squares fitting 
%   fminsmodelfit             - fminsearch subroutine for modelfit.m 
%
% Regression :
%
%   regression                - parameters estimation in a linear regression model 
%   kernelregression          - prediction using a Gaussian kernel regression method 
%   designmatrix              - design matrix in a linear regression model 
%
% Some classical probability distribution functions :
%
%   exponentialpdf            - exponential probability distribution function 
%   exponentialcdf            - exponential cumulative distribution function 
%   exponentialinv            - inverse exponential cumulative distribution function 
%   exponentialstat           - Mean, variance and median of the exponential distribution 
%
%   gausspdf                  - Gaussian probability distribution function 
%   gausscdf                  - Gaussian cumulative distribution function 
%   gaussinv                  - inverse Gaussian cumulative distribution function 
%   gaussstat                 - Mean, variance and median of the Gaussian distribution 
%
%   gaussbipdf                - bivariate Gaussian probability distribution function 
%   gaussbicdf                - bivariate Gaussian cumulative distribution function 
%
%   gaussmultpdf              - multivariate Gaussian probability distribution function
%   gaussmultcdf              - multivariate Gaussian cumulative distribution function 
%
%   loggausspdf               - log-Gaussian probability distribution function 
%   loggausscdf               - log-Gaussian cumulative distribution function 
%   loggaussinv               - inverse log-Gaussian cumulative distribution function 
%   loggaussstat              - Mean, variance and median of the log-Gaussian distribution 
%
%   triangularpdf             - triangular probability distribution function 
%   triangularcdf             - triangular cumulative distribution function 
%   triangularinv             - inverse triangular cumulative distribution function 
%   triangularstat            - Mean, variance and median of the triangular distribution 
%   fmintriangularinv         - fminbnd subroutine for triangularinv.m 
%
%   uniformcdf                - uniform cumulative distribution function 
%   uniforminv                - inverse uniform cumulative distribution function 
%   uniformpdf                - uniform probability distribution function 
%   uniformstat               - Mean, variance and median of the uniform distribution 
