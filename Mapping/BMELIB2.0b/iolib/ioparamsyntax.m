% ioparamsyntax             - Syntaxical help for BME parameter files (Jan 1, 2001)
%
% Provides a help file explaining the syntax of BME parameter files.
%
% SYNTAX :
%
% help ioparamsyntax;
%
% or, more simply, you may just type:
%
% ioparamsyntax;
%
% CONTENT OF HELP FILE:
%
% A BME parameter file is used to store all the BME parameters necessary
% to perform a BME estimation. This file contains the name of the 
% BME method to use to calculate the BME estimate (either
% the mode or the mean estimate), the name of the file containing 
% the hard data, the name of the file containing the soft data, 
% the name of the file containing the coordinates of the estimation
% points, and all the BME parameters necessary for the BME estimation.
%
% Following is a sample BME parameter file, and below are explanations
% for each of the items in the parameter file:
%
%BeginningOfFile-------------------------------------------------------
%         Parameters for BMEproba 
%            ---------------------
%START OF PARAMETERS:
%BMEprobaMode                     \BME method
%BMEhard1.dat                     \file for hard data
%1                                \nv, number of variables primary + others
%1 2 0 0 3                        \column for x,y,z,t and variable(s)
%BMEproba1.dat                    \file for probabilistic soft data
%0 1.000000e+030                  \data trimming limits
%BMEest1.dat                      \file with the estimation points
%BME1.out                         \output file for BME results
%BME1.dbg                         \output file for BME messages
%2                                \number of covariance models
%nuggetC                          \covariance model 1
%0.05                             \variance for model 1
%                                 \additional parameters for model 1
%exponentialC                     \covariance model 2
%1                                \variance for model 2
%3                                \additional parameters for model 2
%4                                \max number of hard data
%4                                \max number of soft data
%100                              \max search radius
%-1                               \order of mean trend (use -1 for no mean trend)
%1000000 0.02                     \maxpts and rEps for numerical integration
%0.0001 100                       \relative tolerance and maxeval for fminbnd
%50 0.001                         \options(6) and options(7) to plot BME posterior pdf
%3                                \options(8) for BMEprobaMoments, number of BME moments to calculate
%0.25 0.5 0.75 0.9                \Confidence Interval Probabilities, up to 10 values
%0                                \dotransform, 0 for no transformation, 1 otherwise
%EndOfFile--------------------------------------------------------------
%
%
% The first 3 lines of the parameter file are used for comments only.
% Each of the next lines provide input variables, followed by comments (if any). 
% The comments should be separated from the input variable using a few blank spaces, 
% as shown above. The input variables on each of the lines is explained in the
% following
% 
% BME method            string       string indicating the BME method to use. The string must be one 
%                                    of the following: 'BMEprobaMode', 'BMEprobaTMode',
%                                    'BMEprobaMoment', or 'BMEprobaTMoment'. Use the help for
%                                    any of these methods to learn more about it.
% hard data file        string       filename of the hard data file. This is a file in GeoEAS format
%                                    see readGeoEAS for more explanation on that format.
% nv                    scalar       number of variables primary + others
% column number         1 by 4+nv    column number for the x, y, z, t and the variables.
%                                    x, y,and z are the spatial coordinates, and t is the time
%                                    for the hard data given in the hard data file above.
%                                    Use 0 for those coordinates that do not exist.
%                                    If the number of variables nv is greater than 1.
%                                    then BME performs a vector estimation.
% proba data file       string       filename of the probabilistic soft data files. 
%                                    see readProba for more explanation on that file format.
% trimming limits       1 by 2       The lower and upper trimming value for the data. Any values
%                                    below or above these cutoffs are taken as NaN (not a number)
% est coord file        string       filename of the file specifying the coordinates of the estimation 
%                                    points. This is a file in GeoEAS format (see GeoEASsyntax)
% BME output file       string       filename for the output file for the BME estimates.
% BME message file      string       filename for the file containing any messages (if any) generated
%                                    during the BME estimation.
% ncov                  scalar       number of nested covariance models. 
% covname1              string       name of the first covariance model (see modelsyntax)
% variance1(s)          nv by nv     matrix of variance(s) (and crossvariances if nv>1), for the nv
%                                    variables
% param1                1 by k       additional parameters for first covariance model
% covnamei              string           These three input variables are repeated 
% variancei(s)          nv by nv         ncov times (i.e. i=2,...,ncov) in order to 
% parami                1 by k           specify each of the nested covariance models
% nhmax                 scalar       maximum number of hard data values to use in the estimation
% nsmax                 scalar       maximum number of soft data points to use in the estimation
% dmax                  1 by q       maximum distance between hard/soft data and the estimation point
%                                    for purely spatial case (no time), then q=1 and dmax 
%                                    is a spatial distance
%                                    for space/time case (time exitst), then q=3 and 
%                                    dmax(1) is the max spatial distance between data and est pt
%                                    dmax(2) is the max temporal difference between data and est pt
%                                    dmax(3) is the space/time metric (see BMEprobaMode)
% order                scalar        order of the polynomial drift along the spatial axes at the
%                                    estimation locations. For the zero-mean case, NaN (Not-a-Number)
%                                    is used.
%                                    order=-1 for a zero mean (i.e. no mean trend)
%                                         =0 for a constant mean
%                                         =1 for a constant+linear mean along each axis
%                                         =2 for a constant+linear+quadratic mean along each axis, etc.
% maxpts rEps          1 by 2        parameter for numerical integration (see BMEprobaMode)
% options(6:7)         1 by 2        options(6) and option(7) for BME moments (see BMEprobaMoment)
% options(8)           scalar        options(8) for BME moment estimation (see BMEprobaMode)
% PCI                  1 by np       up to 10 values of Confidence Interval Probabilities (not used)
% dotransform,         scalar        0 for no transformation, BME method = BMEprobaMode or BMEprobaMoment
%                                    1 for transformation, BME method = BMEprobaTMode or BMEprobaTMoment
%                                    when dotransform=1, then you need to specify the vectors
%                                    yfile and cdfyfile (see BMEprobaTMode) for each variable. This 
%                                    is done by specifying the number ny of elements in the yfile vector,
%                                    then , by listing on the following ny lines the elements 
%                                    yfile(i) cdfyfile(i). This procedure is repeated for each nv variables. 
%                                    
% NOTE :
%
% Run the program IOLIBtest.m in the testslib directory to experiment with BME parameter files.

help ioparamsyntax
