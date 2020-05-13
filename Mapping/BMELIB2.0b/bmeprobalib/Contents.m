% BMELIB spatiotemporal estimation, using hard data and probabilistic soft data (Jan 1, 2001)
%  
% BME estimation :
%   BMEprobaMode               - BME mode prediction with probabilistic data 
%   BMEprobaMoments            - BME Moments prediction with probabilistic data 
%   BMEprobaMomentsXvalidation - Cross validation for BMEprobaMoments
%   BMEprobaPdf                - BME pdf prediction with proba data 
%   BMEprobaCI                 - BME Confidence Interval and posterior pdf 
%   BMEprobaTMode              - BME mode prediction with transform and probabilistic data 
%   BMEprobaTPdf               - BME pdf prediction with transform and proba data 
%   BMEprobaTCI                - BME Confidence Interval and posterior pdf with transform 
%   BMEoptions                 - Default parameters used in the BME functions 
%
% Manipulating probabilistic soft data :
%   proba2interval            - Transform probabilistic soft data to interval soft data 
%   proba2probdens            - Normalizes the probability density function 
%   proba2stat                - Transform probabilistic soft data to mean/variance data  
%   proba2val                 - Calculates the value of probabilistic soft pdf 
%   probacat                  - Concatenates two sets of probabilistic soft data 
%   probaneighbours           - Radial neighbourhood selection 
%   probaoffset               - Offset the probabilistic soft data 
%   probaother2gauss          - Transform soft data for arbitrary RV Y to soft data for Gaussian RV Z 
%   probasplit                - Splits probabilitic data in two sets 
%   probasyntax               - Syntaxical help for soft probabilistic data 
%   softpdftypeCheckArgs      - Check soft pdf arguments 
%   probaUniform              - creates Uniform soft probabilistic data  at a set of data points
%   probaGaussian             - creates Gaussian soft probabilistic data at a set of data points
%   probaStudentT             - creates student-t soft probabilistic data  at a set of data points
%   probacombinedupli         - combine duplicate soft data of probabilistic type
%   see also combinedupli in genlib
%
% Lower level functions :
%   BMEprobaCheckArgs         - Check arguments for the BMEproba functions 
%   coord2Kinterface          - Interface for coord2K 
%   fixempty                  - Fix empty variables 
%   fminBMEprobaDupliMode     - fminbnd subroutine for BMEprobaMode.m when there is soft data at est pt 
%   fminBMEprobaDupliTMode    - fminbnd subroutine for BMEprobaTMode.m when there is soft data at est pt 
%   fminBMEprobaMode          - fminbnd subroutine for BMEprobaMode.m 
%   fminBMEprobaTMode         - fminbnd subroutine for BMEprobaTMode.m 
%   FyTransformCheckArgs      - Check arguments for the BMEprobaT functions 
%   momentsDupFun             - Moments of the BME posterior pdf with soft data at est pt 
%   momentsFun                - Moments of the BME posterior pdf 
%   pdfgauss2other            - Transform pdf for Gaussian RV Z to pdf for arbitrary RV Y 

