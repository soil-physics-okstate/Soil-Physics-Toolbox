% BMELIB miscellaneous, general functions (Jan 1, 2001)
%
%   aniso2iso                 - convert anisotropic to isotropic coordinates 
%   col2mat                   - columnwise defined gridded variable toward matrix 
%   coord2dist                - distance matrix from coordinates 
%   coord2K                   - covariance/variogram matrix from coordinates 
%   creategrid                - create a regular grid in 1D, 2D or 3D 
%   findpairs                 - identify identical coordinates 
%   index2noindex             - split grouped coordinates and values for several variables 
%   invdist                   - prediction using an inverse distance weighting 
%   iso2aniso                 - convert isotropic to anisotropic coordinates 
%   isspacetime               - syntaxical analysis of the string 'model' 
%   kernelsmoothing           - prediction using a Gaussian kernel smoothing 
%   neighbours                - radial neighbourhood selection 
%   noindex2index             - separated coordinates and variables toward grouped 
%   split                     - split coordinates in two sets 
%   trapezint                 - integration using trapezoidal formula 
% 
%   Duplicate data
%   isdupli                   - True if there are duplicate coordinates             
%   finddupli                 - find the duplicate coordinates             
%   avedupli                  - averages duplicate values
%   avedupli2D                - averages duplicate values for 2D coordinate x-y
%   combinedupli              - combine duplicate hard and soft data
%   see also  probacombinedupli  in bmeprobalib
