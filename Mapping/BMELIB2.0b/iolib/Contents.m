% BMELIB exploratory data analysis, file input/output (Jan 1, 2001)
%
% Reading and writing hard data from files in GEO EAS format 
%   readGeoEAS                - Reads data files in GEO EAS format 
%   writeGeoEAS               - Writes data files in Geo EAS format 
%
% Reading and writing probabilistic soft data from files 
%   writeProba                - Writes probabilistic soft data to files 
%   readProba                 - Reads probabilistic soft data from files 
%
% Reading and writing files for BME estimation with hard and probabilistic data
%   ioparamsyntax             - Syntaxical help for BME parameter files
%   readBMEproba              - Reads parameter and data files for BME estimation
%   writeBMEproba             - Writes parameter and data files for BME estimation
%   BMEprobaFromFile          - BME estimation using input from a file
%
% Reading lines from miscellaneous GIS file format
%   readARCe00                - reads ARCINFO files in .e00 format 
%   readAtlasGIS              - reads GIS files in Atlas GIS format 
%
% Comparison of files
%   fcmp                      - Compare two files 
%
% Conversion between space/time grid  and space/time vector formats
%   stgridsyntax              - Syntaxical help for s/t grid format
%   valstg2stv                - Converts values from s/t grid coord to s/t vector coord 
%   valstv2stg                - Converts values from s/t vector coord to s/t grid coord 
%   probastv2stg              - Converts probabilistic data from s/t vector coord to s/t grid coord 
%   probastg2stv              - Converts probabilistic data from s/t grid coord to s/t vector coord 
