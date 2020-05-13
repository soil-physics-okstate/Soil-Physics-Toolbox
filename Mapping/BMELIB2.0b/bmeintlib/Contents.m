% BMELIB spatiotemporal estimation, using hard data and interval soft data (Jan 1, 2001)
%  
% BME estimation :
%   BMEintervalMode           - BME mode prediction with interval data 
%   BMEintervalPdf            - BME pdf prediction with interval data 
%   BMEintervalTMode          - BME mode prediction with transform and interval data 
%   BMEintervalTPdf           - BME mode prediction with transform and interval data  
%
% Lower level functions :
%   fminBMEintervalMode       - fminbnd subroutine for BMEintervalMode.m 
%   fminBMEintervalTModeBME   - fminbnd subroutine for BMEintervalTMode.m, BME case 
%   fminBMEintervalTModeSK    - fminbnd subroutine for BMEintervalTMode.m, SK case  
%   gausscondpdf              - conditional gaussian pdf for intervals 
%   intBMEinterval            - integration subroutine for BMEinterval*.m 
%   localmeanBME              - local mean estimation for BME 
%   transformderiv            - derivative of a transform function 

