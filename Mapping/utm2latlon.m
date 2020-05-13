function [Lat,Lon] = utm2latlon(X,Y,utm_zone)
%
% AP - 22 Apr 2015
% z1=utmzone(p1); % Alternative to autodetect UTM_Zone.
if iscell(utm_zone) && numel(utm_zone)==1
    utm_zone = char(utm_zone);
end
[ellipsoid] = utmgeoid(utm_zone); %Obtain a suggested ellipsoid

%% Set the UTM system based on the information obtained
utmstruct = defaultm('utm'); % Initializes a map projection structure
utmstruct.zone = utm_zone;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct); % The empty latitude limits will be set properly by defaultm.

% Calculate the grid coordinates
[Lat,Lon] = minvtran(utmstruct,X,Y); % Project geographic features to map coordinates
% Use mfwdtran function to convert point locations, lines, and 
% polygon vertices given in latitudes and longitudes to a planar, projected 
% map coordinate system.
format long