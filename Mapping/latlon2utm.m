function [X,Y,zone] = latlon2utm(Lat,Lon,zone,ellipsoidname)
%%LATLON2UTM converts latitude and longitude to planar projected UTM map coorinates.
%
% Syntax
%
%       [X,Y] = latlon2utm(Lat,Lon) converts latitude and longitude to 
%               projected UTM coorinates. UTM_Zone is determined
%               automatically. Lat and Lon must be column vectors. 
%               Since the UTM zone is not specified in this case, Lat and Lon must
%               have their corresponding sign to indicate the correct hemisphere.
%
%       [X,Y] = latlon2utm(Lat,Lon,UTM_Zone) converts latitude and longitude to 
%               projected UTM coorinates. UTM_Zone is a character array
%               such as '14S'. In case you don't know the UTM zone try typing
%               'utmzoneui' in the command window.
%
%       [X,Y,UTM_Zone] = latlon2utm(Lat,Lon,UTM_Zone,ellipsoidname) allows
%                        to specify a custom ellipsoid model. ellipsoidname 
%                        should be character or a EPSG numeric code.
%                        For instance: 
%
%                        LONG NAME                       SHORT NAME   CODE        
%                        World Geodetic System 1984        'wgs84'    7030
%                        Geodetic ReferenceSystem 1980     'grs80'    7019
%
%                        *Note: by default the function selects a geoid
%                        that best fits the current utm zone. You must use
%                        the ellipsoidname argument to set a custom
%                        ellipsoid.
%
%
%       [X,Y,UTM_Zone] = latlon2utm(Lat,Lon) In addition, this syntax retrieves the 
%       predominant UTM zone of the input geo-coordinates.
%
% Note: "The Transverse Mercator and Polar Stereographic projections are 
%       used to organize a worldwide coordinate grid. This system of 
%       projections is generally called Universal Transverse Mercator (UTM)."
%       For more info search for "Working with the UTM System" in the Matlab help.
%
% Adapted from the Matlab help example: "Working with the UTM System" by
% the soil physics research group at Oklahoma State University. Written by AP.
%
% Andres Patrignani - 18 Jul 2015 - Added ellipsoidname and documentation.


%% Check user inputs
if isrow(Lat) 
    Lat = Lat';
end

if isrow(Lon)
    Lon = Lon';
end

if nargin == 2
    zone = utmzone([Lat Lon]); % Alternative to autodetect UTM_Zone.
    [ellipsoid] = utmgeoid(zone); 
    
elseif nargin == 3 
    if ~ischar(zone)
        error('UTM:Format','UTM_zone must be char e.g. ''14S''')
    end
    [ellipsoid] = utmgeoid(zone); 
    
elseif nargin == 4 
    ellipsoid = referenceEllipsoid(ellipsoidname); 
end

%% Build structure array to set the UTM system based on the information obtained
utmstruct = defaultm('utm'); % Initializes a map projection structure
utmstruct.zone = zone;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct); % The empty latitude limits will be set properly by defaultm.

%% Calculate the grid coordinates
[X,Y] = mfwdtran(utmstruct,Lat,Lon); % Project geographic features to map coordinates
format long % Display long -fixed-decimal format.