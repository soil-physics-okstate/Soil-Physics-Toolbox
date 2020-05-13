function [easting,northing]=transverseMercatorProj(la,lo,param);

% transverseMercatorProj    - Transverse mercator projection of latitude-longitude (Jan 1, 2001)
% 
% Perform a Transverse Mercator projection to obtain Easting-Northing (x-y)
% coordinates from longitude-latitude coordinates. This projection requires to
% specify the tangent point between the projection plane and the earth (this
% point should be the longitute-latitude in the middle of the area of interest).
% You must also specify the value of Easting-Northing you want to assign to
% that tangent (origin) point.
% 
% SYNTAX :
%  
% [easting,northing]=transverseMercatorProj(la,lo,param);
%
% INPUT :
%
% la         n by 1       vector of latitude (in degree)
% lo         n by 1       vector of longitude (in degree)
% param      1 by 5 or 4  an optional vector of parameter where
%            param(1): latitude of tangent point for the projection (default value=37 deg)
%            param(2): longitude of tangent point for the projection (default value=-120 deg)
%            param(3): value of easting at tangent point (default value=2500 Km)
%            param(4): value of northing at tangent point (default value=2000 Km)
%            param(5): optional value, radius of the earth (default value=6367 Km)
%
% OUTPUT :
%
% easting    n by 1   vector of eastings, in Km 
% northing   n by 1   vector of northings, in Km 

if nargin<3
  la0 = 37;           % latitude at origin (CS paper of PM10 in NC: la0=37 deg)
  lo0 = -120;         % longitude at origin (CS paper of PM10 in NC: la0=-120 deg)
  disp('Transverse Mercator projection using a tangent point at 37N, 120W');
  x0 = 2500;          % easting at origin (CS paper of PM10 in NC: x0=2500 Km)
  y0 = 2000;          % northing at origin (CS paper of PM10 in NC: y0=2000 Km)
  R = 6367;           % earth's radius (in Km)
else
  la0 = param(1);     % latitude at origin (in degree)
  lo0 = param(2);     % longitude at origin (in degree)
  x0  = param(3);     % easting at origin (in Km)
  y0  = param(4);     % northing at origin (in Km)
  if length(param)>=5
    R = param(5);     % earth's radius (in Km)
  else
    R = 6367;         % earth's radius (in Km)
  end
end;

dtor = pi/180;        % degrees to radians conversion factor
k0 = 1;               % scaling factor
%disp(sprintf('Please note: Using proj. Origin at lat-long=(%.0f,%.0f), east-nor=(%.0f,%.0f)',la0,lo0,x0,y0));
la0 = la0*dtor;       % convert from degrees to radians
lo0 = lo0*dtor;       % convert from degrees to radians

la  = la*dtor;        % convert from degrees to radians
lo  = lo*dtor;        % convert from degrees to radians
B = cos(la).*sin(lo-lo0);
term1 = log((1.0+B)./(1.0-B));
easting   = 0.5 * R * k0 * term1 + x0;
term2 = atan( tan(la) ./ cos(lo-lo0));
northing   = R * k0 * ( term2 - la0 ) + y0;

easting=reshape(easting,size(lo));
northing=reshape(northing,size(la));
