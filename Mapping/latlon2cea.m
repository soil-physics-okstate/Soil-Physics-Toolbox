function [x, y]=latlon2cea(lat, lon, lat_ts, lon_0, x_0, y_0)

  %%% --------------------------------------------------------------------------
  %%% Function to convert from geographic coordinates (i.e. longitude/latitude)
  %%% to Cylindrical Equal Area coordinates on Earth in units of meters.
  %%% 
  %%% Default parameter values are set to apply to the state of Oklahoma,
  %%% specifically the candidate grid to be used in mapping soil moisture.
  %%%
  %%% Parameters:
  %%% lat - latitude [degrees], required
  %%% lon - longitude [degrees], required
  %%% lat_ts - latitude of true scale [degrees], optional
  %%%          default: 35.350
  %%% lon_0 - longitude of origin [degrees], optional
  %%%         default: -97.674
  %%% x_0 - false easting [m], optional
  %%%       default: approx 0
  %%% y_0 - false northing [m], optional
  %%%       default: approx -4,600,000
  %%%
  %%% Output:
  %%% [x, y] - x and y values
  %%%
  %%% Jason Patton
  %%% Last updated 2015-07-07
  %%% Edit: Changed 'endfunction' to 'end' (line 65)
  %%% --------------------------------------------------------------------------

  %%% check variables and set to default values

  % latitude of true scale
  if ~exist('lat_ts', 'var')
     lat_ts = 35.350; % middle latitude of Oklahoma
  end

  % longitude of origin
  if ~exist('lon_0', 'var')
     lon_0 = -97.674; % corner of Kingfisher/Logan/Canadian/Oklahoma counties
  end

  % false easting
  if ~exist('x_0', 'var')
     x_0 = -2.310249342e-9; % essentially 0 since lon_0 already set
  end

  % false northing
  if ~exist('y_0', 'var')
     y_0 = -4566092.606432601; % shifts origin from Equator to corner of counties
  end

  %%% set radius of Earth

  radius = 6378137; % [m]

  %%% compute x and y, see
  % http://en.wikipedia.org/wiki/Cylindrical_equal-area_projection#Formulae

  x = ((lon - lon_0)*pi/180)*cosd(lat_ts)*radius - x_0;
  y = sind(lat)/cosd(lat_ts)*radius + y_0;

  % The results from these equations match up to within a millimeter
  % of the results from Proj.4's projection function.

end
