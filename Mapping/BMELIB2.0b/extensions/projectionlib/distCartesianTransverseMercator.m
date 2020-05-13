function [d]=distCartesianTransverseMercator(c1,c2,param);

% distCartesianTransverseMercator - Cartesian distance for transverse mercator projection of latitude-longitude (Jan 1, 2001)
%
% Calculates the cartesian distance between two set of points in longitude-latitude
% by using a Transverse Mercator projection to obtain Easting-Northing (x-y)
% coordinates from longitude-latitude coordinates. This projection allows to
% specify the tangent point between the projection plane and the earth (this
% point should be the longitute-latitude in the middle of the area of interest).
% If you do not specify the tangent point for the projection, this program 
% uses the average of the longitute-latitude of the two sets of points, which 
% minimizes the error in the calculated distance for any point on the earth, 
% provided that the set of points are not far from one another.
%
% SYNTAX :
%  
% [d]=distCartesianTransverseMercator(c1,c2,param);
%  
% INPUT :
%
% c1     n1 by 2       matrix of latitude-longitude coordinates, in degrees
% c2     n2 by 2       matrix of latitude-longitude coordinates, in degrees
% param  1 by 2 or 5   an optional vector of parameter where
%                      param(1): latitude of tangent point for the projection
%                      default value = mean([c1(:,1);c2(:,1)]) 
%                      param(2): longitude of tangent point for the projection 
%                      default value = mean([c1(:,2);c2(:,2)]) 
%                      param(3): not used
%                      param(4): not used
%                      param(5): optional value, radius of the earth (default value=6367 Km)
%
% OUTPUT
%
% d    n1 by n2        matrix of distances (in Km) between the points c1 and c2

%
% Verify the input
%
n1=size(c1,1);
if size(c1,2)~=2, error('c1 must have 2 columns'); end;
n2=size(c2,1);
if size(c2,2)~=2, error('c2 must have 2 columns'); end;

%
% set the parameters
%
if nargin<3
  param(1)=mean([c1(:,1);c2(:,1)]);  
  param(2)=mean([c1(:,2);c2(:,2)]);  
end;
param(3)=0;
param(4)=0;

%
% Project the latitude-longitute coordinates to obtain easting and northing
%
[x1,y1]=transverseMercatorProj(c1(:,1),c1(:,2),param);
[x2,y2]=transverseMercatorProj(c2(:,1),c2(:,2),param);

%
% Calculate the cartesian distance between the easting and northing projections
%
dx=kron(x1,ones(1,n2))-kron(x2',ones(n1,1));
dy=kron(y1,ones(1,n2))-kron(y2',ones(n1,1));
d=sqrt(dx.^2+dy.^2);
