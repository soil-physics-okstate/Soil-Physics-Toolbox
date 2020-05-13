function [d]=distSpheroidal(c1,c2,param);

% distSpheroidal            - Spheroidal distance for latitude longitude (Jan 1, 2001)
%
% Calculates the spheroidal distance between two set of points in longitude-latitude.
%
% SYNTAX :
%  
% [d]=distSpheroidal(c1,c2,param);
%
% INPUT :
%
% c1     n1 by 2       matrix of latitude-longitude coordinates, in degrees
% c2     n2 by 2       matrix of latitude-longitude coordinates, in degrees
% param  1 by 5        an optional vector of parameters where
%                      param(3): not used
%                      param(4): not used
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

if nargin>=3 & length(param)>=5
  R = param(5);     % earth's radius (in Km)
else
  R = 6367;         % earth's radius (in Km)
end;

d2r=pi/180;

la1=c1(:,1)*d2r;
lo1=c1(:,2)*d2r;
la2=c2(:,1)*d2r;
lo2=c2(:,2)*d2r;
x1=cos(lo1).*cos(la1)*R;
y1=sin(lo1).*cos(la1)*R;
z1=sin(la1)*R;
x2=cos(lo2).*cos(la2)*R;
y2=sin(lo2).*cos(la2)*R;
z2=sin(la2)*R;
dx=kron(x1,ones(1,n2))-kron(x2',ones(n1,1));
dy=kron(y1,ones(1,n2))-kron(y2',ones(n1,1));
dz=kron(z1,ones(1,n2))-kron(z2',ones(n1,1));
t=sqrt(dx.^2+dy.^2+dz.^2);
d=(2*R)*asin((0.5/R)*t);

