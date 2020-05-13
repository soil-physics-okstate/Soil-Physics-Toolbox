function [cx] = trans(cx,model,im)
%%TRANS returns rotated and reduced coordinates.
%
% TRANS is called from COKRI2. It takes as input original coordinates and
% returns the rotated and reduced cooridnates following specifications
% described in model (im,:)
%
% Coded by Andres Patrignani on 13-Jun-2015 from 
% code published by Denis Marcotte, Cokriging with matlab, Computers & 
% Geosciences, Volume 17, Issue 9, 1991, Pages 1265-1280.

% Some constants are defined

[n,d] = size(cx);
[m,p] = size(model);

% Check for 1-D or isotropic model

if p-1>d
    
    % Perform rotation counterclockwise
    
    if d==2
        ang = model(im,4);
        cang = cos(ang/180*pi);
        sang = sin(ang/180*pi);
        rot = [cang,-sang; sang,cang];
    else
        
        % Rotation matrix in 3-D (the first three among d coordinates) is
        % computed around x, y, and z in that order.
        
        rot = eye(3);
        for i=1:3
            ang = model(im,4+i);
            cang = cos(ang/180*pi);
            sang = sin(ang/180*pi);
            axe = ones(3,1);
            axe(i) = 0;
            rot(axe,axe) = rot(axe,axe)*rot2;
        end
    end

% Rotation is performed aroundf x, y, and zin that order, the other
% coordinates are left unchanged.

dm = min(3,d);
cx(:,1:dm) = cx(:,1:dm)*rot;
t = diag([model(im,2:1+dm),ones(d-dm,1)]);
else
    t = eye(d)*model(im,2);
end

% Perform contractions or dilatations (reduced h)
t = max(t,1e-10);
cx = cx/t;

    

