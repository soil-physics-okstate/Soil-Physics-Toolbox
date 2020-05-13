function J = jet2(m)
%JET2    Variant of JET.
%   JET2(M), a variant of JET(M), is the colormap used with the
%   NCSA fluid jet image, except red is replaced by pink
%   JET2, by itself, is the same length as the current colormap.
%   Use COLORMAP(JET2).
%
%   See also JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   C. B. Moler, 5-10-91, 8-19-92.
%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 5.5 $  $Date: 2000/06/01 02:53:46 $
%   modified by M. L. SERRE, 04-20-02.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = max(round(m/4),1);
x = (1:n)'/n;
y = (n/2:n)'/n;
e = ones(length(x),1);
r = [0*y; 0*e; x; e; flipud(y)];
g = [0*y; x; e; flipud(x); 0*y];
%b = [y; e; flipud(x); 0*e; 0*y];        % old line for jet
b = [y; e; flipud(x); x; flipud(y)];     % new line for jet2
J = [r g b];
while size(J,1) > m
   J(1,:) = [];
   if size(J,1) > m, J(size(J,1),:) = []; end
end
