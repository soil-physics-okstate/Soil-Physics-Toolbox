function [xu,yu,za,i] = avedupli2D(x,y,z)
% avedupli2D               - averages duplicate values for 2D coordinate x-y
%
%  all the duplicate values (with same x-y coordinates) are averaged.
%  The returned coordinates xu and yu are unique, while za is the corresponding
%  averaged value.  The index i is such that xu=x(i) and yu=y(i)
%
%  SYNTAX :
% 
%  [xu,yu,za,i] = avedupli2D(x,y,z);
% 
%  INPUT :
%
%  x   vector      vector of x coordinates
%  y   vector      vector of y coordinates
%  z   vector      vector of z values
%
%  OUTPUT :
%
%  xu and yu  vectors    vector of x-y coordinates such that (xu,yu) are a set of unique points
%  za         vector     vector of corresponding z-values obtained by averaging the duplicate points
%  i          vector     index such that xu=x(i) and yu=y(i)

[msg,x,y,z] = xyzchk(x,y,z);
if ~isempty(msg), error(msg); end

% Need x,y and z to be column vectors
sz = prod(size(x));
x = reshape(x,sz,1);
y = reshape(y,sz,1);
z = reshape(z,sz,1);
i = (1:sz)';
sxyzi = sortrows([x y z i],[2 1]);
xs = sxyzi(:,1);
ys = sxyzi(:,2);
zs = sxyzi(:,3);
is = sxyzi(:,4);
ind = [0; ys(2:end) == ys(1:end-1) & xs(2:end) == xs(1:end-1); 0];
if sum(ind) > 0
  %disp('Duplicate x-y data points detected: using average of the z values');
  fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
  fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);
  for i = 1 : length(fs)
    % averaging z values
    zs(fe(i)) = mean(zs(fs(i):fe(i)));
  end
  xs = xs(~ind(2:end));
  ys = ys(~ind(2:end));
  zs = zs(~ind(2:end));
  is = is(~ind(2:end));
end
[xyzu]=sortrows([xs ys zs is],4);
xu=xyzu(:,1);
yu=xyzu(:,2);
za=xyzu(:,3);
i=sort(is);

