function [iu,ir]=finddupli(p)
%  finddupli             - find the duplicate coordinates             
%
%  Find in a set of points those that are in duplicates (i.e. having the same coordinates)
%
%  SYNTAX :
% 
%  [iu,ir]=finddupli(p);
% 
%  INPUT :
%
%  p          n by d     matrix of the coordinates of n points in a space of dimension d
% 
%  OUTPUT :
%
%  iu        nu by 1    vector with the indices of the points that are not in duplicates
%                       When there are no duplicates, then p(iu,:) is equal to p
%                       otherwise nu<n and p(iu,:) is the subset of p that do not have duplicate
%                       coordinates
%  ir        1 by nr    cell array of the duplicate points.  
%                       When there are no duplicates, then ir is empty
%                       Otherwise ir{k} is a k=th cluster of duplicate points, so that
%                       p(ir{k},:) is the k-th subset of p having all the same coordinates

iu=[1:size(p,1)]';
if ~isdupli(p)
  ir=[];
else
  [n d]=size(p);
  i = (1:n)';
  spi = sortrows([p i],1:d);
  ps = spi(:,1:d);
  is = spi(:,d+1);
  ind = ps(2:end,1) == ps(1:end-1,1);
  for k=2:d
    ind = ind & ps(2:end,k) == ps(1:end-1,k);
  end
  ind = [0; ind; 0];
  fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
  fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);
  irv=[];
  for k = 1 : length(fs)
    ir{k}=is(fs(k):fe(k));
    irv=[irv;ir{k}];
  end
  iu=setdiff(iu,irv);
  if isempty(iu), iu=[]; end
end;
