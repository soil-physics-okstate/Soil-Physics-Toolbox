function [pu,za,i,dummyout] = avedupli(p,z,dummyin)
% avedupli               - averages duplicate values
%
%  all the duplicate values (with same s/t p coordinates) are averaged.
%  The returned coordinates pu are unique, while za is the corresponding
%  averaged value.  The index i is such that pu=p(i,:)
%
%  SYNTAX :
% 
%  [pu,za,i] = avedupli(p,z);
% 
%  INPUT :
%
%  p          n by d     matrix of the coordinates of n points in a s/t space of dimension d
%  z          n by 1     vector of z values
%
%  OUTPUT :
%
%  pu         nu by d    matrix of the coordinates of nu unique points in a space of dimension d, 
%                           note that nu<=n
%  za         nu by 1    vector of corresponding z-values obtained by averaging the duplicate points
%  i          nu by 1    index such that pu=p(i,:)

%  The following is for backward compatibility, it will be removed in later releases
if nargin >2 | nargout>3
  warning('Mismatch of input/output arguments. Using avedupli2D.m instead of avedupli.m');
  [pu,za,i,dummyout] = avedupli2D(p,z,dummyin);
  return
end

[n,d]=size(p);
if size(z,1)~=n, error('z must have the same number of rows as p'); end
if n==0,  pu=p;  za=z;  i=[]; return; end
if size(z,2)~=1, error('z must have one column'); end

if ~isdupli(p)
  pu=p;
  za=z;
  i=[1:length(z(:))]';
else
  [iu,ir]=finddupli(p);
  pu=p(iu,:);
  za=z(iu);
  i=iu;
  nu=length(iu);
  for k=1:length(ir)
    nu=nu+1;
    i(nu,1)=ir{k}(1);
    za(nu,1)=mean(z(ir{k},1));
  end  
  pu=p(i,:);
end

