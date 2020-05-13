function [index]=pairsplot(c1,c2,cl,options);

% pairsplot                 - display pairs of points separated by a given distance interval (Jan 1,2001)
%
% Plot line segments for pairs of locations in a two dimensional
% space that are separated by a distance belonging to an interval
% of values. The pairs of locations that meet the distance
% requirements are displayed as connected by a line segment. This
% function is especially useful as an exploratory tool when used
% with the vario.m and covario.m functions.
%
% SYNTAX :
%
% [index]=pairsplot(c1,c2,cl,options);
%
% INPUT : 
%
% c1        n1 by 2   matrix of coordinates for the locations in the first
%                     set. Each line in c1 corresponds to the vector of
%                     coordinates at a location, so the number of columns is
%                     equal to two. The circles symbols plotted at the extremity
%                     of the line segments refer to coordinates in the c1 matrix.
% c2        n2 by 2   matrix of coordinates for the locations in the second set,
%                     using the same conventions as for c1. The star symbols
%                     plotted at the extremity of the line segments refer to
%                     coordinates in the c2 matrix.
% cl        1 by 2    vector that contains the minimum and maximum values for the
%                     distances to be considered. The distance interval is open on
%                     the left and closed on the right. 
% options   1 by 2    optional vector that contains the minimum and maximum values
%                     in degrees for the angles to be considered, with
%                     -90<=options(1) and (2)<=90. The angular class is open on the
%                     left and closed on the right. Angles are counted counterclockwise
%                     from the horizontal axis. Two cases must be distinguished :
%                     i) if options(1)<options(2), then the angles a that are considered
%                     are those for which options(1)<a<=options(2) ;
%                     ii) if options(1)>=options(2), then the angles a that are considered
%                     are those for which a>options(1) or a<=options(2).
%
% OUTPUT :
%
% index     k by 2    matrix of index values for the coordinates in c1 and c2 that satisfy
%                     the distances and angles requirements. Each line of index refers to
%                     a couple of indexes for the coordinates in c1 and c2. The first and
%                     second column of index refer to coordinates belonging to c1 and c2,
%                     respectively, so the number of columns is equal to two.

%%%%%% Initialize the parameters

if size(c1,2)~=2,
  error('Coordinates must be given in a 2D space');
end;

cl=sort(cl);
if cl(1)<0,
  error('Minimum class distance must be >=0');
end;
n1=size(c1,1);
n2=size(c2,1);

%%%%%% Select pairs according to the distances

unit1=ones(n1,1);
unit2=ones(n2,1);
index=[kron((1:n1)',unit2) kron(unit1,(1:n2)')];
dc=kron(unit1,c2)-kron(c1,unit2);
dist=sqrt(sum((dc.^2)')');

isOK=find(dist>max([0 cl(1)]) & dist<=cl(2));
index=index(isOK,:);
dc=dc(isOK,:);

%%%%%% Select pairs according to the angles

if nargin==4,
  a=options(1)*2*pi/360;
  b=options(2)*2*pi/360;
  if (a==b)|(min([a,b])<-pi/2)|(max([a,b])>pi/2),
    error('Angle limits must be different and between or equal to -90 and 90');
  end;
  finddc1null=find(dc(:,1)==0);
  finddc1notnull=find(dc(:,1)~=0);
  ang=zeros(size(dc,1),1);
  ang(finddc1null)=(pi/2)*sign(dc(finddc1null,2));
  ang(finddc1notnull)=atan(dc(finddc1notnull,2)./dc(finddc1notnull,1));
  conda=(ang>a);
  condb=(ang<=b);
  if a<b,
    cond=(conda & condb);
  else
    cond=(conda | condb);
  end;
  index=index(cond,:);
end;

%%%%%% Display the couples

test=(ishold==1);
plot(c1(:,1),c1(:,2),'o');hold on;
plot(c2(:,1),c2(:,2),'*');
nindex=size(index,1);
for i=1:nindex,
  indexi1=index(i,1);
  indexi2=index(i,2);
  plot([c1(indexi1,1) c2(indexi2,1)],[c1(indexi1,2) c2(indexi2,2)]);
end;
if test==0,
  hold off;
end;

