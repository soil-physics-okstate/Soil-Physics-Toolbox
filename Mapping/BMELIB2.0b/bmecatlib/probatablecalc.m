function [d,P,o]=probatablecalc(c,z,cl,method,options);

% probatablecalc         - multicategory bivariate probability tables estimation
%                          (December 1,2003)
% 
% Estimate the bivariate probability tables for a set of categories which
% are known at a set of coordinates. 
%
% SYNTAX :
%
% [d,P,o]=probatablecalc(c,z,cl,method,options);
%
% INPUT : 
%
% c         n by d       matrix of coordinates for the locations where the
%                        categories are known. A line corresponds to the vector
%                        of coordinates at a location, so the number of columns
%                        is equal to the dimension of the space. There is no
%                        restriction on the dimension of the space.
% z         n by 1       vector of codes for the categories at the coordinates
%                        specified in c. Categories are coded as integers ranging
%                        from 1 to nc, where nc is the number of categories.
% cl        ncl by 1     vector giving the limits of the distance classes that
%                        are used for estimating the probability tables. Distance
%                        classes are open on the left and closed on the right. The
%                        lower limit for the first class is >=0.
% method    string       that contains the name of the method used for computing
%                        the distances between pairs of locations. method='kron'
%                        uses a Kronecker product, whereas method='loop' uses a loop
%                        over the locations. Using the Kronecker product is faster
%                        for a small number of locations but may suffer from memory
%                        size limitations depending on the memory available, as it
%                        requires the storage of a distance matrix. The loop method
%                        may be used whatever the number of data locations is and
%                        must be used if an Out of Memory error message is generated.
%                        Both  methods yield exactly the same estimates.
% options   1 by 1 or 3  vector of optional parameters that can be used if default
%                        values are not satisfactory (otherwise this vector can simply
%                        be omitted from the input list of variables), where :
%                        options(1) displays the estimated probability tables as functions
%                        of distance separating locations if the value is set to one
%                        (default value is 0),
%                        options(2) and options(3) are the minimum and maximum values
%                        for the angles to be considered, using the same conventions as
%                        for the pairsplot.m function. Angles can only be specified for
%                        planar coordinates, i.e. when the number of columns in c is
%                        equal to two.
%
% OUTPUT :
%
% d         ncl by 1     vector giving the sorted values of the mean distance separating
%                        the pairs of points that belong to the same distance class. The
%                        first value corresponds to a null distance. 
% P         nc by nc     symmetric array of cells that contains the bivariate probability
%                        tables estimates between the nc categories for the distance classes
%                        specified in d. Diagonal cells contain the ncl by 1 vector of
%                        probability estimates for the same category, whereas off-diagonal
%                        cells contain the ncl by 1 vector of cross-category probability
%                        estimates.For diagonal cells, the first value of each vector is
%                        the estimated proportion of the corresponding category, whereas
%                        this first value is always equal to 0 for off-diagonal cells.
% o         ncl by 1     vector giving the number of pairs of points that belong to the 
%                        corresponding distance classes. The first value corresponds to the
%                        number of locations specified in c. 
%
% NOTE :
%
% 1- The d, P and o output variables can be used without modification as input
% for the probatableplot.m and probatablefit.m functions.
%
% 2- When a distance class do not contain any pairs of points, the function
% output a warning message. The d and P elements for the corresponding distance
% class are thus coded as NaN's, whereas the corresponding o element is equal to 0.

%%%%%% Initialize the parameters

if ~ischar(method),
  error('method should be a char string');
end;

cl=sort(cl);
if cl(1)<0,
  error('minimum class distance must be >=0');
end;
n=size(c,1);
nc=length(cl)-1;
minim=cl(1);
maxim=cl(nc+1);
nv=max(z);
P=cell(nv,nv);

if nargin==4,
  options(1)=0;
  noptions=1;
else
  noptions=length(options);
end;

if noptions==3,
  a=options(2)*2*pi/360;
  b=options(3)*2*pi/360;
  if size(c,2)~=2,
    error('angle limits are specified only for planar coordinates');
  end;
  if (a==b)|(min([a,b])<-pi/2)|(max([a,b])>pi/2),
    error('angle limits must be different and between or equal to -90 and 90');
  end;
end;

if strcmp(method,'kron')==1,   %%%%% Uses a Kronecker product for computing distances

  %%% Compute the distances

  unit=ones(n,1);
  dc=kron(unit,c)-kron(c,unit);
  if size(dc,2)==1,
    dist=abs(dc);
  else
    dist=sqrt(sum((dc.^2)')');
  end;

  %%% Compute the angles

  if noptions==3,
    finddc1null=find(dc(:,1)==0);
    finddc1notnull=find(dc(:,1)~=0);
    ang=zeros(size(dc,1),1);
    ang(finddc1null)=(pi/2)*sign(dc(finddc1null,2));
    ang(finddc1notnull)=atan(dc(finddc1notnull,2)./dc(finddc1notnull,1));
  end;

  %%% Select couples for appropriate distances and angles

  cond=(dist>max([0,minim]))&(dist<=maxim);
  if noptions==3,
    conda=(ang>a);
    condb=(ang<=b);
    if a<b,
      cond=cond & (conda & condb);
    else
      cond=cond & (conda | condb);
    end;
  end;
  dist=dist(cond);
  m=length(dist);
  if m==0,
    error('no couples of values within the specified classes');
  end;

  %%% Loop over the number of variables and compute probability tables

  isclass=cell(nc);
  d=zeros(nc,1)*NaN;
  o=zeros(nc,1);
  for k=1:nc,
    isclass{k}=find((dist>cl(k))&(dist<=cl(k+1)));
    o(k)=length(isclass{k})/2;
    if o(k)~=0,
      d(k)=sum(dist(isclass{k}))/(2*o(k));
    end;
  end;

  for i=1:nv,
    for j=i:nv,
      zi=kron(unit,z);
      zj=kron(z,unit);
      iscatij=((zi==i)&(zj==j));
      iscatij=iscatij(cond);
      p=zeros(nc,1)*NaN;
      for k=1:nc,
        if o(k)~=0,
           p(k)=sum(iscatij(isclass{k}))/(2*o(k));
        end;
      end;
      if i==j,
        p=[sum(z==i)/length(z);p];
      else
        p=[0;p];
      end;
      P{i,j}=p;
      if i~=j,
        P{j,i}=p;
      end;
    end;
  end;
  d=[0;d];
  o=[length(z);o];
  
else                      %%%%% Uses a loop over the data for computing distances

  d=zeros(nc,1);
  o=zeros(nc,1);
  for i=1:nv,
    for j=1:nv,
      P{i,j}=zeros(nc,1);
    end;
  end;

  for i=1:n,
    for j=i+1:n,
      dist=sqrt(sum((c(i,:)-c(j,:)).^2));
      cond=(dist>max([0 minim]))&(dist<=maxim);
      if noptions==3,
        dc=c(i,1:2)-c(j,1:2);
        if dc(1)==0,
          ang=(pi/2)*sign(dc(2));
        else
          ang=atan(dc(2)/dc(1));
        end;
        conda=(ang>a);
        condb=(ang<=b);
        if a<b,
          cond=cond & (conda & condb);
        else
          cond=cond & (conda | condb);
        end;
      end;
      if cond==1,
        index=sum(dist>cl);
        if (index>=1) & (index<=nc),
          d(index)=d(index)+dist;
          o(index)=o(index)+1;
          for k=1:nv,
            for l=k:nv,
              if k==l,
                P{k,l}(index)=P{k,l}(index)+((z(i)==k)&(z(j)==k));
              else
                P{k,l}(index)=P{k,l}(index)+(((z(i)==k)&(z(j)==l))|((z(i)==l)&(z(j)==k)))/2;
              end;
            end;
          end;
        end;
      end;
    end;
  end;

  for i=1:nc,
    if o(i)==0,
      d(i)=NaN;
      for j=1:nv,
        for k=j:nv,
          P{j,k}(i)=NaN;
          P{k,j}(i)=NaN;
        end;
      end;
    else
      d(i)=d(i)/o(i);
      for j=1:nv,
        for k=j:nv,
          P{j,k}(i)=P{j,k}(i)/o(i);
          P{k,j}(i)=P{j,k}(i);
        end;
      end;
    end;
  end;
  for k=1:nv,
    for l=1:nv,
      if k==l,
        p0=sum(z==l)/length(z);
      else
        p0=0;
      end;
      P{k,l}=[p0;P{k,l}];
    end;
  end;
  d=[0;d];
  o=[length(z);o];
end;

%%%%%% display the computed tables if options(1)=1

if options(1)==1,
  test=(ishold==1);
  for i=1:nv,
    for j=i:nv,
      minPij=min(P{i,j});
      maxPij=max(P{i,j});
      subplot(nv,nv,(i-1)*nv+j);
      plot(d,P{i,j},'.');hold on;
      set(gca,'FontSize',8);
      axis([0 max(d) min([0;-1.1*sign(minPij)*minPij]) max([0;1.1*sign(maxPij)*maxPij])]);
      plot([0 max(d)],[0 0],':');
      if i==j,
        xlabel('Distance','FontSize',8);
        ylabel('Probability','FontSize',8);
      end;
      title(['Categories ',num2str(i),'-',num2str(j)],'FontSize',8);
      if test==0,
        hold off;
      end;
    end;
  end;
end;

%%%%%% Check if there are NaN in output

if length(find(isnan(d)))~=0,
  disp('Warning : some distance classes do not contain pairs of points');
end;

