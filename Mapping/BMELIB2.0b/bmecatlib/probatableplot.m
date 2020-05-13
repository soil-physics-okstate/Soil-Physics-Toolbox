function probatableplot(d,P,Property,Value);

% probatableplot         - plot of the multicategory bivariate probability tables
%                          (December 1,2003)
% 
% Plot the bivariate probability tables for a set of categories
% as a function of the distance separating locations. 
%
% SYNTAX :
%
% probatableplot(d,P,Property,Value);
%
% INPUT :
%
% d         ncl by 1     vector giving the sorted values of the mean distance separating
%                        the pairs of points that belong to the same distance class. 
% P         nc by nc     cell array that contains the bivariate probability table estimates
%                        between the nc categories for the distance classes specified in d.
%                        Diagonal cells contain the ncl by 1 vector of probability estimates
%                        for the same category, whereas off-diagonal cells contain the ncl by
%                        1 vector of cross-category probability estimates.
% Property  1 by p       cell array where each cell cell is a string that contains
%                        a legal name of a plot object property. This variable is
%                        optional, as default values are used if Property is missing
%                        from the input list of variables. Execute get(H), where H is
%                        a plot handle, to see a list of plot object properties and
%                        their current values. Execute set(H) to see a list of plot
%                        object properties and legal property values. See also the help
%                        for plot.m.
% Value     1 by p       cell array where each cell is a legal value for the corresponding
%                        plot object property as specified in Property.
%
% NOTE :
%
% For example, when Property={'Color','Linewidth'} and Value={'[0.5 0.5 1]',1},
% the values will be displayed as purple broken lines with a width of 1 pixel.
% By default, probatableplot.m will use the default properties for plot.m

if nargin>2,
  if ~iscell(Property),
    Property={Property};
    Value={Value};
    noptions=1;
  else
    noptions=length(Property);
  end;
else
  noptions=0;
end;

nc=size(P,1);
test=(ishold==1);
for i=1:nc,
  for j=i:nc,
    minPij=min(P{i,j});
    maxPij=max(P{i,j});
    subplot(nc,nc,(i-1)*nc+j);
    a=plot(d,P{i,j});
    for k=1:noptions,
      set(a,Property{k},Value{k});
    end
    hold on;
    plot([0 max(d)],[0 0],':');
    set(gca,'FontSize',8);
    axis([0 max(d) min([0;-1.1*sign(minPij)*minPij]) max([0;1.1*sign(maxPij)*maxPij])]);
    if i==j,
      xlabel('Distance','FontSize',8);
      ylabel('Probability','FontSize',8);
    end;
    title(['Categories ',num2str(i),'-',num2str(j)],'FontSize',8);
    title(['Categories ',num2str(i),'-',num2str(j)],'FontSize',8);
    if test==0,
      hold off;
    end;
  end;
end;
