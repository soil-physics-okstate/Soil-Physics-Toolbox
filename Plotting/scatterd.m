function [xb,yb,n1] = scatterd(X,binrange,outliers)


n = hist3(X,binrange); % default is to 10x10 bins
n1 = n';
n1(size(n,1) + 1, size(n,2) + 1) = 0;
xb = linspace(min(X(:,1)),max(X(:,1)),size(n,1)+1);
yb = linspace(min(X(:,2)),max(X(:,2)),size(n,1)+1);
n1(n1==0) = NaN;
if exist('outliers','var')
    n1(n1<outliers) = NaN;
end
% h = pcolor(xb,yb,n1);
% colormap(parula);
% set(h, 'EdgeColor', 'none');
% colorbar
