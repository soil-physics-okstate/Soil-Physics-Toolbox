function [Pfit]=probatablefit(dfit,d,P,o,kstd,options);

% probatablefit          - fitting of the multicategory bivariate probability tables
%                          (December 1,2003)
%
% SYNTAX : [Pfit]=probatablefit(dfit,d,P,o,kstd,options);
%
% INPUT :
%
% dfit      n by 1       vector of values that specify the distances for which fitted
%                        bivariate probability values are sought. Maximum distance in
%                        dfit cannot exceed maximum distance in d.
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
%                        corresponding distance classes. 
% kstd      scalar       standard deviation of the Gaussian kernel function which is used
%                        in the kernel regression smoothing (the higher the value is, the
%                        stronger the smoothing will be). Suggested starting value is
%                        about 0.5*max(d)/length(d).
% options   1 by 3       optional vector of parameters that can be used if default values are
%                        not satisfactory (otherwise this vector can simply be omitted from the
%                        input list of variables), where :
%                        options(1) displays the fitted probability tables if the value is set
%                        to one (default value is 0).
%                        options(2) is the order of the polynomial used in the kernel regression
%                        smoothing (default value is 0).
%                        options(3) uses the number of obervations as weights for the fitting
%                        if the value is set to one (default value is 1).
%
% OUTPUT :
%
% Pfit      nc by nc     cell array that contains the fitted bivariate probability tables
%                        between the nc categories for the distance classes specified in dfit.
%                        Diagonal cells contain the n by 1 vector of probability estimates
%                        for the same category, whereas off-diagonal cells contain the n by
%                        1 vector of cross-category probability estimates.

%%% Initialize the parameters

if max(dfit)>max(d),
  error('maximum distance in dfit cannot exceed maximum distance in d');
end;

if nargin<6,
  options(1)=0;
  options(2)=0;
  options(3)=1;
end;

%%% Compute table of bivariate probabilities at null distance

P0=probamodel2bitable(0,d,P);

%%% Create symmetric vectors for distance and number of 
%%% observations around null distance

nv=size(P,1);
nd=size(d,1);
ndfit=size(dfit,1);
dsym=[flipud(-d);d(2:nd)];
osym=[flipud(o);o(2:nd)];

%%% Interpolate bivariate probabilities at dfit using a kernel 
%%% regression method. Values are mirrored first around null distance ;
%%% this makes possible to avoid the border effect when smoothing for
%%% dfit values close to zero

Pfit=cell(nv,nv);
Ptot=zeros(ndfit,1);
for i=1:nv,
  for j=i:nv, 
    psym=[flipud(-P{i,j})+2*P0(i,j);P{i,j}(2:nd)];
    Pfit{i,j}=smooth(dfit,dsym,psym,osym,kstd,options(2),options(3));
    Pfit{j,i}=Pfit{i,j};
    minPfitij=min(Pfit{i,j});
    if (minPfitij<0)&(minPfitij>-eps),
      index=(Pfit{i,j}<0);
      Pfit{i,j}(index)=0;
    end;
    if max(Pfit{i,j})<0,
      disp(['Warning: Non valid model - probabilities below 0 for categories ',num2str(i),'-',num2str(j)]);
      disp(['Try to change polynomial order or smoothing parameter']);
    end;
    if max(Pfit{i,j})>1,
      disp(['Warning: Non valid model - probabilities above 1 for categories ',num2str(i),'-',num2str(j)]);
      disp(['Try to change polynomial order or smoothing parameter']);
    end;
    if isnan(sum(Pfit{i,j})),
      disp(['Warning : some fitted values are NaN''s for categories ',num2str(i),'-',num2str(j)]);
      disp(['Increase smoothing parameter or use other distance classes in probatablecalc.m']);
    end;
    if i~=j,
      Ptot=Ptot+Pfit{i,j};
    end;
  end;
end;

%%%%%% display the fitted probability tables if options(1)=1

if options(1)==1,
  test=(ishold==1);
  for i=1:nv,
    for j=i:nv,
      minPij=min(P{i,j});
      maxPij=max(P{i,j});
      subplot(nv,nv,(i-1)*nv+j);
      plot(d,P{i,j},'.');hold on;
      plot(dfit,Pfit{i,j});
      plot([0 max(d)],[0 0],':');
      set(gca,'FontSize',8);
      axis([0 max(d) min([0;-1.1*sign(minPij)*minPij]) max([0;1.1*sign(maxPij)*maxPij])]);
      if i==j,
        xlabel('Distance','FontSize',8);
        ylabel('Probability','FontSize',8);
      end;
      title(['Couple ',num2str(i),'-',num2str(j)'],'FontSize',8);
    end;
    if test==0,
      hold off;
    end;
  end;
end;

