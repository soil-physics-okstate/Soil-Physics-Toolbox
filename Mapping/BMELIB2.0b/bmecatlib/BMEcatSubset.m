function [pk]=BMEcatSubset(ck,cs,zs,dmodel,Pmodel,nsmax,dmax,options);

% BMEcatSubset           - categorical data prediction using BME with subset data
%                          (December 1, 2003)
%
% Prediction of the conditional categorical probability distribution
% function at a set of locations using the knowledge of possible 
% categories at neighbouring locations. The conditional distribution
% is obtained by conditioning a probability table that has been built
% using a maximum entropy algorithm that incorporates the knowledge of
% the bivariate probabilities between two categories as a function of
% the distance between the corresponding locations.
%
% SYNTAX : 
%
% [pk]=BMEcatSubset(ck,cs,zs,dmodel,Pmodel,nsmax,dmax,options);
%
% INPUT :
%
% ck        nk by d    matrix of coordinates for the locations where the
%                      probabilities for each category have to be estimated.
%                      A line corresponds to the vector of coordinates at a 
%                      location, so the number of columns is equal to the 
%                      dimension of the space. There is no restriction on the
%                      dimension of the space.
% cs        ns by d    matrix of coordinates for the locations where possible
%                      categories are known, with the same convention as for ck.
% zs        ns by nc   matrix of binary values for the subset of possible
%                      categories at the coordinates specified in cs among the
%                      set of nc categories. A value coded as 1 means that the
%                      corresponding category is possible, whereas a value coded
%                      as 0 means that the corresponding category is not possible.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of
%                      bivariate probability values between two categories at
%                      distances specified in dmodel.
% nsmax     scalar     maximum number of locations in the neighbourhood used
%                      for the estimation at the locations specified in ck.
% dmax      scalar     maximum distance between an estimation location and
%                      locations where categories are known. All locations
%                      separated by a distance smaller than dmax from an
%                      estimation location will be included in the estimation
%                      process for that location, whereas other locations are 
%                      neglected.
% options   1 by 2     optional vector of parameters that can be used if the
%                      default values are not satisfactory (otherwise it can simply
%                      be omitted from the input list of variables). 
%                      options(1) is taking the value 1 or 0 if the user wants or
%                      does not want to display the order number of the location
%                      which is currently processed, respectively. Default value
%                      is equal to 0.
%                      options(2) is the stopping criterion when fitting the
%                      maximum entropy table. Default value is equal to 1e-3
%                      (see maxentropytable.m).
%
% OUTPUT :
%
% pk        nk by nc   matrix of conditional probability values, where each
%                      column refers to a category and each line refers to an
%                      estimation location. Each line sum up to one. Values coded
%                      as NaN mean that no estimation has been performed at that 
%                      location due to the lack of available data. 

%%%%%% Remove locations where all categories are possible
%%%%%% (zero specific information content)

isnotall1=find(~all(zs'));
cs=cs(isnotall1,:);
zs=zs(isnotall1,:);

%%%%%% Initialize the parameters

if nargin<8,
  options(1)=0;
  options(2)=1e-3;
end;

nk=size(ck,1);
ncat=size(Pmodel,1);

if options(1)==1,
  num2strnk=num2str(nk);
end;

pk=ones(nk,ncat)*NaN;

%%%%%% Main loop starts here

for i=1:nk,

  ck0=ck(i,:);
  [cslocal,zslocal,ds,sumnslocal]=neighbours(ck0,cs,zs,nsmax,dmax);

  if sumnslocal>0,
      
    isdszero=find(ds==0);
    if isempty(isdszero),
      test=0;
    else
      zs0=zslocal(isdszero(1),:);
      zslocal(isdszero(1),:)=[];
      cslocal(isdszero(1),:)=[];
      sumnslocal=sumnslocal-1;
      test=1;
    end;

    Pkh=maxentropytable([ck0;cslocal],dmodel,Pmodel,options(2));
    Pdenom=sumvaluesatindex(Pkh,[ones(1,ncat);zslocal]);
    for j=1:ncat,
      Pnum=sumvaluesatindex(Pkh,[(1:ncat)==j;zslocal]);
      pk(i,j)=Pnum/Pdenom;
    end;

    if test==1,
      pk(i,:)=zs0.*pk(i,:)/sum(zs0.*pk(i,:));
    end;
  end;

  if options(1)==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;

