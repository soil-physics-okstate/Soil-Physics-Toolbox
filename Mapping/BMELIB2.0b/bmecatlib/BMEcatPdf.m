function [pk]=BMEcatPdf(ck,cs,ps,dmodel,Pmodel,nsmax,dmax,options);

% BMEcatPdf              - categorical data prediction using BME with pdf data
%                          (December 1, 2003)
%
% Prediction of the conditional categorical probability distribution
% function at a set of locations using the knowledge of the probability
% distributions of categories at neighbouring locations. The conditional
% distribution is obtained by conditioning a probability table that has
% been built using a maximum entropy algorithm that incorporates the
% knowledge of the bivariate probabilities between two categories as a
% function of the distance between the corresponding locations.
%
% For prediction using both hard and soft information, see the note at
% the end of this help.
%
% SYNTAX : 
%
% [pk]=BMEcatPdf(ck,cs,ps,dmodel,Pmodel,nsmax,dmax,options);
%
% INPUT :
%
% ck        nk by d    matrix of coordinates for the locations where the
%                      probabilities for each category have to be estimated.
%                      A line corresponds to the vector of coordinates at a 
%                      location, so the number of columns is equal to the 
%                      dimension of the space. There is no restriction on the
%                      dimension of the space.
% cs        ns by d    matrix of coordinates for locations where distributions
%                      of categories are known, with same convention as for ck.
% ps        ns by nc   matrix of [0,1] real values for the probabilities of the
%                      categories at the coordinates specified in cs. Each line
%                      correspond to a location and each column to a category,
%                      so that each line must sum up to 1.
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
%
% NOTE : This program can process any kind of mixture of hard and soft information
%        using the following conventions :
%
%        (1) for a hard data (the jth category is observed at the ith location),
%            the corresponding ith line ps(i,:) has its jth probability value
%            equal to one and other probability values equal to 0.
%        (2) for a subset of n possible categories at the ith location, the
%            corresponding ith line ps(i,:) will have probability values
%            equal to 1/n for the elements corresponding to these categories
%            whereas other elements have a zero probability value.
%       
%        Programs BMEcatHard.m and BMEcatSubset.m which are both specific cases
%        of BMEcatPdf.m have been provided for the user's facility. They will
%        provide the same results.

%%%%%% Remove locations where all non null probabilities are equal
%%%%%% (zero specific information content)

isnotOK=sum(abs(diff(ps')))==0;
cs(isnotOK,:)=[];
ps(isnotOK,:)=[];

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
  [cslocal,pslocal,ds,sumnslocal]=neighbours(ck0,cs,ps,nsmax,dmax);

  if sumnslocal>0,
      
    isdszero=find(ds==0);
    if isempty(isdszero),
      test=0;
    else
      ps0=pslocal(isdszero(1),:);
      pslocal(isdszero(1),:)=[];
      cslocal(isdszero(1),:)=[];
      sumnslocal=sumnslocal-1;
      test=1;
    end;

    Pkh=maxentropytable([ck0;cslocal],dmodel,Pmodel,options(2));
    Pkh=producttablepdf(Pkh,[ones(1,ncat);pslocal]);
    Pdenom=sumvaluesatindex(Pkh,ones(sumnslocal+1,ncat));
    for j=1:ncat,
      Pnum=sumvaluesatindex(Pkh,[(1:ncat)==j;ones(sumnslocal,ncat)]);
      pk(i,j)=Pnum/Pdenom;
    end;

    if test==1,
      pk(i,:)=ps0.*pk(i,:)/sum(ps0.*pk(i,:));
    end;
  end;

  if options(1)==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;

