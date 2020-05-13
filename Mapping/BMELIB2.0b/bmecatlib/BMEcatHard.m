function [pk]=BMEcatHard(ck,ch,zh,dmodel,Pmodel,nhmax,dmax,options);

% BMEcatHard             - categorical data prediction using BME with hard data
%                          (December 1, 2003)
%
% Prediction of the conditional categorical probability distribution
% function at a set of locations using the knowledge of observed 
% categories at neighbouring locations. The conditional distribution
% is obtained by conditioning a probability table that has been built
% using a maximum entropy algorithm, that incorporates the knowledge of
% the bivariate probabilities between two categories as a function of
% the distance between the corresponding locations.
%
% SYNTAX : 
%
% [pk]=BMEcatHard(ck,ch,zh,dmodel,Pmodel,nhmax,dmax,options);
%
% INPUT :
%
% ck        nk by d    matrix of coordinates for the locations where the
%                      probabilities for each category have to be estimated.
%                      A line corresponds to the vector of coordinates at a 
%                      location, so the number of columns is equal to the 
%                      dimension of the space. There is no restriction on the
%                      dimension of the space.
% ch        nh by d    matrix of coordinates for the locations where categories
%                      are observed, with the same convention as for ck.
% zh        nh by 1    vector of codes for the categories at the coordinates
%                      specified in ch. Categories are coded as integers ranging
%                      from 1 to nc, where nc is the number of categories.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of
%                      bivariate probability values between two categories at
%                      distances specified in dmodel.
% nhmax     scalar     maximum number of locations in the neighbourhood used
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
  [chlocal,zhlocal,dh,sumnhlocal]=neighbours(ck0,ch,zh,nhmax,dmax);

  if sumnhlocal>0,

    isdhzero=find(dh==0);
    if isempty(isdhzero),
      Pkh=maxentropytable([ck0;chlocal],dmodel,Pmodel,options(2));
      Ph=sumoverone(Pkh,1);
      Pdenom=extractvalueatindex(Ph,[1;zhlocal]);
      for j=1:ncat,
        Pnum=extractvalueatindex(Pkh,[j;zhlocal]);
        pk(i,j)=Pnum/Pdenom;
      end;
    else
      pk(i,:)=zeros(1,ncat);
      pk(i,zhlocal(isdhzero(1)))=1;
    end;
  end;

  if options(1)==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;

