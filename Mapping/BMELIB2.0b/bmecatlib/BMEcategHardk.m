function [pk]=BMEcatHardk(ck,ch,zh,zk,bitable,dmodel,Pmodel,nhmax,dmax,options);

% BMEcatHardk            - categorical data prediction using BME with hard data and table
%                          (December 1, 2003)
%
% Prediction of the discrete categorical probability distribution
% function for a primary variable at a set of locations using the
% knowledge of categories for that variable at neighbouring location
% as well as a the knowledge of a secondary categorical variables at
% each estimation location. Both variables have the same number of
% categories and are linked by a bivariate probability table
% (see notes at the end of this help).
% The conditional probability distribution is obtained by conditioning
% a probability table that has been obtained using a maximum entropy
% algorithm that incorporates the knowledge of the bivariate probabilities
% between two categories as a function of the distance between the
% corresponding locations.
%
% SYNTAX : 
%
% [pk]=BMEcatHardk(ck,ch,zh,zk,bitable,dmodel,Pmodel,nhmax,dmax,options);
%
% INPUT :
%
% ck        nk by d    matrix of coordinates for the locations where the
%                      probabilities for each category of the primary
%                      variable have to be estimated. A line corresponds to
%                      the vector of coordinates at a location, so the number
%                      of columns is equal to the dimension of the space.
%                      There is no restriction on the dimension of the space.
% ch        nh by d    matrix of coordinates for the locations where categories
%                      of primary variables are observed, with the same
%                      convention as for ck.
% zh        nh by 1    vector of codes for the categories of the primary variable
%                      at the coordinates specified in ch. Categories are coded
%                      as integers ranging from 1 to nc, where nc is the number 
%                      of categories.
% zk        nk by 1    vector of codes for the categories of the secondary variable
%                      at the coordinates specified in ck. Categories are coded
%                      as integers ranging from 1 to nc, where nc is the number 
%                      of categories.
% bitable   nc by nc   table of bivariate probabilities between categories
%                      of the primary and secondary variables at a same location
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of
%                      bivariate probability values between two categories of the
%                      primary variable at distances specified in dmodel.
% nhmax     scalar     maximum number of locations in the neighbourhood used
%                      for the estimation at the locations specified in ck.
% dmax      scalar     maximum distance between an estimation location and
%                      locations where categories for primary variable are known.
%                      All locations separated by a distance smaller than dmax
%                      from an estimation location will be included in the estimation
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
%                      estimation location
%
% NOTES :
%
% (1) This program allows the user to incoporate into the maximum entropy procedure
%     the extra information which is brought by the joint probabilities of observing
%     category i for the primary variable and category j for the secondary variable
%     (wihth i,j=1,...nc) at a same location, the category of the secondary variable
%     being known at all the estimation locations.
% (2) If no primary information is available in the selected neighbourhood, the
%     program will update the a priory probability distribution function of the
%     primary variable with the information brought by the knowledge of the
%     observed category for the secondary variable at the estimation location.

%%%%%% Initialize the parameters

if nargin<10,
  options(1)=0;
  options(2)=1e-3;
end;

nk=size(ck,1);
ncat=size(Pmodel,1);

if options(1)==1,
  num2strnk=num2str(nk);
end;

pk=kron(ones(1,ncat),probamodel2bitable(0,dmodel,Pmodel));

%%%%%% Main loop starts here

for i=1:nk,
    
  ck0=ck(i,:);
  [chlocal,zhlocal,dh,sumnhlocal]=neighbours(ck0,ch,zh,nhmax,dmax);
        
  isdhzero=find(dh==0);
  if isempty(isdhzero),
    Pkh=maxentropytablek([ck0;chlocal],dmodel,Pmodel,bitable,options(2));
    if ~isnan(zk(i)),
      Ph=sumoverone(Pkh,1);
      Pdenom=extractvalueatindex(Ph,[1;zhlocal;zk(i)]);
      for j=1:ncat,
        Pnum=extractvalueatindex(Pkh,[j;zhlocal;zk(i)]);
        pk(i,j)=Pnum/Pdenom;
      end;
    else
      ndim=ndims(Pkh);
      Pkh=sumoverone(Pkh,ndim);
      Ph=sumoverone(Pkh,1);
      Pdenom=extractvalueatindex(Ph,[1;zhlocal]);
      for j=1:ncat,
        Pnum=extractvalueatindex(Pkh,[j;zhlocal]);
        pk(i,j)=Pnum/Pdenom;
      end;
    end;
  else
    pk(i,:)=zeros(1,ncat);
    pk(i,zhlocal(isdhzero(1)))=1;
  end;

  if options(1)==1,
    disp([num2str(i),'/',num2strnk]);
  end;
end;

