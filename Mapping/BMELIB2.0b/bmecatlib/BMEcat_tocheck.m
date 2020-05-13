function [pk]=BMEcatHardk(ck,ch,zh,dmodel,Pmodel,nhmax,dmax,options);

% BMEcatHardk            - categorical data prediction using BME with hard data (May 1, 2001)
%
% Prediction of the discrete categorical probability distribution
% function at a set of locations using the knowledge of categories
% at neighbouring locations (see BMEcategHard.m).
% This function allows the user to specify an additional index
% for specifying which probability table should be used at the
% prediction and at the hard data points.
%
% SYNTAX : 
%
% [pk]=BMEcatHardk(ck,ch,zh,dmodel,Pmodel,nhmax,dmax,options);
%
% INPUT :
%
% ck        1 by 2     cell array. The first cell is the nk by d matrix
%                      of coordinates for the locations where the probabilities
%                      for each category have to be estimated. A line corresponds
%                      to the vector of coordinates at a location, so the number
%                      of columns is equal to the dimension of the space. There is
%                      no restriction on the dimension of the space. The second
%                      cell is a nk by 1 vector of integer values, ranging from 1
%                      to nt, specifying the probability tables that must be used.
% ch        nh by d    matrix of coordinates for the locations where the
%                      categories are known, with the same convention as for ck.
% zh        nh by 1    vector of values for the categories at the coordinates
%                      specified in ch. Categories are coded as integers ranging
%                      from 1 to nc, where nc is the number of categories.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nt by nt   cell array. Each cell is composed by a nc by nc cell array,
%                      where each cell is a nd by 1 vector of bivariate
%                      probability values between two categories at distances
%                      specified in d.
% nhmax     scalar     maximum number of locationss that are considered
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
%                      options(1) is taking the value 1 or 0 if the user wants or does
%                      not want to display the order number of the location which is
%                      currently processed, respectively.
%                      options(2) is the stopping criterion when fitting the maximum
%                      entropy table. Default value is equal to 1e-3 (see maxentropytablek.m).
%
% OUTPUT :
%
% pk        nk by nc   matrix of conditional probability values, where each
%                      column refers to a category and each line refers to an
%                      estimation location

%%%%%% Initialize the parameters

if nargin<8,
  options(1)=0;
  options(2)=1e-3;
end;

if iscell(ck),
  typek=ck{2};
  typeh=ch{2};
  ck=ck{1};
  ch=ch{1};
  nk=size(ck,1);
  ncat=size(Pmodel{1,1},1);
  ntype=max([typek;typeh]);
  pk=zeros(nk,ncat);
  ptype=cell(ntype,1);
  for i=1:ntype,
    ptype{i}=diag(probamodel2bitable(0,dmodel,Pmodel{i,i}))';
  end;
  for i=1:nk,
    pk(i,:)=ptype{typek(i)};
  end;
else
  nk=size(ck,1);
  nh=size(ch,1);
  typek=ones(nk,1);
  typeh=ones(nh,1);
  ntype=1;
  ncat=size(Pmodel,1);
  pk=kron(ones(nk,1),diag(probamodel2bitable(0,dmodel,Pmodel))');
  Pmodel={Pmodel};
end;

if options(1)==1,
  num2strnk=num2str(nk);
end;

%%%%%% Main loop starts here

for i=1:nk,

  ck0=ck(i,:);
  typek0=typek(i);
  [chlocal,ztypehlocal,dh,sumnhlocal]=neighbours(ck0,ch,[zh typeh],nhmax,dmax);
  zhlocal=ztypehlocal(:,1);
  typehlocal=ztypehlocal(:,2);

  if sumnhlocal>0,

    isdhzero=find(dh==0);
    if isempty(isdhzero),
      Pkh=maxentropytablek([ck0;chlocal],[typek0;typehlocal],dmodel,Pmodel,options(2));
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

