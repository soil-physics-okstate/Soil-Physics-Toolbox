function [zh]=simucatcondPdf(ch,cs0,ps0,dmodel,Pmodel,nsmax,dmax,options);

% simucatcondPdf         - conditional sequential simulation of categorical variables with pdf data
%                          (December 1, 2003);
%
% Implementation of the sequential conditional simulation
% method for a categorical variable based on the iterated use
% of conditional distributions, so that values are simulated
% one by one. The conditioning information are a set of 
% probability distribution functions. The simucatcondPdf.m
% function is based on the use of the BMEcatPdf.m function.
%
% For conditional simulation using both hard and soft information,
% see the note at the end of this help.
%
% SYNTAX :
%
% [zh]=simucatcondPdf(ch,cs0,ps0,dmodel,Pmodel,nsmax,dmax,options);
%
% INPUT :
%
% ch        nh by d    matrix of coordinates for the locations where categories
%                      have to be simulated. A line corresponds to the vector of 
%                      coordinates at a simulation location, so the number of
%                      columns corresponds to the dimension of the space. There
%                      is no restriction on the dimension of the space.
% cs0       ns0 by d   matrix of coordinates for locations where distributions
%                      of categories are known, with same convention as for ck.
% ps0       ns0 by nc  matrix of [0,1] real values for the probabilities of the
%                      categories at the coordinates specified in cs0. Each line
%                      correspond to a location and each column to a category,
%                      so that each line must sum up to 1.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of
%                      bivariate probability values between two categories at
%                      distances specified in dmodel.
% nsmax     scalar     maximum number of locations in the neighbourhood used
%                      for the simulation at the locations specified in ch.
% dmax      scalar     maximum distance between an estimation location and
%                      locations where categories have already been simulated.
%                      All locations separated by a distance smaller than dmax
%                      from a imulation location will be included in the simulation
%                      process for that location, whereas other locations are 
%                      neglected.
% options   1 by 3     vector of optional parameters that can be used if default
%                      values are not satisfactory (otherwise this vector can
%                      simply be omitted from the input list of variables), where :
%                      options(1)=1 for displaying the simulation location
%                      currently processed (default value is 0),
%                      options(2)=1 for a random selection of the visited simulation
%                      locations in the sequence (default value is 1).
%                      options(3) is the stopping criterion when fitting the
%                      maximum entropy table. Default value is equal to 1e-3
%                      (see maxentropytable.m).
%
% OUTPUT :
%
% zh        nh by 1    vector of codes for the simulated categories at the coordinates
%                      specified in ch. Categories are coded as integers ranging
%                      from 1 to nc, where nc is the number of categories.
%
% NOTE : This program can process any kind of mixture of hard and soft information
%        using the following conventions :
%
%        (1) for a hard data (the jth category is observed at the ith location),
%            the corresponding ith line ps0(i,:) has its jth probability value
%            equal to one and other probability values equal to 0.
%        (2) for a subset of n possible categories at the ith location, the
%            corresponding ith line ps0(i,:) will have probability values
%            equal to 1/n for the elements corresponding to these categories
%            whereas other elements have a zero probability value.

%%%%%% Remove locations where all non null probabilities are equal
%%%%%% (zero specific information content)

isnotOK=sum(abs(diff(ps0')))==0;
cs0(isnotOK,:)=[];
ps0(isnotOK,:)=[];

%%% Initialize the parameters

if nargin<8,
  options(1)=0;
  options(2)=1;
  options(3)=1e-3;
end;

nh=size(ch,1);
ns0=size(cs0,1);
ncat=size(Pmodel,1);

if options(1)==1,
  num2strnh=num2str(nh);
end;

if options(2)==1,
  index=randperm(nh);
  ch=ch(index,:);
end;

cs=[cs0;ch];
ps=[ps0;ones(nh,ncat)*NaN];

%%%%%% Main loop starts here

for i=ns0+1:ns0+nh,
    
  ck=cs(i,:);
  dmaxiter=0;
  pk=ones(1,ncat)*NaN;
  while isnan(sum(pk)),
    dmaxiter=dmaxiter+dmax;    
    [pk]=BMEcatPdf(ck,cs(1:i-1,:),ps(1:i-1,:),dmodel,Pmodel,nsmax,dmaxiter,[0 options(3)]);
  end;
  r=rand;
  catk=sum(r>[0 cumsum(pk)]); 
  ps(i,:)=zeros(1,ncat);
  ps(i,catk)=1;

  if options(1)==1,
    disp([num2str(i-ns0),'/',num2strnh]);
  end;
end;

[indexi,indexj]=find(ps(ns0+1:ns0+nh,:)==1);
indexij=sortrows([indexi indexj],1);
zh=indexij(:,2);

if options(2)==1,
  [ignore,index]=sort(index);
  zh=zh(index);
end;
