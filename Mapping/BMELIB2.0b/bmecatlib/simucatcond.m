function [zh]=simucatcond(ch,ch0,zh0,dmodel,Pmodel,nhmax,dmax,options);

% simucatcond            - conditional sequential simulation of categorical variables with hard data
%                          (December 1, 2003);
%
% Implementation of the sequential conditional simulation
% method for a categorical variable based on the iterated use
% of conditional distributions, so that values are simulated
% one by one.The conditioning information are a set of 
% observed categories.The simucatcond.m function is based on
% the use of the BMEcatHard.m function. For conditional
% simulations using both hard and soft information, see
% the note at the end of this help
%
% SYNTAX :
%
% [zh]=simucatcond(ch,ch0,zh0,dmodel,Pmodel,nhmax,dmax,options);
%
% INPUT :
%
% ch        nh by d    matrix of coordinates for the locations where categories
%                      have to be simulated. A line corresponds to the vector of 
%                      coordinates at a simulation location, so the number of
%                      columns corresponds to the dimension of the space. There
%                      is no restriction on the dimension of the space.
% ch0       nh0 by d   matrix of coordinates for the hard conditioning
%                      categories, with the same convention as for ch.
% zh0       nh0 by 1   vector of code for the hard conditioning categories at the
%                      coordinates specified in ch0.
% dmodel    nd by 1    vector of values for distances for which the bivariate
%                      probabilities between categories have been modeled. 
% Pmodel    nc by nc   cell array, where each cell is a nd by 1 vector of
%                      bivariate probability values between two categories at
%                      distances specified in dmodel.
% nhmax     scalar     maximum number of locations in the neighbourhood used
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

%%% Error message

index=findpairs(ch,ch0);
if ~isempty(index),
  error('ch and ch0 cannot contain identical coordinates');
end;

%%% Initialize the parameters

if nargin<8,
  options(1)=0;
  options(2)=1;
  options(3)=1e-3;
end;

nh=size(ch,1);
nh0=size(ch0,1);
ncat=size(Pmodel,1);

if options(1)==1,
  num2strnh=num2str(nh);
end;

if options(2)==1,
  index=randperm(nh);
  ch=ch(index,:);
end;

ch=[ch0;ch];
zh=[zh0;ones(nh,1)*NaN];

%%%%%% Main loop starts here

for i=nh0+1:nh0+nh,
    
  ck=ch(i,:);
  dmaxiter=0;
  pk=ones(1,ncat)*NaN;
  while isnan(sum(pk)),
    dmaxiter=dmaxiter+dmax;    
    [pk]=BMEcatHard(ck,ch(1:i-1,:),zh(1:i-1),dmodel,Pmodel,nhmax,dmaxiter,[0 options(3)]);
  end;

  zh(i)=sum(rand>[0 cumsum(pk)]);

  if options(1)==1,
    disp([num2str(i-nh0),'/',num2strnh]);
  end;
end;

zh=zh(nh0+1:nh0+nh);
if options(2)==1,
  [ignore,index]=sort(index);
  zh=zh(index);
end;