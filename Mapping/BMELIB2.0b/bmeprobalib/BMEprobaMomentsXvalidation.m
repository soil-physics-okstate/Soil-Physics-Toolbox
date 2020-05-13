function [momentsXval,info,MSE,MAE,ME]=BMEprobaMomentsXvalidation(ckOption,ch,cs,zh,softpdftype,...
        nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
      
% BMEprobaMomentsXvalidation  - Cross validation for BMEprobaMoments
%
% This program performs a cross-validation using BMEprobaMoments.
% The cross-validation estimates are obtained by removing each datum
% in turn from the dataset, and estimating that value on the basis 
% of the remaining data.  The data for which cross validation estimates
% are calculated may be only the hard data, or the hard and soft data.
% 
% 
% SYNTAX :
%
% [momentsXval,info,MSE,MAE,ME]=BMEprobaMomentsXvalidation(ckOption,ch,cs,zh,softpdftype,...
%        nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
%      
% INPUT :
%
% ckOption   scalar       Determine or which data the cross validation estimates are calculated
%                         1 to calculate cross-validation estimation at the hard data only, i.e. ck=ch
%                         2 to calculate cross-validation estimation at the hard and soft data, i.e. ck=[ch;cs];
% 
% See BMEprobaMoments.m for explanations for all the other input variables
%
% OUTPUT :
%
% momentsXval  nk by 3    matrix of the moments of the cross-validation estimate at each estimation point.
%                         The first column are the BME mean, 
%                         The second column has the BME estimation error variance (if options(8)>=2)
%                         The third column has the coeff of skewness (if options(8)>=3)
% info         nk by 1    vector of information about the computation in the local 
%                         neighbourhood around the estimation point :
%                         info=NaN : no computation, no hard or soft data
%                         info=0   : computation using BME with soft data 
%                         info=1   : dubious results, integration error above tolerance
%                                    when integrating over the soft pdf. Try to increase 
%                                    the value of options(3)
%                         info=3   : computation using only hard data
%                         info=4   : computation provides hard data value at estimation point
%                         info=10  : dubious results from integration routine. Try reducing
%                                    the number of soft data points ns, or option(3) or option(4)
% MSE          scalar     Mean Square Error between the cross validation estimate and the data
%                         where the mean value of the soft pdf for the soft data are used when ckOption=2
% MAE          scalar     Mean Absolute Error between the cross validation estimate and the data
% ME           scalar     Mean Error, i.e. the mean value of the cross validation estimate minus the data

%%%%%% Set some parameters

if isempty(ch)
  momentsXval=[];
  info=[];
  MSE=NaN;
  MAE=NaN;
  ME=NaN;
  return;
end

%%%%%% Error messages

[ch,zh,cs,nl,limi,probdens]=fixempty(ch,ch,zh,cs,nl,limi,probdens);
BMEprobaCheckArgs(ch,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax);

%%%%%% Set some parameters

noindex=~iscell(ch);       % test if there is an index for the variables

switch ckOption
case 1, 
  ck=ch;
case 2, 
  if noindex==1
    ck=[ch;cs];
  else  
    ck={[ch{1};cs{1}],[ch{2};cs{2}]};
  end
  otherwise, error('Bad value for ckOptions');
end

if isempty(ck)
  momentsXval=[];
  info=[];
  MSE=NaN;
  MAE=NaN;
  ME=NaN;
  return;
end

if noindex==1,
  nk=size(ck,1);           % nk is the number of estimation points
  nh=size(ch,1);           % nh is the number of hard data
else
  nk=size(ck{1},1);
  nh=size(ch{1},1);
end;

        
%%%%%% Main loop starts here

for ik=1:nk
  if ik <= nh
    chX=ch;
    if noindex==1
      chX(ik,:)=[];
    else
      chX{1}(ik,:)=[];
      chX{2}(ik,:)=[];
    end
    zhX=zh; zhX(ik)=[];
    csX=cs; nlX=nl; limiX=limi; probdensX=probdens;
  else
    chX=ch;
    zhX=zh;
    [dummy,csX,dummy,dummy,dummy,nlX,limiX,probdensX]=probasplit(cs,softpdftype,nl,...
      limi,probdens,ik-nh);
  end
  if noindex==1
    ckX=ck(ik,:);
  else  
    ckX={ck{1}(ik,:),ck{2}(ik,:)};
  end
  [momentsXval(ik,1:3),info(ik,1:3)]=BMEprobaMoments(ckX,chX,csX,zhX,softpdftype,nlX,limiX,probdensX,...
     covmodel,covparam,nhmax,nsmax,dmax,order,options);
end

%%%%%% Set some parameters

switch ckOption
case 1, zdata=zh;
case 2, zdata=[zh;proba2stat(softpdftype,nl,limi,probdens)];
end

MSE=mean((momentsXval(:,1)-zdata).^2);
MAE=mean(abs(momentsXval(:,1)-zdata));
ME=mean(momentsXval(:,1)-zdata);
