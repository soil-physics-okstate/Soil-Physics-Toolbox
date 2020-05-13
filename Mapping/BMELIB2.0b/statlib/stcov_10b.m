function [C,np]=stcov_10b(Zi,cMSi,tMEi,Zj,cMSj,tMEj,rLag,rLagTol,tLag,tLagTol,dist);

% stcov_10b                 - Estimate space/time cross covariance values from data, version 1.0b (Jan 1,2001)
% 
% Estimates the space/time cross covariance between variable Zi and
% variable Zj with measurements at fixed monitoring stations. The 
% monitoring stations for Zi are located at coordinates cMSi, and 
% measuring events are at times tMEi, and similarly cMSj and tMEj are
% the coordinates of monitoring sites and times of measuring events
% for Zj.
%
% SYNTAX :
%
% [C np]=stcov_10b(Zi,cMSi,tMEi,Zj,cMSj,tMEj,rLag,rLagTol,tLag,tLagTol,dist)
%
% INPUT :
%
% Zi      nMSi by nMEi matrix with measurements of Zi at the nMSi monitoring
%                      sites and nMEi measuring events. Zi may have NaN values.
% cMSi    nMSi by 2    matrix of spatial x-y coordinates for the nMSi monitoring
%                      sites
% tMEi    1 by nMEi    vector with the time of the measuring events
% Zj      nMSj by nMEj matrix with measurements of Zj at the nMSj monitoring
%                      sites and nMEj measuring events. Zj may have NaN values.
% cMSj    nMSj by 2    matrix of spatial x-y coordinates for the nMSj monitoring
%                      sites
% tMEj    1 by nMEj    vector with the time of the measuring events
% rLag    nr by 1      vector with the r lags
% rLagTol nr by 1      vector with the tolerance for the r lags
% tLag    1 by nt      vector with the t lags
% tLagTol nt by 1      vector with the tolerance for the t lags
% dist   cell          optional parameter giving the function to calculate distances.
%                      dist{1} is the name of the distance function, such that the
%                      distances between the points c1 and c2 is calculated using 
%                      [d]=eval([dist{1},'(c1,c2)']), where c1 and c2 are n1 by 2 and 
%                      n2 by 2 matrices, respectively, and d is a n1 by n2 matrix of 
%                      distances between c1 and c2.
%                      dist{2} is an optional parameter. Hence when length(dist)>1, then the
%                      distance matrix is calculated using [d]=eval([dist{1},'(c1,c2,dist{2})'])
%                      Note: The default is dist={'coord2dist'}
% 
% OUTPUT :
%
% C      nr by nt     matrix with the covariance values at lags rLag
%                     and tlags
% np     nr by nt     matrix with the number of pairs
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format


%
% Verify input
%
nMSi=size(Zi,1);
nMEi=size(Zi,2);
if size(cMSi,1)~=nMSi | size(cMSi,2)~=2
  error('cMSi must be a nMSi by 2 matrix'); 
end;
if size(tMEi,1)~=1 | size(tMEi,2)~=nMEi
  error('tMEi must be a 1 by nMEi vector'); 
end;
nMSj=size(Zj,1);
nMEj=size(Zj,2);
if size(cMSj,1)~=nMSj | size(cMSj,2)~=2
  error('cMSj must be a nMSj by 2 matrix'); 
end;
if size(tMEj,1)~=1 | size(tMEj,2)~=nMEj
  error('tMEj must be a 1 by nMEj vector'); 
end;
nr=length(rLag);
if length(rLagTol)~=nr
  error('rLag and rLagTol must have same length');
end;
nt=length(tLag);
if length(tLagTol)~=nt
  error('tLag and tLagTol must have same length');
end;

if nargin<11
  %disp('using default distances function: coord2dis');
  dist{1}='coord2dist';
end

%
% Calculate the distance between points
%
if length(dist)==1
  rMat=eval([dist{1},'(cMSi,cMSj)']);
else
  rMat=eval([dist{1},'(cMSi,cMSj,dist{2})']);
end;
iMSmat=kron( (1:nMSi)',ones(1,nMSj) );
jMSmat=kron( (1:nMSj) ,ones(nMSi,1) );
tMat=kron(tMEj,ones(nMEi,1))-kron(tMEi',ones(1,nMEj));
iMEmat=kron( (1:nMEi)',ones(1,nMEj) );
jMEmat=kron( (1:nMEj) ,ones(nMEi,1) );

for ir=1:nr
  ijMSmat=( (rLag(ir)-rLagTol(ir)<=rMat) & (rMat<=rLag(ir)+rLagTol(ir)) );
  nMSpairs=sum(sum(ijMSmat));
  if nMSpairs==0,
    C(ir,1:nt)=NaN;
    np(ir,1:nt)=0;  
  else
    iMSmat_r=iMSmat(ijMSmat);
    jMSmat_r=jMSmat(ijMSmat);
    for it=1:nt 
      ijMEmat=( (tLag(it)-tLagTol(it)<=tMat) & (tMat<=tLag(it)+tLagTol(it)) );
      nMEpairs=sum(sum(ijMEmat));
      if nMEpairs==0
        C(ir,it)=NaN;
        np(ir,it)=0;  
      else
        iMEmat_t=iMEmat(ijMEmat);
        jMEmat_t=jMEmat(ijMEmat);
        Xhead=Zi(iMSmat_r,iMEmat_t);
        Xhead=Xhead(:);
        Xtail=Zj(jMSmat_r,jMEmat_t);
        Xtail=Xtail(:);
        idxValid=(~isnan(Xhead) & ~isnan(Xtail));
        np(ir,it)=sum(idxValid);
        if np(ir,it)==0
          C(ir,it)=NaN;
        else
          Xhead=Xhead(idxValid);
          Xtail=Xtail(idxValid);
          headMean(ir,it)=mean(Xhead);  
          tailMean(ir,it)=mean(Xtail);
          headTailMean=mean( Xhead.*Xtail );
          C(ir,it)=headTailMean-headMean(ir,it)*tailMean(ir,it);    
          %disp(sprintf(...
          %  'r=%f, t=%f, np=%d, headMean=%.3g, tailMean=%.3g, C=%.3g',...
          %  rLag(ir),tLag(it),np(ir,it),headMean(ir,it),tailMean(ir,it),...
          %  C(ir,it)));
        end  
      end
    end
  end
end

