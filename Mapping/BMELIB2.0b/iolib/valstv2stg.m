function [Z,cMS,tME,nanratio]=valstv2stg(p,z,cMS,tME);

% valstv2stg                - Converts values from s/t vector coord to s/t grid coord (Jan 1, 2001)
%
% Converts the values of a space/time variable from a s/t vector 
% format (i.e. the variable z is listed as a vector of n values)
% to a grid format (i.e. the variable Z is given as a nMS by nME matrix 
% corresponding to nMS Monitoring Sites and nME Measuring Events).
% Use help stgridsyntax for information on the s/t grid format.
%
% SYNTAX :
%
% [Z,cMS,tME,nanratio]=valstv2stg(ch,z,cMS,tME);
%
% INPUT :
%
% ch         n by d+1     matrix of space/time coordinates for spatial domain of dimension d
% z          n by 1       vector of field value at coordinate ch
% cMS        nMS by 2     optional matrix of spatial coordinates for the nMS Measuring Sites
% tME        1 by nME     optional vector of times of the tME Measuring Events
%
% OUTPUT :
%
% Z          nMS by nME   matrix of values for the variable Z corresponding to 
%                         nMS Monitoring Sites and nME Measuring Event
% cMS        nMS by d     matrix of spatial coordinates for the nMS Measuring Sites
% tME        1 by nME     vector of times of the tME Measuring Events
% nanratio   scalar       ratio of the NaNs in Z (0<=nanratio<=1) 
%
% NOTE : 
%
% cMS and tME can be provided as input if they are both known.  In that case ch 
% must be a nMS*nME by 3 matrix of the points corresponding to nMS Monitoring Sites
% and nME Measuring Events, listed with space cycling quicker then time.
% 
% See also stgridsyntax for help on s/t grid format


if nargin>=4

  nMS=size(cMS,1);
  nME=length(tME);
  nh=size(p,1);
  if size(z,1)~=nh,
    error('z and p must have the same number of rows');
  end;
  if length(z)~=nMS*nME,
    error('length(z) must be equal to nMS*nME');
  end;
  if sum(p(1:nMS,1)==cMS(:,1))~=nMS
    error('p(1:nMS,1) and cMS(:,1) do not match');
  end;
  if sum(p(1:nMS,2)==cMS(:,2))~=nMS
    error('p(1:nMS,1) and cMS(:,1) do not match');
  end;
  if sum(p(1:nMS:nMS*nME,3)==tME')~=nME
    error('p(1:nMS:nMS*nME,3) and tME do not match');
  end;
  Z=reshape(z,nMS,nME);
  
else 
  
  [n,m]=size(p);
  if m<2, error('p must have at least two columns'); end
  d=m-1;
  if size(z,1)~=n, error('p and z must have the same number of rows'); end
  if size(z,2)~=1, error('z must have only one column'); end
 
  if isdupli(p)
    disp('warning in valstv2stg: duplicated data were averaged');
    [p,z] = avedupli(p,z);
  end
  cMS=unique(p(:,1:d),'rows');
  nMS=size(cMS,1);
  tME=unique(p(:,end))';
  nME=length(tME);  
  
  Z(1:nMS,1:nME)=NaN;
  for iMS=1:nMS
    idxMS= cMS(iMS,1)==p(:,1);
    for id=2:d
      idxMS = idxMS & cMS(iMS,id)==p(:,id);
    end
    idxMS=find(idxMS);
    if length(idxMS)>0
      ti=p(idxMS,end);                    % ti are measurement times for station iMS
      Zi=z(idxMS);                        % Zi are the measurements at ti for station iMS
      [ti2,idx1,idx2]=intersect(tME',ti); %   ti2 == tME(idx1) == ti(idx2)
      Z(iMS,idx1)=Zi(idx2)';              % hence  Z(iMS,idx1) =  Zi(idx2)
    end
  end
end

nanratio=sum(isnan(Z(:)))/(nMS*nME);
