function [Nl,Limi,Probdens]=probastv2stg(cs,softpdftype,nl,limi,probdens,cMS,tME);

% probastv2stg              - Converts probabilistic data from s/t vector coord to s/t grid coord (Jan 1, 2001)
%
% Converts the probabilisitic data of a space/time variable from a s/t vector
% format (i.e. the soft data is listed as a vector of nMS*nME soft pdf,
% corresponding to points with space/time coordinates, where the spatial coordinate
% cycle quicker than the time coordinates),
% to a s/t grid format (i.e. the soft data is given as a nMS by nME matrix  
% of soft pdf corresponding to nMS Measuring Sites and nME Measuring Events).
%
% SYNTAX :
%
% [Nl,Limi,Probdens]=probastv2stg(cs,softpdftype,nl,limi,probdens,cMS,tME);
%
% INPUT :
%
% cs          nMS*nME by 3    matrix of space time coordinates, listing the space/time
%                             locations of the points corresponding to nMS Monitoring Sites
%                             and nME Measuring Event (space cycles quicker then time) 
% softpdftype scalar          type of soft pdf (see probasyntax) representing the probabilitic  
%                             soft data at the nMS Monitoring Sites and nME Measusing Events.
% nl          nMS*nME by 1    vector of number of limits  nl (see probasyntax)
%                             corresponding to the s/t points cs
% limi        nMS*nME by l    matrix representing the limits  values (see probasyntax)
%                             corresponding to the s/t points cs
% probdens    nMS*nME by p    matrix representing the probability densities (see probasyntax)
%                             corresponding to the s/t points cs
% cMS         nMS by 2        matrix of 2D spatial coordinates for the nMS Measuring Sites
% tME         1 by nME        vector of times of the tME Measuring Events
%
% OUTPUT :
%
% Nl          nMS by nME      matrix representing the number of limits  nl (see probasyntax)
%                             at the nMS Monitoring Sites and nME Measusing Events.
% Limi        nMS by nME by l array reprenting the nMS by nME matrix of vectors limi
%                             which have the limits of intervals (see probasyntax)
% Probdens    nMS by nME by p array reprenting the nMS by nME matrix of vectors probdens
%                             which have the probability density values (see probasyntax)
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format



nMS=size(cMS,1);
nME=length(tME);
ns=size(cs,1);
if size(nl,1)~=ns,
  error('cs and nl must have the same number of rows');
end;
if length(nl)~=nMS*nME,
  error('length(nl) must be equal to nMS*nME');
end;
if sum(cs(1:nMS,1)==cMS(:,1))~=nMS | sum(cs(1:nMS,2)==cMS(:,2))~=nMS
  if max(cs(1:nMS,1)-cMS(:,1))>1e-5*max(cMS(:,1)) ...
      & max(cs(1:nMS,2)-cMS(:,2))>1e-5*max(cMS(:,2))
    warning('cs(1:nMS,1:2) and cMS(:,1:2) do not match');
  end;  
end;
if sum(cs(1:nMS:nMS*nME,3)==tME')~=nME
  warning('cs(1:nMS:nMS*nME,3) and tME do not match');
end;

n3Limi=size(limi,2);
n3Probdens=size(probdens,2);
Nl=NaN*ones(nMS,nME,1);
Limi=NaN*ones(nMS,nME,n3Limi);
Probdens=NaN*ones(nMS,nME,n3Probdens);
for iME=1:nME
  for iMS=1:nMS
    Nl(iMS,iME)=nl((iME-1)*nMS+iMS,1);
    Limi(iMS,iME,1:n3Limi)=limi((iME-1)*nMS+iMS,:);
    Probdens(iMS,iME,1:n3Probdens)=probdens((iME-1)*nMS+iMS,:);
  end;
end;

