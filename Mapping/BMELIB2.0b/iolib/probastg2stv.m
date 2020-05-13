function [cs,nl,limi,probdens]=probastg2stv(softpdftype,Nl,Limi,Probdens,cMS,tME);

% probastg2stv              - Converts probabilistic data from s/t grid coord to s/t vector coord (Jan 1, 2001)
%
% Converts the coordinates and probabilisitic data of a space/time variable 
% from a grid format (i.e. the soft data is given as a nMS by nME matrix  
% of soft pdf corresponding to nMS Measuring Sites and nME Measuring Events),
% to a s/t vector format (i.e. the soft data is listed as a vector of nMS*nME soft pdf,
% corresponding to points with space/time coordinates, where the spatial coordinate
% cycle quicker than the time coordinates).
%
% SYNTAX :
%
% [cs,nl,limi,probdens]=probastg2stv(softpdftype,Nl,Limi,Probdens,cMS,tME);
%
% INPUT :
%
% softpdftype scalar          type of soft pdf (see probasyntax) representing the probabilitic  
%                             soft data at the nMS Monitoring Sites and nME Measusing Events.
% Nl          nMS by nME      matrix representing the number of limits  nl (see probasyntax)
%                             at the nMS Monitoring Sites and nME Measusing Events.
% Limi        nMS by nME by l array reprenting the nMS by nME matrix of vectors limi
%                             which have the limits of intervals (see probasyntax)
% Probdens    nMS by nME by p array reprenting the nMS by nME matrix of vectors probdens
%                             which have the probability density values (see probasyntax)
% cMS         nMS by 2        matrix of 2D spatial coordinates for the nMS Measuring Sites
% tME         1 by nME        vector of times of the tME Measuring Events
%
% OUTPUT :
%
% cs          nMS*nME by 3    matrix of space time coordinates, listing the space/time
%                             locations of the points corresponding to nMS Monitoring Sites
%                             and nME Measuring Event (space cycles quicker then time) 
% nl          nMS*nME by 1    vector of number of limits  nl (see probasyntax)
%                             corresponding to the s/t points cs
% limi        nMS*nME by l    matrix representing the limits  values (see probasyntax)
%                             corresponding to the s/t points cs
% probdens    nMS*nME by p    matrix representing the probability densities (see probasyntax)
%                             corresponding to the s/t points cs
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format

nMS=size(cMS,1);
nME=length(tME);


nMS=size(cMS,1);
nME=length(tME);
cs=[kron(ones(nME,1),cMS(:,1)) kron(ones(nME,1),cMS(:,2))...
    kron(tME',ones(nMS,1))];
n2limi=size(Limi,3);
n2probdens=size(Probdens,3);
nl=NaN*ones(nMS*nME,1);
limi=NaN*ones(nMS*nME,n2limi);
probdens=NaN*ones(nMS*nME,n2probdens);
for iME=1:nME
  for iMS=1:nMS
    nl((iME-1)*nMS+iMS,1)=Nl(iMS,iME);
    limi((iME-1)*nMS+iMS,1:n2limi)=Limi(iMS,iME,:);
    probdens((iME-1)*nMS+iMS,1:n2probdens)=Probdens(iMS,iME,:);
  end;
end;