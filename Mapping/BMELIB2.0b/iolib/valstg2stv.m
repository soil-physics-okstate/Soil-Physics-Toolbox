function [ch,z]=valstg2stv(Z,cMS,tME);

% valstg2stv                - Converts values from s/t grid coord to s/t vector coord (Jan 1, 2001)
%
% Converts the coordinates and values of a space/time variable 
% from a grid format (i.e. the variable Z is given as a nMS by nME matrix  
% corresponding to nMS Measuring Sites and nME Measuring Events),
% to a s/t vector format (i.e. the variable z is listed as a vector of nMS*nME values,
% corresponding to points with space/time coordinates, where the spatial coordinate
% cycle quicker than the time coordinates).
%
% SYNTAX :
%
% [ch,z]=valstg2stv(Z,cMS,tME);
%
% INPUT :
%
% Z        nMS by nME       matrix of values for the variable Z corresponding to 
%                           nMS Monitoring Sites and nME Measuring Event
% cMS      nMS by 2         matrix of 2D spatial coordinates for the nMS Measuring Sites
% tME      1 by nME         vector of times of the tME Measuring Events
%
% OUTPUT :
%
% ch       nMS*nME by 3     matrix of space time coordinates, listing the space/time
%                           locations of the points corresponding to nMS Monitoring Sites
%                           and nME Measuring Event (space cycles quicker then time) 
% z        nMS*nME by 1     vector of values for the variable Z corresponding to the
%                           s/t points ch
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format


nMS=size(cMS,1);
nME=length(tME);
if size(cMS,2)~=2,
  error('cMS must have 2 columns');
end;
if size(tME,1)~=1,
  error('tME must be a row vector');
end;
if size(Z,1)~=nMS | size(Z,2)~=nME
  error('Z must be a nMS by nME matrix');
end;
z=reshape(Z,nMS*nME,1);
ch=[kron(ones(nME,1),cMS(:,1)) kron(ones(nME,1),cMS(:,2))...
    kron(tME',ones(nMS,1))];

