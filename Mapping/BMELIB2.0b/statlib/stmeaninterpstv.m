function [mstI]=stmeaninterpstv(cMS,tME,ms,mt,pI);

% stmeaninterpstv           - Interpolate s/t mean onto s/t vector data
%
% Calculates space/time interpolation of mean trend values using 
% an additive separable s/t mean trend model.
% Using the spatial smoothed mean component ms  and the temporal
% smoothed mean component mt obtained by the function mst from
% from measurements at monitoring sites coordinates cMS and measuring
% events times tME, the s/t mean trend at cooridnate cI and time cI
% is given by
%       mst(s,t)= ms(cI) + mt(tI) - mean(mt),
% where ms(cI) is interpolated from ms at the monitoirng sites, and 
% mt(tI) is interpolated from the mt at the measuring events tME.
%
% SYNTAX :
%
% [mstI]=stmeaninterpstv(cMS,tME,ms,mt,pI)
%
% INPUT :
%
% cMS    nMS by 2   matrix of spatial coordinates for the Measuring Sites
% tME    1 by nME   vector of the times of the Measuring Events
% ms     nMS by 1   vector with the smoothed spatial mean at the Measuring  
%                   Sites ms is usually calculated with the mst function
% mt     1 by nME   vector with the smoothed temporal mean at the Measuring Events
%                   mt is usually calculated with the mst function
% pI     ncI by 3   matrix of space/time coordinates for the Interpolation points
%
% OUTPUT :
%
% mstI   ncI by 1   vector of interpolated values of the s/t mean of Z at the 
%                   interpolations. 
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format


%
% Check input data
%
nMS=size(cMS,1);
if size(cMS,2)~=2, error('cMS must have two columns'); end;
if size(tME,1)~=1, error('tME must have one row'); end;
nME=size(tME,2);
if size(ms,1)~=nMS, error('ms must the have same number of rows as cMS'); end;
if size(ms,2)~=1, error('ms must have one column'); end;
if size(mt,1)~=1, error('mt must have one row'); end;
if size(mt,2)~=nME, error('mt must the have same length as tME'); end;
if isempty(pI),
  mstI=[];
  return
end
npI=size(pI,1);
if size(pI,2)~=3, error('pI must have three columns'); end;

%
% Interpolate the spatial data
%
xMS=cMS(:,1);
yMS=cMS(:,2);
xI=pI(:,1);
yI=pI(:,2);

addPadding=0;                  % Note: I have to pad xI because griddata doe not
if length(xI)==1;              % work when xI is a scalar
  addPadding=1;
  xI=[xI;0];   
  yI=[yI;0];     
end;

msI=griddata(xMS,yMS,ms,xI,yI);  % Do the linear interpolation of ms values

if addPadding                       % Note: remove the padding from msI
  msI=msI(1);   
end;

msI(isnan(msI))=mean(ms);        % set the nan values (outside of linear
                                   % interpolation area) to the mean of ms

%
% Add the temporal component:   m(s,t)= ms(s) + mt(t) - mean(mt)
%
tI=pI(:,3);
mtI=interp1([min([tI;tME'])-1;tME';max([tI;tME'])+1],[mt(1);mt';mt(end)],tI);
mstI=msI+mtI-mean(mt);

