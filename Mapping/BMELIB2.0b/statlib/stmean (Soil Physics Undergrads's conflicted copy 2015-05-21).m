function [ms,mss,mt,mts]=stmean(Z,cMS,idMS,tME,p);

% stmean                    - Estimates space/time mean values from measurements (Jan 1,2001)
%
% Assuming a separable additive space/time mean trend model, this 
% function calculates the spatial mean component ms and temporal mean 
% component mt of space/time random field Z, using measurements at 
% fixed measuring  sites cMS and fixed measuring events tME. 
% The spatial mean component ms is obtained by averaging the measurements
% at each measuring sites. Then a smoothed spatial mean component mss
% is obtained by applying an exponential spatial filter to ms.
% Similarly mt is obtained by averaging the measurement for each
% measuring event, and a smoothed temporal mean component mts is obtained
% by applying an exponential temporal filter to mt.
% Then the space/time mean trend is simply given by 
% mst(s,t)=mss(s)+mts(t)-mean(mts)
%
% SYNTAX :
%
% [ms,mss,mt,mts]=stmean(Z,cMS,idMS,tME,p);
%
% INPUTS :
%  
%  Z     nMS by nME matrix of measurements for Z at the nMS monitoring
%                   sites and nME measuring events. Z may have NaN values.
%  cMS   nMS by 2   matrix of spatial x-y coordinates for the nMS monitoring
%                   sites
%  idMS  nMS by 1   vector with a unique id number for each monitoring site
%  tME   1 by nME   vector with the time of the measuring events
%  p     1 by 5     parameters to smooth the spatial and temporal average
%                   p(1)=dNeib  distance (radius) of spatial neighborhood
%                   p(2)=ar     spatial range of exponential smoothing function
%                   p(3)=tNeib  time (radius) of temporal neighborhood
%                   p(4)=at     temporal range of exponential smoothing function
%                   p(5)=tloop  is an optional input used for temporal smoothing.
%                        When tloop>0, the measuring events are looped in a
%                    cycle of duration tloop.
%
% OUTPUT :
%
%  ms    nMS by 1   vector of spatial average
%  mss   nMS by 1   vector of smoothed spatial average
%  mt    1 by nMS   vector of temporal average
%  mts   nMS by 1   vector of smoothed temporal average
%
% NOTE : 
% 
% Use help stgridsyntax for help on s/t grid format



%
% Check the input arguments
%
nMS=size(Z,1);
nME=size(Z,2);
if size(cMS,1)~=nMS | size(cMS,2)~=2
  error('cMS must be a nMS by 2 matrix'); 
end;
if size(idMS,1)~=nMS | size(idMS,2)~=1
  error('idMS must be a nMS by 1 vector'); 
end;
if size(tME,1)~=1 | size(tME,2)~=nME
  error('tME must be a 1 by nME vector'); 
end;

%
% set parameters and variables
%
dNeib=p(1);      % Radius use to select local neighborhood to average Z
ar=p(2);         % spatial range of exponential smoothing function
tNeib=p(3);         % Radius use to select local neighborhood to average PM10
at=p(4);           % temporal range of exponential smoothing function
if length(p)>4
  tloop=p(5);
else
  tloop=0;
end; 
xMS=cMS(:,1);
yMS=cMS(:,2);

%
%  Calculate the spatial averages ms
%
for iMS=1:nMS
  ms(iMS)=mean(Z(iMS,~isnan(Z(iMS,:))));
end;

%
%  smooth the spatial average with an exponential filter to get mss
%
for iMS=1:nMS
  d=sqrt((xMS-xMS(iMS)).^2+(yMS-yMS(iMS)).^2);
  idxMSloc=find(d<dNeib);
  mss(iMS)=sum(ms(idxMSloc).*exp(-d(idxMSloc)/ar)');
  mss(iMS)=mss(iMS)/(sum(exp(-d(idxMSloc)/ar)));
end;

%
%  Calculate the temporal averages mt
%
for iME=1:nME
  mt(iME)=mean(Z(~isnan(Z(:,iME)),iME));
end;

%
%  smooth the temporal average with an exponential filter to get mts
%
for iME=1:nME
  t=abs(tME-tME(iME));
  if tloop>0
    t=min(t,tloop-t); 
  end;
  idxMEloc=find(t<tNeib);
  mts(iME)=sum(mt(idxMEloc).*exp(-t(idxMEloc)/at));
  mts(iME)=mts(iME)/(sum(exp(-t(idxMEloc)/at)));
end;

ms=ms';
mss=mss';