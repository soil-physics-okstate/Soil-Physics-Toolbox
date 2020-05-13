function [softpdftype,nl,limi,probdens]=probacat(softpdftype1,nl1,limi1,probdens1,...
   softpdftype2,nl2,limi2,probdens2);

% probacat                  - Concatenates two sets of probabilistic soft data (Jan 1, 2001)
%
% Concatenates the two sets of probabilistic soft data
% with the same soft probabilitic type
%
% SYNTAX :
%
% [softpdftype,nl,limi,probdens]=probacat(softpdftype1,nl1,limi1,probdens1,...
%   softpdftype2,nl2,limi2,probdens2);
%
% INPUT:
%
% softpdftype1,nl1,limi1,probdens1    first set of probabilistic soft data with ns1 points
%                                     for more info of soft probabilistic data see help BMEprobaMode.
% softpdftype2,nl2,limi2,probdens2    second set set of probabilistic soft data with ns2 points
%  
%
% OUTPUT :
%
% softpdftype,nl,limi,probdens        union of the two sets of soft probabilistic data (ns1+ns2 points)

if softpdftype1~=softpdftype2,
  error('softpdftype1 must be equal to softpdftype2'); 
end;
softpdftype=softpdftype1;

ns1=length(nl1);
ns2=length(nl2);
nl=[nl1;nl2];
ml1=size(limi1,2);
ml2=size(limi2,2);
limi(1:ns1,1:ml1)=limi1;
limi(ns1+1:ns1+ns2,1:ml2)=limi2;
mp1=size(probdens1,2);
mp2=size(probdens2,2);
probdens(1:ns1,1:mp1)=probdens1;
probdens(ns1+1:ns1+ns2,1:mp2)=probdens2;

