function [a,b]=simuinterval(Zs,I,Vmethod,Vvalue);

% simuinterval              - Simulate soft interval data (Jan 1,2000)
% 
% Simulate the lower and upper bounds of soft interval data using
% a set of known Z values at the soft data points
%
% SYNTAX :
%
% [a,b]=simuinterval(Zs,I,Vmethod,Vvalue)
%
% INPUT :
%
% Zs         nh by nSim   matrix of Z values at the soft data points
% I          vector       value(s) for the interval length
% Vmethod    integer      integer between 0 and 3 indicating the method
%                         used to generate intervals as follow
%    Vmethod=0, 1 or 2: I is a scalar, a=Zs+V-I/2, b=Zs+V+I/2
%      Vmethod=0: V=0
%      Vmethod=1: V is uniformly distributed between -I/2 and I/2
%      Vmethod=2: V=Vvalue
%    Vmethod=3: I is a vector of increasing values I(i), a=I(i)<Zs<b=I(i+1) 
%
% OUTPUT :
%
% a          ns by nSim   vector of values for the lower bound of the intervals
% b          ns by nSim   vector of values for the upper bound of the intervals
%                         a and b are such that a<=Zh<=b

[ns nSim]=size(Zs);
if Vmethod==0
  V=0;
elseif Vmethod==1
  V=(rand(ns,nSim)-0.5).*I;
elseif Vmethod==2
  if (Vvalue<-I/2) | (Vvalue>I/2)
    error('Vvalue must be in interval [-I/2,I/2]');
  end
  V=Vvalue;
elseif Vmethod==3
  mInf=-5;pInf=5;
  ZkGrid=I;
  ZkGridInf=[mInf,ZkGrid,pInf];
  for iSim=1:nSim
    for i=1:ns
      ia=sum(ZkGrid<Zs(i,iSim));  
      a(i,iSim)=ZkGridInf(ia+1);
      b(i,iSim)=ZkGridInf(ia+2);      
    end
  end  
else
  error(sprintf('Vmethod=%n is not a valid option',Vmethod));
end

if (Vmethod==0) | (Vmethod==1) | (Vmethod==2)
  a=Zs+V-I/2;
  b=Zs+V+I/2;
end;
