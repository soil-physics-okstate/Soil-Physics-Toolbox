function  [nl,limi,probdens]=simuprobabilistic(Zs,softpdftype,NV,WidthV,probdensV);

% simuprobabilistic         - simulate probabilistic data (Jan 1,2000)
%
% Given a set of hard data Zs and a template for the ns-dimensional
% pdf fs, provides the actual ns-dimensional pdf fs (i.e. the soft 
% probabilistic data) at the soft data points
% 
% SYNTAX :
%
% [nl,limi,probdens]=simuprobabilistic(Zs,softpdftype,NV,WidthV,probdensV);
%
% INPUT :
%
% Zs         ns by nSim   matrix of simulated hard data 
%softpdftype scalar       type of the ns-dimensional pdf fs: 
%                           1-Histogram, 2-Linear, 3-Grid Histogram, 
%                           4-Grid Linear (see helpsoftpdftype)
% NV         ns by 1      number of intervals in each dimension
% WidthV     matrix       interval widths
%                         If softpdftype=1 or 2, WidthV is a ns by 
%                           max(NV) matrix of interval widths. 
%                         If softpdftype=3 or 4 WidthV is a ns by 1 vector
%                           with the constant interval width for each dim.
% probdensV  matrix       values of fs along each dimension 
%                         If softpdftype=1 or 3 probdensV is a ns by max(NV)  
%                            matrix of fs value for the intervals along each dim.
%                         If softpdftype=2 or 4 probdensV is a ns by max(NV)+1
%                            matrix of fs value for the interval limits along
%                            each dimension
%
% OUTPUT :
%
% nl         ns by 1      nl=NV+1
% limi       array        limits of intervals for which the soft pdf is
%                         specified (see helpsoftpdftype)
%                         If softpdftype=1 or 2, limi is a ns by max(nl)  
%                           by nSim array of interval limits 
%                         If softpdftype=3 or 4 limi is a ns by 3 by nSim
%                           array with the lower limits, increment and upper 
%                           limit of intervals
% probdens   array        value of the soft pdf (probability density values)
%                           at the intervals limits given by limi
%                         If softpdftype=1 or 3, probdens is a ns by  
%                           max(nl)-1 by nSim array of interval probabilities
%                         If softpdftype=2 or 4 probdens is a ns by max(nl) 
%                           by nSim array of proba value at the interval limits. 

%
%Check inputs
%
[ns nSim]=size(Zs);

if (ns<1), error('must have at least one soft data point'); end;
if (nSim<1), error('must generate at one least simulation'); end;
if size(NV,1)~=ns
  error('NV must be a vector of length ns');
end;
NVmax=max(NV);

%
%  Check size of WidthV
%
if (softpdftype==1) | (softpdftype==2),
  if (size(WidthV,1)~=ns) | (size(WidthV,2)~=NVmax)
    error('WidthV must be a ns by max(NV) matrix');
  end
elseif (softpdftype==3) | (softpdftype==4),
  if (size(WidthV,1)~=ns) | (size(WidthV,2)~=1)
    error('WidthV must be a ns by 1 matrix');
  end;
else
  error('Unacceptable value for softpdftype');
end;  

%
% Check size of probdensV
%
if (softpdftype==1) | (softpdftype==3),
  if (size(probdensV,1)~=ns) | (size(probdensV,2)~=NVmax)
    error('probdensV must be a ns by max(NV) matrix');
  end
elseif (softpdftype==2) | (softpdftype==4),
  if (size(probdensV,1)~=ns) | (size(probdensV,2)~=NVmax+1)
    error('probdensV must be a ns by max(NV)+1 matrix');
  end
end;

nl=NV+1;

%
% Calculate limiV, the limits of the intervals for V starting at 0
% and UV, the limits of intervals for the uniformly distributed U
%

%
% Do histogram case first
%
if (softpdftype==1) | (softpdftype==3),
  if (softpdftype==3),
    WidthV=kron(WidthV,ones(1,max(NV)));
  end;
  for is=1:ns
    limiV(is,1)=0;
    UV(is,1)=0;
    for iv=2:nl(is),
        limiV(is,iv)=limiV(is,iv-1)+WidthV(is,iv-1);
        UV(is,iv)=UV(is,iv-1)+WidthV(is,iv-1)*probdensV(is,iv-1);
    end;
  end;

%
%  Check that WidthV and probdensV are properly normalized
%
  Eps=1e-10;
  for is=1:ns,
    if (UV(is,nl(is))<1.0-Eps) | (UV(is,nl(is))>1.0+Eps),
      error('WidthV and probdensV not normalized probability');
    end;
  end;

%
%  Generate U uniformly between 0 and 1
%  then transform U to V so that V has a distribution described
%  by WidthV and probdensV
%
  U=rand(ns,nSim);
  for is=1:ns
    UVgrid=UV(is,1:nl(is));
    for iSim=1:nSim
      iv=sum(UVgrid<U(is,iSim));
      V(is,iSim)=limiV(is,iv)+(U(is,iSim)-UVgrid(iv))/probdensV(is,iv);
    end;
  end;

  Y=Zs-V;
  if (softpdftype==1),
    limi=zeros(ns,max(nl),nSim)*NaN;
    for is=1:ns,
      for iv=1:nl(is),
        limi(is,iv,:)=Y(is,:)+limiV(is,iv);
      end;
    end;
  elseif(softpdftype==3),
    limi=zeros(ns,3,nSim)*NaN;
    for is=1:ns,
      limi(is,1,:)=Y(is,:);
      limi(is,2,:)=WidthV(is,1);
      limi(is,3,:)=Y(is,:)+NV(is)*WidthV(is,1)+Eps;
    end;
  end;
  probdens=zeros(ns,max(nl)-1,nSim)*NaN;
  for is=1:ns,
    for iv=1:NV(is),
      probdens(is,iv,:)=probdensV(is,iv);
    end;
  end;
elseif (softpdftype==2) | (softpdftype==4),
  if (softpdftype==4),
    WidthV=kron(WidthV,ones(1,max(NV)));
  end;
  for is=1:ns
    limiV(is,1)=0;
    UV(is,1)=0;
    for iv=2:nl(is),
     limiV(is,iv)=limiV(is,iv-1)+WidthV(is,iv-1);
     UV(is,iv)=UV(is,iv-1)+WidthV(is,iv-1)*0.5*(probdensV(is,iv-1)+probdensV(is,iv));
    end;
  end;

%
%  Check that WidthV and probdensV are properly normalized
%
  Eps=1e-10;
  for is=1:ns,
    if (UV(is,nl(is))<1.0-Eps) | (UV(is,nl(is))>1.0+Eps),
      error('WidthV and probdensV not normalized probability');
    end;
  end;

%
%  Generate U uniformly between 0 and 1
%  then transform U to V so that V has a distribution described
%  by WidthV and probdensV
%
  U=rand(ns,nSim);
  for is=1:ns
    UVgrid=UV(is,1:nl(is));
    for iSim=1:nSim      
      iv=sum(UVgrid<U(is,iSim));
      fsp=(probdensV(is,iv+1)-probdensV(is,iv))/(limiV(is,iv+1)-limiV(is,iv));
      if (fsp~=0),
        fso=probdensV(is,iv)-limiV(is,iv)*fsp;
        V(is,iSim)=1/fsp*(-fso+sqrt(probdensV(is,iv)^2+2*fsp*(U(is,iSim)-UVgrid(iv))));
      elseif (probdensV(is,iv)~=0)
        V(is,iSim)=limiV(is,iv)+(U(is,iSim)-UVgrid(iv))/probdensV(is,iv);
      else
        V(is,iSim)=limiV(is,iv);;
      end;
    end;
  end;

  Y=Zs-V;
  if (softpdftype==2),
    limi=zeros(ns,max(nl),nSim)*NaN;
    for is=1:ns,
      for iv=1:nl(is),
        limi(is,iv,:)=Y(is,:)+limiV(is,iv);
      end;
    end;
  elseif(softpdftype==4),
    limi=zeros(ns,3,nSim)*NaN;
    for is=1:ns,
      limi(is,1,:)=Y(is,:);
      limi(is,2,:)=WidthV(is,1);
      limi(is,3,:)=Y(is,:)+NV(is)*WidthV(is,1)+Eps;
    end;
  end;
  probdens=zeros(ns,max(nl),nSim)*NaN;
  for is=1:ns,
    for iv=1:nl(is),
      probdens(is,iv,:)=probdensV(is,iv);
    end;
  end;
end;
