function MVNLIBtest;

% MVNLIBtest                - Test for the mvnlib directory (Jan 1, 2001)
%
% Performs a serie of tests of the functions in the mvnlib directory. 
%
% SYNTAX :
%
% MVNLIBtest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the mvnAG1
%

disp(' ');
disp('---------------------------- ');
disp('Testing the mvnAG1 function');
disp(' ');

%
% Generating the covariance matrix C
% 

N=4;
rand('state',1);
C=rand(N,N);C=C'*C;
disp('Covariance matrix C=');
disp(C);

%
% Generating the lower and upper bound vectors a and b
%
a=zeros(N,1);
softInterv=1;
b=softInterv*ones(N,1);
disp('Lower bound of integration domain a=');
disp(a');
disp('Upper bound of integration domain b=');
disp(b');

%
% Testing the mvnAG1 function.
%
maxpts=1000000;
aEps=0;rEps=0.001;
disp(sprintf('method     name           value          error      Info CPUtime(sec)'));
t0=clock;[P Err Info]=mvnAG1(a,b,C,maxpts,aEps,rEps);ts=etime(clock,t0);
disp(sprintf('mvnAG1     P=E[1]  %15.10f %15.10f %3d %10.5f',P,Err,Info,ts));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the uvProNR and uvMomVecNR functions
%

disp(' ');
disp('---------------------------- ');
disp('Testing the uvProNR and uvMomVecNR functions');
disp(' ');

%
% setting variance C, mean, bounds a and b, and soft pdf
%
C=.01;
Mean=0.5;
disp(['Variance C=' num2str(C)]);
disp(['Mean=' num2str(Mean)]);
a=0;b=1;
disp(['Lower bound of integration domain a=' num2str(a)]);
disp(['Upper bound of integration domain b=' num2str(b)]);
softPdfType=1;Nl=2;Limi=[a b];Prob=1/(b-a);

%
% Setting some integration parameters
%
aEps=0;rEps=1e-8;

%
% testing uvProNR
%
[Val Err Info]=mvPro(softPdfType,Nl,Limi,Prob,Mean,C,maxpts,aEps,rEps);
disp(sprintf('method     name                 value          error     info'));
disp(sprintf('mvPro      P=E[1]        %15.10f %15.10g %3d',Val,Err,Info));

%
% testing uvMomVecNR
%
As=[0 1 1 1]';
Bs=[1 0 0 0];
P=[1 1 2 3];
nMom=length(P);
[Val Err Info]=mvMomVec(softPdfType,Nl,Limi,Prob,Mean,C,nMom,As,Bs,P,maxpts,aEps,rEps);
disp(sprintf('mvMomVec   P=E[1]        %15.10f %15.10g %3d',Val(1),Err(1),Info(1)));
disp(sprintf('mvMomVec   E[X]          %15.10f %15.10g %3d',Val(2),Err(2),Info(2)));
for i=3:4,
   disp(sprintf('mvMomVec   E[X^(%d)] %15.10f %15.10g %3d',P(i),Val(i),Err(i),Info(i)));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the mvProAG2 and mvMomVecAG2 functions
%

disp(' ');
disp('----------------------------------------------');
disp('Testing the mvProAG2 and mvMomVecAG2 functions');
disp(' ');

clear;
N=4;
rand('state',1);
C=rand(N,N);C=C'*C;

  % Seems to be a bug in MATLAB 6 !!
  if sum(sum(isnan(C)))>0           % If there is a NaN in C
    rand('state',1);                %   try to get rid of it by   
    C=rand(N,N);C=C'*C;             %   regenerating it
    if sum(sum(isnan(C)))>0         %   and if it is still there
      disp('Covariance matrix C='); %     give an error message and stop
      disp(C);
      error('MATLAB bug, there is a NaN in C');
    end
  end;
  
disp('Covariance matrix C=');
disp(C);

%
% Generating the lower and upper bound vectors a and b
%
a=zeros(N,1);
softInterv=1;
b=softInterv*ones(N,1);
disp('Lower bound of integration domain a=');
disp(a');
disp('Upper bound of integration domain b=');
disp(b');

%
% Setting some integration parameters
%
maxpts=1000000;
aEps=0;rEps=0.001;

%
% testing mvPro
%
softPdfType=1;Nl=2*ones(N,1);Limi=[a b];Prob=(1/softInterv)*ones(N,1);Mean=zeros(N,1);
t0=clock;
[Val Err Info]=mvPro(softPdfType,Nl,Limi,Prob,Mean,C,maxpts,aEps,rEps);
ts=etime(clock,t0);
disp(sprintf('method     name           value          error      Info CPUtime(sec)'));
disp(sprintf('mvPro    P=E[1]  %15.10f %15.10f %3d %10.5f',Val,Err,Info,ts));

%
% testing mvMomVec
%

nMom=1;
As=zeros(N,1);bs=1;p=1;
t0=clock;
[Val Err Info]=mvMomVec(softPdfType,Nl,Limi,Prob,Mean,C,nMom,As,bs,p,maxpts,aEps,rEps);
ts=etime(clock,t0);
disp(sprintf('mvMomVec P=E[1]  %15.10f %15.10f %3d %10.5f',Val,Err,Info,ts));

nMom=N;
As=eye(N);
bs=zeros(1,N);
p=ones(1,N);
t0=clock;
[vV vE Info]=mvMomVec(softPdfType,Nl,Limi,Prob,Mean,C,nMom,As,bs,p,maxpts,aEps,rEps);
ts=etime(clock,t0);
for i=1:N
  disp(sprintf('mvMomVec E[X(%d)] %15.10f %15.10f %3d %10.5f',i,vV(i),vE(i),Info,ts/N));
end;

disp(' ');
disp('test complete');

