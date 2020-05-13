% example05_2               - Example VI.3.9 : Part 2
%
% There are 50 categorical-type soft data, randomly located over
% a square of unit size. The categories are ordered and are built
% from classes ranging from -3 to 3 by changing steps.
% The objective is to reconstruct at best the values at the nodes
% of the grid from the intervals. Two approaches are compared : 
%
% 1- taking the middle of the interval at the estimation points
% 2- BME with interval data at the estimation points
%
% The MSE are compared for various choices of the step (from 0.25 to 2)
%
% This example is computationally very heavy due to the large
% number of soft data used in each local neighbourhood.

clear;
rand('seed',1);
randn('seed',1);
disp('Creating data');
covmodel='sphericalC';
param=[1 0.5];
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
nk=size(ck,1);
nhmax=5;
dmax=1;
order=NaN;

[zk]=simuchol(ck,covmodel,param,1);

step=[0.25,0.5,0.75,1,1.5,2,3]';
nstep=length(step);
zkMI=zeros(400,nstep);
zkBME=zeros(400,nstep);
MEMI=zeros(nstep,1);
MEBME=zeros(nstep,1);
MSEMI=zeros(nstep,1);
MSEBME=zeros(nstep,1);
for i=1:nstep,
  disp(['Estimation for step ',num2str(step(i))]);
  q=(-3:step(i):3)';
  a=zeros(nk,1);
  b=zeros(nk,1);
  for j=1:nk,
    index=sum(zk(j)>q);
    a(j)=q(index);
    b(j)=q(index+1);
  end;
  [zkBME]=BMEintervalMode(ck,[],ck,[],a,b,covmodel,param,nhmax,nhmax,dmax,order);
  ZKBME(:,i)=zkBME;
  zkMI=(a+b)/2;
  ZKMI(:,i)=zkMI;
  MEMI(i)=mean(zk-zkMI);
  MEBME(i)=mean(zk-zkBME);
  [MEMI MEBME]
  MSEMI(i)=(zk-zkMI)'*(zk-zkMI)/400;
  MSEBME(i)=(zk-zkBME)'*(zk-zkBME)/400;
  [MSEMI MSEBME]
end;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,14,10]);
set(gcf,'Position',[2,2,14,10]);
plot([0;step],[0;sqrt(MSEMI)],'*');hold on;
plot([0;step],[0;sqrt(MSEMI)]);
plot([0;step],[0;sqrt(MSEBME)],'o');
plot([0;step],[0;sqrt(MSEBME)]);
set(gca,'XTick',(0:0.5:3));
set(gca,'FontSize',10);
xlabel('Interval length','FontSize',12);
ylabel('Square root of MSE','FontSize',12);

print -deps example05_2;

save example05_2

