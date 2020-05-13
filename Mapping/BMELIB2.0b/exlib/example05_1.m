% example05_1               - Example VI.3.9 : Part 1 
%
% There are 50 categorical-type soft data, randomly located over
% a square of unit size. The categories are ordered and are built
% from classes ranging from -3 to 3 by unit steps.
% The objective is to reconstruct at best the values at the nodes
% of the grid from the intervals. Two approaches are compared : 
%
% 1- taking the middle of the interval at the estimation points
% 2- BME with interval data at the estimation points
%
% This example is computationally very heavy due to the large
% number of soft data used in each local neighbourhood.

clear;
rand('seed',1);
randn('seed',1);
disp('Creating data');
covmodel='sphericalC';
param=[1 0.5];
q=(-3:1:3)';
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
nk=size(ck,1);
nhmax=5;
dmax=1;
order=NaN;

[zk]=simuchol(ck,covmodel,param,1);
a=zeros(nk,1);
b=zeros(nk,1);
for i=1:nk,
  index=sum(zk(i)>q);
  a(i)=q(index);
  b(i)=q(index+1);
end;

disp('Estimation with middle of intervals');
zkMI=(a+b)/2; 

disp('BME with interval data');
[zkBME]=BMEintervalMode(ck,[],ck,[],a,b,covmodel,param,nhmax,nhmax,dmax,order);

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,14]);
set(gcf,'Position',[0,0,17,14]);

subplot(2,2,1);
pcolor((0:1/19:1),(0:1/19:1),reshape(zk,20,20));
caxis([min(zk),max(zk)]);
axis('square');
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
title('(a)');

subplot(2,2,2);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkMI,20,20));
caxis([min(zk),max(zk)]);
axis('square');
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
title('(b)');

subplot(2,2,3);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkBME,20,20));
caxis([min(zk),max(zk)]);
axis('square');
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
title('(c)');

subplot(2,2,4);
[pdfMI]=kerneldensity(zk-zkMI,(-2:0.1:2)',0.01);
[pdfBME]=kerneldensity(zk-zkBME,(-2:0.1:2)',0.01);
plot((-2:0.1:2),pdfMI,'--');hold on;
plot((-2:0.1:2),pdfBME,'-');hold off;
axis('square');
set(gca,'XTick',(-2:1:2));
set(gca,'FontSize',8);
xlabel('estimation error');
ylabel('pdf');
title('(d)');

print -deps example05_1

save example05_1

