% example04                 - Example VI.3.8 : Varying size of soft interval
%
% This example illustrates the use of BMEintervalMode
% using a single covariance model in the univariate
% case. The data used are simulated using a Choleski
% decomposition, so they are Gaussian-distributed.
% The estimation is performed for a single point 
% surrounded by 5 hard data and 1 soft data, with
% varying interval lengths
% Full disks are hard data, open disk is soft data
% and cross is the estimation point on (0.5,0.5).

%%%%%% Clear memory and initialize random generators

clear
rand('seed',6);
randn('seed',6);

%%%%%% Create the matrix of coordinates for
%%%%%% the estimation point, the 5 hard data
%%%%%% and the 5+1 soft data

nh=5;
ck=[0.5 0.5];
ch=rand(nh,2);
cs=[0.45 0.55];

covmodel='exponentialC';
covparam=[1 1];
nhmax=5;
nsmax=1;
dmax=Inf;
order=NaN;

%%%%%% Simulate a realization of the values
%%%%%% for an exponential model

nSim=1000;
intlength=(0.1:2.1:4.3)';
for i=1:nSim,
  disp([num2str(i),'/' num2str(nSim)]);
  [zall]=simuchol([ck;ch;cs],covmodel,covparam,1);
  zk(i)=zall(1);
  zh=zall(2:nh+1);
  zs=zall(nh+2);

  [zkkrig1]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
  [zkkrig2]=kriging(ck,[ch;cs],[zh;zs],covmodel,covparam,nhmax,dmax,order);
  Zkkrig1(i)=zkkrig1;
  Zkkrig2(i)=zkkrig2;
  for j=1:length(intlength),
    a=zs-rand*intlength(j);
    b=a+intlength(j);
    [zkBME]=BMEintervalMode(ck,ch,cs,zh,a,b,covmodel,covparam,nhmax,nsmax,dmax,order);
    ZkBME(i,j)=zkBME;
  end;
end;

%%%%%% Display the results

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 18.5 7.5]);
set(gcf,'Position',[0.5 0.5 18.5 7.5]);
subplot(1,2,1);
plot(ck(1),ck(2),'xb');hold on;
plot(ch(:,1),ch(:,2),'.b','MarkerSize',20);
plot(cs(:,1),cs(:,2),'ob','MarkerSize',6);
hold off;
axis('square');
axis([0 1 0 1]);
set(gca,'FontSize',8);
set(gca,'XTick',0:0.2:1);
set(gca,'YTick',0:0.2:1);
ax=axis;
text(mean(ax(1:2)),ax(1)-0.11*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',10);
ylabel('s_2','FontSize',10);
title('(a)');
hold off;

subplot(1,2,2);
linetypes={':','--','-'};
for i=1:length(intlength),
  [pdf]=kerneldensity(zk'-ZkBME(:,i),(-3.5:0.1:3.5)',0.05);hold on;
  plot((-3.5:0.1:3.5)',pdf,linetypes{i});
end;
axis('square');
axis([-2 2 0 0.7]);
set(gca,'FontSize',8);
set(gca,'XTick',-3:1:3);
set(gca,'YTick',0:0.1:0.8);
xlabel('estimation errors','FontSize',10);
ylabel('pdf','FontSize',10);
title('(b)');
box on;
hold off;

print -deps example04

biaskriging1=mean(zk-Zkkrig1)
MSEkriging1=(zk-Zkkrig1)*(zk-Zkkrig1)'/nSim
biaskriging2=mean(zk-Zkkrig2)
MSEkriging2=(zk-Zkkrig1)*(zk-Zkkrig2)'/nSim
for i=1:length(intlength),
  biasBME(i)=mean(zk'-ZkBME(:,i));
  MSEBME(i)=(zk'-ZkBME(:,i))'*(zk'-ZkBME(:,i))/nSim;
end;
biasBME
MSEBME

save example04





