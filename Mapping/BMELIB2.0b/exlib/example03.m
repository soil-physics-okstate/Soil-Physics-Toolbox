% example03                 - Example VI.3.7 : Changing location of soft interval
%
% This example illustrates the use of BMEintervalMode
% using a single covariance model in the univariate
% case. The data used are simulated using a Choleski
% decomposition, so they are Gaussian-distributed.
% The estimation is performed for a single point 
% surrounded by five hard data and a single soft
% data that is maker closer and closer to the 
% estimation point
% Blue circles are hard data, red crosses are soft
% data and the estimation point is on (0.5,0.5).

%%%%%% Clear memory and initialize random generators

clear
rand('seed',6);
randn('seed',6);

%%%%%% Create the matrix of coordinates for
%%%%%% the estimation point, the 5 hard data
%%%%%% and the 3 soft data

nh=5;
ck=[0.5 0.5];
ch=rand(nh,2);
cs=[[0.46;0.42;0.38],[0.54;0.58;0.62]];

covmodel='exponentialC';
covparam=[1 1];
nhmax=5;
nsmax=1;
dmax=Inf;
order=NaN;

%%%%%% Simulate a realization of the values
%%%%%% for an exponential model

for i=1:400,
  disp([num2str(i),'/400']);
  [zall]=simuchol([ck;ch;cs],covmodel,covparam,1);
  zk(i)=zall(1);
  zh=zall(2:nh+1);
  zs=zall(nh+2:nh+4);
  a=zs-rand(3,1)/2;
  b=a+0.5;

  [zkkrig]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
  Zkkrig(i)=zkkrig;
  for j=1:3,
    csj=cs(j,:);
    aj=a(j,:);
    bj=b(j,:);
    [zkBME]=BMEintervalMode(ck,ch,csj,zh,aj,bj,covmodel,covparam,nhmax,nsmax,dmax,order);
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

subplot(1,2,2);
[pdf]=kerneldensity(zk-Zkkrig,(-3.5:0.1:3.5)',0.05);
a=plot((-3.5:0.1:3.5)',pdf,'-b');hold on;
set(a,'LineWidth',1.5);
linetypes={':','--','-'};
for i=1:3,
  [pdf]=kerneldensity(zk'-ZkBME(:,i),(-3.5:0.1:3.5)',0.05);
  plot((-3.5:0.1:3.5)',pdf,linetypes{i});
  axis([-3 3 0 0.8]);
end;
hold off;
axis('square');
axis([-2 2 0 0.8]);
set(gca,'FontSize',8);
set(gca,'XTick',-3:1:3);
set(gca,'YTick',0:0.1:0.8);
xlabel('estimation errors','FontSize',10);
ylabel('pdf','FontSize',10);
title('(b)');

print -deps example03

biaskriging=mean(zk-Zkkrig)
MSEkriging=(zk-Zkkrig)*(zk-Zkkrig)'/400
for i=1:3,
  biasBME(i)=mean(zk'-ZkBME(:,i));
  MSEBME(i)=(zk'-ZkBME(:,i))'*(zk'-ZkBME(:,i))/400;
end;
biasBME
MSEBME

save example03





