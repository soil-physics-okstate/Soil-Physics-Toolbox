% example11_1               - Example VI.5.5 : Part 1
%
% There are 50 interval-type soft data, randomly located over
% a square of unit size. The limits of the intervals
% are taken as the 0.001, 0.1, 0.2,...,0.9 and 0.999 quantiles
% of the zero-mean unit variance Gaussian distribution.
% The objective is to predict at best probabilities at the nodes
% of the grid. Three approaches are compared : 
%
% 1- kriging with exact 50 simulated hard data
% 2- BME with 50 interval data
% 3- indicator kriging with interval data
  
clear;
rand('seed',1);
randn('seed',1);
disp('Creating data');
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
nh=50;
ch=rand(nh,2);
cs=ch;
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([ch;ck],covmodel,covparam,ns);
prob=[0.001;(0.1:0.1:0.9)';0.999];
classe=gaussinv(prob,[0 1]);
zh=Z(1:nh);
zk=Z(nh+1:nh+400);
a=zeros(nh,1);
b=zeros(nh,1);
for i=1:nh,
  a(i)=sum(zh(i)>classe);
  b(i)=a(i)+1;
end;
a=classe(a);
b=classe(b);
z=(-4:0.05:4)';
nhmax=3;
nsmax=3;
dmax=1;
order=NaN;

disp('Kriging with hard data...');
[mKrig,vKrig]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
pdfKrig=zeros(length(z),400);
for i=1:400,
  pdfKrig(:,i)=gausspdf(z,[mKrig(i),vKrig(i)]);
end;

disp('BME with interval data...');
pdfBME=zeros(length(z),400);
mBME=zeros(400,1);
vBME=zeros(400,1);
for i=1:400,
  [z,pdf,info]=BMEintervalPdf(z,ck(i,:),[],cs,[],a,b,covmodel,covparam,nhmax,nsmax,dmax,order);
  [m,v]=pdfstat(z,pdf);
  pdfBME(:,i)=pdf;
  mBME(i)=m;
  vBME(i)=v;
end;

disp('Indicator kriging...');
cdfInd=zeros(400,length(classe));
for i=1:length(classe),
  [pInd,vInd]=krigingIndGauss(ck,ch,zh,covmodel,covparam,prob(i),nhmax,dmax);
  cdfInd(:,i)=pInd;
end;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 18.5 7.5]);
set(gcf,'Position',[0.5 0.5 18.5 7.5]);
subplot(1,2,1);
plot(z,gausspdf(z,[0 1]));hold on;
for i=1:max(size(classe)),
  plot([classe(i), classe(i)],[0 gausspdf(classe(i),[0 1])])
end;
xlabel('X','FontSize',10);
ylabel('pdf','FontSize',10);
axis('square');
axis([-4 4 0 0.45]);
set(gca,'XTick',(-4:1:4));
set(gca,'YTick',(0:0.05:0.45));
set(gca,'FontSize',8);
title('(a)');
hold off;

subplot(1,2,2);
for i=2:length(classe)-1,
  Cind=model2indic((0:0.01:1)',covmodel,covparam,prob(i));
  plot((0:0.01:1)',Cind,'b-');hold on
end;
xlabel('distance','FontSize',10);
ylabel('indicator covariance','FontSize',10);
axis('square');
axis([0 1 0 0.25]);
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.05:0.25));
set(gca,'FontSize',8);
title('(b)');
hold off;

print -deps example11_1_1

lpos=4;
upos=8;
lb=classe(lpos);
ub=classe(upos);
lp=prob(lpos);
up=prob(upos);

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
Sizelim=[4 4+eps];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 0],[1 1 1]};
markerplot(ch,zh,Sizelim,Property,Value);
axis([0 1 0 1]);
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.14*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(a)');
hold off;

subplot(3,2,2);
cond=(zk>lb)&(zk<ub);
pcolor((0:1/19:1),(0:1/19:1),reshape(cond,20,20));
caxis([0 1]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.14*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(b)');
hold off;

subplot(3,2,3);
pKrig=zeros(400,1);
for i=1:400,
  pKrig(i)=gausscdf(ub,[mKrig(i),vKrig(i)])-gausscdf(lb,[mKrig(i),vKrig(i)]);
end;
pcolor((0:1/19:1),(0:1/19:1),reshape(pKrig,20,20));
caxis([0 1]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.14*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(c)');
hold off;

subplot(3,2,4);
pBME=zeros(400,1);
for i=1:400,
  pBME(i)=trapezint(z,pdfBME(:,i),lb,ub);
end;
pcolor((0:1/19:1),(0:1/19:1),reshape(pBME,20,20));
caxis([0 1]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.14*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(d)');
hold off;

subplot(3,2,5);
pInd=cdfInd(:,upos)-cdfInd(:,lpos);
pcolor((0:1/19:1),(0:1/19:1),reshape(pInd,20,20));
caxis([0 1]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.14*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(e)');
hold off;

print -deps example11_1_2

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[1 1 15 15]);
set(gcf,'Position',[1 1 15 15]);

subplot(2,2,1);
map=[1 1 1;0.8 0.8 0.8];
colormap(map);
contourf((0:1/19:1),(0:1/19:1),reshape(pBME,20,20),[min(pBME) 0.4]);
axis('square');
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
ax=axis;
hold on;
plot([ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'k');
title('(a)');
hold off;

subplot(2,2,2);
map=[1 1 1;0.8 0.8 0.8];
colormap(map);
contourf((0:1/19:1),(0:1/19:1),reshape(pBME,20,20),[min(pBME) 0.6]);
axis('square');
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
ax=axis;
hold on;
plot([ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'k');
title('(b)');
hold off;

subplot(2,2,3);
map=[1 1 1;0.8 0.8 0.8];
colormap(map);
contourf((0:1/19:1),(0:1/19:1),reshape(pInd,20,20),[min(pInd) 0.4]);
axis('square');
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
ax=axis;
hold on;
plot([ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'k');
title('(c)');
hold off;

subplot(2,2,4);
map=[1 1 1;0.8 0.8 0.8];
colormap(map);
contourf((0:1/19:1),(0:1/19:1),reshape(pInd,20,20),[min(pInd) 0.6]);
axis('square');
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
xlabel('s_1');
ylabel('s_2');
ax=axis;
hold on;
plot([ax(1) ax(2) ax(2) ax(1) ax(1)],[ax(3) ax(3) ax(4) ax(4) ax(3)],'k');
title('(d)');
hold off;

print -deps example11_1_3

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 18.5 18.5]);
set(gcf,'Position',[0.5 0.5 18.5 18.5]);
threshold=(0.01:0.01:0.99);
%subplot(2,2,1);
%for i=1:length(threshold),
%  ClassifyKrig(i)=sum((pKrig>threshold(i))&cond) + sum((pKrig<threshold(i))&~cond);
%  ClassifyBME(i)=sum((pBME>threshold(i))&cond) + sum((pBME<threshold(i))&~cond);
%  ClassifyInd(i)=sum((pInd>threshold(i))&cond) + sum((pInd<threshold(i))&~cond);
%end;
%plot(threshold,ClassifyKrig/400,'-b');hold on;
%plot(threshold,ClassifyBME/400,'-.b');
%plot(threshold,ClassifyInd/400,'--b');
%axis('square');
%axis([0 1 0.45 0.80])
%set(gca,'XTick',(0:0.2:1));
%set(gca,'YTick',(0.45:0.05:0.80));
%set(gca,'FontSize',10);
%xlabel('probability threshold','FontSize',10);
%ylabel('Correct classification rate','FontSize',10);
%title('(a)');
%hold off;

subplot(2,1,1);
for i=1:length(threshold),
  CorrectKrig(i)=sum((pKrig>threshold(i))&cond)/sum(pKrig>threshold(i));
  CorrectBME(i)=sum((pBME>threshold(i))&cond)/sum(pBME>threshold(i));
  CorrectInd(i)=sum((pInd>threshold(i))&cond)/sum(pInd>threshold(i));
end;
plot(threshold,CorrectKrig,'-b');hold on;
plot(threshold,CorrectBME,'-.b');
plot(threshold,CorrectInd,'--b');
axis('square');
axis([0 1 0.5 1])
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0.5:0.1:1));
set(gca,'FontSize',10);
xlabel('probability threshold','FontSize',10);
ylabel('Correct selection rate','FontSize',10);
title('(a)');
hold off;

map=[0.8 0.8 0.8;1 1 1];
colormap(map);
subplot(2,2,3);
contourf((0:1/19:1),(0:1/19:1),reshape(zk,20,20),[lb ub]);hold on;
condKrig=pBME>0.8;
plot(ck(condKrig,1),ck(condKrig,2),'.','MarkerSize',20);
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',10);
xlabel('s_1');
ylabel('s_2');
title('(b)');
hold off;

subplot(2,2,4);
contourf((0:1/19:1),(0:1/19:1),reshape(zk,20,20),[lb ub]);hold on;
condInd=pInd>0.8;
plot(ck(condInd,1),ck(condInd,2),'.','MarkerSize',20);
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',10);
xlabel('s_1');
ylabel('s_2');
title('(c)');
hold off;

print -deps example11_1_4

save example11_1



