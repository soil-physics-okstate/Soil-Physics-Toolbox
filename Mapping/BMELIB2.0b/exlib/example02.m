% example02                 - Example VI.3.5 :  Using probabilistic data from a secondary variable
%
% There are 50 hard data randomly located and 400 pdf-type
% soft data at the nodes of a 20 x 20 grid. The pdf are
% obtained by conditioning with respect to another a
% secondary variable correlated with the first one (rho=0.9).
% The objective is to predict at best the values at the nodes
% of the grid. Two approaches are proposed : 
%
% 1- kriging with only 50 hard data
% 2- kriging with 50 hard data and the conditional soft pdf at
%    prediction point taken as "hardened" value
% 3- bme with rigorous incorporation of hard data and the
%    soft pdf at the estimation point
%
% The first approach is deficient because it does not make use of
% the 400 soft data. For the second approach, because kriging is
% an exact predictor, using as hard data the middle of the intervals
% does not allow to take into account the hard data at all. Only BME
% is doing a correct job in this case.
  
clear;
rand('seed',2);
randn('seed',2);
disp('Creating data');

z=(0:0.05:5)';
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=rand(25,2);
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([ch;ck],covmodel,covparam,ns);
zh=Z(1:25)+2.5;
zk=Z(26:425)+2.5;
yh=zh;
yk=zk;
wk=randn(400,1);
mzk=yk+(yk/4).*wk;
vzk=(yk.^2)/16;
nhmax=10;
nsmax=3;
dmax=1;
order=0;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[5,5,12,9]);
set(gcf,'Position',[5,5,12,9]);
plot([0 5],[0 5]);hold on;
plot(yk,mzk,'r.');
for i=1:length(z),
  ub(i)=gaussinv(0.95,[z(i),(z(i)^2)/12]);
  lb(i)=gaussinv(0.05,[z(i),(z(i)^2)/12]);
end;
plot(z,ub','--');
plot(z,lb','--');
axis([0 5 0 8]);
set(gca,'FontSize',10);
set(gca,'XTick',0:0.5:5);
set(gca,'YTick',0:1:8);
ax=axis;
text(mean(ax(1:2)),ax(1)-0.12*diff(ax(3:4)),'X_2',...
  'HorizontalAlignment','center','FontSize',10);
text(ax(1)-0.07*diff(ax(1:2)),mean(ax(3:4)),'X_1',...
  'HorizontalAlignment','center','FontSize',10,'Rotation',90);
set(gca,'Position',[0.13 0.13 0.775 0.815])

print -deps example02_1;

disp('Kriging with hard data...');
[zkKRIGH,vk]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
disp('BME with hard+ pdf at estimation point');
softpdftype=2;
options=BMEoptions;
options(8)=1;
for i=1:400,
  limi=(-5:0.1:5)+mzk(i);
  nl=length(limi);
  [probdens]=gausspdf(limi,[mzk(i),vzk(i)]);
  [moments,info]=BMEprobaMoments(ck(i,:),ch,ck(i,:),zh,softpdftype,...
    nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
  zkBME(i)=moments(1);  
end;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
pcolor((0:1/19:1),(0:1/19:1),reshape(zk,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
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
title('(a)');
hold off;

subplot(3,2,2);
Sizelim=[4 4+eps];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 0],[0 0 0]};
markerplot(ch,zh,Sizelim,Property,Value)
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
title('(b)');
hold off;

subplot(3,2,3);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGH,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
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
pcolor((0:1/19:1),(0:1/19:1),reshape(mzk,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
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
pcolor((0:1/19:1),(0:1/19:1),reshape(zkBME,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
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

subplot(3,2,6);
[pdfmzk]=kerneldensity(zk-mzk,(-3:0.1:3)',0.1);
[pdfzkKRIGH]=kerneldensity(zk-zkKRIGH,(-3:0.1:3)',0.1);
[pdfzkBME]=kerneldensity(zk-zkBME',(-3:0.1:3)',0.1);
plot((-3:0.1:3)',pdfmzk,':');hold on
plot((-3:0.1:3)',pdfzkKRIGH,'--');hold on;
plot((-3:0.1:3)',pdfzkBME,'-');
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(-3:1:3));
set(gca,'FontSize',8);
axis([-2 2 0 .8]);
xlabel('estimation error');
ylabel('pdf');
title('(f)');
hold off;

print -deps example02_2

biasExpect=mean(zk-mzk)
biasKriging=mean(zk-zkKRIGH)
biasBME=mean(zk-zkBME')
MSEExpect=(zk-mzk)'*(zk-mzk)/400
MSEkriging=(zk-zkKRIGH)'*(zk-zkKRIGH)/400
MSEBME=(zk-zkBME')'*(zk-zkBME')/400

save example02


