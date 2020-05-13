% example13                 - Example VI.6.2 : Merging Maps
%
% There are 50 hard data randomly located.
% The objective is to predict at best the values at the nodes
% of the grid. Three approaches are proposed : 
%
% 1- kriging with only 25 hard data above 0
% 1- kriging with only 25 hard data below 0
% 2- kriging with the 50 hard data
% 3- combining the two kriging  maps

clear;
rand('seed',1);
randn('seed',1);
disp('Creating data');

z=(-5:0.1:5)';
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=rand(50,2);
c=[ck;ch];
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol(c,covmodel,covparam,ns);
zh=Z(401:450);
index1=find(zh<0);
index2=find(zh>0);
ch1=ch(index1,:);
ch2=ch(index2,:);
zh1=zh(index1);
zh2=zh(index2);
zk=Z(1:400);
nhmax=10;
dmax=1;
order=NaN;

disp('Kriging set 1...');
[zkKRIGH1,vkKRIGH1]=kriging(ck,ch1,zh1,covmodel,covparam,nhmax,dmax,order);
disp('Kriging set 2...');
[zkKRIGH2,vkKRIGH2]=kriging(ck,ch2,zh2,covmodel,covparam,nhmax,dmax,order);
disp('Kriging set 1+2...');
[zkKRIGH12,vkKRIGH12]=kriging(ck,[ch1;ch2],[zh1;zh2],covmodel,covparam,nhmax,dmax,order);
disp('Combining maps...');
for i=1:400,
  pdf1=gausspdf(z,[zkKRIGH1(i) vkKRIGH1(i)]);
  pdf2=gausspdf(z,[zkKRIGH2(i) vkKRIGH2(i)]);
  pdf=pdfprod(z,pdf1,pdf2);
  zkKRIGUPD(i)=pdfstat(z,pdf);
end;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[3.5 0.4 13 21.5]);
set(gcf,'PaperPosition',[3.5 0.4 13 21.5]);

subplot(4,2,1);
Sizelim=[4 4+eps];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 0],[0 0 0]};
markerplot(ch1,zh1,Sizelim,Property,Value);
hold on;
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 0],[1 1 1]};
markerplot(ch2,zh2,Sizelim,Property,Value);
axis([0 1 0 1]);
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ylabel('s_2');
title('(a)');
hold off;

subplot(4,2,2);
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
ylabel('s_2');
title('(b)');
hold off;

subplot(4,2,3);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGH1,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ylabel('s_2');
title('(c)');
hold off;

subplot(4,2,4);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGH2,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ylabel('s_2');
title('(d)');
hold off;

subplot(4,2,5);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGH12,20,20));
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


subplot(4,2,6);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGUPD,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ylabel('s_2');
title('(f)');
hold off;

subplot(4,2,8);
[pdfzkKRIGH1]=kerneldensity(zk-zkKRIGH1,(-3:0.1:3)',0.1);
[pdfzkKRIGH2]=kerneldensity(zk-zkKRIGH2,(-3:0.1:3)',0.1);
[pdfzkKRIGH12]=kerneldensity(zk-zkKRIGH12,(-3:0.1:3)',0.1);
[pdfzkKRIGUPD]=kerneldensity(zk-zkKRIGUPD',(-3:0.1:3)',0.1);
plot((-3:0.1:3)',pdfzkKRIGH1,':');hold on;
plot((-3:0.1:3)',pdfzkKRIGH2,':');
plot((-3:0.1:3)',pdfzkKRIGH12,'--');
plot((-3:0.1:3)',pdfzkKRIGUPD,'-');
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(-3:1:3));
set(gca,'FontSize',8);
xlabel('estimation error');
ylabel('pdf');
title('(g)');
axis([-2 2 0 .8])
hold off;

print -deps example13

biasKrigingH1=mean(zk-zkKRIGH1)
biasKrigingH2=mean(zk-zkKRIGH2)
biasKrigingH12=mean(zk-zkKRIGH12)
biasKRIGUPD=mean(zk-zkKRIGUPD')
MSEkrigingH1=(zk-zkKRIGH1)'*(zk-zkKRIGH1)/400
MSEkrigingH2=(zk-zkKRIGH2)'*(zk-zkKRIGH2)/400
MSEkrigingH12=(zk-zkKRIGH12)'*(zk-zkKRIGH12)/400
MSEkrigingUPD=(zk-zkKRIGUPD')'*(zk-zkKRIGUPD')/400

save example13


