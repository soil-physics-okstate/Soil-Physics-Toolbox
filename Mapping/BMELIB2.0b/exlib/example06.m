% example06                 - Example VI.4.1 : Mapping efficiency of soft interval data
%
% There are 25 hard data and 100 interval-type soft data, randomly
% located over a square of unit size. The width of the intervals is
% equal to 1.5, such that the simulated values are uniformly
% distributed over the intervals.
% The objective is to predict at best the values at the nodes of
% the grid. Three approaches are proposed : 
%
% 1- kriging with only 25 hard data
% 2- kriging with 50 hard data
% 3- BME with 25 hard and 25 interval data
%
% A map of the reduction of variance can be computed between
% on one side the two kriging situations and on the other side
% the first kriging situation and BME
  
clear;
rand('seed',1);
randn('seed',1);
disp('Creating data');
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=rand(25,2);
cs=rand(25,2);
z=(-4:0.1:4)';
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([ch;cs;ck],covmodel,covparam,ns);
zh=Z(1:25);
zs=Z(26:50);
a=zs-rand(25,1)*1.5;
b=a+1.5;
nhmax=10;
nsmax=3;
dmax=1;
order=NaN;

disp('Kriging with 25 hard data...');
[zkKRIGHA,vkA]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
disp('Kriging with 50 hard data...');
[zkKRIGHB,vkB]=kriging(ck,[ch;cs],[zh;zs],covmodel,covparam,nhmax,dmax,order);
disp('BME with 25 hard data + 25 interval data');
zkBME=zeros(400,1);
vkBME=zeros(400,1);
for i=1:400,
  [z,pdf,info]=BMEintervalPdf(z,ck(i,:),ch,cs,zh,a,b,covmodel,covparam,nhmax,nsmax,dmax,order);
  [m,v]=pdfstat(z,pdf);
  zkBME(i)=m;
  vkBME(i)=v;
end;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
Sizelim=[4 4+eps];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 0],[0 0 0]};
markerplot(ch,zh,Sizelim,Property,Value);
hold on;
Sizelim=[4 4+eps];
Value={'o',[0 0 0],[1 1 1]};
markerplot(cs,(a+b)/2,Sizelim,Property,Value);
axis([0 1 0 1]);
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(a)');
hold off;

subplot(3,2,2);
pcolor((0:1/19:1),(0:1/19:1),reshape(vkA,20,20));
caxis([0 0.8]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(b)');
hold off;

subplot(3,2,3);
pcolor((0:1/19:1),(0:1/19:1),reshape(vkB,20,20));
caxis([0 0.8]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(c)');
hold off;

subplot(3,2,4);
pcolor((0:1/19:1),(0:1/19:1),reshape(vkBME,20,20));
caxis([0 0.8]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
ax=axis;
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(d)');
hold off;

subplot(3,2,5);
index=(vkA-vkBME)./(vkA-vkB);
pcolor((0:1/19:1),(0:1/19:1),reshape(index,20,20));
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
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(e)');
hold off;

subplot(3,2,6);
[nf,nb]=hist(index,0:0.05:1);
[obj]=bar(nb,nf);
set(obj,'FaceColor',[1 1 1]);
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(0:0.1:1));
set(gca,'FontSize',8);
axis([0 1 0 80]);
xlabel('efficiency index');
ylabel('pdf');
title('(f)');

print -deps example06
save example06



