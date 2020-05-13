% example10                 - Example VI.5.4 : Using 50 hard data, 100 intervals and 400 grid node intervals
%
% This show how BME deals with soft data when they become more
% and more informative. BME starts with 25 hard data and 25
% soft data that are weekly informative. The maps is compare to
% the kriging map to assert this. The width of the intervals is
% progressively reduced, so that and the end the 25 soft data
% are close to hard data. This is also shown by comparison with
% a kriging map with 50 hard data.
% The result are also compared to kriging with 25 hard data and
% 25 hardened soft data, in order to show the deficiency of the
% approach
  
clear;
rand('seed',10);
randn('seed',10);
disp('Creating data');
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=rand(25,2);
cs=rand(100,2);
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([ch;cs;ck],covmodel,covparam,ns);

classe1=gaussinv([0.001 0.5 0.999]',[0 1]);
classe2=gaussinv([0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.999]',[0 1]);
zh=Z(1:25);
zs=Z(26:125);
a=zeros(100,1);
b=zeros(100,1);
for i=1:100,
 a1(i)=sum(Z(25+i)>classe1);
 b1(i)=a1(i)+1;
 a2(i)=sum(Z(25+i)>classe2);
 b2(i)=a2(i)+1;
end;
a1=classe1(a1);
b1=classe1(b1);
a2=classe2(a2);
b2=classe2(b2);
nhmax=10;
nsmax=3;
dmax=1;
order=NaN;

disp('Kriging with hard data/25 hard data...');
[zkKRIGH1,vk]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
disp('Kriging with hard data/125 hard data...');
[zkKRIGH2,vk]=kriging(ck,[ch;cs],[zh;zs],covmodel,covparam,nhmax,dmax,order);
disp('Kriging with hard+"hardened" soft data/classe1');
[zkKRIGHS1,vk]=kriging(ck,[ch;cs],[zh;(a1+b1)/2],covmodel,covparam,nhmax,dmax,order);
disp('Kriging with hard+"hardened" soft data/classe2');
[zkKRIGHS2,vk]=kriging(ck,[ch;cs],[zh;(a2+b2)/2],covmodel,covparam,nhmax,dmax,order);
disp('BME with hard+soft data/classe1');
[zkBME1,info]=BMEintervalMode(ck,ch,cs,zh,a1,b1,covmodel,covparam,nhmax,nsmax,dmax,order);
disp('BME with hard+soft data/classe2');
[zkBME2,info]=BMEintervalMode(ck,ch,cs,zh,a2,b2,covmodel,covparam,nhmax,nsmax,dmax,order);

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[0.5 0.5 18.5 9]);
set(gcf,'PaperPosition',[0.5 0.5 18.5 9]);

subplot(1,2,1);
z=(-3.5:0.1:3.5)';
plot(z,gausspdf(z,[0 1]));hold on;
for i=1:max(size(classe2)),
  plot([classe2(i), classe2(i)],[0 gausspdf(classe2(i),[0 1])])
end;
xlabel('X','FontSize',12);
ylabel('pdf','FontSize',12);
axis([-4 4 0 0.45]);
set(gca,'XTick',(-4:1:4));
set(gca,'YTick',(0:0.05:0.45));
set(gca,'FontSize',10);
title('(a)');

subplot(1,2,2);
z=(-3.5:0.1:3.5)';
plot(z,gausspdf(z,[0 1]));hold on;
for i=1:max(size(classe1)),
  plot([classe1(i), classe1(i)],[0 gausspdf(classe1(i),[0 1])])
end;
xlabel('X','FontSize',12);
ylabel('pdf','FontSize',12);
axis([-4 4 0 0.45]);
set(gca,'XTick',(-4:1:4));
set(gca,'YTick',(0:0.05:0.45));
set(gca,'FontSize',10);
title('(b)');

print -deps example10_1

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'Position',[3.5 0.4 13 21.5]);
set(gcf,'PaperPosition',[3.5 0.4 13 21.5]);

subplot(4,2,1);
Sizelim=[4 4+eps];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 0],[0 0 0]};
markerplot(ch,zh,Sizelim,Property,Value);
hold on;
Sizelim=[4 4+eps];
Value={'o',[0 0 0],[1 1 1]};
markerplot(cs,(a1+b1)/2,Sizelim,Property,Value);
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
pcolor((0:1/19:1),(0:1/19:1),reshape(Z(126:525),20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(b)');
hold off;

subplot(4,2,3);
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
title('(c)');
hold off;

subplot(4,2,4);
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
title('(d)');
hold off;

subplot(4,2,5);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGHS2,20,20));
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
title('(e)');
hold off;

subplot(4,2,6);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGHS1,20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(f)');
hold off;

subplot(4,2,7);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkBME2,20,20));
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
text(mean(ax(1:2)),ax(3)-0.18*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(g)');
hold off;

subplot(4,2,8);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkBME1,20,20));
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
text(mean(ax(1:2)),ax(3)-0.18*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
title('(h)');
hold off;

print -deps example10_2

biasKrigingH1=mean(Z(126:525)-zkKRIGH1);
biasKrigingHS1=mean(Z(126:525)-zkKRIGHS1);
biasBME1=mean(Z(126:525)-zkBME1);
MSEkrigingH1=(Z(126:525)-zkKRIGH1)'*(Z(126:525)-zkKRIGH1)/400;
MSEkrigingHS1=(Z(126:525)-zkKRIGHS1)'*(Z(126:525)-zkKRIGHS1)/400;
MSEBME1=(Z(126:525)-zkBME1)'*(Z(126:525)-zkBME1)/400;

biasKrigingH2=mean(Z(126:525)-zkKRIGH2);
biasKrigingHS2=mean(Z(126:525)-zkKRIGHS2);
biasBME2=mean(Z(126:525)-zkBME2);
MSEkrigingH2=(Z(126:525)-zkKRIGH2)'*(Z(126:525)-zkKRIGH2)/400;
MSEkrigingHS2=(Z(126:525)-zkKRIGHS2)'*(Z(126:525)-zkKRIGHS2)/400;
MSEBME2=(Z(126:525)-zkBME2)'*(Z(126:525)-zkBME2)/400;

[biasKrigingH1 biasKrigingH2;biasKrigingHS1 biasKrigingHS2;biasBME1 biasBME2]
[MSEkrigingH1 MSEkrigingH2;MSEkrigingHS1 MSEkrigingHS2;MSEBME1 MSEBME2]

save example10



