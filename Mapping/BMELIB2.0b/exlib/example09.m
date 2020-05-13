% example09                 - Example VI.5.3 : Using 50 hard data and 400 grid node intervals
%
% There are 50 hard data randomly located and 400 interval-type
% soft data at  the nodes of a 20 x 20 grid. The limits of the
% intervals are taken as the 0.001, 0.25, 0.5 and 0.999 quantiles
% of the zero-mean unit variance Gaussian distribution.
% The objective is to predict at best the values at the nodes of
% the grid. Three approaches are proposed : 
%
% 1- kriging with only 50 hard data
% 2- kriging with hard 50 hard data and "hardened" soft data
%    taken as the middle of the intervals
% 3- BME with hard and soft data
%
% The first approach is deficient because it does not make use of
% the 400 soft data. For the second approach, because kriging is
% an exact predictor, using the middle of the intervals does not
% allow to take into account the hard data at all and just restitute
% the middle of the class at the prediction point. Only BME is doing
% a correct job in this case.
  
clear;
rand('seed',2);
randn('seed',2);
disp('Creating data');
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=rand(50,2);
cs=ck;
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([ch;cs],covmodel,covparam,ns);
classe=gaussinv([0.001 0.25 0.5 0.999]',[0 1]);
zh=Z(1:50);
a=zeros(400,1);
b=zeros(400,1);
for i=1:400,
 a(i)=sum(Z(50+i)>classe);
 b(i)=a(i)+1;
end;
a=classe(a);
b=classe(b);
nhmax=10;
nsmax=3;
dmax=1;
order=NaN;

disp('Kriging with hard data...');
[zkKRIGH,vk]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
disp('Kriging with hard+"hardened" soft data');
[zkKRIGHS,vk]=kriging(ck,[ch;cs],[zh;(a+b)/2],covmodel,covparam,nhmax,dmax,order);
disp('BME with hard+soft data');
[zkBME,info]=BMEintervalMode(ck,ch,cs,zh,a,b,covmodel,covparam,nhmax,nsmax,dmax,order);

figure;
z=(-3.5:0.1:3.5)';
plot(z,gausspdf(z,[0 1]));hold on;
for i=1:max(size(classe)),
  plot([classe(i), classe(i)],[0 gausspdf(classe(i),[0 1])])
end;
xlabel('values','FontSize',12);
ylabel('density','FontSize',12);
axis([-4 4 0 0.45]);
set(gca,'XTick',(-4:1:4));
set(gca,'YTick',(0:0.05:0.45));
set(gca,'FontSize',10);

print -deps example09_1

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
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
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(a)');
hold off;

subplot(3,2,2);
pcolor((0:1/19:1),(0:1/19:1),reshape(Z(51:450),20,20));
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
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
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
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(c)');
hold off;

subplot(3,2,4);
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGHS,20,20));
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
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
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
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
title('(e)');
hold off;

subplot(3,2,6);
[pdfzkKRIGH]=kerneldensity(Z(51:450)-zkKRIGH,(-3:0.1:3)',0.1);
[pdfzkKRIGHS]=kerneldensity(Z(51:450)-zkKRIGHS,(-3:0.1:3)',0.1);
[pdfzkBME]=kerneldensity(Z(51:450)-zkBME,(-3:0.1:3)',0.1);
plot((-3:0.1:3)',pdfzkKRIGH,':');hold on
plot((-3:0.1:3)',pdfzkKRIGHS,'--');
plot((-3:0.1:3)',pdfzkBME,'-');
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(-3:1:3));
set(gca,'FontSize',8);
xlabel('estimation error');
ylabel('pdf');
title('(f)');
hold off;

print -deps example09_2

biasKrigingH=mean(Z(51:450)-zkKRIGH)
biasKrigingHS=mean(Z(51:450)-zkKRIGHS)
biasBME=mean(Z(51:450)-zkBME)
MSEkrigingH=(Z(51:450)-zkKRIGH)'*(Z(51:450)-zkKRIGH)/400
MSEkrigingHS=(Z(51:450)-zkKRIGHS)'*(Z(51:450)-zkKRIGHS)/400
MSEBME=(Z(51:450)-zkBME)'*(Z(51:450)-zkBME)/400

save example09


