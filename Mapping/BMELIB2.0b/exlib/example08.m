% example08                 - Example VI.5.2 : Using 100 soft interval data
%
% There are only 100 interval-type soft data, randomly
% located over a square oF unit size. The limits of the intervals
% are taken as the 0.001, 0.25, 0.5 and 0.999 quantiles of the 
% zero-mean unit variance Gaussian distribution.
% The objective is to predict at best the values at the nodes of
% the grid. Three approaches are proposed : 
%
% 1- kriging with 100 "hardened" soft data taken as the middle
%    of the intervals
% 2- BME with soft data
%
% The first approach is deficient because it uses the middle of
% the intervals and yields a map that gives too much weight to the 
% middle of the interval and is unable to account for the uncertainty
% that varies from one class to another. Only BME is doing a correct
% job in this case.
  
clear;
rand('seed',1);
randn('seed',1);
disp('Creating data');
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=[];
cs=rand(100,2);
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([cs;ck],covmodel,covparam,ns);
classe=gaussinv([0.001 0.25 0.5 0.999]',[0 1]);
zh=[];
zk=Z(101:500);
a=zeros(100,1);
b=zeros(100,1);
for i=1:100,
 a(i)=sum(Z(i)>classe);
 b(i)=a(i)+1;
end;
a=classe(a);
b=classe(b);
nhmax=10;
nsmax=3;
dmax=1;
order=NaN;

disp('Kriging with hard+"hardened" soft data');
[zkKRIGHS,vk]=kriging(ck,[ch;cs],[zh;(a+b)/2],covmodel,covparam,nhmax,dmax,order);
disp('BME with hard+soft data');
[zkBME,info]=BMEintervalMode(ck,ch,cs,zh,a,b,covmodel,covparam,nhmax,nsmax,dmax,order);

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[2 2 14 8]);
set(gcf,'Position',[2 2 14 8]);
z=(-3.5:0.1:3.5)';
plot(z,gausspdf(z,[0 1]));hold on;
for i=1:max(size(classe)),
  plot([classe(i), classe(i)],[0 gausspdf(classe(i),[0 1])])
end;
ax=axis;
text(mean(ax(1:2)),ax(3)-0.13*diff(ax(3:4)),'X',...
  'HorizontalAlignment','center','FontSize',10);
ylabel('pdf','FontSize',10);
axis([-4 4 0 0.45]);
set(gca,'XTick',(-4:1:4));
set(gca,'YTick',(0:0.05:0.45));
set(gca,'FontSize',10);
set(gca,'Position',[0.13 0.13 0.775 0.815]);

print -deps example08_1

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
Sizelim=[4 4+eps];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
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
pcolor((0:1/19:1),(0:1/19:1),reshape(Z(101:500),20,20));
caxis([min(zk)+0.5 max(zk)-0.5]);
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
pcolor((0:1/19:1),(0:1/19:1),reshape(zkKRIGHS,20,20));
caxis([min(zk)+0.5 max(zk)-0.5]);
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
pcolor((0:1/19:1),(0:1/19:1),reshape(zkBME,20,20));
caxis([min(zk)+0.5 max(zk)-0.5]);
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

subplot(3,2,6);
[pdfzkKRIGHS]=kerneldensity(Z(101:500)-zkKRIGHS,(-3:0.1:3)',0.1);
[pdfzkBME]=kerneldensity(Z(101:500)-zkBME,(-3:0.1:3)',0.1);
plot((-3:0.1:3)',pdfzkKRIGHS,'--');hold on;
plot((-3:0.1:3)',pdfzkBME,'-');
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(-3:1:3));
set(gca,'FontSize',8);
xlabel('estimation error');
ylabel('pdf');
title('(e)');
hold off;

print -deps example08_2

biasKriging=mean(Z(101:500)-zkKRIGHS)
biasBME=mean(Z(101:500)-zkBME)
MSEkriging=(Z(101:500)-zkKRIGHS)'*(Z(101:500)-zkKRIGHS)/400
MSEBME=(Z(101:500)-zkBME)'*(Z(101:500)-zkBME)/400

save example08
