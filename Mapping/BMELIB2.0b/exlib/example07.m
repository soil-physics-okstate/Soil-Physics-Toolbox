% example07                 - Example VI.5.1 :  Using 25 hard data and 100 intervals 
%
% There are 25 hard data and 100 interval-type soft data, randomly
% located over a square oF unit size. The limits of the intervals
% are taken as the 0.001, 0.25, 0.5 and 0.999 quantiles of the 
% zero-mean unit variance Gaussian distribution.
% The objective is to predict at best the values at the nodes of
% the grid. Three approaches are proposed : 
%
% 1- kriging with only 25 hard data
% 2- kriging with hard 25 hard data and 100 "hardened" soft data
%    taken as the middle of the intervals
% 3- BME with hard and soft data
%
% The first approach is deficient because it does not make use of
% the 100 soft data. For the second approach, using the middle of
% the intervals yields a map that gives too much weight to the 
% middle of the interval and is unable to account for the uncertainty
% that varies from one class to another. Only BME is doing a correct
% job in this case.
  
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
classe=gaussinv([0.001 0.25 0.5 0.999]',[0 1]);
zh=Z(1:25);
a=zeros(100,1);
b=zeros(100,1);
for i=1:100,
 a(i)=sum(Z(25+i)>classe);
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
print -deps example07_1

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
[pdfzkKRIGH]=kerneldensity(Z(126:525)-zkKRIGH,(-3:0.1:3)',0.1);
[pdfzkKRIGHS]=kerneldensity(Z(126:525)-zkKRIGHS,(-3:0.1:3)',0.1);
[pdfzkBME]=kerneldensity(Z(126:525)-zkBME,(-3:0.1:3)',0.1);
plot((-3:0.1:3)',pdfzkKRIGH,':');hold on;
plot((-3:0.1:3)',pdfzkKRIGHS,'--');
plot((-3:0.1:3)',pdfzkBME,'-');
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(-3:1:3));
set(gca,'FontSize',8);
xlabel('X');
ylabel('pdf');
title('(f)');
hold off;

print -deps example07_2

biasKrigingH=mean(Z(126:525)-zkKRIGH)
biasKrigingHS=mean(Z(126:525)-zkKRIGHS)
biasBME=mean(Z(126:525)-zkBME)
MSEkrigingH=(Z(126:525)-zkKRIGH)'*(Z(126:525)-zkKRIGH)/400
MSEkrigingHS=(Z(126:525)-zkKRIGHS)'*(Z(126:525)-zkKRIGHS)/400
MSEBME=(Z(126:525)-zkBME)'*(Z(126:525)-zkBME)/400

save example07



