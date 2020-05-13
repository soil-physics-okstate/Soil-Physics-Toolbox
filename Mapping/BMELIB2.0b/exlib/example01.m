% example01                 - Example VI.3.4 : Using probabilistic data from a secondary variable
%
% There are 25 hard data randomly located and 400 soft 
% probabilistic data at the nodes of a 20 x 20 grid. The soft 
% pdf are obtained by conditioning with respect to two classes 
% (value<0 or value>0) the variables of interest, using a
% secondary variable correlated with the first one (rho=0.9).
% The objective is to predict at best the values at the nodes
% of the grid. Three approaches are proposed : 
%
% 1- kriging with only 25 hard data
% 2- kriging with 25 hard data and the conditional soft pdf at
%    prediction point taken as "hardened" value
% 3- bme with rigorous incorporation of hard data and the
%    soft pdf at the estimation point
%
% The first approach is deficient because it does not make use of
% the 400 soft data. For the second approach, because kriging is
% an exact predictor, using as hard data the expectation of the pdf
% does not allow to take into account the hard data at all. Only BME
% is doing a correct job in this case.
  
clear;
rand('seed',3);
randn('seed',3);
disp('Creating data');

classe=[-5,0,5];
z=(-5:0.1:5)';
K=[1,0.9;0.9,1];
for i=1:length(classe)-1,
  [pdfclasse{i}]=gausscondpdf(z,classe(i),classe(i+1),[0;0],K);
  [Epdf(i),Vpdf(i)]=pdfstat(z,pdfclasse{i});
end;

[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
ch=rand(25,2);
covmodel='exponentialC';
covparam=[1 1];
ns=1;
[Z,L]=simuchol([ch;ck],covmodel,covparam,ns);
zh=Z(1:25);
index=zeros(400,1);
for i=1:400,
 index(i)=sum(Z(25+i)>classe);
end;
Epdf=Epdf(index);Epdf=Epdf';
nhmax=10;
nsmax=3;
dmax=1;
order=NaN;

disp('Kriging with hard data...');
[zkKRIGH,vk]=kriging(ck,ch,zh,covmodel,covparam,nhmax,dmax,order);
disp('BME with hard data + pdf...');
softpdftype=2;
nl=length(z);
options=BMEoptions;
options(8)=1;
for i=1:400,
  [moments,info]=BMEprobaMoments(ck(i,:),ch,ck(i,:),zh,softpdftype,...
    nl,z',pdfclasse{index(i)}',covmodel,covparam,nhmax,nsmax,dmax,order,options);
  zkBME(i)=moments(1);  
end;

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0.5 0.5 18.5 7.5]);
set(gcf,'Position',[0.5 0.5 18.5 7.5]);
subplot(1,2,1);
[x,y]=meshgrid(-3:0.2:3,-3:0.2:3);
f=gaussbipdf(x,y,[0 0 1 1 0.9]);
contour((-3:0.2:3)',(-3:0.2:3)',f,(0:0.05:0.3));hold on;
plot([-3,3],[0 0],'--');
axis([-3 3 -3 3]);
axis('square');
set(gca,'FontSize',8);
set(gca,'XTick',-3:3);
set(gca,'YTick',-3:3);
title('(a)');
ax=axis;
text(mean(ax(1:2)),ax(3)-0.11*diff(ax(3:4)),'X_1',...
  'HorizontalAlignment','center','FontSize',10);
ylabel('X_2','FontSize',10);

subplot(1,2,2);
for i=1:length(classe)-1,
 plot(z,pdfclasse{i});hold on;
end;
axis([-4 4 0 0.7]);
axis('square');
set(gca,'FontSize',8);
set(gca,'XTick',-4:4);
set(gca,'YTick',0:0.1:0.7);
title('(b)');
ax=axis;
text(mean(ax(1:2)),ax(3)-0.11*diff(ax(3:4)),'X_1',...
  'HorizontalAlignment','center','FontSize',10);
ylabel('conditional pdf','FontSize',10);

print -deps example01_1

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
pcolor((0:1/19:1),(0:1/19:1),reshape(Z(26:425),20,20));
caxis([min(zh)+0.5 max(zh)-0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(a)');
ax=axis;
text(mean(ax(1:2)),ax(3)-0.14*diff(ax(3:4)),'s_1',...
  'HorizontalAlignment','center','FontSize',8);
ylabel('s_2');
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
pcolor((0:1/19:1),(0:1/19:1),reshape(Epdf,20,20));
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
[pdfzkKRIGE]=kerneldensity(Z(26:425)-Epdf,(-3:0.1:3)',0.1);
[pdfzkKRIGH]=kerneldensity(Z(26:425)-zkKRIGH,(-3:0.1:3)',0.1);
[pdfzkBME]=kerneldensity(Z(26:425)-zkBME',(-3:0.1:3)',0.1);
plot((-3:0.1:3)',pdfzkKRIGH,':');hold on;
plot((-3:0.1:3)',pdfzkKRIGE,'--');
plot((-3:0.1:3)',pdfzkBME,'-');
Position=get(gca,'Position');
Position(3)=Position(3)*0.85;
set(gca,'Position',Position);
set(gca,'XTick',(-3:1:3));
set(gca,'FontSize',8);
axis([-2 2 0 0.9]);
xlabel('estimation error');
ylabel('pdf');
title('(f)');
hold off;

print -deps example01_2

biasExpect=mean(Z(26:425)-Epdf)
biasKriging=mean(Z(26:425)-zkKRIGH)
biasBME=mean(Z(26:425)-zkBME')
MSEinterval=(Z(26:425)-Epdf)'*(Z(26:425)-Epdf)/400
MSEkriging=(Z(26:425)-zkKRIGH)'*(Z(26:425)-zkKRIGH)/400
MSEBme=(Z(26:425)-zkBME')'*(Z(26:425)-zkBME')/400

save example01


