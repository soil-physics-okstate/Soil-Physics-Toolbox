% example00                 - Example II.4.3 : Using 1-norm metric
%
% Simulation of a space random field using a 1-norm metric,
% with an exponential covariance model having a sill of 1
% and a range equal to 3.
% There are 20 simulated values on the field. Simple kriging
% is performed using the correct model and the Euclidean metric
% on a 20 x 20 grid over a square of unit size.

clear;
rand('seed',5);
randn('seed',3);
[ck]=creategrid([0 0],[1/19 1/19],[20 20]);
nh=20;
ch=rand(nh,2);
unit=ones(400+nh,1);
covmodel='exponentialC';
covparam=[1 3];

%%%%%% Simulation with 1-norm metric

disp('Simulating data...');
dif=kron(unit,[ck;ch])-kron([ck;ch],unit);
dist=sum(abs(dif'))';
[K]=exponentialC(dist,covparam);
K=reshape(K,400+nh,400+nh);
L=chol(K);
Z=L'*randn(400+nh,1);
zk=Z(1:400);
zh=Z(401:400+nh);

%%%%%% Kriging with 1-norm metric

disp('Kriging with 1norm metric...');
zk1norm=zeros(400,1);
Kh=K(401:400+nh,401:400+nh);
for i=1:400,
  dif=[ch(:,1)-ck(i,1) ch(:,2)-ck(i,2)];
  dist=sum(abs(dif'))';
  [kh]=exponentialC(dist,covparam);
  lam=inv(Kh)*kh;
  zk1norm(i)=zh'*lam;
end;

%%%%%% Kriging with Euclidean metric

disp('Kriging with Euclidean metric...');
[zkEuclid]=kriging(ck,ch,zh,covmodel,covparam,Inf,Inf,NaN);

%%%%%% Figures

disp('Creating figures...');

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);

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
title('(a)');
hold off

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
title('(b)');
hold off;

plotmore=0;       % Set this variable to 1 in order to create additional subplots
if plotmore==1
  subplot(3,2,3);
  pcolor((0:1/19:1),(0:1/19:1),reshape(zk1norm,20,20));
  caxis([min(zh)+0.5 max(zh)-0.5]);
  cb=colorbar('vert');
  set(cb,'FontSize',8);
  hold on;
  colormap(gray);
  shading flat;
  set(gca,'XTick',(0:0.2:1));
  set(gca,'YTick',(0:0.2:1));
  set(gca,'FontSize',8);
  title('(c)');
  hold off
  
  subplot(3,2,4);
  pcolor((0:1/19:1),(0:1/19:1),reshape(zkEuclid,20,20));
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
  hold off
  
  subplot(3,2,6);
  [pdfzk1norm]=kerneldensity(zk-zk1norm,(-3:0.1:3)',0.1);
  [pdfzkEuclid]=kerneldensity(zk-zkEuclid,(-3:0.1:3)',0.1);
  plot((-3:0.1:3)',pdfzk1norm,'-');hold on;
  plot((-3:0.1:3)',pdfzkEuclid,'--');
  Position=get(gca,'Position');
  Position(3)=Position(3)*0.85;
  set(gca,'Position',Position);
  set(gca,'XTick',(-3:1:3));
  set(gca,'FontSize',8);
  title('(e)');
  axis([-1.3 1 0 0.9]);
  hold off;
  
end
  
bias1norm=mean(zk-zk1norm)
biasEuclid=mean(zk-zkEuclid)
MSE1norm=(zk-zk1norm)'*(zk-zk1norm)/400
MSEEuclid=(zk-zkEuclid)'*(zk-zkEuclid)/400

print -deps example00_1

[Gxzk,Gyzk]=gradient(reshape(zk,20,20));
[Gxzk1norm,Gyzk1norm]=gradient(reshape(zk1norm,20,20));
[GxzkEuclid,GyzkEuclid]=gradient(reshape(zkEuclid,20,20));

figure;
set(gcf,'Units','centimeters');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,17,21]);
%set(gcf,'Position',[0,0,17,21]);

subplot(3,2,1);
pcolor((0:1/19:1),(0:1/19:1),reshape(Gxzk,20,20));
caxis([-0.5 0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(a)');
hold off

subplot(3,2,2);
pcolor((0:1/19:1),(0:1/19:1),reshape(Gyzk,20,20));
caxis([-0.5 0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(b)');
hold off

subplot(3,2,3);
pcolor((0:1/19:1),(0:1/19:1),reshape(Gxzk1norm,20,20));
caxis([-0.5 0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(c)');
hold off

subplot(3,2,4);
pcolor((0:1/19:1),(0:1/19:1),reshape(Gyzk1norm,20,20));
caxis([-0.5 0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(d)');
hold off

subplot(3,2,5);
pcolor((0:1/19:1),(0:1/19:1),reshape(GxzkEuclid,20,20));
caxis([-0.5 0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(e)');
hold off

subplot(3,2,6);
pcolor((0:1/19:1),(0:1/19:1),reshape(GyzkEuclid,20,20));
caxis([-0.5 0.5]);
cb=colorbar('vert');
set(cb,'FontSize',8);
hold on;
colormap(gray);
shading flat;
set(gca,'XTick',(0:0.2:1));
set(gca,'YTick',(0:0.2:1));
set(gca,'FontSize',8);
title('(f)');
hold off

print -deps example00_2

save example00
