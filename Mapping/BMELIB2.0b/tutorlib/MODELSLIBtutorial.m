% MODELSLIBtutorial         - tutorial for the modelslib directory (Jan 1,2001)

%%% Clear memory content and set echo to on.

clear;
clc;
echo on;

%%% Plot the graphs of the various variogram models
%%% with unit sill (or sill equivalent) parameter

d=(0:0.02:1.5)';

figure;
set(gcf,'Units','centimeters');
set(gcf,'Position',[3 3 14 17]);

subplot(4,2,1);
  plot(d,nuggetV(d,[1]));
  title('nuggetV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,2);
  plot(d,exponentialV(d,[1 1]));
  title('exponentialV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,3);
  plot(d,sphericalV(d,[1 1]));
  title('sphericalV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,4);
  plot(d,gaussianV(d,[1 1]));
  title('gaussianV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,5);
  plot(d,holesinV(d,[1 0.5]));
  title('holesinV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.5]);
subplot(4,2,6);
  plot(d,holecosV(d,[1 0.5]));
  title('holecosV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 2]);
subplot(4,2,7);
  plot(d,linearV(d,1));
  title('linearV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 2]);
  xlabel('distance','FontSize',10);
  ylabel('variogram','FontSize',10);
subplot(4,2,8);
  plot(d,powerV(d,[2 1]));
  title('powerV','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 2.5]);

%%% Type any key for continuing...

pause;
clc;

%%% Plot the graphs of the various covariance models
%%% with unit sill parameter

figure;
set(gcf,'Units','centimeters');
set(gcf,'Position',[4 4 14 17]);

subplot(4,2,1);
  plot(d,nuggetC(d,[1]));
  title('nuggetC','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,2);
  plot(d,exponentialC(d,[1 1]));
  title('exponentialC','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,3);
  plot(d,sphericalC(d,[1 1]));
  title('sphericalC','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,4);
  plot(d,gaussianC(d,[1 1]));
  title('gaussianC','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 0 1.2]);
subplot(4,2,5);
  plot(d,holesinC(d,[1 0.5]));
  title('holesinC','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 -0.5 1.2]);
  xlabel('distance','FontSize',10);
  ylabel('covariance','FontSize',10);
subplot(4,2,6);
  plot(d,holecosC(d,[1 0.5]));
  title('holecosC','FontSize',10);
  set(gca,'FontSize',8);
  axis([0 1.5 -1 1]);

%%% Type any key for continuing...

pause;
clc;

%%% Plot a nested variogram models that includes a nugget
%%% effect with sill equal to 0.2, and two spherical models
%%% with sill equal to 0.4 and range equal to 0.3 and 1.2,
%%% respectively.

figure;
model={'nuggetV','sphericalV','sphericalV'};
param={[0.2],[0.4 0.3],[0.4 1.2]};

modelplot(d,model,param);

%%% End of the MODELSLIB tutorial

echo off;






