% SIMULIBtutorial           - tutorial for the simulib directory (Jan 1,2001)

%%% Clear memory content and set echo to on.

clear;
clc;
echo off;
if exist('BMEINTLIBtutorial.mat')~=2,
  error('You must first run BMEINTLIBtutorial once before running this tutorial');
end;
echo on;

%%% Load the MATLAB binary data file called 'BMEINTLIBtutorial.mat'
%%% (see the BMEINTLIBtutorial.m script file).

load BMEINTLIBtutorial

%%% Type any key for continuing...

pause;
clc;

%%% Simulate sequentially the sand content values at the nodes
%%% of a two dimensional grid with regular spacing. The grid spacing
%%% is 1000 meters along both axes. Simulation is performed without
%%% taking into account the silt and clay content values. 
%%% At each simulation location, there is a maximum of 10 hard sand
%%% content values used for the conditional simulation. The mean is
%%% supposed to be constant over the area.

nhmax=10;
dmax=5000;
order=0;

[sandSIM1]=simuseqcond(ck,ch,sand,modelsand,paramsand,nhmax,dmax,order);
[sandSIM2]=simuseqcond(ck,ch,sand,modelsand,paramsand,nhmax,dmax,order);
[sandSIM3]=simuseqcond(ck,ch,sand,modelsand,paramsand,nhmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the three
%%% simulations as well as the kriging estimates using pcolor.m.

figure;
subplot(2,2,1),
  pcolor(ck1,ck2,SandKRIGU);
  colormap(hot);
  colorbar;
  axis equal;
  title('Kriging estimate','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,2),
  pcolor(ck1,ck2,reshape(sandSIM1,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #1','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,3),
  pcolor(ck1,ck2,reshape(sandSIM2,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #2','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,4),
  pcolor(ck1,ck2,reshape(sandSIM3,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #3','FontSize',10);
  set(gca,'FontSize',9);

%%% Type any key for continuing...

pause;
clc;

%%% Simulate sequentially the sand content values at the nodes
%%% of a two dimensional grid with regular spacing as above, but
%%% now taking into account the additional soft information at the
%%% simulation locations : the simulated sand content values must
%%% lie in the [0%,100%] interval.

ak=zeros(size(ck(:,1)));
bk=ones(size(ck(:,1)))*100;
nsmax=1;

[sandSIMLIM1]=simuseqcondInt(ck,ch,ck,sand,ak,bk,modelCsand,paramsand,nhmax,nsmax,dmax,order);
[sandSIMLIM2]=simuseqcondInt(ck,ch,ck,sand,ak,bk,modelCsand,paramsand,nhmax,nsmax,dmax,order);
[sandSIMLIM3]=simuseqcondInt(ck,ch,ck,sand,ak,bk,modelCsand,paramsand,nhmax,nsmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the three
%%% simulations as well as the kriging estimates using pcolor.m.

figure;
subplot(2,2,1),
  pcolor(ck1,ck2,SandKRIGU);
  colormap(hot);
  colorbar;
  axis equal;
  title('Kriging estimate','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,2),
  pcolor(ck1,ck2,reshape(sandSIMLIM1,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #1','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,3),
  pcolor(ck1,ck2,reshape(sandSIMLIM2,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #2','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,4),
  pcolor(ck1,ck2,reshape(sandSIMLIM3,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #3','FontSize',10);
  set(gca,'FontSize',9);

%%% Type any key for continuing...

pause;
clc;

%%% Simulate sequentially the sand content values at the nodes
%%% of a two dimensional grid with regular spacing as above, by
%%% taking into account the soft information at the simulation
%%% locations as above, plus the soft information from the
%%% IntervalData.txt file. The interval information are merged 
%%% together. There is one additional interval used for the
%%% conditional simulation, which is not located at the simulation
%%% location.

cs_all=[cs;ck];
a_all=[a_sand;ak];
b_all=[b_sand;bk];
nsmax=2;

[sandSIMLIMSOFT1]=simuseqcondInt(ck,ch,cs_all,sand,a_all,b_all,modelCsand,paramsand,nhmax,nsmax,dmax,order);
[sandSIMLIMSOFT2]=simuseqcondInt(ck,ch,cs_all,sand,a_all,b_all,modelCsand,paramsand,nhmax,nsmax,dmax,order);
[sandSIMLIMSOFT3]=simuseqcondInt(ck,ch,cs_all,sand,a_all,b_all,modelCsand,paramsand,nhmax,nsmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the three
%%% simulations as well as the BME estimates using pcolor.m.

figure;
subplot(2,2,1),
  pcolor(ck1,ck2,SandBMEU);
  colormap(hot);
  colorbar;
  axis equal;
  title('BME estimate','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,2),
  pcolor(ck1,ck2,reshape(sandSIMLIMSOFT1,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #1','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,3),
  pcolor(ck1,ck2,reshape(sandSIMLIMSOFT2,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #2','FontSize',10);
  set(gca,'FontSize',9);
subplot(2,2,4),
  pcolor(ck1,ck2,reshape(sandSIMLIMSOFT3,19,17));
  colormap(hot);
  colorbar;
  axis equal;
  title('Simulation #3','FontSize',10);
  set(gca,'FontSize',9);

%%% Type any key for continuing...

pause;
clc;

%%% Save all the variables to the 'SIMULIBtutorial.mat' binary data
%%% file for subsequent use

save SIMULIBtutorial

%%% End of the SIMULIB tutorial

echo off;

