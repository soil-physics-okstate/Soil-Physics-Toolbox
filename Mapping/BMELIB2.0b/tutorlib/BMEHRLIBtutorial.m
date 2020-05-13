% BMEHRLIBtutorial          - tutorial for the bmehrlib directory (Jan 1,2001)

%%% Clear memory content and set echo to on.

clear;
clc;
if exist('STATLIBtutorial.mat')~=2,
  error('You must run STATLIBtutorial once before running this tutorial');
end;
echo on;

%%% Load the MATLAB binary data file called 'STATLIBtutorial.mat'
%%% (see the STATLIBtutorial.m script file).

load STATLIBtutorial

%%% Type any key for continuing...

pause;
clc;

%%% Compute the kriging estimate of the sand content at the nodes
%%% of a two dimensional grid with regular spacing. The grid spacing
%%% is 1000 meters along both axes. Kriging is performed on the raw
%%% sand content values without taking into account the silt and
%%% clay content values. 
%%% At each estimation location, there is a maximum of 5 sand content
%%% values used for the estimation. These values are lying not further
%%% than 10 000 meters from the estimation location. The mean is supposed
%%% to be constant over the area.

minc=[178000 90000];
dc=[1000 1000];
nc=[17 19];
[ck]=creategrid(minc,dc,nc);

nhmax=5;
dmax=10000;
order=0;

[sandKRIGU,vsandKRIGU]=kriging(ck,ch,sand,modelsand,paramsand,nhmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% using pcolor.m, after a conversion for the data format has been
%%% performed with col2mat.m

figure;
[ck1,ck2,SandKRIGU]=col2mat(ck,sandKRIGU);
pcolor(ck1,ck2,SandKRIGU);
colormap(hot);
colorbar;
axis equal;
title('Sand estimate - univariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Compute the cokriging estimate of the sand content at the nodes
%%% of the same two dimensional grid. Kriging is realized on the
%%% raw values for sand, silt and clay content values. Silt and clay
%%% content values are used for improving the quality of the estimation
%%% for the sand content.
%%% At each estimation location, there is a maximum of 5 values
%%% for the sand, silt and clay contents used for the estimation.
%%% These value are lying not further than 10000 meters from the
%%% estimation location. As for the sand content, the mean is supposed
%%% to be constant over the area for the silt and clay contents.
%%% For input in kriging.m, the data and coordinates for the sand, silt
%%% and clay contents as well as the coordinates for the estimation
%%% locations are converted to the appropriate cell arrays.

ckindex={ck,1};
[chindex,all_hard]=noindex2index({ch,ch,ch},{sand,silt,clay});

nhmax=[5;5;5];
dmax=10000;
order=[0;0;0];

[sandKRIGM,vsandKRIGM]=kriging(ckindex,chindex,all_hard,model,param,nhmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% using pcolor.m, after a conversion for the data format has been
%%% performed with col2mat.m

figure;
[ck1,ck2,SandKRIGM]=col2mat(ck,sandKRIGM);
pcolor(ck1,ck2,SandKRIGM);
colormap(hot);
colorbar;
axis equal;
title('Sand estimate - multivariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Compute the cokriging estimate of the sand content at the nodes
%%% of the same two dimensional grid, but Kriging is realized on the
%%% transformed sand, silt and clay content values. 

yfile={sandfileK,siltfileK,clayfileK};
cdfyfile={cdfsandfileK,cdfsiltfileK,cdfclayfileK};

[sandKRIGMT,vsandKRIGMT]=krigingT(ckindex,chindex,all_hard,modelT,paramT,nhmax,dmax,yfile,cdfyfile);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% using pcolor.m, after a conversion for the data format has been
%%% performed with col2mat.m

figure;
[ck1,ck2,SandKRIGMT]=col2mat(ck,sandKRIGMT);
pcolor(ck1,ck2,SandKRIGMT);
colormap(hot);
colorbar;
axis equal;
title('Sand estimate with transformation - multivariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Save all the variables to the 'BMEHRLIBtutorial.mat' binary data
%%% file for subsequent use

save BMEHRLIBtutorial

%%% End of the BMEHRLIB tutorial

echo off;




