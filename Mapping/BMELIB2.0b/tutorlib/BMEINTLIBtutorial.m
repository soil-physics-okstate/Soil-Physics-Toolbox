% BMEINTLIBtutorial         - tutorial for the bmeintlib directory (Jan 1,2001)

%%% Clear memory content and set echo to on.

clear;
clc;
echo off;
if exist('BMEHRLIBtutorial.mat')~=2,
  error('You must run BMEHRLIBtutorial once before running this tutorial');
end;
echo on;

%%% Load the MATLAB binary data file called 'BMEHRLIBtutorial.mat'
%%% (see the BMEHRLIBtutorial.m script file). This file contains
%%% the values for many of the input variables that are used here.


load BMEHRLIBtutorial

%%% Type any key for continuing...

pause;
clc;

%%% Load the text data file called 'IntervalData.txt' that contains
%%% soft information for the sand, silt and clay contents values.
%%% Set the variables from the content of the file. The two first
%%% columns are the spatial coordinates (in meters). The next six
%%% columns are the lower and upper limits of the interval of values
%%% for the sand, silt and clay contents (in %) at the same additional
%%% coordinates, respectively.

load IntervalData.txt;

cs=IntervalData(:,1:2);
a_sand=IntervalData(:,3);
b_sand=IntervalData(:,4);
a_silt=IntervalData(:,5);
b_silt=IntervalData(:,6);
a_clay=IntervalData(:,7);
b_clay=IntervalData(:,8);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the values for
%%% sand content using colorplot.m, with the 'hot' color
%%% map. Hard data are represented by squares, whereas soft
%%% data are represented by black disks having a size which
%%% is proportional to the class order.

figure;

Sizelim=[5 20];
Property={'Marker','MarkerFaceColor','MarkerEdgeColor'};
Value ={'o','none',[0 0 0]};
markerplot(cs,a_sand,Sizelim,Property,Value);
hold on;

Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',12,[0 0 0]};
colorplot(ch,sand,'hot',Property,Value);
axis equal

%%% Type any key for continuing...

pause;
clc;

%%% Compute the BME estimate of the sand content at the nodes
%%% of a two dimensional grid with regular spacing. The grid spacing
%%% is 1000 meters along both axes. BME is realized on the raw
%%% sand content values without taking into account the silt and
%%% clay content values. 
%%% At each estimation location, there is a maximum of 5 hard and 3
%%% soft sand content values used for the estimation. These values
%%% are lying not further than 10 000 meters from the estimation location.
%%% The mean is supposed to be constant over the area.

modelCsand={'nuggetC','exponentialC'};
nhmax=5;
nsmax=3;
dmax=10000;
order=0;

[sandBMEU]=BMEintervalMode(ck,ch,cs,sand,a_sand,b_sand,modelCsand,paramsand,nhmax,nsmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% for BME with interval data using pcolor.m, after a conversion for
%%% the data format has been performed with col2mat.m

figure;
[ck1,ck2,SandBMEU]=col2mat(ck,sandBMEU);
pcolor(ck1,ck2,SandBMEU);
colormap(hot);
colorbar;
axis equal;
title('BME sand estimate - univariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% for Kriging with only hard data, as the output obtained from the
%%% BMEHRLIBtutorial.m script file.

figure;
pcolor(ck1,ck2,SandKRIGU);
colormap(hot);
colorbar;
axis equal;
title('Kriging sand estimate - univariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Compute the multivariate BME estimate of the sand content at
%%% the nodes of the same two dimensional grid. Silt and clay
%%% content values are used for improving the quality of the estimation
%%% for the sand content.
%%% At each estimation location, there is a maximum of 5 hard and 1 soft
%%% values for the sand, silt and clay contents used for the estimation.
%%% These value are lying not further than 10000 meters from the
%%% estimation location. As for the sand content, the mean is supposed
%%% to be constant over the area for the silt and clay contents.
%%% For input in the function, the data and coordinates for the sand, silt
%%% and clay intervals are converted to the appropriate cell arrays.

[csindex,all_a]=noindex2index({cs,cs,cs},{a_sand,a_silt,a_clay});
[csindex,all_b]=noindex2index({cs,cs,cs},{b_sand,b_silt,b_clay});

modelC={'nuggetC','exponentialC'};
nhmax=[5;5;5];
nsmax=[1;1;1];
dmax=10000;
order=[0;0;0];

[sandBMEM]=BMEintervalMode(ckindex,chindex,csindex,all_hard,all_a,all_b,modelC,param,nhmax,nsmax,dmax,order);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% using pcolor.m, after a conversion for the data format has been
%%% performed with col2mat.m

figure;
[ck1,ck2,SandBMEM]=col2mat(ck,sandBMEM);
pcolor(ck1,ck2,SandBMEM);
colormap(hot);
colorbar;
axis equal;
title('BME sand estimate - multivariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% for Kriging with only hard data, as the output obtained from the
%%% BMEHRLIBtutorial.m script file.

figure;
pcolor(ck1,ck2,SandKRIGM);
colormap(hot);
colorbar;
axis equal;
title('Kriging Sand estimate - multivariate case');

%%% Type any key for continuing...

pause;
clc;

%%% Save all the variables to the 'BMEINTLIBtutorial.mat' binary data
%%% file for subsequent use

save BMEINTLIBtutorial

%%% End of the BMEINTLIB tutorial

echo off;




