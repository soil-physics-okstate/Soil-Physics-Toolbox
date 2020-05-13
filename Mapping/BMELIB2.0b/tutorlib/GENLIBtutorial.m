% GENLIBtutorial            - tutorial for the genlib directory (Jan 1,2001)

%%% Clear memory content and set echo to on.

clear;
clc;
echo on;

%%% Load the data file called 'Falmagne.txt' and display 
%%% its content (see the Falmagne.txt file for the meaning
%%% of these variables).

[val,valname,filetitle]=readGeoEAS('Falmagne.txt');
ch=val(:,1:2);
sand=val(:,3);
silt=val(:,4);
clay=val(:,5);
code=val(:,6);
whos

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display as circles the locations
%%% given in the file with correct axes proportions for the graph,
%%% then freeze the graphic.

figure;
plot(ch(:,1),ch(:,2),'o');
axis equal
hold on;

%%% Type any key for continuing...

pause;
clc;

%%% Create a two dimensional grid with regular spacing along the
%%% axes and superimpose this grid over the area. The grid spacing
%%% is 1000 meters along both axes.

minc=[178000 90000];
dc=[1000 1000];
nc=[17 19];
[ck]=creategrid(minc,dc,nc);

plot(ck(:,1),ck(:,2),'+r');

%%% Type any key for continuing...

pause;
clc;

%%% Perform an estimation of the sand content at the nodes of the
%%% grid using two methods :
%%%
%%% i)  a nearest neighbour estimation ;
%%%
%%% ii) a kernel smoothing with a variance of the Gaussian kernel
%%%     equal to 1e6.
%%%
%%% In both cases, the values taken into account for the estimations
%%% are the 20 closest values, as long as their distance from the
%%% estimation location does not exceed 50000 meters. Only the values
%%% of the sand content over the area are taken into account for the
%%% estimation.

nhmax=20;
dmax=50000;

[zkID]=invdist(ck,ch,sand,Inf,nhmax,dmax);
[zkKS]=kernelsmoothing(ck,ch,sand,1e6,nhmax,dmax);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the result of the estimation
%%% for both methods using pcolor.m, after a conversion for the data
%%% format has been performed with col2mat.m

figure;
[ck1,ck2,ZkID]=col2mat(ck,zkID);
pcolor(ck1,ck2,ZkID);
colormap(hot);
colorbar;
axis equal;
title('Nearest neighbour');

figure;
[ck1,ck2,ZkKS]=col2mat(ck,zkKS);
pcolor(ck1,ck2,ZkKS);
colormap(hot);
colorbar;
axis equal
title('Kernel smoothing');

%%% End of the GENLIB tutorial

echo off;




