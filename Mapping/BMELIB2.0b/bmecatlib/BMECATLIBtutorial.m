% BMECATLIBtutorial      - tutorial for the bmecatlib directory
%                          (December 1, 2003)

%%% Clear memory and set echo to on.

clear;
clc;
echo on;

%%% Load the text data file called 'SoilsData.txt' that contains
%%% the spatial coordinates of the locations where soil types
%%% have been observed.

load SoilsData.txt;
c=SoilsData(:,1:2);
zh=SoilsData(:,3);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the categories using
%%% colorplot.m, where categories are coded as 1 (black), 2 (yellow),
%%% 3 (red) and 4 (white).

figure;
colorplot(c,zh);
axis square;
xlabel('Easting');
ylabel('Northing');
title('Observed categories');
drawnow;

%%% Type any key for continuing...

pause;
clc;

%%% Compute the omnidirectional bivariate probability tables
%%% and display theses probabilities as a function of the
%%% distance separating locations.Distance classes are 0.05
%%% units wide and are ranging from 0 to 0.5
%%% The proportion of each category is also displayed.

cl=(0:5:50)';
method='kron';
[d,P,o]=probatablecalc(c,zh,cl,method);

[P{1,1}(1),P{2,2}(1),P{3,3}(1),P{4,4}(1)]

figure;
probatableplot(d,P);
drawnow;

%%% Type any key for continuing...

pause;
clc;

%%% Fits the bivariate probability tables as functions of the
%%% distance separating a pair of locations. Fitting is realized
%%% at distances ranging from 0 to 0.45 with 0.01 step.
%%% The smoothing factor is equal to twice the suggested initial
%%% value, for a polynomial of order 1


dfit=(0:1:45)';
kstd=1*max(d)/length(d);
options=[1 1 1];

figure;
[Pfit]=probatablefit(dfit,d,P,o,kstd,options);
drawnow;

%%% Type any key for continuing...

pause;
clc;

%%% Compute the conditional probability estimates at the nodes
%%% of a two dimensional grid with regular spacing. The grid spacing
%%% is 5 units along both axes.
%%% At each estimation location, there is a maximum of 3 neighbouring
%%% values used for the estimation. These values are lying not further
%%% than 20 distance units from the estimation location.
%%% (the chosen values for nmax and options(2) have been selected
%%% for faster computations in this tutorial).

minc=[0 0];
dc=[5 5];
nc=[21 21];
[ck]=creategrid(minc,dc,nc);

nmax=3;
dmax=20;
options=[1 1e-2];

%%% Type any key for executing the estimation (takes time)...

pause;
clc;

[pkh]=BMEcatHard(ck,c,zh,dfit,Pfit,nmax,dmax,options);

%%% Open a new graphic window and display the categories having
%%% a maximum probability of occurence at the grid locations. 
%%% Data for observed categories are superimposed on the graph.

[trash,zkhmax]=max(pkh');
[ck1,ck2,zkhmax]=col2mat(ck,zkhmax');

figure;
pcolor(ck1,ck2,zkhmax);
colormap hot
hold on;
shading flat
colorplot(c,zh);
axis square;
xlabel('Easting');
ylabel('Northing');
title('Maximum probability estimates (hard data)');
axis square;
drawnow;

%%% Type any key for continuing...

pause;
clc;

%%% Redo the same estimation after data recoding. It is assumed
%%% now that we only know if category 1 or 2 has been observed
%%% or if category 3 or 4 has been observed.

n=size(c,1);
zs=zeros(n,4);
is1or2=find(zh<=2);
is3or4=find(zh>2);
zs(is1or2,1:2)=1;
zs(is3or4,3:4)=1;

%%% Type any key for executing the estimation (takes time)...

pause;
clc;

[pks]=BMEcatSubset(ck,c,zs,dfit,Pfit,nmax,dmax,options);

%%% Open a new graphic window and display the categories having
%%% a maximum probability of occurence at the grid locations. 
%%% Data for observed categories are superimposed on the graph.

[trash,zksmax]=max(pks');
[ck1,ck2,zksmax]=col2mat(ck,zksmax');

figure;
pcolor(ck1,ck2,zksmax);
colormap hot
hold on;
shading flat
colorplot(c,zh);
axis square;
xlabel('Easting');
ylabel('Northing');
title('Maximum probability estimates (categories subsets)');
axis square;
drawnow;

%%% Type any key for continuing...

pause;
clc;

%%% Conditional simulation of categories at the data
%%% locations using a sequential procedure.Simulation
%%% is performed on the grid with observed categories
%%% as conditioning data

%%% Type any key for executing the simulation (takes time)...

pause;
clc;

options=[1 1 1e-3];
[zhsimul]=simucatcond(ck,c,zh,dfit,Pfit,nmax,dmax,options);

%%% Open a new graphic window and display the simulated
%%% categories at the grid locations. 

[ck1,ck2,zhsimul]=col2mat(ck,zhsimul);

figure;
pcolor(ck1,ck2,zhsimul);
colormap hot
shading flat
hold on;
colorplot(c,zh);
axis square;
xlabel('Easting');
ylabel('Northing');
title('Simulated categories');
drawnow;

%%% Type any key for continuing...

pause;
clc;

%%% Save all the variables to the 'BMECATLIBtutorial.mat' binary data
%%% file for possible subsequent use

save BMECATLIBtutorial

%%% End of the BMECATLIB tutorial

echo off;
