% GRAPHLIBtutorial          - tutorial for the graphlib directory (Jan 1,2001)

%%% Clear memory content and set echo to on.

clear;
clc;
echo on;

%%% Load the data file called 'Falmagne.txt' and display 
%%% its content (see the Falmagne.txt file for the meaning
%%% of these variables). ch has two columns which are the  
%%% the spatial coordinates (in meters). The next three
%%% variables have the values for the sand, silt and clay contents (in %).
%%% The last variable is a numeric code for the corresponding soil type.


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

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the histograms
%%% and the scatter plots for the sand; silt and clay
%%% contents

figure;
histscatterplot([sand silt clay]);

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the values for
%%% sand content using colorplot.m, with the 'hot' color
%%% map and squares having a black edge color. Axes are
%%% set to correct proportions

figure;
Property={'Marker','MarkerSize','MarkerEdgeColor'};
Value ={'s',12,[0 0 0]};
colorplot(ch,sand,'hot',Property,Value);
axis equal

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display the values for
%%% silt content using markerplot.m, with circles having
%%% a blue edge color and a transparent face color. Axes
%%% are set to correct proportions.

figure;
sizelim=[3 20];
Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
Value={'o',[0 0 1],'none'};
markerplot(ch,silt,sizelim,Property,Value);
axis equal

%%% Type any key for continuing...

pause;
clc;

%%% Open a new graphic window and display :
%%%
%%% i)  in the left hand side the values for clay content using
%%%     poleplot.m, with aquamarine circles having a black
%%%     border and the default size ;
%%%
%%% ii) in the right hand side the value 0 for clay content below
%%%     15 % and the value 1 for clay content above 15 % using a
%%%     blue Times font with letter size equal to 9 points
%%%
%%% The graphic window can be manually resized for improving the
%%% vizualization

figure;
subplot(1,2,1);
  Property={'Marker','MarkerEdgeColor','MarkerFaceColor'};
  Value={'o',[0 0 0],[127/255 1 212/255]};
  poleplot(ch,clay,Property,Value);
subplot(1,2,2);
  Property={'Color','FontName','FontSize'};
  Value={[0.5 0 0.9],'Times',9};
  valplot(ch,clay>15,Property,Value);

%%% Type any key for continuing...

pause;
clc;

%%% Save for further use the variables from the Falmagne.txt file
%%% into a MATLAB binary data file called Falmagne.mat

save Falmagne ch sand silt clay code

%%% End of the GRAPHLIB tutorial

echo off;

