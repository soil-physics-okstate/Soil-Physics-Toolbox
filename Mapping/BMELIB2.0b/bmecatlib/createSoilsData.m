clear;
rand('seed',1);
randn('seed',1);
c=rand(200,2)*100;
z=simuchol(c,'exponentialC',[1 50]);
z=(z>-0.6745)+(z>0)+(z>0.6745)+1;
Data=[c,z];
save SoilsData.txt Data -ascii
clear;