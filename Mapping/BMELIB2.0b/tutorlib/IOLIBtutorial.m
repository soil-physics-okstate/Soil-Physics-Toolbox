% IOLIBtutorial             - tutorial for the iolib directory (Jan 1,2001)

clear;
clc;
echo on;

%%% Tutorial for the iolib directory. In this tutorial we learn
%%% how to write hard data to a file, and how to read it from that file
%%% how to write soft data to a file, and how to read it from that file
%%% how to write all the parameters and data needed for a BME estimation
%%% to a file, and how to read that file to do the BME estimation.
%%%
%%% This tutorial does not generate any graphic, so make the window as 
%%% wide as possible or you won't be able to read all the lines (some
%%% commands are very long).

%%% Type any key for continuing...
pause;
clc;

%%% Generate the coordinates and values of the hard data,
%%% write this hard data to a file using the writeGeoEAS.m command,
%%% read that file back with the readGeoEAS.m command, and verify
%%% that the value read back are the same as the original.

ch=[0 0;1 3;2 4];
zh=[0.4;1.2;1.3];

writeGeoEAS([ch zh],{'x','y','Variable 1'},'File with hard data','hard.dat');

[val,valname,filetitle]=readGeoEAS('hard.dat');
chF=val(:,1:2);
zhF=val(:,3);

disp(chF)
disp(zhF)

%%% Type any key for continuing...
pause;
clc;

%%% Generate the coordinates and values of the soft data,
%%% write this soft data to a file using the writeProba.m command,
%%% read that file back with the readProba.m command, and verify
%%% that the value read back for limi is the same as the original.

cs=[1 4;3 3];
indexs=[1;1];
softpdftype=2;
nl=[4;3];
limi=[0.1 0.2 0.3 1.0;0.3 0.6 1.3 NaN];
probdens=[0 2.0 2.0 0;0 2.0 0 NaN];

writeProba({cs,indexs},0,softpdftype,nl,limi,probdens,'File with proba data','proba.dat');

[csF,isST,softpdftypeF,nlF,limiF,probdensF,filetitle]=readProba('proba.dat');
csF=csF{1};

disp(limiF)

%%% Type any key for continuing...
pause;
clc;

%%% Now that we have the hard and soft data, set the 
%%% coordinate of the estimation point, set the parameters
%%% needed for the BME estimation, and calculate the 
%%% BME mode estimate at the estimation point.

ck=[2 2];

BMEmethod='BMEprobaMode';
cNugget=0.05;               
cc=1;aa=1.0;
covmodel={'nuggetC','exponentialC'};
covparam={cNugget,[cc 3*aa]};
nhmax=4; 
nsmax=4;
dmax=100;
order=0;
options=BMEoptions;

zk1=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
disp(['zk1=' num2str(zk1)]);

%%% Type any key for continuing...
pause;

clc;
disp(['zk1=' num2str(zk1)]);

%%% Please note the value of zk1.

%%% Type any key for continuing...
pause;

%%% We can write  the BME parameters, the hard data, the soft data,
%%% and the coordinates of the estimation points in four separate 
%%% files with a single command using writeBMEproba.m. To try this 
%%% out, write all the BME parameters and data to files using
%%% writeBMEproba.m, clear all variables from MATLAB memory,
%%% read back all the BME parameters and data using readBMEproba.m,
%%% and using this BME data re-calculate the BME mode estimate.

fileNames={'BME.par';'hard.dat';'proba.dat';'est.dat';'BME.out';'BME.dbg'};

writeBMEproba(BMEmethod,fileNames,ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);

clear;

[BMEmethod,outputfiles,ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options]=readBMEproba('BME.par');

[zk2,info]=BMEprobaMode(ck,ch,cs,zh,softpdftype,nl,limi,probdens,covmodel,covparam,nhmax,nsmax,dmax,order,options);
disp(['zk2=' num2str(zk2)]);

%%% Note that the value of zk2 is equal to zk1 calculated earlier, which 
%%% suggests that all the paramters and data was correctly stored in the files.

%%% Type any key for continuing...
pause;

clc;
disp(['zk2=' num2str(zk2)]);

%%% Note the value of zk2.

%%% Type any key for continuing...
pause;

%%% BMEprobaFromFile.m is a completely automated command 
%%% for doing BME estimation using BME files. This command
%%% reads the data from files, perform the BME estimation
%%% calculation, and return the estimated value in 
%%% an output file. Let's try this out to re-calculate
%%% the BME mode estimate.

clear;
[zk3,info]=BMEprobaFromFile('BME.par');
disp(['zk3=' num2str(zk3)]);

%%% Note that the value of zk3 above is equal to zk2 and zk1 calculated 
%%% earlier, which suggest that BMEprobaFromFile.m is working correctly.

%%% Type any key for continuing...
pause;
clc;

%%% Hence we have learned that BME stores data in different files.
%%% To learn more about the format of these files, check out the
%%% following files that were generated in this tutorial:
%%%
%%% hard.dat        GeoEAS file with the hard data,
%%% soft.dat        Proba file with the probabilistic data,
%%% BME.par         File with the BME paramters,
%%% est.dat         GeoEAS file with the coordinates of estimation pts,
%%% BME.out         GeoEAS output file with the BME estimates.
%%%
%%% To learn more about the syntax of BME parameter file, see help ioparamsyntax.

echo off;
