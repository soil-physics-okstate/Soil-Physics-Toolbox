% HYDRUS_main.m
% Created by Derek Groenendyk
% 5/4/2012
% Class that calls and runs HYRDUS in Matlab

clear all;
clc;

exp = 'Central_2';
% exp = 'Cyclic';
mainDirectory = 'C:\Users\Tyson Ochsner\Dropbox (OSU Soil Physics)\Maples\Conventional\';
%expDirectory = [mainDirectory 'Projects\' exp];
expDirectory = [mainDirectory exp];

% order of the soils included in usdaSoils
% sand, loamy sand, sandy loam,	loam, silt,	silt loam, sandy clay loam, clay loam,
% silty clay loam, sandy clay loam,	silty clay,	clay
% usdaSoils = csvread('USDARosettaSoils.csv');

%%
method = 1;

% Single Run
if method == 1

         atmosphIN = ATMOSPHIN(expDirectory);
         fluxestotest = [0.5 1.0 5.0 10.0];
         numfluxes = length(fluxestotest);
         
         for jj=1:numfluxes
         atmosphIN.setData('Prec',fluxestotest(jj));
       atmosphIN.update();
        soilstotest = [1 2 3 4];
    numsoils = length(soilstotest);
        trial = 1;
        Hydrus = HYDRUS_Class(expDirectory);
      %       create Selector.IN object
        selectIN = SELECTORIN(expDirectory);

        noCMDWindow = 1;
    %     set soil parameters
        paramList = {'thr','ths','Alfa','n','Ks'};   
           for ii=1:numsoils;
              ind = soilstotest(ii);
              %ind = 1;
            paramValues = usdaSoils;

            for j=1:length(paramList)
                if j == 5
                    paramValue = roundn(paramValues(ind,j),-2);
                elseif j == 3
                    paramValue = roundn(paramValues(ind,j),-4);
                else
                    paramValue = roundn(paramValues(ind,j),-4);
                end
                selectIN.setData(paramList(j),paramValue,1) % last arg selects material
            end
            %       set other selector.in parameters
                selectIN.setData('TPrint(1),TPrint(2),...,TPrint(MPL)',[1,2,3,4,5,6,7,8,9,10,11,12,13,14])
                selectIN.setData('MPL',14)

        %       commits/writes changes to the file (overwrites previous file)
                selectIN.update(); 

            Hydrus.run_hydrus(noCMDWindow)
            Hydrus.outputResults(exp,[num2str(ii),num2str(jj)])
           end; 
end;
    % Batch Run or Multiple Simulations
elseif method == 2
        
    Hydrus = HYDRUS_Class(expDirectory);
    noCMDWindow = 1;
    
    numSims = 1; % number of simulations to run
    
    paramList = {'thr','ths','Alfa','n','Ks'};
    
    load('newSoilProps.txt');
    paramValues = newSoilProps;
    for i=1:numSims
%       create Selector.IN object
        selectIN = SELECTORIN(expDirectory);
        
%       set soil parameters
        ind = i;
        for j=1:length(paramList)
            if j == 5
                paramValue = roundn(paramValues(ind,j+1)*100.0*(60*60*24),-2);
            elseif j == 3
                paramValue = roundn(paramValues(ind,j+1)/100.0,-4);
            else
                paramValue = roundn(paramValues(ind,j+1),-4);
            end
            selectIN.setData(paramList(j),paramValue,2) % last arg selects material
        end
       
%       set other selector.in parameters
        selectIN.setData('TPrint(1),TPrint(2),...,TPrint(MPL)',[1,2,3,4,5,6,7,8,9,10,11,12,13,14])
        
%       commits/writes changes to the file (overwrites previous file)
        selectIN.update();        
        
%       create ATMOSPH.IN object
        atmosphIN = ATMOSPHIN(expDirectory);
%       set ATMOSPH.IN parameters
        atmosphIN.setData('SinusVar','t')
        atmosphIN.setData('tAtm',[1,6])
        atmosphIN.setData('Prec',[2,0])
        atmosphIN.setData('rSoil',[-1,1])
        atmosphIN.setData('hCritA',[9999,10000])
        atmosphIN.setData('MaxAL',2)
%       commits/writes changes to the file (overwrites previous file)
        atmosphIN.update()        
        
%       create PROFILE.DAT object
        profileDAT = PROFILEDAT(expDirectory);
        
%       values used for setting parameters
        columnLength = 100; % row 4, 2
        numLayers = 101; % row 5, 1
        layerThickness = 1;
        matArray = []; % column 4
        layerArray = []; % column 5
        beta = 0.000000; % column 6
        xz = 1.00000; % column 7, 8, 9

%         set PROFILE.DAT parameters
        profileDAT.setData('columnLength',columnLength)
        profileDAT.setData('numLayers',numLayers) 
        profileDAT.setData('x',0:1:101)
        profileDAT.setData('h',ones([101,1])*-10.0)
        profileDAT.setData('numObs',7)
        profileDAT.setData('obsLoc',[10,35,77])
        
%       commits/writes changes to the file (overwrites previous file)        
        profileDAT.update()

        Hydrus.run_hydrus(noCMDWindow)
        trial = i;
        Hydrus.outputResults(exp,num2str(trial))
    end
end




