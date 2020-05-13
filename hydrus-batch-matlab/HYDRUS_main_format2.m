% HYDRUS_main.m
% Created by Derek Groenendyk
% c: 5/4/2012
% m: 6/1/2015
% Class that calls and runs HYRDUS in Matlab

clear all;
clc;

noCMDWindow = 0; %set to 1 to NOT have to press enter after each simulation
        
exp = 'Heat_Demo_colin'; % directory main name
% exp = 'Cyclic'; % directory main name
% mainDirectory = 'C:\Temp\HYDRUS_Data\';
mainDirectory = 'C:\Derek\ProgrammingFolder\HYDRUS_Data\Projects\';
expDirectory = [mainDirectory exp];

Hydrus = HYDRUS_Class(expDirectory);

% order of the soils included in usdaSoils:
% sand, loamy sand, sandy loam,	loam, silt,	silt loam, sandy clay loam, 
% clay loam, silty clay loam, sandy clay loam, silty clay, clay
%
% all of the van Genuchten soils for the 12 main soil types
% usdaSoils = csvread([mainDirectory 'USDARosettaSoils.csv']); 

% paramValues=usdaSoils;
% paramList = {'thr','ths','Alfa','n','Ks'};


% list of variable properties with min and max values
topflux=[0 5];
topperiod=[1 30];
soiltype=[1 12];
numprinttimes=14;
endtime=140;
infiltcycles=3;
topBCtimes = (1:2*infiltcycles)*endtime/infiltcycles/2;
% rsoilhold=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
rsoilhold = repmat([-1],1,36);
% Prechold=[1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
Prechold = repmat([1 0],1,18);

% number of realizations to run
numrealizs=0;


if false
%     create PROFILE.DAT object
    profileDAT = PROFILEDAT(expDirectory);

%     values used for setting parameters
    columnLength = 100; % location in file: row 4, 2
    numLayers = 101; % always 1 greater then the length, location in file: row 5, 1
    layerThickness = 1; % in simulation length units
    matArray = []; % location in file: column 4
    layerArray = []; % location in file: column 5
    beta = 0.000000; % location in file: column 6
    xz = 1.00000; % location in file: column 7, 8, 9

%     set PROFILE.DAT parameters
    profileDAT.setData('columnLength',columnLength)
    profileDAT.setData('numLayers',numLayers)
    profileDAT.setData('x',0:1:100)
    profileDAT.setData('h',ones([numLayers,1])*-50.)
    profileDAT.setData('obsLoc',[10,20,30,40,50,60,70])

%     commits/writes changes to the file (overwrites previous file)        
    profileDAT.update()
end

    
paramList = {'Qn','Qo','Dispr.','B1','B2','B3','Cn', 'Co','Cw'};
paramList = {'B3','Cn', 'Co','Cw'};
param_heat = [0.57 0 5 1.56728e+016 2.53474e+016 9.89388e+016 1.43327e+014 1.8737e+014 3.12035e+014 ]; %loam
param_heat = [9.89388e+016 1.43327e+014 1.8737e+014 3.12035e+014];

selectIN = SELECTORIN(expDirectory);

ind = 1; %Soil type number
paramValues = param_heat;
for ii = 1:2
    if ii == 5
        paramValue = roundn(paramValues(ind,ii),-2);
    elseif ii == 3
        paramValue = roundn(paramValues(ind,ii),-4);
    else
        paramValue = roundn(paramValues(ind,ii),-4);
    end
    paramList(ii)
    selectIN.setData(paramList(ii),paramValue,ind) % last arg selects material
    selectIN.update();
end



for ii=1:numrealizs
%     later apply Latin hypercube sampling here
%     topfluxval = rand(1)*(max(topflux)-min(topflux))+min(topflux);
%     topperiodval = rand(1)*(max(topperiod)-min(topperiod))+min(topperiod);
%     soiltypeval = floor(0.5+rand(1)*(max(soiltype)-min(soiltype))+min(soiltype));

    topfluxval = 1;
    topperiodval = 30;
    soiltypeval = 4;

%   create ATMOSPH.IN object
%     atmosphIN = ATMOSPHIN(expDirectory);
    
%   updates the parameter info
%     atmosphIN.setData('Prec',[topfluxval]);
%     atmosphIN.setData('tAtm',topBCtimes)
%     atmosphIN.setData('Prec',-Prechold(1:infiltcycles*2))
%     atmosphIN.setData('rSoil',rsoilhold(1:infiltcycles*2))
%     atmosphIN.setData('hCritA',rsoilhold(1:infiltcycles*2)*100000)
%     atmosphIN.setData('MaxAL',infiltcycles*2)
    
%   sets up sinusoidal inputs
%     atmosphIN.setData('SinusVar','f')
    
%   commits/writes changes to the file (overwrites previous file)
%     atmosphIN.update()

    selectIN = SELECTORIN(expDirectory);
    selectIN.setData('B3',10.89388e+016);
    selectIN.setData('Alfa',0.054);
    selectIN.update();


%   set soil parameters
%     selectIN = SELECTORIN(expDirectory);
%     for j=1:length(paramList)
%         if j == 5
%             paramValue = roundn(paramValues(soiltypeval,j),-2); 
%         elseif j == 3
%             paramValue = roundn(paramValues(soiltypeval,j),-4);
%         else
%             paramValue = roundn(paramValues(soiltypeval,j),-4);
%         end
%         selectIN.setData(paramList(j),paramValue,1) % last arg selects material
%     end
    
%   set other selector.in parameters
%     selectIN.setData('TPrint(1),TPrint(2),...,TPrint(MPL)',[endtime/numprinttimes:endtime/numprinttimes:endtime])
%     selectIN.setData('MPL',numprinttimes)
%     selectIN.setData('tMax',endtime)
%   commits/writes changes to the file (overwrites previous file)
%     selectIN.update();        

%   runs the simulation(s)
%     Hydrus.run_hydrus(noCMDWindow)
%   moves the input & output files for saving
%     Hydrus.outputResults(exp,num2str(ii))

%   collect results 
%     nodinf = NODINF(expDirectory);
%     data = nodinf.getAllData();
%     qvals(ii,:,:)=data(:,:,5); % flux = index 5
% 
%     holdtopfluxval(ii)=topfluxval;
%     holdtopperiodval(ii)=topperiodval;
%     holdsoiltypeval(ii)=soiltypeval;

end;

% plot flux profiles in time
% for ii=1:numrealizs
%     toplotz=-1:-1:-size(qvals,3);    
%     toplotq=squeeze(qvals(ii,:,:));
%     plot(toplotq,toplotz);
%     hold on;
% end;



