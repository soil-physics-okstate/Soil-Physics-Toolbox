% processData.m
% 
%

clear all;
clc;

exp = 'Test';
% exp = 'Cyclic';
mainDirectory = 'C:\Temp\HYDRUS_Data\';
expDirectory = [mainDirectory 'Projects\' exp];

trialNums = [11,12,13,14,21,22];
% trialNums = [11];
fluxDepth = 50;
time = 10;

fluxData = zeros(length(trialNums),1);

for ii=1:length(trialNums)
    trialNum = trialNums(ii);
    expDirectory = [mainDirectory 'Results\' exp '\Trial= ' num2str(trialNum)]
    nodinf = NODINF(expDirectory);
    data = nodinf.getAllData();
    fluxData(ii) = data(time,fluxDepth,5);
end

% nodinf = NODINF(expDirectory);
obsnode = OBSNODE(expDirectory);

% nodinf.getWCData();
% nodinf.getHeadData();
% nodinf.getKData();
% nodinf.getFluxData();
% nodinf.getCapData();

% data = nodinf.getAllData();
data = obsnode.getObsData();
