function [mesoNum MP VWC Se HC]=soilDrainage(mesonetfile,rosettafile,outputfile,siteID)
%Reads in Mesonet Data and converts to soil drainage
%using van Genuchten soil water retention parameters
%*************************************************************************%

% K:\Mesonet Data\vgparameters.xlsx  ==> use for rosettafile

%Read in ROSETTA Data
[rosettaNum,rosettaTxt,rosettaRaw]=xlsread(rosettafile,2); 
%K:\Mesonet Data\vgparameters.xlsx, sheet2

%*************************************************************************%
%find correct ROSETTA data
depth='55';                         %depth 55-65 cm for 60 cm HC
searchParam=strcat(siteID,depth);   %take 4 letter siteID, add depth

[row,col]=size(rosettaTxt);
found=false;
col=2;
while(~found&&row>0)
    if(strcmp(rosettaTxt(row,col),searchParam))
        found=true;
    else
        row=row-1;
    end
end

%*************************************************************************%

%Grab ROSETTA parameters from vgparameters file at specified excel row
%Otherwise, set all parameters to NaN to indicate error 

%Account for the fact that we searched rosettaTxt, which has 3 lines of 
%headers that rosettaNum does not.
row=row-3;

if(~(row==0))
    ThetaR=rosettaNum(row,9);    %cm^3/cm^3
    ThetaS=rosettaNum(row,10);   %cm^3/cm^3
    a=rosettaNum(row,11);        %1/kPa
    n=rosettaNum(row,12);        %unitless
    Ko=rosettaNum(row,14);       %cm/day
    L=rosettaNum(row,15);        %unitless
else
    ThetaR=nan;                 
    ThetaS=nan;   
    a=nan;        
    n=nan;        
    Ko=nan;       
    L=nan;
end

%Write ROSETTA parameters to output file, sheet 2 for QA
sheet=2;
xlswrite(outputfile,{'ThetaR','ThetaS','Alpha','n','Ko','L'},sheet,'a1:f1');
xlswrite(outputfile,[ThetaR,ThetaS,a,n,Ko,L],sheet,'a2:f2');


%*************************************************************************%
%Read in Mesonet *.csv file and create new excel file
%row1,col5 if downloading Temperature Reference AND Rainfall data
[mesoNum,mesoTxt,mesoRaw] = xlsread(mesonetfile);

%*************************************************************************%

%Read in Mesonet Temperature Reference values at 60cm and convert 
%to hydraulic conductivity
%when reading in mesonet csv data, TR60 is at column 8
c=8;
%TR60=mesoNum(:,c);
TR60=mesoNum(:,c);
    
%find Mesonet error messages (all negative) and convert to NaN
i=TR60<0;
TR60(i)=nan;
    
%Convert Mesonet Temperature Reference Value to Matric Potential (MP)
MP=-.717*exp(1.788*TR60);%MP=kPa
%c=1.788 degC^-1
%a=.717kPa
        
%Convert Matric Potential to Volumetric Water Content (VWC)
VWC=ThetaR+((ThetaS-ThetaR)./((1+(-a*MP).^n).^(1-1/n))); %cm^3/cm^3

%Convert Volumetric Water Content to Effective Saturation (Se)
Se=(VWC-ThetaR)/(ThetaS-ThetaR); %cm^3/cm^3
  
%Convert Effective Saturation to Hydraulic Conductivity (HC)
HC=Ko*Se.^L.*(1-(1-Se.^(n/(n-1))).^(1-1/n)).^2; %cm/d

%*************************************************************************%

%Write headers into new excel file
xlswrite(outputfile,{'Year','Month','Day','SiteID','RAIN','TR05','TR25','TR60','60cm MP','60cm VWC','60cm Se','60cm HC'},'a1:l1');
xlswrite(outputfile,{'    ','     ','   ','      ',' in ','degC','degC','degC','kPa    ','cm3/cm3 ','cm3/cm3','cm/d   '},'a2:l2');

%How long is new input data?
[row,col]=size(mesoNum);
%create range to input generated data, accounting for headers
row=row+2;
range='a3:l';
row=num2str(row);
range=strcat(range,row);

%Write generated data into output file
xlswrite(outputfile,[mesoNum MP VWC Se HC],range);

%Write site identifier into output file
nameRange='d3:d';
nameRange=strcat(nameRange,row);
xlswrite(outputfile,{siteID},nameRange);

%Working original raw data output writer
%xlswrite(outputfile,[mesoNum MP VWC Se HC]);

%Call separate summary function to get a sum of HC for each year
S=summarizeDrainage(HC,mesoNum);
[rs,cs]=size(S);
xlswrite(outputfile,{'Year','HC Sum','Missing'},3,'a1:c1');
xlswrite(outputfile,{'    ','cm/yr ','Data   '},3,'a2:c2');
rangeS='a3:c';
rs=num2str(rs+2);
rangeS=strcat(rangeS,rs);
xlswrite(outputfile,[S],3,rangeS);
