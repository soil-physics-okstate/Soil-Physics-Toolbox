% STATLIBtutorial           - tutorial for the statlib directory (Jan 1,2001)

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

%%% Compute and plot the density scaled histograms for the
%%% sand, silt and clay contents. 

nbins=20;
bounds=[0 100];

figure;
subplot(2,2,1);
  histscaled(sand,nbins,bounds);
  title('Sand');
  hold on;
subplot(2,2,2);
  histscaled(silt,nbins,bounds);
  title('Silt');
  hold on;
subplot(2,2,3);
  histscaled(clay,nbins,bounds);
  title('Clay');
  hold on;

%%% Type any key for continuing...

pause;
clc;

%%% Compute the basic statistics for the sand, silt and clay
%%% contents, respectively.

m=mean([sand silt clay])
v=std([sand silt clay]).^2
s=skewness([sand silt clay])
k=kurtosis([sand silt clay])

%%% Type any key for continuing...

pause;
clc;

%%% Compute and plot the kernel
%%% density estimates evaluated from 0 to 100 % by step of 1%.
%%% The variance of the Gaussian kernel is equal to 20.

sandfileK=(-15:0.5:55)';
siltfileK=(25:0.5:100)';
clayfileK=(-10:1:55)';
v=20;

pdfsandfileK=kerneldensity(sand,sandfileK,v);
pdfsiltfileK=kerneldensity(silt,siltfileK,v);
pdfclayfileK=kerneldensity(clay,clayfileK,v);

subplot(2,2,1);
  plot(sandfileK,pdfsandfileK,'r');
subplot(2,2,2);
  plot(siltfileK,pdfsiltfileK,'r');
subplot(2,2,3);
  plot(clayfileK,pdfclayfileK,'r');

%%% Type any key for continuing...

pause;
clc;

%%% Estimate and plot the cumulative distribution function for
%%% the sand, silt and clay contents. Superimpose on the graph
%%% the cumulative distribution function obtained by integrating
%%% the kernel density estimates

[sandfileE,cdfsandfileE]=cdfest(sand);
[siltfileE,cdfsiltfileE]=cdfest(silt);
[clayfileE,cdfclayfileE]=cdfest(clay);

[cdfsandfileK]=pdf2cdf(sandfileK,pdfsandfileK);
[cdfsiltfileK]=pdf2cdf(siltfileK,pdfsiltfileK);
[cdfclayfileK]=pdf2cdf(clayfileK,pdfclayfileK);

figure;

subplot(2,2,1);
  plot(sandfileE,cdfsandfileE);
  hold on;
  plot(sandfileK,cdfsandfileK,'r');
  title('Sand');
  axis([0 100 0 1]);
subplot(2,2,2);
  plot(siltfileE,cdfsiltfileE);
  hold on;
  plot(siltfileK,cdfsiltfileK,'r');
  title('Silt');
  axis([0 100 0 1]);
subplot(2,2,3);
  plot(clayfileE,cdfclayfileE);
  hold on;
  plot(clayfileK,cdfclayfileK,'r');
  title('Clay');
  axis([0 100 0 1]);

%%% Type any key for continuing...

pause;
clc;

%%% Compute the directional variograms for the sand content along
%%% the N-S (stars) and the W-E (circles) directions, using an
%%% angular tolerance of 90 degrees. Distances classes are 500
%%% meters wide and are ranging from 0 to 10000 meters.

cl=(0:500:10000)';
figure;

optionsNS=[0,45,-45];
[dsandNS,vsandNS,osandNS]=vario(ch,sand,cl,'kron',optionsNS);
plot(dsandNS,vsandNS,'*');

hold on;

optionsWE=[0,-45,45];
[dsandWE,vsandWE,osandWE]=vario(ch,sand,cl,'kron',optionsWE);
plot(dsandWE,vsandWE,'o');

%%% Type any key for continuing...

pause;
clc;

%%% Compute the omnidirectional variograms for the sand content.
%%% Distances classes are the same than above.

cl=(0:500:10000)';
[dsand,vsand,osand]=vario(ch,sand,cl,'kron');
plot(dsand,vsand);

%%% Type any key for continuing...

pause;
clc;

%%% Plot a nested variogram models, composed of a nugget effect
%%% with a sill equal to 30, and an exponential model with a sill
%%% equal to 60 and a range equal to 3000.

modelsand={'nuggetV','exponentialV'};
paramsand0={[30],[60 3000]};
modelplot(dsand,modelsand,paramsand0);

%%% Type any key for continuing...

pause;
clc;

%%% Fit by a weighted least squares method the estimated
%%% variogram using a nested model that includes a nugget
%%% effect and an exponential model, and display the results.
%%% The initial values for the parameters are as above.

paramsand=modelfit(dsand,vsand,osand,modelsand,paramsand0);
modelplot(dsand,modelsand,paramsand,'Color',[1 0 0]);
paramsand{1}
paramsand{2}

%%% Type any key for continuing...

pause;
clc;

%%% Compute the whole set of omnidirectional variograms and
%%% cross variograms for the raw values of the sand, silt
%%% and clay contents. Distance classes are the same than above.

options=1;
figure;
[d,V,o]=vario(ch,[sand silt clay],cl,'kron',options);

%%% Type any key for continuing...

pause;
clc;

%%% Fit by a weighted least squares method the whole set of
%%% variograms and cross variograms, using a nested model that
%%% includes a nugget effect and an exponential model with range
%%% equal to 3000. Display the estimated sill coefficients for
%%% the nugget effect and the exponential models

model={'nuggetV','exponentialV'};
param0={[30],[30 3000]};
options=1;
[param]=coregfit(d,V,o,model,param0,options);

param{1}{1}
param{2}{1}

%%% Type any key for continuing...

pause;
clc;

%%% Transform the sand, silt and clay content values toward
%%% zero mean unit variance distributed values and display the
%%% histograms and scatter plots for the transformed values.

sandT=other2gauss(sand,sandfileE,cdfsandfileE);
siltT=other2gauss(silt,siltfileE,cdfsiltfileE);
clayT=other2gauss(clay,clayfileE,cdfclayfileE);

figure;
histscatterplot([sandT siltT clayT]);

%%% Type any key for continuing...

pause;
clc;

%%% Compute the whole set of omnidirectional variograms and
%%% cross variograms for the Gaussian transformed sand, silt
%%% and clay contents. Distance classes are the same than above.

options=1;
figure;
[dT,VT,oT]=vario(ch,[sandT siltT clayT],cl,'kron',options);

%%% Type any key for continuing...

pause;
clc;

%%% Fit by a weighted least squares method the whole set of
%%% variograms and cross variograms, using a nested model that
%%% includes a nugget effect and an exponential model with range
%%% equal to 3000. Display the estimated sill coefficients for
%%% the nugget effect and the exponential models

modelT={'nuggetV','exponentialV'};
paramT0={[0.5],[0.5 3000]};
options=1;
[paramT]=coregfit(dT,VT,oT,modelT,paramT0,options);

paramT{1}{1}
paramT{2}{1}

%%% Type any key for continuing...

pause;
clc;

%%% Save all the variables to the 'STATLIBtutorial.mat' binary data
%%% file for subsequent use

save STATLIBtutorial

%%% End of the STATLIB tutorial

echo off;




