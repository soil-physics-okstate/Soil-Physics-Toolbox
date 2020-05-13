function probaGenerationTest
% probaGenerationTest       - Test functions to generate probabilistic soft data
%
% Test the probaGaussian, probaUniform and probaStudentT functions.
%
% SYNTAX :
%
% probaGaussianTest
%

rand('state',0);
ns=2;

% Consider ns soft data points.  At each soft data point we have a set of
% observed values.  We want to know the soft probabilistic distribution of the 
% mean of the set of observed values at each data point.  To do that, we 
% calculate the arithmetic average zave and the estimated variance zvar
% at each data point.  Then we can genrerate probabilistic soft data in 
% each of three fashion:
% Using a Gaussian pdf with mean zvar and variance = zvar / nobvs, where
%    nobvs is the number of observation values at each data point,
% Using a Stundent-t pdf so that t = (z-zave) / sqrt(zvar/nobs) is 
%    student-t distributed with nobvs-1 degree of freedom
% Using a uniform pdf spaning 2 times the standard deviation sig=sqrt(zvar)
%    on each side of the arithmetic average zvar 

%  For each data point generate a set of observed values
nobvs=ceil(2+rand(ns,1)*5);
for is=1:ns
  Xi=rand(nobvs(is),1);
  zave(is,1)=mean(Xi);
  zvar(is,1)=var(Xi);
end
zlow=zave-2*sqrt(zvar./nobvs);
zup=zave+2*sqrt(zvar./nobvs);

% test probaGaussian
[softpdftype,nl,limi,probdens]=probaGaussian(zave,zvar./nobvs);
figure;
[h]=probaplot(softpdftype,nl,limi,probdens);
ax=axis;
title('Gaussian probabilistic soft data');
xlabel('x_s');
ylabel('f_S(x_s)');
ax=axis;

% test probaStudentT if the tinv function exists
if exist('tinv')==2
  [softpdftype,nl,limi,probdens]=probaStudentT(zave,zvar,nobvs);
  figure;
  [h]=probaplot(softpdftype,nl,limi,probdens);
  axis(ax);
  title('Student-T probabilistic soft data');
  xlabel('x_s');
  ylabel('f_S(x_s)');
  axis(ax);
end
  
% test probaUniform
[softpdftype,nl,limi,probdens]=probaUniform(zlow,zup);
figure;
[h]=probaplot(softpdftype,nl,limi,probdens);
title('Uniform probabilistic soft data');
xlabel('x_s');
ylabel('f_S(x_s)');
axis(ax);

disp('test complete');

