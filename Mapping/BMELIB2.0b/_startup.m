% 
% Comment to install BMELIB:
%
% When you install BMELIB on your machine, edit this file so
% it correctly reflects the name of the directory in which you 
% installed the BMELIB packages
%

dirBME=pwd;
addpath([dirBME '\iolib'],'-end');
addpath([dirBME '\graphlib'],'-end');
addpath([dirBME '\modelslib'],'-end');
addpath([dirBME '\statlib'],'-end');
addpath([dirBME '\bmeprobalib'],'-end');
addpath([dirBME '\bmeintlib'],'-end');
addpath([dirBME '\bmehrlib'],'-end');
addpath([dirBME '\simulib'],'-end');
addpath([dirBME '\genlib'],'-end');
addpath([dirBME '\mvnlib'],'-end');
addpath([dirBME '\tutorlib'],'-end');
addpath([dirBME '\exlib'],'-end');
addpath([dirBME '\testslib'],'-end');
addpath([dirBME '\extensions\stmapping'],'-end');
addpath([dirBME '\extensions\projectionlib'],'-end');
addpath([dirBME '\bmecatlib'],'-end');

disp('Search path set for BMELIB version 2.0b on MATLAB 6');