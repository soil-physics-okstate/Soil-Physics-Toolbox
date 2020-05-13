
% mvnlibcompile             - Compile the functions in the mvnlib directory (Jan 1, 2001)
%
% Several functions in the mvnlib directory are written in FORTRAN and 
% they need to be compiled before they can be used in BMELIB. This 
% program will use the MATLAB mex compiler to compile all the necessary  
% functions in the mvnlib directory. In order for this program to work,
% the MATLAB mex compiler needs to be working with a FORTRAN
% compiler. If you have problems running this program, check
% that the mex compiler is working correctly with the FORTRAN
% compiler on your platform (by for example running the MATLAB
% tutorial for mex)
%
% IMPORTANT NOTE: You need to first change to the mvnlib directory
% before running this program!
%
% SYNTAX :
%
% (first make sure you are in the mvnlib directory, then type:)
% mvnlibcompile;
%
% NOTES :
%
% 1- Once you have used this program, you can check that the 
% FORTRAN functions in the mvnlib directory were properly compiled 
% by running the MVNLIBtest function (see the testslib directory).
%
% 2- It is worth reenforcing upon you the fact that you need to 
% first change to the mvnlib directory before running this program
%
% 3- Edit this file if you want to change the compilation options.
% It may be possible to optimize further the mvnlib functions by compiling
% the dcuhre.f file using a FORTRAN compiler for your platform, and then linking
% the object file created (say dcuhre.obj), with the different mex functions.
% This two step approach allows to optimize aggressively dcuhre.f, where
% most of the computation happens (using for example a -O4 switch). 
% Then to link the object file dcuhre.obj use a mex command similar to:
%    mex -O mvMomVecAG2.f mvMomVecAG2g.f dcuhre.obj symInv.f
%                                               ^^^

disp('Compiling the FORTRAN functions in the mvnlib directory');

disp('Compiling uvProNR');
mex -O uvProNR.f uvProNRg.f qromb3.f trapzd3.f polint3.f
disp('Compiling uvMomNR');
mex -O uvMomNR.f uvMomNRg.f qromb3.f trapzd3.f polint3.f
disp('Compiling mvMomVecAG2');
mex -O mvMomVecAG2.f mvMomVecAG2g.f dcuhre.f symInv.f
disp('Compiling mvProAG2');
mex -O mvProAG2.f mvProAG2g.f dcuhre.f symInv.f
disp('Compiling mvnAG1');
mex -O mvnPack.f mvnAG1g.f -output mvnAG1

disp('Compilation finished. Use MVNLIBtest to test the newly compiled functions');
