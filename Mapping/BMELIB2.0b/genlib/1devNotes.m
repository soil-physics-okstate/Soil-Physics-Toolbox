

*****************************************
**** Notes for BMElib developers
****
**** This file contains a log of the modifications
**** for each version of BMElib
****
*****************************************
  

******   Creating BMELIB zip  ******************
1) remove the extension directory
2) remove *.mat and *.eps from exlib
3) remove *.mat and *.dat from tutorlib
4) remove WS_FTP.LOG from all directories
6) remove *.dat, *.out, *.par, *.txt from testslib


******   Version 1.0b  ******************

9/18/01 - Updated files: (version 1.0b)

AK provided the compiled mex files for MATLAB6-linux86, 
as well as a fix for compiling FORTRAN mex on linux

AK provided the compiled mex files for MATLAB6-linux86, 
as well as a fix for compiling FORTRAN mex on linux

JK provided the fix to make uvProNR.f and uvMomNR.f
compile and work on SGI.  Fix involved correcting the following
line in qromb3.f:

OLD line does not work on SGI:

call polint3(h(j-KM),s(j-KM),K,0.,ss,dss)

NEW line does work on SGI:

call polint3(h(j-KM),s(j-KM),K,0.0d0,ss,dss)

9/20/01 - Updated files: (version 1.0b)

tutorlib: (apparently this did not get done in version 1.0)
STATLIBtutorial  :  histscatterplot.m name change
GRAPHLIBtutorial :  histscatterplot.m name change

graphlib:
removed files scatterplot.m (now using histscatterplot.m)


******   Version 1.0c  ******************

7/16/03 - Corrections for Matlab v.6.5

Correction of issue with prod function: Does not accept arguments of type logical

C:\package\amendedBMELIB2.0_MATLAB6.5\bmeprobalib\softpdftypeCheckArgs.m
             l.63: if ~prod(double(sl==[ns nlMax])), error('limi must be a ns by max(nl) matrix'); end;
        WAS: l.59: if ~prod(sl==[ns nlMax]), error('limi must be a ns by max(nl) matrix'); end;
	          l.65: if ~prod(double(sl==[ns 3])), error('limi must be a ns by 3 matrix'); end;
        WAS: l.61: if ~prod(sl==[ns 3]), error('limi must be a ns by 3 matrix'); end;
             l.76: if ~prod(double(sp==[ns nlMax-1])), error('probadens must be a ns by max(nl)-1 matrix'); end;
        WAS: l.72: if ~prod(sp==[ns nlMax-1]), error('probadens must be a ns by max(nl)-1 matrix'); end;
             l.78: if ~prod(double(sp==[ns nlMax])), error('probadens must be a ns by max(nl) matrix'); end;
        WAS: l.74: if ~prod(sp==[ns nlMax]), error('probadens must be a ns by max(nl) matrix'); end;

C:\package\amendedBMELIB2.0_MATLAB6.5\genlib\coord2K.m
	     l.126:  issymm=prod(prod(double(c1{1}==c2{1})));    % issymm to 0 or 1 
	WAS: l.122:  issymm=prod(prod(c1{1}==c2{1}));    % issymm to 0 or 1 

C:\package\amendedBMELIB2.0_MATLAB6.5\iolib\writeBMEproba.m
	     l.52:  if ~prod(double(size(fileNames)==[6 1])), error('fileNames must be a 6 by 1 cell'); end;
	WAS: l.48:  if ~prod(size(fileNames)==[6 1]), error('fileNames must be a 6 by 1 cell'); end;

7/16/03 - Corrections for Minor bug fixes

C:\package\amendedBMELIB2.0_MATLAB6.5\simulib\simuinterval.m
	l.34: V = rand(ns,nSim)-0.5).*I;
     WAS: V = rand(ns,nSim)-0.5)*I;

C:\package\amendedBMELIB2.0_MATLAB6.5\statlib\crosscovario.m
	l.99: if (size(c1,2)~=2 & size(c2,2)~=2),
     WAS: if size(c,2) ~= 2,

C:\package\amendedBMELIB2.0_MATLAB6.5\statlib\crossvario.m
    l.107: if (size(c1,2)~=2 & size(c2,2)~=2),
    WAS: if size(c,2) ~= 2,

      
******   Version 2.0  ******************

10/8/01
Added the parameter zvalue in colorplot.m and markerplot.m

1/25/02
Changed writing format in writeGeoEAS (from %12g to %12.10g)

3/6/02
Corrected bug in BMEprobaFromFile when using BMEmethod='BMEprobaMoments':
Changed line 113 from
  val(:,d+isST+1:d+isST+3)=BMEmoments;
to
  val(:,d+isST+1:d+isST+3)=BMEest;

Also changed the IOLIBtest file so it can now test for the two following cases 
%   testtype=5 to test BMEprobaMoments
%   testtype=6 same as testtype=5, but WITH soft data at estimation pt

3/7/02
Added the function pairsindex.m in statlib which implements methods 'kron', 'kronloop' and 
'superblock' to find pairs of points corresponding to given intervals of distance.

Modified stcov.m so it now uses pairsindex.m

Changed Contents.m of statlib to add pairsindex.m

Added pairsindextest.m in testslib and changed Contents.m of testslib.

Added avedupli.m in genlib to average duplicate values, and added this function in Contents.m of genlib


4/8/02
Added directory extensions\stmapping with the following programs:
estZsBME.m
plotZsBME.m
testestZsBME.m
testplotZsBME.m

Added the following path in startup.m:
addpath C:\package\BMELIB2.0_MATLAB6.0\extensions\stmapping -end

put stmeaninterpstv.m in statlib

added stmeaninterpstv.m in Contents.m of statlib directory 

4/31/02
Added jet2.m in graphlib and modified accordingly contents.m in graphlib

Fixed bug in creategrid (changed min(c1) min(c2) and min(c3) with 
minc(1) minc(2) and minc(3)

5/12/02
Corrected krigingIndGauss.m (from Patrick) 
Le bug concernait le cas ou deux valeurs sont fournies dans la variable p
(calcul pour la probabilite d'un intervalle) au lieu d'une seule (calcul
d'une probabilité cumulée)

7/15/02
Changed test_pairsindex to pairsindextest in testslib.

7/17/02
Added probaUniform.m, probaGaussian.m and probaStudentT.m in bmeprobalib for 
generation of probabilistic soft data, and modified Content.m accordingly

Added probaGenerationTest.min in testslib to test probaUniform.m, probaGaussian.m 
and probaStudentT.m, and modified Content.m accordingly


6/20/03
Fixed bug in BMEprobaCI when only one confidence level is requested. Change
the following
  for iCI=1:nCI,
    zlCI(ik,iCI)=min(zlCIik{iCI});
    zuCI(ik,iCI)=max(zuCIik{iCI});
  end;
to this
  if nCI==1, 
    zlCI(ik,1)=min(zlCIik);
    zuCI(ik,1)=max(zuCIik);
  else
    for iCI=1:nCI,
      zlCI(ik,iCI)=min(zlCIik{iCI});
      zuCI(ik,iCI)=max(zuCIik{iCI});
    end;
  end;

8/12/2003
Added functions isdupli.m and finddupli.m in genlib, and
updated Contents.m
  
9/2/2003
Changed function avedupli.m so it works in a s/t space of any dimension p.
Added avedupli2D.m which averages duplicates in a 2D space
updated Contents.m

Modified valstv2stg.m so that it can handle case where cMS and tME are not known


******   Version 2.0a  ******************

9/10/2003:
modified markerplot.m and colorplot.m so they can handle NaNs in z


******   Version 2.0b  ******************

1/20/2004: 
created a STATLIBtest.

1/25/2004: 
Changed fmin to fminbnd in the code of modelfit and Triangularinv

1/25/2004: 
Changed fmin to fminbnd in the text of BMEprobaMode, BMEprobaTMode, BMEIntervalMode, 
BMEIntervalTMode, Ioparamsyntax, WriteBMEproba, BMEoptions, BMEprobaCI, BMEprobaMode
BMEprobaMoments, BMEprobaPdf, BMEprobaTCI, BMEprobaTMode, BMEprobaTPdf, BMEprobalib/contents
BMEstatlib/contents, BMEIntervalMode, BMEIntervalTMode, BMEintlib/Contents, FminBMEIntervalMode
FminBMEIntervalTModeBME, FminBMEIntervalTModeSK

1/25/2004: 
Changed the normalization at the end of probaGaussian and
probaStudentT so it is quicker.

1/25/2004: 
Error in  BMEPROBALIBtutorial.m  using colorplot.
Line 296 replace 
colorplot(ch((th==2),:),zh(th==2),'hot',Property,Value);
with 
colorplot(ch((th==2),1:2),zh(th==2),'hot',Property,Value);

Replace line 321 with
colorplot(ch((th==3),1:2),zh(th==3),'hot',Property,Value);

Replace line 346 with
colorplot(ch((th==4),1:2),zh(th==4),'hot',Property,Value);

Replace line 301 with
markerplot(cs((ts==2),1:2),ms(ts==2),[5 20],Property,Value);

Replace line 326 with
markerplot(cs((ts==3),1:2),ms(ts==3),[5 20],Property,Value);

Replace line 351 with
markerplot(cs((ts==4),1:2),ms(ts==4),[5 20],Property,Value);

1/28/2004
Added quantest, and updated stalib/Contents.m accordingly

1/30/2004
fixed bug in histline and histscaled by replace the following line:
n(nbins)=n(nbins)+sum(w(z>class(nbins)));
with
n(nbins)=n(nbins)+sum(w(z>class(nbins+1)));

2/3/2004
Modified proba2val so it uses interp1 for the case of softpdftype=2
Modified probaStudentT so it uses tinv to calculate xnorm
Added combinedupli in genlib and modified Contents accordingly
Added probacombinedupli in bmeproblib and modified Contents accordingly
Added combinedupliTest in BMEPROBALIBtest and modified Contents accordingly
Modified BMEprobaCheckArgs so it uses isdupli to check for duplicate data points

4/15/2004
Added the bmecatlib directory provided by Patrick Bogaert

4/30/2004
Added BMEprobaMomentsXvalidation.m in and updated bmeprobalib/Contents.m accordingly

