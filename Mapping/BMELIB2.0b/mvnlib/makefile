#------------------------------------------------------------------
#------------------------------------------------------------------
#  Set compiler and option        
#


#OPTS          = -g +z
OPTS         = +O2 -DNDEBUG +z
F77          =  f77


# Note: the +z option is needed to make object files that can be 
# used with matlab mex compiler. 
#------------------------------------------------------------------

#  General rules
#

.f.o:
	$(F77) $(OPTS) -c $*.f

#------------------------------------------------------------------
# mvMomVecAG2.mexhp7 
# in unix: make dcuhre.o and symInv.o 
# in matlab: mex mvMomVecAG2.f mvMomVecAG2g.f dcuhre.o symInv.o

#
# mvMomVecAG2Test.f: tests the subroutine in mvMomVecAG2.f
#
mvMomVecAG2Test: mvMomVecAG2Test.o mvMomVecAG2.o dcuhre.o symInv.o
	$(F77) $(OPTS) mvMomVecAG2Test.o mvMomVecAG2.o dcuhre.o symInv.o -o mvMomVecAG2Test -lm

mvMomVecAG2Test: mvMomVecAG2Test.f
mvMomVecAG2.o: mvMomVecAG2.f

#------------------------------------------------------------------
# mvProAG2.mexhp7 
# in unix: make dcuhre.o and symInv.o 
# in matlab: mex mvProAG2.f mvProAG2g.f dcuhre.o symInv.o

#
# mvProAG2Test.f: tests the subroutine in mvProAG2.f
#
mvProAG2Test: mvProAG2Test.o mvProAG2.o dcuhre.o symInv.o
	$(F77) $(OPTS) mvProAG2Test.o mvProAG2.o dcuhre.o symInv.o -o mvProAG2Test -lm

mvProAG2Test: mvProAG2Test.f
mvProAG2.o: mvProAG2.f


#------------------------------------------------------------------
# dcuhreTest.f: tests the subroutine in dcuhre.h program 
#               for multiple integrals
#
dcuhreTest: dcuhre.o dcuhreTest.o
	$(F77) $(OPTS) dcuhre.o dcuhreTest.o -o dcuhreTest

#------------------------------------------------------------------
# mvnAG1.mexhp7
# in unix: make mvnPack.o
# in matlab: mex mvnPack.o mvnAG1g.f -output 'mvnAG1.mexhp7' 

#  
# mvnAG1Test.f: tests the subroutine in mvnAG1.f
#

mvnAG1Test: mvnAG1Test.o mvnAG1.o
	$(F77) $(OPTS) mvnAG1Test.o mvnPack.o -o mvnAG1Test

#
# cleaning up stuff
#

clean:
	rm -f *.o *~



