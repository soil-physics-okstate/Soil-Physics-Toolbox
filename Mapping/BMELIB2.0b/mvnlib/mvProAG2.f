      SUBROUTINE mvProAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C,
     & maxpts, aEps, rEps, key, Value, Error, ifail)
*
*     A subroutine for computing the nDim multiple integral
*      b
*     I dXs fs(Xs) mvnpdf(Xs,Mean,Cov),
*      a                 
*     where a and b define the lower and upper limit of a
*     hyperrectangle, a=Limi(:,1), b=Limi(:,lastColumn), NL a vector
*     of length nDim, and the function fs is defined by softPdfType 
*     and the matrices Limi and Prob. 
*
* From:
* "An Adaptive Multidimensional Integration Routine for a Vector of Integrals"
*          Alan Genz 
*          Department of Mathematics
*          Washington State University 
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
*  Input Parameters
*
*     nDim    INTEGER, the number of dimension.
*     softPdfType: Type of pdf soft data; 1-Histogram, 2-linear, 
*                  3-Grid histogram, 4-Grid Linear
*     Nl      INTEGER, array with number of intervals for each dimension.
*     Limi    DOUBLE, array of interval limits in each dimension
*     Prob    DOUBLE, array of proba for each interval.
*     Mean    DOUBLE, array of mean in each dimension
*     C   DOUBLE, array for the covariance matrix C, in lower diag by colum 
*            The coefficient in row I column J of the covariance matrix
*            C should be stored in C( J + (I*(I-1))/2 ), for J <= I.
*     maxpts INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with maxpts = 1000*nDim, and then
*            increase maxpts if ERROR is too large.
*     aEps   DOUBLE absolute error tolerance.
*     rEps   DOUBLE relative error tolerance.
C     key    Integer to control integration rule by dcuhre program
C            key to selected local integration rule.
C            key = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            key = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            key = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            key = 3 gives the user the degree 9 integration rule.
C            key = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
*
*  Output Parameters
*
*     Value  DOUBLE, estimated value for the integral
*     Error  DOUBLE, 
*     ifail INTEGER, termination status parameter:
C            ifail = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with maxpts or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            ifail = 1 if maxpts was too small for DCUHRE
C              to obtain the required accuracy. In this case DCUHRE
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            ifail = 2 if key is less than 0 or key greater than 4.
C            ifail = 3 if NDIM is less than 2 or NDIM greater than 15.
C            ifail = 4 if key = 1 and NDIM not equal to 2.
C            ifail = 5 if key = 2 and NDIM not equal to 3.
C            ifail = 6 if NUMFUN is less than 1.
C            ifail = 7 if volume of region of integration is zero.
C            ifail = 8 if maxpts is less than 3*NUM.
C            ifail = 9 if maxpts is less than MINPTS.
C            ifail = 10 if EPSABS < 0 and EPSREL < 0.
C            ifail = 11 if NW is too small.
C            ifail = 12 if softPdfType not equal to 1, 2, 3, or 4
C            ifail = 13 if max(Nl) is greater than 101.
C            ifail = 14 if C is not spd.
*
      EXTERNAL mvf
C
      INTEGER nDim, softPdfType, Nl(*), maxpts, key, ifail
      INTEGER minpts, nFun, NWORK
      PARAMETER (NWORK=100000,MAXDIM=20)
      DOUBLE PRECISION Limi(*), Prob(*), C(*), Mean(*), aEps, rEps, 
     &     Value(*), Error(*)
      DOUBLE PRECISION Work(NWORK)
      DOUBLE PRECISION A(MAXDIM), B(MAXDIM)
C
C
      ifail=0
      if ((softPdfType.EQ.1).OR.(softPdfType.EQ.2)) then
         do is=1,nDim
            A(is)=Limi(is)
            B(is)=Limi(is+(nl(is)-1)*nDim)
         enddo
      elseif ((softPdfType.EQ.3).OR.(softPdfType.EQ.4)) then
         do is=1,nDim
            A(is)=Limi(is)
            B(is)=Limi(is+(3-1)*nDim)
         enddo
      else
C         print '(''Error in mvProAG2: Unauthorized softPdfType'')'
          ifail=13
          RETURN
      endif
C
      nFun=1
      minpts=0
C
      CALL mvfInit(nDim,softPdfType,Nl,Limi,Prob,C,Mean,ifail)
      if (ifail.GT.0) then
         Value(1)=0
         Error(1)=0
         RETURN
      endif	
C      
      CALL DCUHRE(nDim, nFun, A, B, minpts, maxpts, mvf, aEps,
     *  rEps, key, NWORK, 0, Value, Error, nEval, ifail, Work)
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         mvf
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE mvf(nDim, X, nFun, F)
C
C     mulvariate function
C
      INTEGER MAXDIM, MAXNL 
      PARAMETER ( MAXDIM = 20, MAXNL=101)
      INTEGER nDim, nFun, ifail
      INTEGER softPdfType, Nl(*)
      INTEGER nlMax, i, j, ij, is, il, nullty, ifault, isil
      DOUBLE PRECISION X(*), F(*)
      DOUBLE PRECISION Limi(*), Prob(*), C(*), Mean(*)
      DOUBLE PRECISION sum, detC, pi, Work(2*MAXDIM), fs
      DOUBLE PRECISION EPSI
      PARAMETER ( EPSI = 1e-10 )
      INTEGER softPdfTypeS, NlS(MAXDIM) 
      DOUBLE PRECISION LimiS(MAXDIM*MAXNL), 
     &  ProbS(MAXDIM*MAXNL), MeanS(MAXDIM),  
     &  COVinv((MAXDIM*(MAXDIM+1))/2), NormVal,
     &  fsoS(MAXDIM*MAXNL), fspS(MAXDIM*MAXNL)
      SAVE softPdfTypeS, NlS
      SAVE LimiS, ProbS, fsoS, fspS, MeanS, COVinv, NormVal
C
      sum=0
      ij=0
      do i=1,nDim
         do j=1,i-1
            ij=ij+1
            sum=sum+2.0*(X(i)-MeanS(i))*COVinv(ij)*(X(j)-MeanS(j))
         enddo
         ij=ij+1
         sum=sum+(X(i)-MeanS(i))*COVinv(ij)*(X(i)-MeanS(i))
      enddo
      F(1)=normVal*exp(-sum/2)
C
      fs=1.0

C Case softPdfTypeS=4: Grid linear pdf
      if (softPdfTypeS.EQ.4) then
         do is=1,nDim
            il=(X(is)+EPSI-LimiS(is))/LimiS(is+nDim)+1
            isil=is+(il-1)*nDim
            fs=fs*(fsoS(isil)+fspS(isil)*X(is))
         enddo
C Case softPdfTypeS=3: Grid histogram pdf
      elseif (softPdfTypeS.EQ.3) then
         do is=1,nDim
            il=(X(is)+EPSI-LimiS(is))/LimiS(is+nDim)+1
            fs=fs*ProbS(is+(il-1)*nDim)
         enddo
C Case softPdfTypeS=2: Linear pdf
      elseif (softPdfTypeS.EQ.2) then
         do is=1,nDim
            il=1
         do while ((LimiS(is+il*nDim).LT.X(is)).AND.(il+1.LT.NlS(is)))
               il=il+1
            enddo
            isil=is+(il-1)*nDim
            fs=fs*(fsoS(isil)+fspS(isil)*X(is))
         enddo
C Case softPdfTypeS=1: Histogram pdf
      elseif (softPdfTypeS.EQ.1) then
        do is=1,nDim
          il=1
          do while ((LimiS(is+il*nDim).LT.X(is)).AND.(il+1.LT.NlS(is)))
             il=il+1
          enddo
          fs=fs*ProbS(is+(il-1)*nDim)
        enddo
      endif
C
      F(1)=fs*F(1)
      RETURN
C
C  entry point for initialization of mvnPdf
C                  
      ENTRY mvfInit(nDim,softPdfType,Nl,Limi,Prob,C,Mean,ifail)
C
C  Check that LimiS and ProbS arrays are big enough
C
      nlMax=-1
      do is=1,nDim
         if (nlMax.LT.Nl(is)) then
            nlMax=Nl(is)
         endif
      enddo
      if (nlMax.GT.MAXNL) then
C         print '(''mvProAG2.f error: max(Nl) can't be greater than 101'')'
          ifail=13
          RETURN
      endif
C
C  save arguments:     softPdfType,  Nl,  Limi,  Prob,  Mean
C  in local variables: softPdfTypeS, NlS, LimiS, ProbS, MeanS,
C  and calculate fsoS and fspS

      softPdfTypeS=softPdfType
C Case softPdfType=4: Grid linear pdf
      if (softPdfType.EQ.4) then
         do is=1,nDim
            NlS(is)=Nl(is)
            MeanS(is)=Mean(is)
            do il=1,3
               isil=is+(il-1)*nDim
               LimiS(isil)=Limi(isil)
            enddo
            do il=1,Nl(is)-1
               isil=is+(il-1)*nDim
               fspS(isil)=(Prob(isil+nDim)-Prob(isil))/Limi(is+nDim)
               fsoS(isil)=Prob(isil)
     &                   -(Limi(is)+(il-1)*Limi(is+nDim))*fspS(isil)
            enddo
            il=Nl(is)
            isil=is+(il-1)*nDim
            fsoS(isil)=Prob(isil)
            fspS(isil)=0
         enddo
C Case softPdfType=3: Grid histogram pdf
      elseif (softPdfType.EQ.3) then
         do is=1,nDim
            NlS(is)=Nl(is)
            MeanS(is)=Mean(is)
            do il=1,3
               isil=is+(il-1)*nDim
               LimiS(isil)=Limi(isil)
            enddo
            do il=1,Nl(is)-1
               isil=is+(il-1)*nDim
               ProbS(isil)=Prob(isil)
            enddo
         enddo
C Case softPdfType=2: Linear pdf
      elseif (softPdfType.EQ.2) then
         do is=1,nDim
            NlS(is)=Nl(is)
            MeanS(is)=Mean(is)
            do il=1,Nl(is)-1
               isil=is+(il-1)*nDim
               LimiS(isil)=Limi(isil)
               fspS(isil)=   (Prob(isil+nDim)-Prob(isil))
     &                     / (Limi(isil+nDim)-Limi(isil))
               fsoS(isil)=Prob(isil)-Limi(isil)*fspS(isil)
            enddo
            il=Nl(is)
            isil=is+(il-1)*nDim
            LimiS(isil)=Limi(isil)
            fsoS(isil)=Prob(isil)
            fspS(isil)=0
         enddo
C Case softPdfType=1: Histogram pdf
      elseif (softPdfType.EQ.1) then
         do is=1,nDim
            NlS(is)=Nl(is)
            MeanS(is)=Mean(is)
            do il=1,Nl(is)-1
               isil=is+(il-1)*nDim
               LimiS(isil)=Limi(isil)
               ProbS(isil)=Prob(isil)
            enddo
            il=Nl(is)
            isil=is+(il-1)*nDim
            LimiS(isil)=Limi(isil)
         enddo
      endif
C
C  get the inverse COVinv and the determinant detC of the matrix C
C
      nn=(nDim*(nDim+1))/2
      call symInv(C,nDim,(nDim*(nDim+1))/2,Work,COVinv,detC,
     &            nullty,ifault)
      if (ifault.EQ.2) then
C        print '(''Error in mvProAG2: C is not spd'')'
         ifail=14
      endif
*     normVal = normalization constant = 1/(2*pi)^(nDim/2) * 1/sqrt(det(C))
      pi=2.0*acos(0.0)
      normVal=1.0/( (2.0*pi)**(nDim/2.0) * sqrt(detC) )
      END


