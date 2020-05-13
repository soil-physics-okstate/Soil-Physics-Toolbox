      SUBROUTINE uvMomNR(softPdfType, Nl, Limi, Prob, Mean, C, As,bs,p,
     & aEps, rEps, Value, Error, ifail)
*
*     A subroutine for computing the integral
*      b
*     I dXs fs(Xs) * (As*Xs+bs)**p * mvnpdf(Xs,Mean,Cov),
*      a                 
*     where a and b define the lower and upper limit of integration
*     a=Limi(1), b=Limi(lastColumn), Nl is the number of limits, LIMI
*     is a vector of lenght Nl and the function fs is defined 
*     by softPdfType and the vectors Limi and Prob. 
*
* From:
* "Numerical Recipies in FORTRAN"
*
*  Input Parameters
*
*     softPdfType: Type of pdf soft data; 1-Histogram, 2-linear, 
*                  3-Grid histogram, 4-Grid Linear
*     Nl     INTEGER, number of intervals.
*     Limi   DOUBLE, vector of interval limits
*     Prob   DOUBLE, vector of proba for each interval.
*     Mean   DOUBLE, mean
*     C      DOUBLE, covariance value
*     As     DOUBLE, linear coeffient
*     bs     DOUBLE, constant
*     p      DOUBLE, power coefficient
*     aEps   DOUBLE absolute error tolerance.
*     rEps   DOUBLE relative error tolerance.
*
*  Output Parameters
*
*     Value  DOUBLE, estimated value for the integral
*     Error  DOUBLE, 
C     ifail  Integer.
C            ifail = 0 for normal exit, when Error <=  aEps or
C              Error <=  abs(Value)*rEps 
C            ifail = 1 too many iterations necessary to obtain 
C                         the required accuracy.
C            ifail = 2 if softPdfType is not equal to 1, 2, 3, or 4
*
      EXTERNAL uvf
C
      INTEGER softPdfType, Nl, ifail
      DOUBLE PRECISION Limi(*), Prob(*), C, Mean, As, bs, p, aEps,rEps, 
     &     Value, Error
      DOUBLE PRECISION a, b
C
C
      if ((softPdfType.EQ.1).OR.(softPdfType.EQ.2)) then
        a=Limi(1)
        b=Limi(Nl)
      elseif ((softPdfType.EQ.3).OR.(softPdfType.EQ.4)) then
        a=Limi(1)
        b=Limi(3)
      else
        ifail=2
	  return
      endif
C
C
      dummy=uvfInit(softPdfType,Nl,Limi,Prob,C,Mean,As,bs,p)
C      
      call qromb3(uvf,a,b,aEps,rEps,Value,Error,ifail)
C      call qromb3(uvf,a,b,Value)
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C       FUNCTION uvf (uNIvARIATE fUNCTION)
C       
      FUNCTION uvf(X)
C
C     univariate function
C
      DOUBLE PRECISION uvf
      DOUBLE PRECISION uvfInit
      INTEGER MAXNL 
      PARAMETER ( MAXNL=101)
      INTEGER nFun
      INTEGER softPdfType, Nl
      INTEGER il 
      DOUBLE PRECISION X
      DOUBLE PRECISION Limi(*), Prob(*), C, Mean, As, bs, p
      DOUBLE PRECISION sum, detC, pi, fs
      DOUBLE PRECISION EPSI
      PARAMETER ( EPSI = 1e-10 )
      INTEGER softPdfTypeS, NlS 
      DOUBLE PRECISION LimiS(MAXNL), 
     &  ProbS(MAXNL), MeanS, AsS, bsS, pS, 
     &  COVinv, NormVal,
     &  fsoS(MAXNL), fspS(MAXNL)
      SAVE softPdfTypeS, NlS
      SAVE LimiS, ProbS, fsoS, fspS, MeanS, COVinv, AsS,bsS,pS, NormVal
C
      uvf=normVal*exp(-(X-MeanS)*COVinv*(X-MeanS)/2)
C
C Case softPdfTypeS=4: Grid linear pdf
      if (softPdfTypeS.EQ.4) then
        il=(X+EPSI-LimiS(1))/LimiS(2)+1
        fs=fsoS(il)+fspS(il)*X
C Case softPdfTypeS=3: Grid histogram pdf
      elseif (softPdfTypeS.EQ.3) then
        il=(X+EPSI-LimiS(1))/LimiS(2)+1
        fs=ProbS(il)
C Case softPdfTypeS=2: Linear pdf
      elseif (softPdfTypeS.EQ.2) then
        il=1
        do while ((LimiS(1+il).LT.X).AND.(il+1.LT.NlS))
          il=il+1
        enddo
        fs=fsoS(il)+fspS(il)*X
C Case softPdfTypeS=1: Histogram pdf
      elseif (softPdfTypeS.EQ.1) then
        il=1
        do while ((LimiS(1+il).LT.X).AND.(il+1.LT.NlS))
          il=il+1
        enddo
        fs=ProbS(il)
      endif
C
	uvf=fs*((AsS*X+bsS)**pS)*uvf
      RETURN
C
C  entry point for initialization of uvfInit
C                  
      ENTRY uvfInit(softPdfType,Nl,Limi,Prob,C,Mean,As,bs,p)
C
C  Check that LimiS and ProbS arrays are big enough
C
      if (Nl.GT.MAXNL) then
       print '(''Error in uvProNR.f: max(Nl) must be less than 101'')'
      endif
C
C  save arguments:     softPdfType,  Nl,  Limi,  Prob,  Mean, As,  bs,  p
C  in local variables: softPdfTypeS, NlS, LimiS, ProbS, MeanS, AsS, bsS, pS
C  and calculate fsoS and fspS

      softPdfTypeS=softPdfType
      AsS=As
      bsS=bs
	pS=p
C Case softPdfType=4: Grid linear pdf
      if (softPdfType.EQ.4) then
        NlS=Nl
        MeanS=Mean
        do il=1,3
          LimiS(il)=Limi(il)
        enddo
        do il=1,Nl-1
          fspS(il)=(Prob(il+1)-Prob(il))/Limi(2)
          fsoS(il)=Prob(il)
     &               -(Limi(1)+(il-1)*Limi(2))*fspS(il)
        enddo
        il=Nl
        fsoS(il)=Prob(il)
        fspS(il)=0
C Case softPdfType=3: Grid histogram pdf
      elseif (softPdfType.EQ.3) then
        NlS=Nl
        MeanS=Mean
        do il=1,3
          LimiS(il)=Limi(il)
        enddo
        do il=1,Nl-1
          ProbS(il)=Prob(il)
        enddo
C Case softPdfType=2: Linear pdf
      elseif (softPdfType.EQ.2) then
        NlS=Nl
        MeanS=Mean
        do il=1,Nl-1
          LimiS(il)=Limi(il)
          fspS(il)=   (Prob(il+1)-Prob(il))
     &               / (Limi(il+1)-Limi(il))
          fsoS(il)=Prob(il)-Limi(il)*fspS(il)
        enddo
        il=Nl
        LimiS(il)=Limi(il)
        fsoS(il)=Prob(il)
        fspS(il)=0
C Case softPdfType=1: Histogram pdf
      elseif (softPdfType.EQ.1) then
        NlS=Nl
        MeanS=Mean
        do il=1,Nl-1
          LimiS(il)=Limi(il)
          ProbS(il)=Prob(il)
        enddo
        il=Nl
        LimiS(il)=Limi(il)
      endif
C
C  get the inverse COVinv and the determinant detC of C
C
      COVinv=1/C
	detC=C
      pi=2.0*acos(0.0)
      normVal=1.0/( (2.0*pi)**(1/2.0) * sqrt(detC) )
      END

