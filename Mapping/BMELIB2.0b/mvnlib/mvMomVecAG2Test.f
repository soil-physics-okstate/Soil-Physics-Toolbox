C
C     This is a summary of the routines in this file. 
C
C  mvMomVecAG2Test is a simple test driver for mvMomVecAG2.f
*test for:
*   [3      [0 0.5 1       [1 1    [2.0 0.2 0.2        [0
* Nl=3  Limi=0 0.5 1 , Prob=1 1  C= 0.2 2.0 0.2  , Mean=0 
*    3]      0 0.5 1]       1 1]    0.2 0.2 2.0]        0]
*
*   [1        
* As=0  bs=0  p=1
*    0]        
*
* results on a HP 715
*
*
*Using softPdfType=1
*       Value       Error
*        .008876    .956E-12
*Using softPdfType=2
*       Value       Error
*        .008876    .956E-12
*Using softPdfType=3
*       Value       Error
*        .008876    .956E-12
*Using softPdfType=4
*       Value       Error
*        .008876    .956E-12
      PROGRAM mvMomVecAG2Test
      INTEGER key, i, j, N, nDim, maxpts, ifail, softPdfType
      PARAMETER (nDim = 3)
      INTEGER Nl(nDim), nMom
      DOUBLE PRECISION Limi(nDim*3), Prob(nDim*3), C((nDim*(nDim+1))/2),
     &   Mean(nDim), As(nDim), bs, p 
      DOUBLE PRECISION Error, Value, aEps, rEps
      DO N = 1,nDim
         Nl(N)=3
         Limi(N) = 0
         Limi(N+nDim) = 0.5
         Limi(N+2*nDim) = 1
         Prob(N) = 1
         Prob(N+nDim) = 1
         Prob(N+2*nDim) = 1
         Mean(N)=0
      ENDDO
      nMom=1
      As(1)=1 
      As(2)=0 
      As(3)=0 
      bs=0
      p=1
      ij=0
      do i=1,nDim
         do j=1,i-1
            ij=ij+1
            C(ij)=0.2
         enddo
         ij=ij+1
         C(ij)=2.0
      enddo
      maxpts = 10000
      key = 0
      aEps = 0
      rEps = 1E-6      
C case 1: histogram
      print '(''Using softPdfType=1'')'
      softPdfType=1
      CALL mvMomVecAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * nMom,As,bs,p,maxpts, aEps, rEps, key, Value, Error, ifail)
      print *, '      Value       Error'
      print 9999, Value, Error
 9999 FORMAT (3X, F12.6, G12.3)
C case 2: linear
      print '(''Using softPdfType=2'')'
      softPdfType=2
      CALL mvMomVecAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * nMom,As,bs,p,maxpts, aEps, rEps, key, Value, Error, ifail)
      print *, '      Value       Error'
      print 9999, Value, Error
C case 3: grid histogram
      print '(''Using softPdfType=3'')'
      softPdfType=3
      CALL mvMomVecAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * nMom,As,bs,p,maxpts, aEps, rEps, key, Value, Error, ifail)
      print *, '      Value       Error'
      print 9999, Value, Error
C case 4: grid linear
      print '(''Using softPdfType=4'')'
      softPdfType=4
      CALL mvMomVecAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * nMom,As,bs,p,maxpts, aEps, rEps, key, Value, Error, ifail)
      print *, '      Value       Error'
      print 9999, Value, Error
      END

