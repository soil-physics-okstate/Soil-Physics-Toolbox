C
C     This is a summary of the routines in this file. 
C
C  mvProAG2Test is a simple test driver for mvProAG2.f
*test for:
*   [3      [0 0.5 1       [1 1    [2.0 0.2 0.2        [0
* Nl=3  Limi=0 0.5 1 , Prob=1 1  C= 0.2 2.0 0.2  , Mean=0 
*    3]      0 0.5 1]       1 1]    0.2 0.2 2.0]        0]
*
* results on a HP 715
*
* Using softPdfType=1
* ifail =  0
*     N   ESTIMATED ERROR    INTEGRAL
*     4        .000000    .183833E-01
* Using softPdfType=2
* ifail =  0
*     N   ESTIMATED ERROR    INTEGRAL
*     4        .000000    .183833E-01
* Using softPdfType=3
* ifail =  0
*     N   ESTIMATED ERROR    INTEGRAL
*     4        .000000    .183833E-01
* Using softPdfType=4
* ifail =  0
*     N   ESTIMATED ERROR    INTEGRAL
*     4        .000000    .183833E-01
      PROGRAM mvProAG2Test
      INTEGER key, i, j, N, nDim, maxpts, ifail, softPdfType
      PARAMETER (nDim = 3)
      INTEGER Nl(nDim)
      DOUBLE PRECISION Limi(nDim*3), Prob(nDim*3), C((nDim*(nDim+1))/2),
     &   Mean(nDim)
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
      CALL mvProAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * maxpts, aEps, rEps, key, Value, Error, ifail)
      PRINT 9999, ifail
 9999 FORMAT ( 'ifail = ', I2, /'    N   ESTIMATED ERROR    INTEGRAL')
      PRINT 9998, N, Error, Value
 9998 FORMAT (3X, I2, F15.6, G15.6)
C case 2: linear
      print '(''Using softPdfType=2'')'
      softPdfType=2
      CALL mvProAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * maxpts, aEps, rEps, key, Value, Error, ifail)
      PRINT 9999, ifail
      PRINT 9998, N, Error, Value
C case 3: grid histogram
      print '(''Using softPdfType=3'')'
      softPdfType=3
      CALL mvProAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * maxpts, aEps, rEps, key, Value, Error, ifail)
      PRINT 9999, ifail
      PRINT 9998, N, Error, Value
C case 4: grid linear
      print '(''Using softPdfType=4'')'
      softPdfType=4
      CALL mvProAG2(nDim, softPdfType, Nl, Limi, Prob, Mean, C, 
     * maxpts, aEps, rEps, key, Value, Error, ifail)
      PRINT 9999, ifail
      PRINT 9998, N, Error, Value
      END

