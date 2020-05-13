C  uvProNRtest is a simple test driver for uvProNR.f
*test for:
* Nl=[3]  Limi=[0 5 10] , Prob=[1 1 1]  C= [2.0]  , Mean=0 
*
* results on a HP 715
*
*Using softPdfType=1
*ifail =  0
*    ESTIMATED ERROR    INTEGRAL
*      0.752988E-07   0.499999
*Using softPdfType=2
*ifail =  0
*    ESTIMATED ERROR    INTEGRAL
*      0.752988E-07   0.499999
*Using softPdfType=3
*ifail =  0
*    ESTIMATED ERROR    INTEGRAL
*      0.752988E-07   0.499999
*Using softPdfType=4
*ifail =  0
*    ESTIMATED ERROR    INTEGRAL
*      0.752988E-07   0.499999
*

      PROGRAM uvProNRtest
      INTEGER maxpts, ifail, softPdfType
      INTEGER Nl
      DOUBLE PRECISION Limi(3), Prob(3), C, Mean
      DOUBLE PRECISION Error, Value, aEps, rEps
      
	softPdfType = 1
      Nl=3
      Limi(1) = 0
      Limi(2) = 5
      Limi(3) = 10
      Prob(1) = 1
      Prob(2) = 1
      Prob(3) = 1
      aEps = 0
      rEps = 1E-6      

      C=2.0
      Mean=0.0

C case 1: histogram
      print '(''Using softPdfType=1'')'
      softPdfType=1
      CALL uvProNR(softPdfType, Nl, Limi, Prob, Mean, C, 
     * aEps, rEps, Value, Error, ifail)
      PRINT 9999, ifail
 9999 FORMAT ( 'ifail = ', I2, /'    ESTIMATED ERROR    INTEGRAL')
      PRINT 9998, Error, Value
 9998 FORMAT (3X, G15.6, G15.6)
C case 2: linear
      print '(''Using softPdfType=2'')'
      softPdfType=2
      CALL uvProNR(softPdfType, Nl, Limi, Prob, Mean, C, 
     * aEps, rEps, Value, Error, ifail)
      PRINT 9999, ifail
      PRINT 9998, Error, Value
C case 3: grid histogram
      print '(''Using softPdfType=3'')'
      softPdfType=3
      CALL uvProNR(softPdfType, Nl, Limi, Prob, Mean, C, 
     * aEps, rEps, Value, Error, ifail)
      PRINT 9999, ifail
      PRINT 9998, Error, Value
C case 4: grid linear
      print '(''Using softPdfType=4'')'
      softPdfType=4
      CALL uvProNR(softPdfType, Nl, Limi, Prob, Mean, C, 
     * aEps, rEps, Value, Error, ifail)
      PRINT 9999, ifail
      PRINT 9998, Error, Value
      END

