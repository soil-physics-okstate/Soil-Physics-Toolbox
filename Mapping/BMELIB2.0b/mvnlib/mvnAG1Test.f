*
*
*     test program for the SADMVN program from Alan Genz
*
*     PROGRAM mvnAG1Test
*
*
      DOUBLE PRECISION ABSEPS, RELEPS
      DOUBLE PRECISION VALR, ERRR, VALS, ERRS, VALK, ERRK, VALD, ERRD
      INTEGER N, NN, I, J, K, IJ, MAXPTS, IFTR, IFTS, IFTK, IFTD
      PARAMETER ( N = 5, NN = ( N - 1 )*N/2, MAXPTS = 5000*N*N*N )
      PARAMETER ( ABSEPS = 0.00005, RELEPS = 0 )
C      DOUBLE PRECISION CORREL(NN), LOW(N), UP(N), CONDIT
      DOUBLE PRECISION CORREL(NN), LOW(N), UP(N)
      INTEGER INFIN(N)
*          Chen Problem
      DATA ( UP(I), I=1,N)  /.0, 1.5198, 1.7817, 1.4755, 1.5949/
      DATA (LOW(I), I=1,N)  /.0,  .0   , 1.7817, 1.4755, 1.5949/
      DATA (INFIN(I), I=1,N)/ 1, 2     , 1     , 1     , 0     /
      DATA (CORREL(I),I=1,NN)/-0.707107,0.0,0.5,0.0,2*0.5,0.0,3*0.5/
      PRINT '(''        Test of RANMVN, SADMVN, KROMVN and SPHMVN'')'
      PRINT '(12X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(''           Number of Dimensions is '',I2)', N
      PRINT '(''     Maximum # of Function Values is '',I7)', MAXPTS
*
      DO K = 1,3
         PRINT '(/'' I     Limits'')'
         PRINT '(4X,''Lower  Upper  Lower Left of Correlation Matrix'')'
         IJ = 0
         DO I = 1,N
            IF ( INFIN(I) .LT. 0 ) THEN 
               PRINT '(I2, '' -infin  infin '', 7F9.5)',
     &              I, ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ELSE IF ( INFIN(I) .EQ. 0 ) THEN 
               PRINT '(I2, '' -infin'', F7.4, 1X, 7F9.5)',
     &              I, UP(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ELSE IF ( INFIN(I) .EQ. 1 ) THEN 
               PRINT '(I2, F7.4, ''  infin '', 7F9.5)',
     &              I, LOW(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ELSE 
               PRINT '(I2, 2F7.4, 1X, 7F9.5)', 
     &              I, LOW(I), UP(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ENDIF
            IJ = IJ + I-1
         END DO
C         PRINT '(20X, ''Condition Number is'', F8.1)', CONDIT(N,CORREL) 
         CALL SADMVN( N, LOW, UP, INFIN, CORREL, 
     &        MAXPTS, ABSEPS, RELEPS, ERRS, VALS, IFTS )
C         CALL KROMVN( N, LOW, UP, INFIN, CORREL, 
C     &        MAXPTS, ABSEPS, RELEPS, ERRK, VALK, IFTK )
C         CALL RANMVN( N, LOW, UP, INFIN, CORREL, 
C     &        MAXPTS, ABSEPS, RELEPS, ERRR, VALR, IFTR )
C         CALL SPHMVN( N, LOW, UP, INFIN, CORREL, 
C     &        MAXPTS, ABSEPS, RELEPS, ERRD, VALD, IFTD )
         PRINT '('' Results for:  RANMVN'',9X,''SADMVN'',9X,''KROMVN'', 
     &                      9X, ''SPHMVN'')'
         PRINT '('' Values:   '',4(F11.6,I4))',
     &        VALR, IFTR, VALS, IFTS, VALK, IFTK, VALD, IFTD
         PRINT '(''Errests:   '',4(2X,''('',F8.6'')'',3X))', 
     &        ERRR, ERRS, ERRK, ERRD
         INFIN(1) = INFIN(1) - 1
      END DO
      END
