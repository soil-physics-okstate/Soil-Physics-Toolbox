C
C uvMomNRg - Gateway function for uvMomNR.f
C
C call in matlab:
C [Val Err info]=uvMomNR(softPdfType,Nl,Limi,Prob,Mean,C,As,bs,p,aEps,rEps)
C 
C call in fortran
C      CALL uvMomNR(softPdfType, Nl, Limi, Prob, Mean, C, As, bs, p,
C     &   aEps, rEps, Val, Err, info)
C
C
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
C
      INTEGER PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS
C
      INTEGER mxCreateFull, mxGetPr
      INTEGER mxGetM, mxGetN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C local integer and real variables
      INTEGER N, M, MAXNL
	PARAMETER (MAXNL = 101)
C pointers to input and output matlab variables
      INTEGER Nlp, softPdfTypep, Limip, Probp, Meanp, Cp, Asp,bsp,pp,
     &     aEpsp, rEpsp, Valp, Errp, infop
C integer input and output for matlab function
      INTEGER Nl, softPdfType, info
C above variables in real format (for passing from/to matlab)
      real*8 Nlr, softPdfTyper, infor
C real input and output for matlab
      real*8  Limi(MAXNL), Prob(MAXNL),Mean, C, As,bs,p, aEps, rEps,
     &  Val, Err
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 11) THEN
        CALL MEXERRMSGTXT('uvMomNR requires 11 input arguments')
      ELSEIF (NLHS .GT. 3) THEN
        CALL MEXERRMSGTXT('uvMomNR requires 3 output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF Nl
C
      M = mxGetM(PRHS(2))
      N = mxGetN(PRHS(2))
C
      IF (M .NE. 1) THEN
         CALL MEXERRMSGTXT('Nl must be an integer')
      ENDIF
C
      IF (N .NE. 1) THEN
         CALL MEXERRMSGTXT('Nl must be an integer')
      ENDIF
C
C-------------------------------------------------------------------
C  Before going any further we need to get softPdfType and Nl. 
C  softPdfType and Nl is needed to extract the other RHS arguments
C
C ASSIGN POINTER
      softPdfTypep    = mxGetPr(PRHS(1))
      Nlp    = mxGetPr(PRHS(2))
C
C COPY POINTER TO LOCAL VARIABLE
      CALL mxCopyPtrToReal8(softPdfTypep, softPdfTyper, 1)
      CALL mxCopyPtrToReal8(Nlp, Nlr, 1)
C
C Convert into integer type
      softPdfType=softPdfTyper
      Nl=Nlr
C-------------------------------------------------------------------
C
C
C Make sure that Nl is less that MAXNL=101
      IF (Nl.GT. MAXNL) THEN
         CALL MEXERRMSGTXT('Nl must be less than 101')
      ENDIF
C
C CHECK THE DIMENSIONS OF Limi
C
      M = mxGetM(PRHS(3))
      N = mxGetN(PRHS(3))
C
C case grid linear or grid histogram
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.3)) then
         IF ((M .NE. 1) .OR. (N .NE. 3)) THEN
            CALL MEXERRMSGTXT('Limi must be a 1 by 3 vector')
         ENDIF
C case linear or histogram
      elseif ((softPdfType.EQ.2).OR.(softPdfType.EQ.1)) then
        IF ((M .NE. 1) .OR. (N .NE. Nl)) THEN
          CALL MEXERRMSGTXT('Limi must be a 1 by Nl vector')
        ENDIF
C other cases
      else
         CALL MEXERRMSGTXT('Unacceptable value for softPdfType')
      endif
C
C CHECK THE DIMENSIONS OF Prob
C
      M = mxGetM(PRHS(4))
      N = mxGetN(PRHS(4))
C
C case linear
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.2)) then
        IF ((M .NE. 1) .OR. (N .NE. Nl)) THEN
          CALL MEXERRMSGTXT('Prob must be a 1 by Nl vector')
        ENDIF
C case histogram
      elseif ((softPdfType.EQ.3).OR.(softPdfType.EQ.1)) then
        IF ((M .NE. 1) .OR. (N .NE. Nl-1)) THEN
         CALL MEXERRMSGTXT('Prob must be a 1 by Nl-1 vector')
        ENDIF
      endif
C
C CHECK THE DIMENSIONS OF Mean
C
      M = mxGetM(PRHS(5))
      N = mxGetN(PRHS(5))
C
      IF ((M .NE. 1).OR.(N .NE. 1)) THEN
         CALL MEXERRMSGTXT('Mean must be a scalar')
      ENDIF
C
C CHECK THE DIMENSIONS OF C
C
      M = mxGetM(PRHS(6))
      N = mxGetN(PRHS(6))
C
      IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
        CALL MEXERRMSGTXT('C must a scalar')
      ENDIF
C
C CHECK THE DIMENSIONS OF As
C
      M = mxGetM(PRHS(7))
      N = mxGetN(PRHS(7))
C
      IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
        CALL MEXERRMSGTXT('As must a scalar')
      ENDIF
C
C CHECK THE DIMENSIONS OF bs
C
      M = mxGetM(PRHS(8))
      N = mxGetN(PRHS(8))
C
      IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
        CALL MEXERRMSGTXT('bs must a scalar')
      ENDIF
C
C CHECK THE DIMENSIONS OF p
C
      M = mxGetM(PRHS(9))
      N = mxGetN(PRHS(9))
C
      IF ((M .NE. 1) .OR. (N .NE. 1)) THEN
        CALL MEXERRMSGTXT('p must a scalar')
      ENDIF
C
C CREATE MATRICES FOR RETURN ARGUMENTS
C
      PLHS(1) = mxCreateFull(1,1,0)
      PLHS(2) = mxCreateFull(1,1,0)
      PLHS(3) = mxCreateFull(1,1,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      Valp = mxGetPr(PLHS(1))
      Errp = mxGetPr(PLHS(2))
      infop = mxGetPr(PLHS(3))
C
      Limip  = mxGetPr(PRHS(3))
      Probp  = mxGetPr(PRHS(4))
      Meanp  = mxGetPr(PRHS(5))
      Cp     = mxGetPr(PRHS(6))
      Asp    = mxGetPr(PRHS(7))
      bsp    = mxGetPr(PRHS(8))
      pp     = mxGetPr(PRHS(9))
      aEpsp  = mxGetPr(PRHS(10))
      rEpsp  = mxGetPr(PRHS(11))
C
C COPY OTHER RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
C
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.3)) then
         CALL mxCopyPtrToReal8(Limip, Limi, 3)
      elseif ((softPdfType.EQ.2).OR.(softPdfType.EQ.1)) then
         CALL mxCopyPtrToReal8(Limip, Limi, Nl)
      endif
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.2)) then
         CALL mxCopyPtrToReal8(Probp, Prob, Nl)
      elseif ((softPdfType.EQ.3).OR.(softPdfType.EQ.1)) then
         CALL mxCopyPtrToReal8(Probp, Prob, Nl-1)
      endif
      CALL mxCopyPtrToReal8(Meanp, Mean, 1)
      CALL mxCopyPtrToReal8(Cp, C, 1)
      CALL mxCopyPtrToReal8(Asp, As, 1)
      CALL mxCopyPtrToReal8(bsp, bs, 1)
      CALL mxCopyPtrToReal8(pp, p, 1)
      CALL mxCopyPtrToReal8(aEpsp, aEps, 1)
      CALL mxCopyPtrToReal8(rEpsp, rEps, 1)
C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
      CALL uvMomNR(softPdfType, Nl, Limi, Prob, Mean, C, As,bs,p, 
     &    aEps, rEps, Val, Err, info)
C
C Convert integer output into real before passing it to matlab
C
      infor=info
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
C

      CALL mxCopyReal8ToPtr(Val, Valp, 1)
      CALL mxCopyReal8ToPtr(Err, Errp, 1)
      CALL mxCopyReal8ToPtr(infor, infop, 1)

      RETURN
      END
