C
C mvMomVecAG2g - Gateway function for mvMomVecAG2.f
C
C call in matlab:
C [Val Err info]=mvMomVecAG2(softPdfType,Nl,Limi,Prob,Mean,C,nMom,As,bs,p,maxpts,aEps,rEps)
C 
C call in fortran
C      CALL mvMomVecAG2(ndim, softPdfType, Nl, Limi, Prob, Mean, CORREL, 
C          nMom,As,bs,p, maxpts, aEps, rEps, key, Val, Err, info)
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
C parameters
      INTEGER MAXNDIM, MAXCDIM, MAXNL, MAXNMOM
      PARAMETER(MAXNDIM   = 20, MAXNL=101, MAXNMOM=10, 
     &          MAXCDIM = (MAXNDIM*(MAXNDIM+1))/2 )
C local integer and real variables
      INTEGER ndim, N, M, nlmax, key 
      real*8 CORREL(MAXCDIM)
C pointers to input and output matlab variables
      INTEGER Nlp, softPdfTypep, Limip, Probp, Meanp, Cp, maxptsp, 
     &     nMomp, Asp, Bsp, Pp, aEpsp, rEpsp, Valp, Errp, infop
C integer input and output for matlab function
      INTEGER Nl(MAXNDIM), softPdfType, nMom, maxpts,
     &     info
C above variables in real format (for passing from/to matlab)
      real*8 Nlr(MAXNDIM), softPdfTyper, nMomr, maxptsr,
     &     infor
C real input and output for matlab
      real*8  Limi(MAXNDIM*MAXNL), Prob(MAXNDIM*(MAXNL-1)),
     &  Mean(MAXNDIM),C(MAXNDIM*MAXNDIM),As(MAXNDIM*MAXNMOM),
     &  Bs(MAXNMOM), P(MAXNMOM), aEps, rEps,
     &  Val(MAXNMOM+1), Err(MAXNMOM+1)
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 13) THEN
        CALL MEXERRMSGTXT('mvMomVecAG2 requires 13 input arguments')
      ELSEIF (NLHS .GT. 3) THEN
        CALL MEXERRMSGTXT('mvMomVecAG2 requires 3 output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF Nl
C
      M = mxGetM(PRHS(2))
      N = mxGetN(PRHS(2))
      ndim = MAX(M,N)
C
      IF (ndim .GT. MAXNDIM) THEN
         CALL MEXERRMSGTXT('mvMomVecAG2 requires that length(Nl) < 20')
      ENDIF
C
      IF (MIN(M,N) .NE. 1) THEN
         CALL MEXERRMSGTXT('mvMomVecAG2 requires Nl to be a N x 1 vec')
      ENDIF
C
C-------------------------------------------------------------------
C  Before going any further we need to get softPdfType, Nl, and nMom. 
C  softPdfType, Nl and nMom are needed to extract the other RHS arguments
C
C ASSIGN POINTER
      softPdfTypep    = mxGetPr(PRHS(1))
      Nlp    = mxGetPr(PRHS(2))
      nMomp  = mxGetPr(PRHS(7))
C
C COPY POINTER TO LOCAL VARIABLE
      CALL mxCopyPtrToReal8(softPdfTypep, softPdfTyper, 1)
      CALL mxCopyPtrToReal8(Nlp, Nlr, ndim)
      CALL mxCopyPtrToReal8(nMomp, nMomr, 1)
C
C Convert into integer type
      softPdfType=softPdfTyper
      do i=1,ndim
         Nl(i)=Nlr(i)
      enddo
      nMom = nMomr
C-------------------------------------------------------------------
C
C
C Make sure that max(Nl) is less that MAXNL=101
      nlmax= -1
      do i=1,ndim
         if (nlmax.LT.Nl(i)) then
            nlmax=Nl(i)
         endif
      enddo
      IF (nlmax.GT. MAXNL) THEN
         CALL MEXERRMSGTXT('max(Nl) must be less than 101')
      ENDIF
C
C CHECK THE DIMENSIONS OF Limi
C
      M = mxGetM(PRHS(3))
      N = mxGetN(PRHS(3))
C
C case grid linear or grid histogram
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.3)) then
         IF ((M .NE. ndim) .OR. (N .NE. 3)) THEN
            CALL MEXERRMSGTXT('Limi must be a length(Nl) by 3 mat')
         ENDIF
C case linear or histogram
      elseif ((softPdfType.EQ.2).OR.(softPdfType.EQ.1)) then
        IF ((M .NE. ndim) .OR. (N .NE. nlmax)) THEN
          CALL MEXERRMSGTXT('Limi must be a length(Nl) by max(Nl) mat')
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
        IF ((M .NE. ndim) .OR. (N .NE. nlmax)) THEN
          CALL MEXERRMSGTXT('Prob must be a length(Nl) by max(Nl) mat')
        ENDIF
C case histogram
      elseif ((softPdfType.EQ.3).OR.(softPdfType.EQ.1)) then
        IF ((M .NE. ndim) .OR. (N .NE. nlmax-1)) THEN
         CALL MEXERRMSGTXT('Prob must be a length(Nl) by max(Nl)-1 mat')
        ENDIF
      endif
C
C CHECK THE DIMENSIONS OF Mean
C
      M = mxGetM(PRHS(5))
      N = mxGetN(PRHS(5))
C
      IF ((MAX(M,N) .NE. ndim).OR.(MIN(M,N) .NE. 1)) THEN
         CALL MEXERRMSGTXT('Mean to be a N x 1 vec')
      ENDIF
C
C CHECK THE DIMENSIONS OF C
C
      M = mxGetM(PRHS(6))
      N = mxGetN(PRHS(6))
C
      IF ((M .NE. ndim) .OR. (N .NE. ndim)) THEN
        CALL MEXERRMSGTXT('mvMomVecAG2 requires C be a N x N matrix')
      ENDIF
C
C CHECK THE DIMENSIONS OF As
C
      M = mxGetM(PRHS(8))
      N = mxGetN(PRHS(8))
C
      IF ((M .NE. ndim).OR.(N .NE. nMom)) THEN
         CALL MEXERRMSGTXT('As to be a N x nMom matrix')
      ENDIF
C
C CHECK THE DIMENSIONS OF Bs
C
      M = mxGetM(PRHS(9))
      N = mxGetN(PRHS(9))
C
      IF ((M.NE.1).OR.(N.NE.nMom)) THEN
         CALL MEXERRMSGTXT('Bs to be a 1 x nMom vec')
      ENDIF
C
C CHECK THE DIMENSIONS OF P
C
      M = mxGetM(PRHS(10))
      N = mxGetN(PRHS(10))
C
      IF ((M.NE.1).OR.(N.NE.nMom)) THEN
         CALL MEXERRMSGTXT('P to be a 1 x nMom vec')
      ENDIF
C
C CREATE MATRICES FOR RETURN ARGUMENTS
C
      PLHS(1) = mxCreateFull(nMom,1,0)
      PLHS(2) = mxCreateFull(nMom,1,0)
      PLHS(3) = mxCreateFull(1,1,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      Valp = mxGetPr(PLHS(1))
      Errp = mxGetPr(PLHS(2))
      infop = mxGetPr(PLHS(3))
C
C      Nlp    = mxGetPr(PRHS(1))
C      softPdfTypep = mxGetPr(PRHS(2))
      Limip  = mxGetPr(PRHS(3))
      Probp  = mxGetPr(PRHS(4))
      Meanp  = mxGetPr(PRHS(5))
      Cp     = mxGetPr(PRHS(6))
C      nMomp    = mxGetPr(PRHS(7))
      Asp    = mxGetPr(PRHS(8))
      bsp    = mxGetPr(PRHS(9))
      pp     = mxGetPr(PRHS(10))
      maxptsp= mxGetPr(PRHS(11))
      aEpsp  = mxGetPr(PRHS(12))
      rEpsp  = mxGetPr(PRHS(13))
C
C COPY OTHER RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
C
C      CALL mxCopyPtrToReal8(Nlp, Nlr, ndim)
C      CALL mxCopyPtrToReal8(softPdfTypep, softPdfTyper, 1)
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.3)) then
         CALL mxCopyPtrToReal8(Limip, Limi, ndim*3)
      elseif ((softPdfType.EQ.2).OR.(softPdfType.EQ.1)) then
         CALL mxCopyPtrToReal8(Limip, Limi, ndim*nlmax)
      endif
      if ((softPdfType.EQ.4).OR.(softPdfType.EQ.2)) then
         CALL mxCopyPtrToReal8(Probp, Prob, ndim*nlmax)
      elseif ((softPdfType.EQ.3).OR.(softPdfType.EQ.1)) then
         CALL mxCopyPtrToReal8(Probp, Prob, ndim*(nlmax-1))
      endif
      CALL mxCopyPtrToReal8(Meanp, Mean, ndim)
      CALL mxCopyPtrToReal8(Cp, C, ndim*ndim)
C      CALL mxCopyPtrToReal8(nMomp, nMomr, 1)
      CALL mxCopyPtrToReal8(Asp, As, ndim*nMom)
      CALL mxCopyPtrToReal8(Bsp, Bs, nMom)
      CALL mxCopyPtrToReal8(Pp, P, nMom)
      CALL mxCopyPtrToReal8(maxptsp, maxptsr, 1)
      CALL mxCopyPtrToReal8(aEpsp, aEps, 1)
      CALL mxCopyPtrToReal8(rEpsp, rEps, 1)
C
C Convert necessary input matlab argument to integer type
C
      maxpts=maxptsr
C
C transfer C array into CORREL array
C
      k=0
      DO i=1,ndim
         DO j=1,i
            k=k+1
            ij= i+(j-1)*ndim
            CORREL(k)=C(ij)
         ENDDO
      ENDDO

C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
      key=0
      CALL mvMomVecAG2(ndim, softPdfType, Nl, Limi, Prob, Mean, CORREL,
     &    nMom,As,Bs,P,maxpts, aEps, rEps, key, Val, Err, info)
C
C Convert integer output into real before passing it to matlab
C
      infor=info
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
C

      CALL mxCopyReal8ToPtr(Val, Valp, nMom)
      CALL mxCopyReal8ToPtr(Err, Errp, nMom)
      CALL mxCopyReal8ToPtr(infor, infop, 1)
C      
C      if (info.EQ.1) then
C         print '(''mvMomVecAG2 warning: info=1, maxpts too small'')'
C      elseif (info.NE.0) then
C         print '(''mvMomVecAG2 warning: info not equal to zero'')'
C      endif

      RETURN
      END


