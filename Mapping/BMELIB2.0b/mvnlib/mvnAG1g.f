C
C mvnAG1g - Gateway function for mvnAG1.f
C

C
C The max dim for vector LOWEWR, UPPER, etc..
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

      INTEGER LOWERp, UPPERp, Cp, MAXPTSp, ABSEPSp, 
     &     RELEPSp, VALUEp, ERRORp, INFORMp

      INTEGER N, M
      INTEGER maxNDIM,maxCDIM
      PARAMETER(maxNDIM   = 20, 
     &          maxCDIM = (maxNDIM*(maxNDIM-1))/2 )

      INTEGER NDIM, INFIN(maxNDIM), MAXPTS, 
     &     INFORM
      real*8 MAXPTSr, 
     &     INFORMr
      real*8  C(maxNDIM*maxNDIM),
     &  LOWER(maxNDIM), UPPER(maxNDIM), CORREL(maxCDIM),ABSEPS, RELEPS,
     &  ERROR, VALUE
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 6) THEN
        CALL MEXERRMSGTXT('mvnAG1 requires 6 input arguments')
      ELSEIF (NLHS .GT. 3) THEN
        CALL MEXERRMSGTXT('mvnAG1 requires 3 output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF LOWER
C
      M = mxGetM(PRHS(1))
      N = mxGetN(PRHS(1))
      NDIM = MAX(M,N)
C
      IF (NDIM .GT. maxNDIM) THEN
         CALL MEXERRMSGTXT('mvnAG1 requires that length(A) < 20')
      ENDIF
C
      IF (MIN(M,N) .NE. 1) THEN
         CALL MEXERRMSGTXT('mvnAG1 requires A to be a N x 1 vec')
      ENDIF
C
C CHECK THE DIMENSIONS OF UPPER
C
      M = mxGetM(PRHS(2))
      N = mxGetN(PRHS(2))
C
      IF ((MAX(M,N) .NE. NDIM) .OR. (MIN(M,N) .NE. 1)) THEN
        CALL MEXERRMSGTXT('mvnAG1 requires B to be a N x 1 vec')
      ENDIF
C
C CHECK THE DIMENSIONS OF CORREL
C
      M = mxGetM(PRHS(3))
      N = mxGetN(PRHS(3))
C
      IF ((M .NE. NDIM) .OR. (N .NE. NDIM)) THEN
        CALL MEXERRMSGTXT('mvnAG1 requires C be a N x N matrix')
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
      VALUEp = mxGetPr(PLHS(1))
      ERRORp = mxGetPr(PLHS(2))
      INFORMp = mxGetPr(PLHS(3))
C
      LOWERp = mxGetPr(PRHS(1))
      UPPERp = mxGetPr(PRHS(2))
      Cp     = mxGetPr(PRHS(3))
      MAXPTSp = mxGetPr(PRHS(4))
      ABSEPSp = mxGetPr(PRHS(5))
      RELEPSp = mxGetPr(PRHS(6))
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
C
      CALL mxCopyPtrToReal8(LOWERp, LOWER, NDIM)
      CALL mxCopyPtrToReal8(UPPERp, UPPER, NDIM)
      CALL mxCopyPtrToReal8(Cp, C, NDIM*NDIM)
      CALL mxCopyPtrToReal8(MAXPTSp, MAXPTSr, 1)
      CALL mxCopyPtrToReal8(ABSEPSp, ABSEPS, 1)
      CALL mxCopyPtrToReal8(RELEPSp, RELEPS, 1)
C
C Get integer arguments from their real value
C
      DO i=1,NDIM
         INFIN(i)=2
      ENDDO
      MAXPTS=MAXPTSr
C
C Scale the vector LOWER and UPPER and the matrix C so that its
C diagonal elements are ones, then put the lower off-diagonal 
C elements of the matrix C, i.e. C(i,j) for j<i, in the vector 
C CORREL by row order.
C

C      PRINT '(''Vector LOWER:'')'
C      PRINT '(7F9.5)', ( LOWER(i),i=1,NDIM )
C      PRINT '(''Vector UPPER:'')'
C      PRINT '(7F9.5)', ( UPPER(i),i=1,NDIM )
C      PRINT '(''C:'')'
C      DO i=1,NDIM
C         PRINT '(7F9.5)', ( C(i+(j-1)*NDIM),j=1,NDIM )
C      ENDDO

      DO i=1,NDIM
         ii = i+(i-1)*NDIM 
         LOWER(i)=LOWER(i)/sqrt(C(ii))
         UPPER(i)=UPPER(i)/sqrt(C(ii))
      EndDO

      IF (NDIM.GT.1) THEN

         DO i=2,NDIM
            ii = i+(i-1)*NDIM
            DO j=1,i-1
               jj= j+(j-1)*NDIM
               ij= i+(j-1)*NDIM
               C(ij)=C(ij)/sqrt(C(ii)*C(jj))
            ENDDO
         ENDDO

         k=0
         DO i=2,NDIM
            DO j=1,i-1
               k=k+1
               ij= i+(j-1)*NDIM
               CORREL(k)=C(ij)
            ENDDO
         ENDDO

      ENDIF

C      PRINT '(''Vector LOWER:'')'
C      PRINT '(7F9.5)', ( LOWER(i),i=1,NDIM )
C      PRINT '(''Vector UPPER:'')'
C      PRINT '(7F9.5)', ( UPPER(i),i=1,NDIM )
C      PRINT '(''Vector INFIN:'')'
C      PRINT '(7I12)', ( INFIN(i),i=1,NDIM )
C      PRINT '(''C:'')'
C      DO i=1,NDIM
C         PRINT '(7F9.5)', ( C(i+(j-1)*NDIM),j=1,NDIM )
C      ENDDO
C      PRINT '(''Vector CORREL:'')'
C      PRINT '(28F9.5)', ( CORREL(i),i=1,(NDIM*(NDIM-1))/2 )
C      PRINT '(''maxpts='',I12,'', aEps='',F9.6,'', rEps='',F9.6)',
C     &      MAXPTS,ABSEPS,RELEPS

C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
      CALL SADMVN(NDIM, LOWER, UPPER, INFIN, CORREL, MAXPTS, 
     &     ABSEPS, RELEPS, ERROR, VALUE, INFORM)
C
C Convert integer output into real before passing it to matrix output
C
      INFORMr=INFORM
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
C
      CALL mxCopyReal8ToPtr(VALUE, VALUEp, 1)
      CALL mxCopyReal8ToPtr(ERROR, ERRORp, 1)
      CALL mxCopyReal8ToPtr(INFORMr, INFORMp, 1)

C      if (INFORM.EQ.1) then
C         print '(''mvnAG1 warning: Inform=1, maxpts too small'')'
C      elseif (INFORM.NE.0) then
C         print '(''mvnAG1 warning: Inform not equal to zero'')'
C      endif

      RETURN
      END
