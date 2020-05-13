C Calculate the inverse and det of positive semi-definite symmetric matrix
C
c This file contains program AS7 and AS6 from statlib 
* (http://lib.stat.cmu.edu/apstat) 
c Modifed by Marc Serre, added determinant, Feb 1998
c
        subroutine symInv(a, n, nn, w, c, detC, nullty, ifault)
c
c       Algorithm AS7, Applied Statistics, vol.17, 1968, p.198.
c
c       Forms in c( ) as lower triangle, a generalised inverse
c       of the positive semi-definite symmetric matrix a( )
c       order n, stored as lower triangle.
c
c       arguments:-
c       a()     = input, the symmetric matrix to be inverted, stored in
c                 lower triangular form
c       n       = input, order of the matrix
c       nn      = input, the size of the a and c arrays     n*(n+1)/2
c       w()     = workspace, dimension at least n.
c       c()     = output, the inverse of a (a generalized inverse if c is
c                 singular), also stored in lower triangular.
c                 c and a may occupy the same locations.
c       detC    = output, det of C
c       nullty  = output, the rank deficiency of a.
c       ifault  = output, error indicator
c                       = 1 if n < 1
c                       = 2 if a is not +ve semi-definite
c                       = 3 if nn < n*(n+1)/2
c                       = 0 otherwise
c
c***************************************************************************
c
        double precision a(nn), c(nn), w(n), x, zero, one, detC
c
        data zero, one /0.0d0, 1.0d0/
c
c       cholesky factorization of a, result in c
c
        call chol(a, n, nn, c, nullty, ifault)
        if(ifault.ne.0) return
C
C  Calculate the determinant of a as the square of the product of diagonal
C  coefficients of the lower diag matrix c, where a=c*c'
C
        detC=1.0
        do i=1,n
           detC=detC*c((i*(i+1))/2)
        enddo
        detC=detC*detC
c
c       invert c & form the product (cinv)'*cinv, where cinv is the inverse
c       of c, row by row starting with the last row.
c       irow = the row number, ndiag = location of last element in the row.
c
        irow=n
        ndiag=nn
   10   l=ndiag
        if (c(ndiag) .eq. zero) goto 60
        do 20 i=irow,n
          w(i)=c(l)
          l=l+i
   20   continue
        icol=n
        jcol=nn
        mdiag=nn
   30   l=jcol
        x=zero
        if(icol.eq.irow) x=one/w(irow)
        k=n
   40   if(k.eq.irow) go to 50
        x=x-w(k)*c(l)
        k=k-1
        l=l-1
        if(l.gt.mdiag) l=l-k+1
        go to 40
   50   c(l)=x/w(irow)
        if(icol.eq.irow) go to 80
        mdiag=mdiag-icol
        icol=icol-1
        jcol=jcol-1
        go to 30
   60   do 70 j=irow,n
          c(l)=zero
          l=l+j
   70   continue
   80   ndiag=ndiag-irow
        irow=irow-1
        if(irow.ne.0) go to 10
        return
        end
C
C
CC This file contains AS6 and the enhanced version ASR44.   See AS7 also.
C 
C 
      SUBROUTINE CHOL (A,N,NN,U,NULLTY,IFAULT)
C 
C       Algorithm AS6, Applied Statistics, vol.17, (1968)
C 
C       Given a symmetric matrix order n as lower triangle in a( )
C       calculates an upper triangle, u( ), such that uprime * u = a.
C       a must be positive semi-definite.  eta is set to multiplying
C       factor determining effective zero for pivot.
C 
C       arguments:-
C       a()     = input, a +ve definite matrix stored in lower-triangula
C                 form.
C       n       = input, the order of a
C       nn      = input, the size of the a and u arrays      n*(n+1)/2
C       u()     = output, a lower triangular matrix such that u*u' = a.
C                 a & u may occupy the same locations.
C       nullty  = output, the rank deficiency of a.
C       ifault  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 3 if nn < n*(n+1)/2
C                       = 0 otherwise
C 
C***********************************************************************
C 
      DOUBLE PRECISION A(NN),U(NN),ETA,ETA2,X,W,ZERO
C 
C       The value of eta will depend on the word-length of the
C       computer being used.  See introductory text.
C 
      DATA ETA,ZERO/1.D-9,0.0D0/
C 
      IFAULT=1
      IF (N.LE.0) RETURN
      IFAULT=3
      IF (NN.LT.N*(N+1)/2) RETURN
      IFAULT=2
      NULLTY=0
      J=1
      K=0
      ETA2=ETA*ETA
      II=0
C 
C       Factorize column by column, icol = column no.
C 
      DO 80 ICOL=1,N
        II=II+ICOL
        X=ETA2*A(II)
        L=0
        KK=0
C 
C       IROW = row number within column ICOL
C 
        DO 40 IROW=1,ICOL
          KK=KK+IROW
          K=K+1
          W=A(K)
          M=J
          DO 10 I=1,IROW
            L=L+1
            IF (I.EQ.IROW) GO TO 20
            W=W-U(L)*U(M)
            M=M+1
 10       CONTINUE
 20       IF (IROW.EQ.ICOL) GO TO 50
          IF (U(L).EQ.ZERO) GO TO 30
          U(K)=W/U(L)
          GO TO 40
 30       IF (W*W.GT.ABS(X*A(KK))) RETURN
          U(K)=ZERO
 40     CONTINUE
 50     IF (ABS(W).LE.ABS(ETA*A(K))) GO TO 60
        IF (W.LT.ZERO) RETURN
        U(K)=SQRT(W)
        GO TO 70
 60     U(K)=ZERO
        NULLTY=NULLTY+1
 70     J=J+ICOL
 80   CONTINUE
      IFAULT=0
      END
