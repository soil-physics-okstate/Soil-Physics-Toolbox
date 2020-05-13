      SUBROUTINE qromb3(func,a,b,aEps,rEps,ss,dss,ifail)
      INTEGER JMAX,JMAXP,K,KM,ifail
      DOUBLE PRECISION a,b,func,aEps,rEps,ss,dss
      EXTERNAL func
      PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      DOUBLE PRECISION h(JMAXP),s(JMAXP)
      ifail = 0
      h(1)=1.
      do 11 j=1,JMAX
        call trapzd3(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint3(h(j-KM),s(j-KM),K,0.,ss,dss)
          dss = abs(dss) 
          if ( (dss.le.aEps) .OR. (dss.le.rEps*abs(ss)) ) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25*h(j)
11    continue
      ifail=1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #4-UZ2.
