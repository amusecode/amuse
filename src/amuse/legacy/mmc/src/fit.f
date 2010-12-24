cAttempt to fit numerical data from numeric-fit-svn to produce ideal model!
      implicit double precision (a-h,o-z)
      parameter (ndatamax = 128)
      dimension data(ndatamax,16),xinc(6),sig(ndatamax)
      integer xinc
      data xinc/1,2,3,4,6,8/
      open (7,file='temp')
      i = 1
 10   continue
      read(7,*,end=20) (data(i,j),j=1,16)
      sig(i) = 1.d0
      i = i+1
      goto 10
 20   continue
      ndata = i - 1
c      print*,ndata
      mwt = 0
      do i = 1,7
         do k = 1,6
            j = xinc(k)
            call fit(data(1,j),data(1,i+9),ndata,sig,mwt,a,b,siga,
     &           sigb,chi2,q)
            write (6,*) i,j,a,siga,b,sigb
         enddo
      enddo
      stop
      end


      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      implicit double precision (a-h,o-z)
      INTEGER mwt,ndata
c      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
      double precision a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
CU    USES gammq
      INTEGER i
c      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      double precision sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
c      print*,x
c      print*,y
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        q=1.
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
         
        do 16 i=1,ndata
c           print*,y(i),a,b,x(i),sig(i)
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
c          print*,'chi2:',chi2
16      continue
c        print*,'fit:',ndata,chi2
c      stop
        q=gammq(0.5*(ndata-2),0.5*chi2)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      double precision FUNCTION gammq(a,x)
      implicit double precision (a-h,o-z)
c      REAL a,gammq,x
      double precision a,gammq,x
CU    USES gcf,gser
c      REAL gammcf,gamser,gln
      double precision gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
c         print*,'gammq:',a
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      double precision FUNCTION gammln(xx)
      implicit double precision (a-h,o-z)
      double precision  gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      SUBROUTINE gser(gamser,a,x,gln)
      implicit double precision (a-h,o-z)
      INTEGER ITMAX
c      REAL a,gamser,gln,x,EPS
      double precision a,gamser,gln,x,EPS
      PARAMETER (ITMAX=1000,EPS=3.e-12)
CU    USES gammln
      INTEGER n
c      REAL ap,del,sum,gammln
      double precision  ap,del,sum,gammln
      gln=gammln(a)
c      print*,'a= ',a
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.

      SUBROUTINE gcf(gammcf,a,x,gln)
      implicit double precision (a-h,o-z)
      INTEGER ITMAX
c      REAL a,gammcf,gln,x,EPS,FPMIN
      double precision a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=1000,EPS=3.e-12,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
c      REAL an,b,c,d,del,h,gammln
      double precision  an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
c        print*,del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
c      print*,a
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.
            
      
