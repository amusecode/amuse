! PRTOFT:
!   Convert pressure and density pair to f and T pair.
!   Uses numerical recipes routine newt2d (to be replaced by something
!   else for license and efficiency reasons - it's only a 2D matrix, we
!   don't need LU decomposition to solve the system of equations).
      SUBROUTINE PRTOFT (logP,logrho,logf,logT)
      use constants
      implicit none
      double precision, intent(out) :: logT, logf
      double precision, intent(in) :: logP, logrho
      double precision :: XA(9), NA(9)
      double precision :: NEO, NIO, NZZ, AVM, NE;
      COMMON /ABUND / XA, NA, NEO, NIO, NZZ, AVM, NE
      double precision :: AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA
      double precision :: SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR(23);
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR
      double precision :: X(2),Y(2), CE
      logical :: check;


! Initial guesses for F, T assuming ionized non-degenerate ideal gas
      NIO = XA(1) + (42.0*XA(2) + 14.0*XA(5) + 12.0*XA(4) + 10.5*XA(3) +
     &      8.4*XA(6) + 7.0*XA(7) + 6.0*XA(8) + 3.0*XA(9))/168.0
      NEO = 0.5*(1.0 + XA(1))

      logT = logP - logrho - LOG(CR*(NIO+NEO))
      logf = (logrho - 1.5*logt + LOG(0.5*CEN*NEO)) - 2.0*(1.0-LOG(2.0))
! Find F, T by Newton-Raphson iteration
      X(1) = logf
      X(2) = logT
      Y(1) = logP
      Y(2) = logrho
      CALL NEWT2D(X,Y,CHECK)
      !IF(CHECK) WRITE(*,*) 'Failed to find global minimum...'
      logf = X(1)
      logT = X(2)
      
      end subroutine
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,y,f,stpmax,check,func)
      IMPLICIT NONE
      INTEGER n
      LOGICAL check
      double precision f,fold,stpmax,g(n),p(n),x(n),y(n),xold(n),func,ALF,TOLX
      PARAMETER (ALF=1.d-4,TOLX=3.d-16)
      EXTERNAL func
C USES func
      INTEGER i
      double precision a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,
     *     temp,test,tmplam
      check=.false.
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.d0
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x,y)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              tmplam=(-b+sqrt(disc))/(3.d0*a)
            endif
            if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1d0*alam)
      goto 1
      END
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      double precision a(np,np),b(n)
      INTEGER i,ii,j,ll
      double precision sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        imax=1
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
      SUBROUTINE newt2d(x,y,check)
      IMPLICIT NONE
      INTEGER MAXITS
      LOGICAL check
      double precision x(2),y(2),TOLF,TOLMIN,TOLX,STPMX
      PARAMETER (MAXITS=200,TOLF=1.d-12,TOLMIN=1.d-12,TOLX=3.d-16,
     &     STPMX=100.d0)
C USES fdjac,fmin,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,indx(2)
      double precision d,den,f,fold,stpmax,sum,temp,test,fvec(2),fjac(2,2),
     &     g(2),p(2),xold(2),fmin
      COMMON /newtv/ fvec,fjac
      EXTERNAL fmin
      f=fmin(x,y)
      test=0.d0
      do i=1,2
        if(abs(fvec(i)) > test)test=abs(fvec(i))
      end do
      if(test < .01d0*TOLF) return
      sum=0.d0
      do i=1,2
        sum=sum+x(i)**2
      end do
      stpmax=STPMX*max(sqrt(sum), dble(2))
      do its=1,MAXITS
        do i=1,2
          sum=0.d0
          do j=1,2
            sum=sum+fjac(j,i)*fvec(j)
          end do
          g(i)=sum
        end do
        do i=1,2
          xold(i)=x(i)
        end do  
        fold=f
        do i=1,2
          p(i)=-fvec(i)
        end do
        call ludcmp(fjac,2,2,indx,d)
        call lubksb(fjac,2,2,indx,p)
        call lnsrch(2,xold,fold,g,p,x,y,f,stpmax,check,fmin)
        test=0.d0
        do i=1,2
          if(abs(fvec(i)) > test)test=abs(fvec(i))
        end do
        if(test < TOLF)then
          check=.false.
          return
        endif
        if(check)then
          test=0.d0
          den=max(f,.5d0*2)
          do i=1,2
            temp=abs(g(i))*max(abs(x(i)),1.d0)/den
            if(temp > test)test=temp
          end do
          if(test < TOLMIN)then
            check=.true.
          else
            check=.false.
          endif
          return
        endif
        test=0.d0
        do i=1,2
          temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
          if(temp > test)test=temp
        end do
        if(test < TOLX) then
           return
        end if
         
      end do
      !pause 'MAXITS exceeded in newt'
      END

      FUNCTION fmin(x,y)
      IMPLICIT NONE
      double precision fmin,x(2),y(2),fvec(2),fjac(2,2)
      COMMON /newtv/ fvec,fjac
C USES funcv
      INTEGER i
      double precision sum
      call funcv(x,y,fvec,fjac)
      sum=0.d0
      do i=1,2
        sum=sum+fvec(i)**2
      end do
      fmin=sum
      return
      END

      SUBROUTINE FUNCV (X,Y,FVEC,FJAC)
      IMPLICIT double precision (A-H, O-Z)
      double precision X(2),Y(2),FJAC(2,2),FVEC(2)
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR

      CALL STATEF(X(1),X(2))
      
      FVEC(1) = AP - Y(1)
      FVEC(2) = ARHO - Y(2)

      FJAC(1,1) = PF
      FJAC(1,2) = PT
      FJAC(2,1) = RF
      FJAC(2,2) = RT
      RETURN
      END
