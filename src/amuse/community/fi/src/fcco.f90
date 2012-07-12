module fccoMod
 implicit none
 private
 public :: fccoInit,fcco,q1c,q1co

 real, parameter :: n0=1520, &    	! cm^-3
                    xifuv=3.      	! FUV absorption constant (2-4)
		    		    
 real, save :: metal
 integer, save :: poly

 contains

subroutine fccoInit
 use ElementsMod
 use MoleculeMod
 
 call InitElements()
 metal=zQ()/solarzQ()
 poly=polyQ()
 if(poly.NE.0.AND.poly.NE.1) call fccoError('poly error',poly)

end subroutine

function fcco(n,G0,Pe,Z)
 real :: fcco,n,G0,Pe,zz
 real, optional :: Z
 
 zz=metal
 if(present(Z)) zz=Z
 if(poly.EQ.0) then
  fcco=fccoco(n,G0,zz,Pe)
 else
  fcco=fccolog(n,G0,zz,Pe)
 endif
end function

function q1c(t)
 real :: q1c,t
 real, parameter :: e1=23.6,e2=62.4 ! K
 integer, parameter :: g0=1,g1=3,g2=5
 q1c=g1*exp(-e1/t)/(g0+g1*exp(-e1/t)+g2*exp(-e2/t)) 
end function

function q1co(t)
 real :: q1co,t
 real, parameter :: e1=5.5 ! K
 integer, parameter :: g1=3
 q1co=g1*exp(-e1/t)/(2*t/e1) 
end function

function e2solvelog(a1,a2) result(xnew)
 real :: a1,a2,x,xnew,en(0:2)
 integer :: i,nmax=100
 xnew=1.0001
 call ENXB(2,alog(xnew)*a1,en)
 if(en(2)-a2*xnew.LT.0) return
 i=0
10 i=i+1
  x=xnew
  call ENXB(2,alog(x)*a1,en)
  xnew=x+(en(2)-a2*x)/(en(1)*a1/x+a2)
  if(xnew.LT.1) then
   xnew=1.
   return
  endif
  if(abs(xnew-x).GT.1.e-5.AND.i.LT.nmax) goto 10  
  call enxb(2,alog(xnew)*a1,en)
  if(abs(en(2)-a2*xnew).GT.1.e-3) then 
   call fccoError('e2solvelog error',i)
  stop
 endif

end function

function e2solve(x)
 real :: x,thnew,th,e2solve
 integer :: i,nmax=100

thnew=1.d-9
if(x.gt.exp(-thnew)-thnew*E1(thnew)) then
 e2solve=thnew
 return
endif 
if(x.LE.1.d-40) then
 e2solve=87.6081992538658
 return
endif  
i=0
10 i=i+1
 th=thnew
 thnew=th+(exp(-th)-th*E1(th)-x)/E1(th)
 if(abs(thnew-th).GT.1.e-5.AND.i.LE.nmax) goto 10
 if((exp(-th)-th*E1(th)-x).GT.1.e-5) then 
  call fccoError(' e2solve error',i)
  stop
 endif
 e2solve=thnew 
end function

function fccolog(n,G0,Z,Pe) result(fcoc)
 real :: n,G0,Z,Pe,a1,a2,x,fcoc

 a1=1.11*Z*n0/1.d3*sqrt(pe/1.d4)*xifuv
 a2=(1+24.8*Z)*8./9.*n/1.d6/G0
 x=e2solvelog(a1,a2)
 fcoc=1/x**2
end function

function fccoco(n,G0,Z,Pe)
 real :: n,G0,Z,Pe,xia,fccoco

 xia=e2solve(REAL(4./3.*n/1.d6/G0*(1+24.8*Z))) 
 fccoco=MAX(0.d0,1-xiA*1.d3/n0/xifuv/1.66/Z/sqrt(Pe/1.d4))**3
end function

        SUBROUTINE ENXB(N,X,EN)
!
!       ===============================================
!       Purpose: Compute exponential integral En(x)
!       Input :  x --- Argument of En(x)
!                n --- Order of En(x)  (n = 0,1,2,...)
!       Output:  EN(n) --- En(x)
!       ===============================================
!
        REAL X,EN,R,S,RP,PS,ENS,T,T0,S0
        INTEGER N,J,K,L,M
        DIMENSION EN(0:N)
        IF (X.EQ.0.0) THEN
           EN(0)=1.0D+300
           EN(1)=1.0D+300
           DO 10 K=2,N
10            EN(K)=1.0D0/(K-1.0)
           RETURN
        ELSE IF (X.LE.1.0) THEN
           EN(0)=EXP(-X)/X
           DO 40 L=1,N
              RP=1.0D0
              DO 15 J=1,L-1
15               RP=-RP*X/J
              PS=-0.5772156649015328D0
              DO 20 M=1,L-1
20               PS=PS+1.0D0/M
              ENS=RP*(-LOG(X)+PS)
              S=0.0D0
              DO 30 M=0,20
                 IF (M.EQ.L-1) GO TO 30
                 R=1.0D0
                 DO 25 J=1,M
25                  R=-R*X/J
                 S=S+R/(M-L+1.0D0)
                 IF (ABS(S-S0).LT.ABS(S)*1.0D-15) GO TO 35
                 S0=S
30            CONTINUE
35            EN(L)=ENS-S
40         CONTINUE
        ELSE
           EN(0)=EXP(-X)/X
           M=15+INT(100.0/X)
           DO 50 L=1,N
              T0=0.0D0
              DO 45 K=M,1,-1
45               T0=(L+K-1.0D0)/(1.0D0+K/(X+T0))
              T=1.0D0/(X+T0)
50            EN(L)=EXP(-X)*T
        ENDIF
        END SUBROUTINE

        FUNCTION E1(X)
!
!       ============================================
!       Purpose: Compute exponential integral E1(x)
!       Input :  x  --- Argument of E1(x) 
!       Output:  E1 --- E1(x) ( x > 0 )
!       ============================================
!
        REAL X,E1,ES1,ES2
        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=-LOG(X)+((((1.07857D-3*X-9.76004D-3)*X+5.519968D-2)*X &
     &        -0.24991055D0)*X+0.99999193D0)*X-0.57721566D0
        ELSE
           ES1=(((X+8.5733287401D0)*X+18.059016973D0)*X &
     &         +8.6347608925D0)*X+0.2677737343D0
           ES2=(((X+9.5733223454D0)*X+25.6329561486D0)*X &
     &         +21.0996530827D0)*X+3.9584969228D0
           E1=EXP(-X)/X*ES1/ES2
        ENDIF
        RETURN
        END FUNCTION

 subroutine fccoError(string,i)
  character(*) :: string
  integer, optional :: i

  print*,' fccoError detected:'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif

  stop
 end subroutine

end module

subroutine fccotest
 use fccoMod
 real :: x,y,n,G0,Z,Pe
 
 call fccoInit
 
 10 read*,n,G0,Z,Pe
 if(n.LT.0) stop
 y=fcco(n,G0,Pe,Z)
 print*,n,G0,Z,Pe,y
 goto 10
end 

