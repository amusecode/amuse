! prtoft:
!   Convert pressure and density pair to f and T pair.
!   Uses numerical recipes routine newt2d (to be replaced by something
!   else for license and efficiency reasons - it's only a 2D matrix, we
!   don't need LU decomposition to solve the system of equations).

subroutine prtoft (logP,logRho,logF,logT, xa)
   use real_kind
   use constants
   
   implicit none
   real(double), intent(out) :: logT, logF
   real(double), intent(in) :: logP, logRho, xa(9)
   real(double) :: neo, nio
   real(double) :: x(2),y(2)
   logical :: check
   
   
   ! Initial guesses for F, T assuming ionized non-degenerate ideal gas
   nio = xa(1) + (42.0*xa(2) + 14.0*xa(5) + 12.0*xa(4) + 10.5*xa(3) +&
        8.4*xa(6) + 7.0*xa(7) + 6.0*xa(8) + 3.0*xa(9))/168.0
   neo = 0.5*(1.0 + xa(1))
   
   logT = logP - logRho - log(cr*(nio+neo))
   logF = (logRho - 1.5*logT + log(0.5*cen*neo)) - 2.0*(1.0-log(2.0))
   ! Find F, T by Newton-Raphson iteration
   x(1) = logF
   x(2) = logT
   y(1) = logP
   y(2) = logRho
   call newt2d(x,y,check, xa)
   !IF(CHECK) WRITE(*,*) 'Failed to find global minimum...'
   logF = x(1)
   logT = x(2)
end subroutine prtoft




subroutine lnsrch(n,xold,fold,g,p,x,y,f,stpmax,check,func, fvec, fjac)
   use real_kind
   
   implicit none
   integer :: n
   logical :: check
   real(double) :: f,fold,stpmax,g(n),p(n),x(n),y(n),xold(n),func,alf,tolx
   real(double) :: fvec(n),fjac(n)
   parameter (alf=1.d-4,tolx=3.d-16)
   external func
   ! USES func
   integer :: i
   real(double) :: a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam
   
   check=.false.
   sum=0.d0
   do i=1,n
      sum=sum+p(i)*p(i)
   end do
   sum=sqrt(sum)
   if(sum.gt.stpmax) then
      do i=1,n
         p(i)=p(i)*stpmax/sum
      end do
   end if
   slope=0.d0
   do i=1,n
      slope=slope+g(i)*p(i)
   end do
   test = 0.d0
   do i=1,n
      temp=abs(p(i))/max(abs(xold(i)),1.d0)
      if(temp.gt.test)test=temp
   end do
   alamin=tolx/test
   alam=1.d0

1  continue  !FIXME Replace with 'do'
   do i=1,n
      x(i)=xold(i)+alam*p(i)
   end do
   f=func(x,y, fvec, fjac)
   if(alam.lt.alamin)then
      do i=1,n
         x(i)=xold(i)
      end do
      check=.true.
      return
   else if(f.le.fold+alf*alam*slope)then
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
         end if
         if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
      end if
   end if
   alam2=alam
   f2=f
   fold2=fold
   alam=max(tmplam,.1d0*alam)
   goto 1  !FIXME Replace with 'end do'
end subroutine lnsrch




subroutine lubksb(a,n,np,indx,b)
   use real_kind
   
   implicit none
   integer :: n,np,indx(n)
   real(double) :: a(np,np),b(n)
   integer :: i,ii,j,ll
   real(double) :: sum
   
   ii=0
   do i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0)then
         do j=ii,i-1
            sum=sum-a(i,j)*b(j)
         end do
      else if (sum.ne.0.) then
         ii=i
      end if
      b(i)=sum
   end do
   do i=n,1,-1
      sum=b(i)
      do j=i+1,n
         sum=sum-a(i,j)*b(j)
      end do
      b(i)=sum/a(i,i)
   end do
end subroutine lubksb




subroutine ludcmp(a,n,np,indx,d)
   use real_kind
   
   implicit none
   integer :: n,np,indx(n),nmax
   real(double) :: d,a(np,np),tiny
   parameter (nmax=500,tiny=1.0e-20)
   integer :: i,imax,j,k
   real(double) :: aamax,dum,sum,vv(nmax)
   
   d=1.
   do i=1,n
      aamax=0.
      do j=1,n
         if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
      end do
      !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
      vv(i)=1./aamax
   end do
   do j=1,n
      do i=1,j-1
         sum=a(i,j)
         do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
         end do
         a(i,j)=sum
      end do
      aamax=0.
      imax=1
      do i=j,n
         sum=a(i,j)
         do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
         end do
         a(i,j)=sum
         dum=vv(i)*abs(sum)
         if (dum.ge.aamax) then
            imax=i
            aamax=dum
         end if
      end do
      if (j.ne.imax)then
         do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
         end do
         d=-d
         vv(imax)=vv(j)
      end if
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=tiny
      if(j.ne.n) then
         dum=1./a(j,j)
         do i=j+1,n
            a(i,j)=a(i,j)*dum
         end do
      end if
   end do
end subroutine ludcmp



function fmin(x,y,fvec, fjac, xa)
   use real_kind
   
   implicit none
   real(double) :: fmin,x(2),y(2),fvec(2),fjac(2,2), xa(9)
   ! USES funcv
   integer :: i
   real(double) :: sum
   
   call funcv(x,y,fvec,fjac, xa)
   sum=0.d0
   do i=1,2
      sum=sum+fvec(i)**2
   end do
   fmin=sum
end function fmin




subroutine funcv (x,y,fvec,fjac, xa)
   use real_kind
   use eostate_types
   
   implicit none
   real(double) :: x(2),y(2),xa(9),fvec(2),fjac(2,2)
   type(eostate) :: eos
   type(abundance) :: abund

   call statef(x(1),x(2), xa(9), abund, eos)
   
   fvec(1) = eos%ap - y(1)
   fvec(2) = eos%arho - y(2)
   
   fjac(1,1) = eos%pf
   fjac(1,2) = eos%pt
   fjac(2,1) = eos%rf
   fjac(2,2) = eos%rt
end subroutine funcv



subroutine newt2d(x,y,check, xa)
   use real_kind
   
   implicit none
   integer :: maxits
   logical :: check
   real(double) :: x(2),y(2),xa(9), tolf,tolmin,tolx,stpmx
   parameter (maxits=200,tolf=1.d-12,tolmin=1.d-12,tolx=3.d-16,stpmx=100.d0)
   ! USES fdjac,fmin,lnsrch,lubksb,ludcmp
   integer :: i,its,j,indx(2)
   real(double) :: d,den,f,fold,stpmax,sum,temp,test,fvec(2),fjac(2,2),g(2),p(2),xold(2),fmin
   
   f=fmin(x,y,fvec,fjac, xa)
   test=0.d0
   do i=1,2
      if(abs(fvec(i)) > test)test=abs(fvec(i))
   end do
   if(test < .01d0*tolf) return
   sum=0.d0
   do i=1,2
      sum=sum+x(i)**2
   end do
   stpmax=stpmx*max(sqrt(sum), dble(2))
   do its=1,maxits
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
      call lnsrch(2,xold,fold,g,p,x,y,f,stpmax,check,fmin,fvec,fjac)
      test=0.d0
      do i=1,2
         if(abs(fvec(i)) > test)test=abs(fvec(i))
      end do
      if(test < tolf)then
         check=.false.
         return
      end if
      if(check)then
         test=0.d0
         den=max(f,.5d0*2)
         do i=1,2
            temp=abs(g(i))*max(abs(x(i)),1.d0)/den
            if(temp > test)test=temp
         end do
         if(test < tolmin)then
            check=.true.
         else
            check=.false.
         end if
         return
      end if
      test=0.d0
      do i=1,2
         temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
         if(temp > test)test=temp
      end do
      if(test < tolx) then
         return
      end if
      
   end do
   !pause 'MAXITS exceeded in newt'
end subroutine newt2d



