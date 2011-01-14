      subroutine potent(ipot,k1,k2,ri,ui)
*
*
*         compute potential 'ui' for radius 'ri'
*         --------------------------------------
*
*
      include 'common.h'
*
      real*8 ri,ui,z,z1,z2,z3,z4,z5
*
      integer i,n,ipot,k1,k2
*
*
      z = 1.d0/ri
*
      if(ipot.eq.0) then
        n = nt
*
        call locate(r,n,k1,k2,ri,i)
*
        if(i.eq.n) then
          ui = -smt*z
          return
        endif
*
        if(i.eq.0) then
          ui = u(1)
          return
        endif
*
        z1 = 1.0d0/r(i)
        z2 = 1.0d0/r(i+1)
        z3 = z1 - z2
        z4 = (u(i+1) - u(i))/z3
        z5 = z1 - z
        ui = u(i) + z5*z4
        return
*
      else
*
        n = nto
*
        call locate(ro,n,k1,k2,ri,i)
*
        if(i.eq.n) then
          ui = -smto*z
          return
        endif
*
        if(i.eq.0) then
          ui = uo(1)
          return
        endif
*
        z1 = 1.0d0/ro(i)
        z2 = 1.0d0/ro(i+1)
        z3 = z1 - z2
        z4 = (uo(i+1) - uo(i))/z3
        z5 = z1 - z
        ui = uo(i) + z5*z4
        return
*
      endif
*        
*
      end
*
*
*
*    
      subroutine locate(xx,n,k1,k2,z,j)
*
      include 'params.h'
*
      real*8 xx(nmax),z
*
      integer n,k1,k2,j,jl,jm,ju
*      
      jl = k1 - 1
      ju = k2 + 1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(z.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl
      return
      end
*
*
*
*

