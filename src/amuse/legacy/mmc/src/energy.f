      subroutine energy(ienerg)
*
*
*       total energy.
*       -------------
*
      include 'common.h'
*
      real*8 zzz
      integer i,n,ienerg,im
*
*
*     ----------------------------------------------------------
*
*         ienerg = 1    energies calculated from x and xdot
*         ienerg = 2    energies calculated from r, vr and vt
*
*     ----------------------------------------------------------
*
      n=nt
      zzz = 0.d0
*
*       calculate the potential and kinetic energy
*
*   
      if(ienerg.eq.1) then
*
        call coepot
*
        zkin = 0.0d0
        pot = 0.0d0
*
        do 20 i=1,n
           im = iname(i)
*
*      sum the potential energy
* 
           pot = pot - 0.5d0*body(im)*u(i)
*
*      sum the kinetic energy (include c.m. bodies but not components)
*
           zkin = zkin + body(im)*(xdot(im,1)**2 + xdot(im,2)**2 + 
     &            xdot(im,3)**2)
*
   20   continue
*
        zkin = 0.5d0*zkin
*
        return
*
      endif
*
*
      if(ienerg.eq.2) then
*
        if(time.eq.0.0d0) call coepot
*
        zkin = 0.0d0
        pot = 0.0d0
*
        do 30 i=1,n
           im = iname(i)
           zzz = zzz + body(im)
*
*      sum the potential energy
* 
           pot = pot - 0.5d0*body(im)*u(i)
*
*      sum the kinetic energy (include c.m. bodies but not components)
*
           zkin = zkin + body(im)*(vr(im)**2 + vt(im)**2)
*
   30   continue
*
      endif
*
      zkin = 0.5*zkin
*
        print*,'energy smt,zkin,pot,n,nt =',zzz,zkin,pot,n,nt
*
      return
*
      end
*
*
*       total energy = zkin - pot + etide + ebin + esub + emerge + ecoll.
*
*
*
*
*
