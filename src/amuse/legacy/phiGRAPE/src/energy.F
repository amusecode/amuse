      SUBROUTINE energy
C=========================================================================
C===
C===  task: compute energies, angular momentum and center of mass
C===  -----
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"
C---  local vars
      INTEGER          i,k,k3
      DOUBLE PRECISION mcm,dttmp,dt2half,dt3over6,
     $                 xx,xy,xz,vx,vy,vz
      LOGICAL          ltest

      ltest = .FALSE.
      dttmp = 0.
C---  predict to current time
C---  (for all particles, dt may be zero in which case current pos
C---  and velocities are just copied)
      DO i=1,N
         
         dttmp    = MAX(ABS(time_cur - t(i)),dttmp)
         dt2half  = 0.5d0 * dttmp * dttmp
         dt3over6 = dttmp * dt2half / 3.d0

         DO k=1,3
            k3 = k + 3
            body_pred_loc(k,i)  = x(k,i) + v(k,i)*dttmp  
     $                                   + a(k,i)*dt2half
     $                                   + adot(k,i)*dt3over6
            body_pred_loc(k3,i) = v(k,i) + a(k,i)*dttmp
     $                                   + adot(k,i)*dt2half

         ENDDO


C---  in case we need to update potential ind_act is needed
         ind_act(i) = ind_loc(i)

      ENDDO

      ltest      = ((dttmp .NE. 0.d0) .OR. ltest)

C---  an extra force calculation to get the predicted potential
      IF (ltest) THEN
         print*,'hellow1',dttmp
         n_act = N

         CALL gravity(time_cur,0)
C             -------

C---  copy potential only 
         DO i=1,N
            pot(i) = body_act_new(ipnew,i)
         ENDDO

      ENDIF


C---  get potential and kinetic energy
      E_pot = 0.d0
      E_kin = 0.d0
      DO i=1,N
         vx    = body_pred_loc(4,i)
         vy    = body_pred_loc(5,i)
         vz    = body_pred_loc(6,i)
         E_pot = E_pot + m(i)*pot(i)
         E_kin = E_kin + m(i)*( vx*vx + vy*vy + vz*vz )
      ENDDO
      E_pot = 0.5d0 * E_pot
      E_kin = 0.5d0 * E_kin
      E_tot = E_kin + E_pot


C---  get center of mass
      mcm = 0.d0
      DO k=1,NDIM
         xcm(k) = 0.d0
         vcm(k) = 0.d0
      ENDDO

      DO i=1,N
         mcm    = mcm + m(i)
         xx     = body_pred_loc(1,i)
         xy     = body_pred_loc(2,i)
         xz     = body_pred_loc(3,i)
         vx     = body_pred_loc(4,i)
         vy     = body_pred_loc(5,i)
         vz     = body_pred_loc(6,i)
         xcm(1) = xcm(1) + m(i)*xx
         xcm(2) = xcm(2) + m(i)*xy
         xcm(3) = xcm(3) + m(i)*xz
         vcm(1) = vcm(1) + m(i)*vx
         vcm(2) = vcm(2) + m(i)*vy
         vcm(3) = vcm(3) + m(i)*vz
      ENDDO

      DO k=1,NDIM
         xcm(k) = xcm(k)/mcm
         vcm(k) = vcm(k)/mcm
      ENDDO
      MASS = mcm

C---  get angular momentum
      DO k=1,NDIM
         mom(k) = 0.d0
      ENDDO
      
      DO i=1,N
         xx     = body_pred_loc(1,i)
         xy     = body_pred_loc(2,i)
         xz     = body_pred_loc(3,i)
         vx     = body_pred_loc(4,i)
         vy     = body_pred_loc(5,i)
         vz     = body_pred_loc(6,i)
         mom(1) = mom(1) + m(i) * ( xy*vz - xz*vy )
         mom(2) = mom(2) + m(i) * ( xz*vx - xx*vz )
         mom(3) = mom(3) + m(i) * ( xx*vy - xy*vx )
      ENDDO
     
      inrg = 1

       
      RETURN
      END
