      SUBROUTINE gravity(tnow,ifirst)
C===  calc the grav for the active particles in first step
      IMPLICIT NONE
#include "general.inc"
#include "stopcond.inc"

      INTEGER is_stopping_condition_enabled
      INTEGER is_collision_detection_enabled
      INTEGER is_pair_detection_enabled
      INTEGER next_index_for_stopping_condition
      INTEGER set_stopping_condition_info
      INTEGER set_stopping_condition_particle_index
      INTEGER stopping_index
      INTEGER reset_stopping_conditions
      INTEGER error

C---  input vars
      INTEGER          ifirst
      DOUBLE PRECISION tnow

C---  GRAPE vars
      INTEGER    NGP
      PARAMETER (NGP=16384)

      INTEGER g6err, g6calc_lasthalf, g6calc_lasthalf2

C---  local vars
      INTEGER ni,nn,i,ii,iii,ig,k,index_i(NGP),nn_i(NGP),id1,id2,j,nj,
     $        ip1,ip2
      DOUBLE PRECISION x_i(NDIM,NGP),v_i(NDIM,NGP),a_i(NDIM,NGP),
     $                 jerk_i(NDIM,NGP),p_i(NGP),h2_i(NGP)
      DOUBLE PRECISION hsum2,r2,r(NDIM),dttmp,dt2half,dt3over6
      DOUBLE PRECISION xp, yp, zp, vxp, vyp, vzp, dtp, dtp2, dtp3
      DOUBLE PRECISION ax, ay, az, jx, jy, jz, potj
      DOUBLE PRECISION xi, yi, zi, vxi, vyi, vzi
      DOUBLE PRECISION dx, dy, dz, dvx, dvy, dvz
      DOUBLE PRECISION rr2, rv, rinv1, rinv2, rinv3, alpha
      DOUBLE PRECISION dv2

C---  set some initial guess for first call of GRAPES
      IF (ifirst.EQ.1) THEN
         DO i=1,n_act

            ig = ind_act(i)

            pot_act_tmp_loc(ig)       = -1.d0
            DO k=1,NDIM
               a_act_tmp_loc(k,ig)    = 1.d0
               adot_act_tmp_loc(k,ig) = 10.d0
            ENDDO
         ENDDO
      ENDIF

C---  set actual time on GRAPE
      CALL g6_set_ti(clusterid, tnow)

C---  n_act is all active particles and on all nodes the same
      ni = n_act

      error = is_stopping_condition_enabled(COLLISION_DETECTION,
     $        is_collision_detection_enabled)
      error = is_stopping_condition_enabled(PAIR_DETECTION,
     $        is_pair_detection_enabled)

C---  define the local phi, a, adot for these active particles
      DO i=0,ni-1,npipe
         nn = npipe                            ! block size 48 or smaller
         IF ((ni - i).LT.npipe) nn = ni - i
         DO ii=1,nn

            iii = i+ii
        
            ig  = ind_act(iii)
            index_i(ii) = ig
            h2_i(ii)    = eps2

            DO k=1,NDIM

               x_i(k,ii)    = body_pred(k,iii)
               v_i(k,ii)    = body_pred(k+3,iii)
               a_i(k,ii)    = a_act_tmp_loc(k,ig) ! some initial guess for a and jerk
               jerk_i(k,ii) = adot_act_tmp_loc(k,ig) ! using partial force from previous
                                                     ! call here
            ENDDO
            p_i(ii)         = pot_act_tmp_loc(ig)
         ENDDO

         IF (ifirst.EQ.1) THEN
C---  get forces etc. based on bad initial guess
           CALL g6calc_firsthalf(clusterid, n_loc-nbh, nn,
     $           index_i,x_i,v_i,
     $           a_i, jerk_i, p_i, eps2, h2_i)
           g6err = g6calc_lasthalf(clusterid, n_loc-nbh,nn,
     $          index_i,x_i,v_i,
     $          eps2, h2_i, a_i, jerk_i, p_i)


           DO ii=1,nn       ! test for bad jerk and reset
              IF (jerk_i(1,ii).EQ.0.d0) jerk_i(1,ii) = 1.d-5
              IF (jerk_i(2,ii).EQ.0.d0) jerk_i(2,ii) = 1.d-5
              IF (jerk_i(3,ii).EQ.0.d0) jerk_i(3,ii) = 1.d-5
           ENDDO

        ENDIF


C---  now get forces based on a (hopefully) good guess from a previous force calculation
        CALL g6calc_firsthalf(clusterid, n_loc-nbh, nn,
     $       index_i, x_i, v_i,
     $       a_i, jerk_i, p_i, eps2, h2_i)
        g6err = g6calc_lasthalf2(clusterid,
     $       n_loc-nbh, nn, index_i, x_i, v_i,
     $       eps2, h2_i, a_i, jerk_i, p_i,
     $       nn_i)


        IF (is_collision_detection_enabled.GT.0
     $      .OR.is_pair_detection_enabled.GT.0) THEN 
            DO ii=1,nn

                iii   = ii+i
                ip1   = index_i(ii)
                ip2   = nn_i(ii)
                id1   = myid(ip1)
                id2   = myid(ip2)
                ip2   = ip2 - ind_loc(1) + 1
                hsum2 = (h_act(iii)+h_loc(ip2))**2
                
C---  get the predicted position of secondary particle which may
C---  or may not be an active particle
                dttmp    = tnow - t_loc(ip2)
                dt2half  = 0.5d0 * dttmp * dttmp
                dt3over6 = dttmp * dt2half / 3.d0
                DO k=1,NDIM
                    r(k) = x_loc(k,ip2) + v_loc(k,ip2)*dttmp
     $                  + a_loc(k,ip2)*dt2half
     $                  + adot_loc(k,ip2)*dt3over6
     $                  - x_i(k,ii)
                ENDDO
                
C---  now compute sqaure of radius
                r2    = (r(1)**2 + r(2)**2 + r(3)**2)

C---  test for collision + bound? (fip)
                IF (hsum2.GE.r2) THEN
                    dv2=0
         
                    icollision  = 1 ! collision detected

C ELLO-- implement this according top stopcond lib C
                    IF (m_act(iii).GE.m_loc(ip2)) THEN
                        idprimary   = id1
                        idsecondary = id2
                    ELSE
                        idprimary   = id2
                        idsecondary = id1
                    ENDIF
            
                    IF (is_collision_detection_enabled.GT.0) THEN
                        call amuse_handle_collisions(id1, id2)
                    ENDIF
                    IF (is_pair_detection_enabled.GT.0) THEN
                        DO k=1,NDIM
                            dv2=dv2+
     $                          (v_loc(k,ip2)+
     $                          a_loc(k,ip2)*dttmp-v_i(k,ii))**2
                        ENDDO
C--                        print*,(m_act(iii)+m_loc(ip2))/sqrt(r2+eps2),
C--     $                         0.5*dv2
                        IF((m_act(iii)+m_loc(ip2))/sqrt(r2+eps2)
     $                                    .GT.0.5*dv2) THEN
                            stopping_index = 
     $                      next_index_for_stopping_condition()

                            error = 
     $                      set_stopping_condition_info(stopping_index,
     $                              PAIR_DETECTION)
                    
                            error = 
     $                      set_stopping_condition_particle_index(
     $                              stopping_index,
     $                              0,
     $                              id1)
                        
                            error = 
     $                      set_stopping_condition_particle_index(
     $                              stopping_index,
     $                              1,
     $                              id2)
                        ENDIF
                    ENDIF
C--                    GOTO 110    
C-- there was a GOTO 110 to prevent double detection of same collision, 
C-- but that is handled by the HL code now!

                ENDIF
            ENDDO
        ENDIF
  110   CONTINUE                ! break after detection of a collision
C---  temporary save forces etc. (will be summed over all PE's later)
        DO ii=1,nn

           iii = ii + i
           body_act_tmp(ipnew,iii)         = p_i(ii)
           DO k=1,NDIM
              body_act_tmp(iaxnew-1+k,iii) = a_i(k,ii)
              body_act_tmp(ijxnew-1+k,iii) = jerk_i(k,ii)
           ENDDO

        ENDDO

      ENDDO            ! finished getting forces etc for all active particles

C---  Store the value of the local partial force etc...
C---  this is needed as an initial guess for the next force calculation
C---  the *_tmp_loc arrays have forces etc excerted only by local particles
C---  and can be accessed by the global particle index
      DO i=1,n_act
         ig = ind_act(i)

         pot_act_tmp_loc(ig)       = body_act_tmp(ipnew,i)
         DO k=1,NDIM
            a_act_tmp_loc(k,ig)    = body_act_tmp(iaxnew-1+k,i)
            adot_act_tmp_loc(k,ig) = body_act_tmp(ijxnew-1+k,i)
         ENDDO

      ENDDO


      do i=1,n_act
         ig = ind_act(i)

         xi  = body_pred(1, i)
         yi  = body_pred(2, i)
         zi  = body_pred(3, i)
         vxi = body_pred(4, i)
         vyi = body_pred(5, i)
         vzi = body_pred(6, i)

         if (ig .le. nbh) then
            nj = n_loc
            body_act_tmp(1, i) = 0
            body_act_tmp(2, i) = 0
            body_act_tmp(3, i) = 0
            body_act_tmp(4, i) = 0
            body_act_tmp(5, i) = 0
            body_act_tmp(6, i) = 0
            body_act_tmp(7, i) = 0
         else
            nj = nbh
         end if

         do j = 1,nj
            dtp = min_t - t_loc(j)
            dtp2 = dtp*dtp/2.0
            dtp3 = dtp2*dtp/3.0

            xp = x_loc(1, j) + v_loc(1, j)*dtp +
     $           a_loc(1,j)*dtp2 + adot_loc(1,j)*dtp3
            yp = x_loc(2, j) + v_loc(2, j)*dtp +
     $           a_loc(2,j)*dtp2 + adot_loc(2,j)*dtp3
            zp = x_loc(3, j) + v_loc(3, j)*dtp +
     $           a_loc(3,j)*dtp2 + adot_loc(3,j)*dtp3

            vxp = v_loc(1, j) + a_loc(1, j)*dtp + adot_loc(1,j)*dtp2
            vyp = v_loc(2, j) + a_loc(2, j)*dtp + adot_loc(2,j)*dtp2
            vzp = v_loc(3, j) + a_loc(3, j)*dtp + adot_loc(3,j)*dtp2

            dx = xp - xi
            dy = yp - yi
            dz = zp - zi
            rr2 = dx*dx + dy*dy + dz*dz
            dvx = vxp - vxi
            dvy = vyp - vyi
            dvz = vzp - vzi

            if (ig .ne. j) then
               rinv1 = 1.0/sqrt(rr2 + bh_eps)
            else
               rinv1 = 0
            end if
            rv    = dx*dvx + dy*dvy + dz*dvz
            rinv2 = rinv1*rinv1
            alpha = rv * rinv2;
            rinv1 = rinv1 * m_loc(j)
            rinv3 = rinv1*rinv2

            ax   = rinv3*dx
            ay   = rinv3*dy
            az   = rinv3*dz
            jx   = rinv3*dvx - 3.0*alpha*ax
            jy   = rinv3*dvy - 3.0*alpha*ay
            jz   = rinv3*dvz - 3.0*alpha*az
            potj = - rinv1

            body_act_tmp(1, i) = body_act_tmp(1,i) + ax
            body_act_tmp(2, i) = body_act_tmp(2,i) + ay
            body_act_tmp(3, i) = body_act_tmp(3,i) + az
            body_act_tmp(4, i) = body_act_tmp(4,i) + jx
            body_act_tmp(5, i) = body_act_tmp(5,i) + jy
            body_act_tmp(6, i) = body_act_tmp(6,i) + jz
            body_act_tmp(7, i) = body_act_tmp(7,i) + potj
         end do

      end do

      RETURN
      END

      subroutine amuse_handle_collisions(id1, id2)
        implicit none
#include "stopcond.inc"
        integer, intent(in) :: id1, id2
        integer :: next_index_for_stopping_condition
        integer :: set_stopping_condition_info
        integer :: set_stopping_condition_particle_index
        integer :: sc_idx, err
        
C       Make sure each collision is reported only once:
        if (id1.GT.id2) return
        
        sc_idx = next_index_for_stopping_condition()
        err = set_stopping_condition_info(sc_idx, COLLISION_DETECTION)
        if (err.LT.0) write(*,*) "Error in amuse_handle_collisions - 1"
        
        err = set_stopping_condition_particle_index(sc_idx, 0, id1)
        if (err.LT.0) write(*,*) "Error in amuse_handle_collisions - 2"
        
        err = set_stopping_condition_particle_index(sc_idx, 1, id2)
        if (err.LT.0) write(*,*) "Error in amuse_handle_collisions - 3"
      end subroutine amuse_handle_collisions


      SUBROUTINE gravity_at_point(tnow, pointx, pointy, pointz,
     $   fx, fy, fz, p, number_of_points)
C===  calc the gravity at a point, cannot use data on the point
C===  for previous timesteps, so is less acurate
      IMPLICIT NONE
#include "general.inc"

C---  input vars
      INTEGER, intent(IN) :: number_of_points
      DOUBLE PRECISION tnow
      DOUBLE PRECISION pointx(number_of_points)
      DOUBLE PRECISION pointy(number_of_points),pointz(number_of_points)
      DOUBLE PRECISION :: fx(number_of_points), fy(number_of_points)
      DOUBLE PRECISION :: fz(number_of_points), p(number_of_points)

C---  GRAPE vars
      INTEGER    NGP
      PARAMETER (NGP=16384)

      INTEGER g6err, g6calc_lasthalf, g6calc_lasthalf2

C---  local vars
      INTEGER k, j, index_i(NGP), nn, i
      DOUBLE PRECISION x_i(NDIM,NGP),v_i(NDIM,NGP),a_i(NDIM,NGP),
     $                 jerk_i(NDIM,NGP),p_i(NGP),h2_i(NGP)

C---  set actual time on GRAPE
      CALL g6_set_ti(clusterid, tnow)

      DO i=0,number_of_points-1,npipe

         nn = npipe                         ! block size 48 or smaller
         IF ((number_of_points - i).LT.npipe) nn = number_of_points - i

C---  define the local phi, a, adot for these active particles
          DO j=1, nn
              index_i(j) = -1
              x_i(1,j) = pointx(j+i)
              x_i(2,j) = pointy(j+i)
              x_i(3,j) = pointz(j+i)
              h2_i(j)  = eps2
              p_i(j) = -1.d0
              DO k=1,NDIM
                a_i(k,j) = 1.d0
                v_i(k,j) = 0.d0
                jerk_i(k,j) = 10.d0
              ENDDO
          ENDDO

C---  get forces etc. based on bad initial guess
      CALL g6calc_firsthalf(clusterid, n_loc-nbh, nn,
     $           index_i,x_i,v_i,
     $           a_i, jerk_i, p_i, eps2, h2_i)
      g6err = g6calc_lasthalf(clusterid, n_loc-nbh, nn,
     $          index_i,x_i,v_i,
     $          eps2, h2_i, a_i, jerk_i, p_i)




C---  now get forces based on a (hopefully) good guess from a previous
C---  force calculation
      CALL g6calc_firsthalf(clusterid, n_loc-nbh, nn,
     $       index_i, x_i, v_i,
     $       a_i, jerk_i, p_i, eps2, h2_i)
      g6err = g6calc_lasthalf(clusterid, n_loc-nbh, nn,
     $          index_i,x_i,v_i,
     $          eps2, h2_i, a_i, jerk_i, p_i)

          DO j=1, nn
            fx(j+i) = a_i(1,j)
            fy(j+i) = a_i(2,j)
            fz(j+i) = a_i(3,j)
            p(j+i) = p_i(j)
          ENDDO
      ENDDO
      RETURN
      END
