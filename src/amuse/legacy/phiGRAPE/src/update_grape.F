      SUBROUTINE update_grape(imode)
C=========================================================================
C===  update particles on GRAPE
C===  imode = 0: active particles only
C===          1: all local particles
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      DOUBLE PRECISION a2by18(NDIM),a1by6(NDIM),aby2(NDIM)
      INTEGER i,j,k,ii,jj,imode

C---  get first index of local active particle
      ii = disp_act_v(myRank+1)

      IF (imode.NE.1) THEN

         DO j=1,n_act_loc

            i  = ind_act_loc(j) ! local index of j-th active particle
            jj = ii+j           ! index of same particle in global body_act

            DO k=1,NDIM

               a2by18(k) = 0.d0 
               a1by6(k)  = adot_act(k,j) * over6 
               aby2(k)   = a_act(k,j)    * over2

            ENDDO

            if (i .gt. nbh) then 

               CALL g6_set_j_particle(clusterid, i-1-nbh, ind_act(j),
     $              t_act(j),
     $              dt_act(j),m_act(j),a2by18,a1by6,aby2,
     $              v_act(1,j), x_act(1,j))
            end if
         
         ENDDO

      ELSE

C---  sending all local particles
         DO j=nbh+1,n_loc
            
            DO k=1,NDIM
               a2by18(k) = 0.d0 
               a1by6(k)  = adot_loc(k,j)*over6 
               aby2(k)   = a_loc(k,j)*over2
            ENDDO

            CALL g6_set_j_particle(clusterid, j-1-nbh, ind_loc(j), 
     $           t_loc(j),
     $                             dt_loc(j), m_loc(j), a2by18, a1by6,
     $                             aby2, v_loc(1,j), x_loc(1,j))
C                -----------------
         ENDDO

      ENDIF


      RETURN
      END
