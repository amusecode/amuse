      SUBROUTINE update_loc_p
      IMPLICIT NONE
#include "general.inc"


C---  local vars
      INTEGER i,j,k

C---  update local particles
      DO i=1,n_act_loc

         j = ind_act_loc(i)           ! get the local index of i-th active particle
         DO k=1,NDIM
            x_loc(k,j)    = x_act(k,i)
            v_loc(k,j)    = v_act(k,i)
            a_loc(k,j)    = a_act(k,i)
            adot_loc(k,j) = adot_act(k,i)
         ENDDO
         pot_loc(j)       = pot_act(i)
         t_loc(j)         = t_act(i)
         dt_loc(j)        = dt_act(i)


      ENDDO


      RETURN
      END


