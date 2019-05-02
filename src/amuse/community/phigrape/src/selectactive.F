      SUBROUTINE selectactive
C=========================================================================
C===
C===  task:   select all local active particles, ie. all local particles
C===  -----   that need a force-update
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      INTEGER i,k



C---  loop over local particles, select and count active particles
      n_act_loc = 0
      DO i=1,n_loc

         IF ((t_loc(i)+dt_loc(i)).EQ.min_t) THEN  ! set particle active

            n_act_loc                = n_act_loc + 1
            ind_act_loc(n_act_loc)   = i              ! remember local index
            ind_act(n_act_loc)       = ind_loc(i)     ! and global index]
            m_act(n_act_loc)         = m_loc(i)
            h_act(n_act_loc)         = h_loc(i)
            DO k=1,NDIM
               x_act(k,n_act_loc)    = x_loc(k,i)
               v_act(k,n_act_loc)    = v_loc(k,i)
               a_act(k,n_act_loc)    = a_loc(k,i)    
               adot_act(k,n_act_loc) = adot_loc(k,i)    
            ENDDO
            pot_act(n_act_loc)       = pot_loc(i)
            t_act(n_act_loc)         = t_loc(i)
            dt_act(n_act_loc)        = dt_loc(i)
           
         ENDIF
      ENDDO

#ifdef DEBUG
      
      WRITE(STDERR,500) myRank,n_act_loc
 500  FORMAT('PE',I2.2,': selected ',I8,' active particles')
#endif

      RETURN
      END
