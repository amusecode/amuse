      SUBROUTINE corrector
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      INTEGER i,k,ii,power,kanew,kjnew
      DOUBLE PRECISION dt_tmp,dt3over6,dt4over24,dt5over120,
     $                 dtinv,dt2inv,dt3inv,
     $                 a2(NDIM),a3(NDIM),a1abs,adot1abs,
     $                 a2dot1(NDIM),a2dot1abs,a3dot1abs,dt_new

C---  functions
      DOUBLE PRECISION dmod


C---  get first index of local active particle
      ii = disp_act_v(myRank+1)

      DO i=1,n_act_loc

         dt_tmp     = dt_act(i)

         dt3over6   = dt_tmp*dt_tmp*dt_tmp/6.d0
         dt4over24  = dt3over6*dt_tmp/4.d0
         dt5over120 = dt4over24*dt_tmp/5.d0

         dtinv      = 1.d0/dt_tmp
         dt2inv     = dtinv*dtinv
         dt3inv     = dt2inv*dtinv


         DO k=1,NDIM
            kanew = iaxnew+k-1 
            kjnew = ijxnew+k-1 
            a2(k) = -  6.d0 * (   a_act(k,i)
     $                          - body_act_new(kanew,ii+i) )*dt2inv
     $              - (  4.d0 * adot_act(k,i)
     $                 + 2.d0 * body_act_new(kjnew,ii+i) )*dtinv
            a3(k) =   12.d0 * (   a_act(k,i)
     $                          - body_act_new(kanew,ii+i))*dt3inv
     $              +  6.d0 * (   adot_act(k,i)
     $                          + body_act_new(kjnew,ii+i))*dt2inv
         ENDDO

         DO k=1,NDIM
            body_pred(k,ii+i)   = body_pred(k,ii+i)   + dt4over24*a2(k)
     $                                                + dt5over120*a3(k)
            body_pred(k+3,ii+i) = body_pred(k+3,ii+i) + dt3over6*a2(k)
     $                                                + dt4over24*a3(k)
         ENDDO

         a1abs    = SQRT(   body_act_new(iaxnew,ii+i)
     $                     *body_act_new(iaxnew,ii+i)
     $                   +  body_act_new(iaynew,ii+i)
     $                     *body_act_new(iaynew,ii+i)
     $                   +  body_act_new(iaznew,ii+i)
     $                     *body_act_new(iaznew,ii+i) )
         adot1abs = SQRT(   body_act_new(ijxnew,ii+i)
     $                     *body_act_new(ijxnew,ii+i)
     $                   +  body_act_new(ijynew,ii+i)
     $                     *body_act_new(ijynew,ii+i)
     $                   +  body_act_new(ijznew,ii+i)
     $                     *body_act_new(ijznew,ii+i))

         DO k=1,NDIM
            a2dot1(k) = a2(k) + dt_tmp*a3(k)
         ENDDO

         a2dot1abs = SQRT( a2dot1(1)*a2dot1(1)
     $                    +a2dot1(2)*a2dot1(2)
     $                    +a2dot1(3)*a2dot1(3))
         a3dot1abs = SQRT( a3(1)*a3(1) + a3(2)*a3(2) + a3(3)*a3(3))

C---  new time step (see Makino&Aarseth, 1992, PASJ, 44, 141)
         dt_new    = SQRT(   eta * (a1abs*a2dot1abs+adot1abs*adot1abs)
     $                     / (adot1abs*a3dot1abs+a2dot1abs*a2dot1abs))


C---  the current time min_t must be an integer multiple of the new time step
C---  that means the new time step can be smaller than or equal to the old time step
C---  for bigger time steps the condition has to be checked and in case it is not 
C---  fulfilled the time step is unchanged (also the new time step is limited to twice
C---  the old time step)

C---  if the new time step is smaller than the old one but not smaller
C---  than dt_min make dt_new a power of two
         IF (dt_new.LT.dt_min) THEN
#ifdef SILENT
            WRITE(0,500) myRank,i+ii,dt_new,dt_min,time_cur
 500        FORMAT('PE',I2.2,': Warning! timestep for particle'
     $              ,I8,' smaller than dt_min',3E14.6)
#endif
            dt_new = dt_min
         ENDIF

         IF ((dt_new.LT.dt_tmp) .AND. (dt_new.GE.dt_min)) THEN
            power  = INT(LOG(dt_new)/LOG(2.d0)) - 1
            dt_tmp = 2.d0**power               
         ENDIF

C---  if the new time step is bigger check integer multiple criteria
C---  also limit to twice the old time step and dt_max
         IF (       (dt_new.GE.2.d0*dt_tmp)
     $        .AND. (dmod(min_t,2.d0*dt_tmp).EQ.0.d0)
     $        .AND. (2.d0*dt_tmp.LE.dt_max)        ) THEN
            dt_tmp = 2.d0*dt_tmp ! twice the old time step is power of two
         ENDIF

C---  set new time step and time for active particle i
         dt_act(i)  = dt_tmp  
         t_act(i)   = min_t


C---  update positions etc of active particle i
         DO k=1,NDIM
            x_act(k,i)    = body_pred(k,ii+i) 
            v_act(k,i)    = body_pred(k+3,ii+i) 
            a_act(k,i)    = body_act_new(iaxnew+k-1,ii+i) 
            adot_act(k,i) = body_act_new(ijxnew+k-1,ii+i)
         ENDDO
         pot_act(i) = body_act_new(ipnew,ii+i)

      ENDDO


      RETURN
      END


C***************************************************************************

      DOUBLE PRECISION FUNCTION dmod(x,y)
      IMPLICIT NONE

      DOUBLE PRECISION x,y

      dmod = (x/y - INT(x/y)) * y

      RETURN
      END
