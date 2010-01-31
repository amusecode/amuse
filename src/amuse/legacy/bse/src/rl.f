***
      REAL*8 FUNCTION RL(Q)
      IMPLICIT NONE
      REAL*8 Q,P
*
* A function to evaluate R_L/a(q), Eggleton 1983.
*
      P = Q**(1.d0/3.d0)
      RL = 0.49d0*P*P/(0.6d0*P*P + LOG(1.d0+P))
*
      RETURN
      END
***
