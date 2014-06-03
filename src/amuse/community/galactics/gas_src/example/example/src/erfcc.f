      FUNCTION erfcc(x)
      REAL erfcc,x
      REAL t,z
      z=abs(x)
      t=1./(1.+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*
     *(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
     *(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.-erfcc
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
