      integer*4 ibuf(15), ibuf1(15)
      integer*4 jbuf(2), jbuf1(2)
      integer*4 kbuf(2), kbuf1(2)
c this allows the gparameters to be passed as a C structure
      equivalence (chalo, ibuf)
      equivalence (psic, jbuf)
      equivalence (fcut_halo, kbuf)
