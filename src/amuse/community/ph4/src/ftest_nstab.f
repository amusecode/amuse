      program ftest_nstab

c     Read the command line and pass the first 7 arguments to nstab.

      real*8 sigma,ei0,eo,relinc,m1,m2,m3
      character*32 arg

      if (iargc() .ge. 7) then
          call getarg(1, arg)
          read(arg, *) sigma
          call getarg(2, arg)
          read(arg, *) ei0
          call getarg(3, arg)
          read(arg, *) eo
          call getarg(4, arg)
          read(arg, *) relinc
          call getarg(5, arg)
          read(arg, *) m1
          call getarg(6, arg)
          read(arg, *) m2
          call getarg(7, arg)
          read(arg, *) m3
          print *, nstab(sigma,ei0,eo,relinc,m1,m2,m3)
      endif

      end
