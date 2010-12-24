      function init_sequence()
      include 'common.h'
      integer init_sequence
      integer iphase
      integer initial_run, res
      print*, 'initwrapper beeing called, now try init'
      call input
      print *,istart
      if (istart.eq.1) then
          call start
          iphase = 1
      else
          call mydump(2)
          call input
          iphase = 2
      endif

      res = initial_run(iphase)

      init_sequence = res
      print*,'leaving amuse code interface'
      end function

      function initial_run(iphase)
      integer evolve_step
      integer iphase
      integer initial_run
      integer res, amuse_output

      print *, iphase

      do
          if (iphase.eq.1) then
c             call output
c             we use our own 'cos the original stops, amuse doesnt
c             like that:
              res = amuse_output()
              print*, res
              if (res.lt.0) exit
          endif
          call zone
          call relaxt
          iphase = 1
      end do

      print*,'leaving init sequence'
      initial_run = res
      end function

      function parameter_test(times)
      implicit none
      integer parameter_test
      double precision times

c      times = TIME
      parameter_test = 0
      end function
