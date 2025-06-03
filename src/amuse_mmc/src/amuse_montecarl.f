      function init_sequence()
      include 'common.h'
      integer init_sequence
      integer run_a_while
      integer iphase
      integer res

c     this is done by amuse_input call because it overwrites
c     AMUSE settings, so only explicit use please..
c     call input

      if (istart.eq.1) then
          write(6,*)   '  calling start..'
          call flush(6)

          call start
          iphase = 1
      else
          write(6,*)   '  not calling start..'
          call flush(6)

          call mydump(2)
c         hmmmm...
c         call input
          iphase = 2
      endif

c     only in evolve...
c     res = run_a_while(iphase)
      init_sequence = res
      end function

      function run_a_while(iphase)
      include 'common.h'
      integer evolve_step
      integer run_a_while
      integer iphase
      integer res, amuse_output

      if (iphase.ne.1) then
          call zone
          call relaxt
          iphase = 1
          res = -20
      endif

      if (iphase.eq.1) then
c     call output
c     we use our own 'cos the original stops, amuse doesnt
c     like that:
          res = amuse_output()
c     if (res.lt.0) exit
c     let python script decide when to stop...
          call zone
          call relaxt
      endif

      run_a_while = res
      end function

      function parameter_test(times)
      include 'common.h'
      integer parameter_test
      double precision times

      times = TIME
      parameter_test = 0
      end function

      function evolve_src(time_end)
      include 'common.h'
      integer evolve_src
      integer iphase
      integer res, dumpres
      double precision tcrit_, time_end

      res = 0

   10 if (time_end.gt.timet) then
          write(6,*) 'call zone'
          call zone
          write(6,*) 'call relaxt', timet
          call relaxt
          goto 10
      endif
      time_end = timet
      evolve_src = 0

      end function

      subroutine amuse_set_mmc_data_directory(datadir_)
        character(len=200) datadir
        common /AMUSE/ datadir
        character(len=200), intent(in) :: datadir_
        datadir=datadir_
      end subroutine

      subroutine amuse_input_src()
        include 'common.h'
        call input()
      end subroutine

