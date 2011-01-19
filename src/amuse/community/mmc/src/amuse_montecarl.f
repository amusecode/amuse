      function init_sequence(res)
      include 'common.h'
      integer init_sequence
      integer initial_run 
      integer res
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

      res = initial_run()

      init_sequence = res
      print*,'leaving amuse code interface'
      end function

      function initial_run(res)
      include 'common.h'
      integer evolve_step
      integer initial_run
      integer res, amuse_output

      if (iphase.ne.1) then
          call zone
          call relaxt
          iphase = 1
      endif


      if (iphase.eq.1) then
c     call output
c     we use our own 'cos the original stops, amuse doesnt
c     like that:
          res = amuse_output()
          print*, res
c     if (res.lt.0) exit
c     let python script decide when to stop...
          call zone
          call relaxt
          initial_run = res
      endif

      initial_run = res
      end function

      function parameter_test(times)
      include 'common.h'
      integer parameter_test
      double precision times

      times = TIME
      parameter_test = 0
      end function

      subroutine amuse_set_mmc_data_directory(datadir_)
        character(len=200) datadir
        common /AMUSE/ datadir
        character(len=200), intent(in) :: datadir_
        datadir=datadir_
      end subroutine

      function total_kinetic_energy(Ek)
      include 'common.h'
      integer total_kinetic_energy
      double precision Ek
      Ek = ZKIN
      total_kinetic_energy = 0
      end function

