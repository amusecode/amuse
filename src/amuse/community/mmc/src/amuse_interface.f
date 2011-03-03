      function commit_particles_src()
      include 'common.h'
      
      call sort2(nt, r, iname)

      end function

      function recommit_particles_src()
      include 'common.h'

      call sort2(nt, r, iname)

      end function

      function get_number_of_particles_src(n_)
      include 'common.h'
      integer get_number_of_particles_src
      integer n_
      
      n_ = NT
      get_number_of_particles_src = 0 
      end function

      function get_time_src(time_)
      include 'common.h'
      integer get_time_src
      double precision timet_
      time_ = time
      get_time_src = 0 
      end function

      function get_timet_src(timet_)
      include 'common.h'
      integer get_timet_src
      double precision timet_
      timet_ = timet
      get_timet_src = 0 
      end function

      function get_crossing_time_src(tcr_)
      include 'common.h'
      integer get_crossing_time_src
      double precision tcr_
      tcr_ = tcr
      get_timet_src = 0 
      end function

      function get_state_src(id, mass_, r_, vr_, vt_)
      include 'common.h'
      integer get_state_src
      integer id
      double precision mass_, r_, vr_, vt_

      mass_ = BODY(id)

      r_ = R(id)
      vr_ = VR(id)
      vt_ = VT(id)
      get_state_src = 0
      end function get_state_src

      function set_state_src(id, mass_, r_, vr_, vt_)
      include 'common.h'
      integer set_state_src
      integer id
      double precision mass_, r_, vr_, vt_

      BODY(id) = mass_

      R(id) = r_ 
      VR(id) = vr_
      VT(id) = vt_
      set_state_src = 0
      end function set_state_src

      function get_total_kinetic_energy_src(T)
      include 'common.h'
      integer get_total_kinetic_energy
      double precision T
      T = ZKIN
      get_total_kinetic_energy = 0
      end function

      function get_total_potential_energy_src(V)
      include 'common.h'
      integer get_total_potential_energy_src
      double precision V
      V = -1
      get_total_potential_energy_src = 0
      end function

      function set_flagr(temp_, N)
      include 'common.h'
      integer, intent(in)::N
      double precision temp_(N)
      integer set_flagr, i
      do i = 1, N
          flagr(i) = temp_(i)
      end do
      set_flagr = 0
      end function

      function get_flagr(id, flagr_)
      include 'common.h'
      integer id
      double precision flagr_
      integer get_flagr
      flagr_ = flagr(id)
      get_flagr = 0
      end function

      function get_nlagra(temp_)
      include 'common.h'
      integer get_nlagra
      integer temp_
      temp_ = nlagra
      get_nlagra = 0
      end function

      function test_sort(n, aa, bb)
      integer n
      real*8 aa
      integer bb
      integer test_sort
      dimension aa(n), bb(n)

      call sort2(n, aa, bb)

      test_sort = 0
      end function 
