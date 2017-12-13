module BarAndSpiralsInterface
  !GALACTIC POTENTIAL. IT HAS BAR AND SPIRAL ARMS
  !THE MODELS ARE IN 2D OR 3D. 
  !THE TOTAL FORCE AND POTENTIAL  IS MEASURED
  !WITH RESPECT TO A COROTATING FRAME WITH THE BAR,
  !THEN POSITIONS AND VELOCITIES HERE ARE ASSUMED TO BE IN THE SAME FRAME
      IMPLICIT NONE
      real*8 :: omegs2, fourpi
      real*8 :: tin,totalmass, MB1
      real*8 :: A2,B2,C2,UA2,UB2,UC2,XK,UA2B2,UA2C2,UB2C2,SUA2C2
      real*8:: V000,V100,V010,V001,V110,V101,V011,V200,V020,V002,V111,V210,V201,V120,V102,V012,V300,V030,V003, V021
      real*8 :: CTE0,CTE1,CTE2
      integer :: niter,idim
      real*8 :: radd,aax,aay,aaz,pr
      real*8 :: rl2,aux,aux2, qq, at, ax,qqbul,qqbar,qqhalo
      real*8 :: phi,d2,axis_a,axis_b,axis_c,gmb, f,e
      real*8 :: baxis,baxis2,barden
      real(8) :: eps_sp, rho0_scaled

      double precision :: G= 1
      double precision :: pi=3.141592653589793d0
      double precision :: FI0=0.D0
      double precision :: us3=1.d0/3.d0 
      double precision :: gravc= 1
      double precision :: time= 0.D0

      ! this flag is for the computation of the circular velocity and the tidal tensor. 
      !Default input positions x,y,z: in inertial frame. Therefore:
      !If xflag = 1 -> x, y, z are defined in a corotating system 
      !If xflag =2  -> x, y, z are defined in an inertial system
      integer :: xflag=2

      !constants for bar (taken from allruns.dat )
      double precision :: bar_phase= 0
      double precision :: barmas= 431.03
      double precision :: aaxis= 3.13 ! in kpc
      double precision :: caxis= 0.000001 ! in kpc 
      double precision :: axirat= 0.32
      double precision :: omegs= 5. ! in 10km/s*kpc^-1. Counterclockwise system
      double precision :: nbt= 0
      double precision :: tgrow

      ! constants for spiral arms
      integer:: spiral_model= 0 ! 0 is TWA, 1 is C&G, 2 is Lepine
      double precision :: spiral_phase= 0
      integer :: m=2 !number of spiral arms
      double precision :: phi0= 0.D0 ! initial phase of the spiral wrt the bar. IN RADIANS!!!!!!
      double precision :: tani= 0.277 !tangent pitch angle
      double precision :: rsigma= 2.5 ! RSigma   kpc
      double precision :: ombra= 2. ! in 10km/s*kpc^-1 
      
      !constant for TWA model 
      double precision :: nsp= 100 ! N
      double precision :: rsp= 1.5 ! in kpc 
      double precision :: asp= 8.5 ! (10km/s)^2/kpc

      ! constants for the Cox and Gomez (2002) model
      real(8) :: rho0_sp= 2.5e+7 ! density amplitude [M_o/kpc^3], for CG02
      real(8) :: r0_sp= 8.0      ! fiducial radius where rho0_sp is measured, for CG02
      real(8) :: h_sp= 0.18       ! scale height [kpc]
    
      !these are for the logaritmic spiral. Unset in the interface
      double precision:: sigma= 1.5 ! width of the spiral arm in kpc
      double precision:: gamma= 0 !phase for the shape factor fm

      !For the Lepine nodel (superposition of SA):
      double precision :: asp2= 6.8 ! Amplitude other pattern. asp2= 0.8*asp
      double precision :: ombra2= 2. ! Pattern speed other SA pattern
      double precision :: tani2= -0.1227845 !pitch angle other pattern (-7 grades)
      double precision :: m2= 2 !number of spiral arms 
      double precision :: phi21= -3.49065850399 !-200 grades. Angle wrt primary SA pattern
      
      
      ! These constants are for transient spirals
      real*8, dimension(:), allocatable:: rand_numbers
      real*8, dimension(:), allocatable:: omegas_sp
      real*8, dimension(:), allocatable:: phi_sp
      integer :: ntot
      double precision:: sigma_s= 1.02269032206 ! Sp timelife. Default: 100 Myr
      double precision:: ts ! time for max phi of spiral arms
      double precision:: aspi= 1.D0 ! amplitude of the non transcient structure
      double precision:: t_sim= 47.0437548146 ! Default: 4600 Myr
      
      !constants for bulge (default values)
      double precision :: B1= 0.3873D0  
      double precision :: MB= 606

      ! constants for disk (default values)
      double precision  :: MD= 3690 
      double precision  :: AD2=5.3178D0 ! kpc
      double precision  :: BD2=0.25D0 !kpc

      !constants for halo (default values)
      double precision  :: MH= 4615 ! 2.32e7MSun
      double precision :: A3=12.D0 !kpc
      double precision :: CS1=2.02D0
      double precision :: CS2=1.02D0
      double precision  :: RMAX=100.D0 !kpc Unset in the interface
      
      double precision:: omega_sys
      double precision:: initial_phase

      !this is for setting the desired galactic model. Default: axi case
      logical:: spiral_contribution= .false.
      logical:: bar_contribution= .false.
      logical:: transient_spiral= .false.
      
    contains
      
      integer function initialize_code()
        ! When the parameters are set in the script, this function is called
        initialize_code = 0
      end function initialize_code
      
      integer function cleanup_code()
        cleanup_code = 0
      end function cleanup_code
      
      integer function commit_parameters()
      ! Putting the initialization of the parameters in here, I can call vel_circ

      ! Initializing parameters of bar . 
      fourpi= 4.d0*pi
      baxis = aaxis * axirat
      baxis2 = baxis ** 2
      barden=barmas * 0.5968310366D0/ (aaxis * baxis2)
      
      axis_a=aaxis
      axis_b=baxis
      !axis_c=baxis-5.d-2
      axis_c= caxis
      a2=axis_a*axis_a
      b2=axis_b*axis_b
      c2=axis_c*axis_c
      ua2=1.d0/a2
      ub2=1.d0/b2
      uc2=1.d0/c2
      gmb=gravc*barmas
      CTE0 = GMB*3.d0/4.d0
      CTE1 = GMB*15.d0/8.d0 ! when n=1
      CTE2 = GMB*105.d0/16.d0 !when n=2. See 2008_manos_phD_thesis
      UA2B2 = 1.d0/(A2 - B2)
      UA2C2 = 1.d0/(A2 - C2)
      UB2C2 = 1.d0/(B2 - C2)
      SUA2C2 = DSQRT(UA2C2)
      PHI = DASIN(DSQRT(1.D0 - C2*UA2))
      XK = DSQRT(UA2C2*(A2-B2))
      CALL ELINT(PHI,XK,F,E)
      D2 = 2.D0*DSQRT(UA2*UB2*UC2)
      
      V000 = 2.D0*F*SUA2C2
      V100 = 2.D0*(F-E)*UA2B2*SUA2C2
      V001 = (D2*B2 - 2.D0*E*SUA2C2)*UB2C2
      V010 = D2 - V100 - V001
      
      V110 = (V010 - V100)*UA2B2
      V101 = (V001 - V100)*UA2C2
      V011 = (V001 - V010)*UB2C2
      
      V200 = (D2*UA2 - V110 - V101)*US3
      V020 = (D2*UB2 - V110 - V011)*US3
      V002 = (D2*UC2 - V011 - V101)*US3
      
      V111 = (V011 - V110)*UA2C2
      V210 = (V110 - V200)*UA2B2
      V201 = (V101 - V200)*UA2C2
      V120 = (V020 - V110)*UA2B2
      V021 = (V011 - V020)*UB2C2
      V102 = (V002 - V101)*UA2C2
      V012 = (V002 - V011)*UB2C2
      
      V300 = (D2*UA2*UA2 - V210 - V201)*.2D0
      V030 = (D2*UB2*UB2 - V120 - V021)*.2D0
      V003 = (D2*UC2*UC2 - V102 - V012)*.2D0
      
      if (bar_contribution .eqv. .false.)then
         OMEGS=0
         BARMAS=0
      endif

      MB1= MB !mass of the bulge
      totalmass=MB

      write(*,*) "total mass", totalmass
      write(*,*)'bulge mass: ', MB1
      write(*,*)'bar mass: ', barmas
 
           
      TIN=0.D0
      IF(OMEGS.NE.0.D0)THEN
         TGROW=TIN+(NBT*2.D0*PI/abs(OMEGS))
      ELSE
         TGROW=0.D0
      ENDIF
      WRITE(*,*)'growth of the bar (Gyr): ',TGROW*9.78462D7/1.d9
      WRITE(*,*)'growth of the bar (ut): ',TGROW

      !!! scaled constants for the CG02 spiral amrs
      if (spiral_model == 1) then
         eps_sp = 1.d0/rsigma
         rho0_scaled = rho0_sp*exp(eps_sp*r0_sp)
!         write(*,*) 'rho', rho0_sp, rho0_scaled*2.32*1d-3
      endif

      
      !Generating different pattern speeds for transcient spirals.
      !In our scheme, the transicents will not be overlapped; so
      ! ts must be 4*sigma_s

      if (transient_spiral .eqv. .True.) then
         ts= 8*sigma_s
         ! generate number of spiral events:
         ntot= nint(t_sim/ts)
         ! generate omega sp and orientations:
         allocate(rand_numbers(ntot))
         allocate(omegas_sp(ntot))
         allocate(phi_sp(ntot))
         call random_number(rand_numbers) !generating ntot random numbers between [0-1]
         omegas_sp= 2.2+ 0.8*rand_numbers ! in  10kms/kpc + is for backw int
         phi_sp= 2*pi*rand_numbers
         ombra= omegas_sp(1)
         write(*,*) "ntot:", ntot
         write(*,*) "tsim, ts:", t_sim, ts
         write(*,*) "omegs_sp for transcient structure:", omegas_sp
         write(*,*) "phi_sp for transcient structure:", phi_sp
       endif
      

      !setting omega of the system
      if((bar_contribution.eqv. .false.).and.(spiral_contribution.eqv. .false.))then
         omega_sys= 0
         initial_phase= 0
      endif

      if((bar_contribution.eqv. .True.).and.(spiral_contribution.eqv. .false.))then
         omega_sys= omegs
         initial_phase= bar_phase
      endif
      
      if((bar_contribution.eqv. .false.).and.(spiral_contribution.eqv. .true.))then
         omega_sys= ombra
         initial_phase= spiral_phase
      endif

      ! In the combined effect, I will be in the ref sys of the bar
      if((bar_contribution.eqv. .true.).and.(spiral_contribution.eqv. .true.))then
         omega_sys= omegs
         initial_phase= bar_phase
         phi0=  spiral_phase- bar_phase
         
      endif
        commit_parameters = 0
      end function commit_parameters
      

      integer function recommit_parameters()
        recommit_parameters = 0
      end function recommit_parameters
      
      ! The following functions change the default values of the constants
      ! these are also legacy functions in interface.py
      
      integer function set_time(value)
        double precision :: value
        time = value
        set_time = 0
      end function set_time
      
      integer function get_time(value)
        double precision :: value
        value = time
        get_time = 0
      end function get_time
      
      integer function get_omega_sys(value)
        double precision :: value
        value = omega_sys
        get_omega_sys = 0
      end function get_omega_sys
      
      integer function get_initial_phase(value)
        double precision :: value
        value = initial_phase
        get_initial_phase = 0
      end function get_initial_phase
      
      ! for computation of tidal tensor
      integer function get_flag(value)
        double precision :: value
        value = xflag
        get_flag = 0
      end function get_flag
      
      integer function set_flag(value)
        double precision :: value
        xflag = value
        set_flag = 0
      end function set_flag
      
      
      ! 4BAR
      integer function set_bar_phase(value)
        double precision :: value
        bar_phase = value
        set_bar_phase = 0
      end function set_bar_phase
      
      integer function get_bar_phase(value)
        double precision :: value
        value = bar_phase
        get_bar_phase = 0
      end function get_bar_phase
      
      integer function set_mass_bar(value)
        double precision :: value
        barmas = value
        set_mass_bar = 0
      end function set_mass_bar
      
      integer function get_mass_bar(value)
        double precision :: value
        value = barmas
        get_mass_bar = 0
      end function get_mass_bar
      
      integer function set_aaxis_bar(value)
        double precision :: value
        aaxis = value
        set_aaxis_bar = 0
      end function set_aaxis_bar
      
      integer function get_aaxis_bar(value)
        double precision :: value
        value = aaxis
        get_aaxis_bar = 0
      end function get_aaxis_bar
      
      
      integer function set_caxis_bar(value)
        real(8) :: value
        caxis = value
        set_caxis_bar = 0
      end function set_caxis_bar
      
      integer function get_caxis_bar(value)
        real(8) :: value
        value = caxis
        get_caxis_bar = 0
      end function get_caxis_bar
      
      
      integer function set_axis_ratio_bar(value)
        double precision :: value
        axirat = value
        set_axis_ratio_bar = 0
      end function set_axis_ratio_bar
      
      integer function get_axis_ratio_bar(value)
        double precision :: value
        value = axirat
        get_axis_ratio_bar = 0
      end function get_axis_ratio_bar
      
      integer function set_omega_bar(value)
        double precision :: value
        omegs = value
        set_omega_bar = 0
      end function set_omega_bar
      
      integer function get_omega_bar(value)
        double precision :: value
        value = omegs
        get_omega_bar = 0
      end function get_omega_bar
      
      integer function set_nbt(value)
        double precision :: value
        nbt = value
        set_nbt = 0
      end function set_nbt
      
      integer function get_nbt(value)
        double precision :: value
        value = nbt
        get_nbt = 0
      end function get_nbt
      
      integer function get_tgrowth(value)
        double precision :: value
        value = tgrow
        get_tgrowth = 0
      end function get_tgrowth
      
      ! SPIRAL ARMS
      
      integer function set_spiral_phase(value)
        double precision :: value
        spiral_phase = value
        set_spiral_phase = 0
      end function set_spiral_phase
      
      integer function get_spiral_phase(value)
        double precision :: value
        value = spiral_phase
        get_spiral_phase = 0
      end function get_spiral_phase
      
      integer function set_m(value)
        double precision :: value
        m = value
        set_m = 0
      end function set_m
      
      integer function get_m(value)
        double precision :: value
        value = m
        get_m = 0
      end function get_m
      
      
      integer function set_tan_pitch_angle(value)
        double precision :: value
        tani = value
        set_tan_pitch_angle = 0
      end function set_tan_pitch_angle
      
      integer function get_tan_pitch_angle(value)
        double precision :: value
        value = tani
        get_tan_pitch_angle = 0
      end function get_tan_pitch_angle
      
      
      integer function set_r_sigma(value)
        double precision :: value
        rsigma = value
        set_r_sigma = 0
      end function set_r_sigma
      
      integer function get_r_sigma(value)
        double precision :: value
        value = rsigma
        get_r_sigma = 0
      end function get_r_sigma
      
      integer function set_omega_spiral(value)
        double precision :: value
        ombra = value
        set_omega_spiral = 0
      end function set_omega_spiral
      
      integer function get_omega_spiral(value)
        double precision :: value
        value = ombra
        get_omega_spiral = 0
      end function get_omega_spiral
      
      ! TWA parameters
      integer function set_rsp(value)
        double precision :: value
        rsp = value
        set_rsp = 0
      end function set_rsp
      
      integer function get_rsp(value)
        double precision :: value
        value = rsp
        get_rsp = 0
      end function get_rsp
      
      
      integer function set_N(value)
        double precision :: value
        nsp = value
        set_N = 0
      end function set_N
      
      integer function get_N(value)
        double precision :: value
        value = nsp
        get_N = 0
      end function get_N
      
      integer function set_amplitude(value)
        double precision :: value
        asp = value
        set_amplitude = 0
      end function set_amplitude
      
      integer function get_amplitude(value)
        double precision :: value
        value = asp
        get_amplitude = 0
      end function get_amplitude
      

      ! CG02 spiral arms
      
      integer function set_spiral_density_amplitude(value)
        real(8) :: value
        rho0_sp = value
        set_spiral_density_amplitude = 0
      end function set_spiral_density_amplitude
      
      integer function get_spiral_density_amplitude(value)
        real(8) :: value
        value = rho0_sp
        get_spiral_density_amplitude = 0
      end function get_spiral_density_amplitude
      
      integer function set_fiducial_radius(value)
        real(8) :: value
        r0_sp = value
        set_fiducial_radius = 0
      end function set_fiducial_radius
      
      integer function get_fiducial_radius(value)
        real(8) :: value
        value = r0_sp
        get_fiducial_radius = 0
      end function get_fiducial_radius
      
      integer function set_scale_height(value)
        real(8) :: value
        h_sp = value
        set_scale_height = 0
      end function set_scale_height
      
      integer function get_scale_height(value)
        real(8) :: value
        value = h_sp
        get_scale_height = 0
      end function get_scale_height
      
      
      !TRANSCIENT STRUCTURE
      
      integer function set_t_sim(value)
        double precision :: value
        t_sim = value
        set_t_sim = 0
      end function set_t_sim
      
      integer function get_t_sim(value)
        double precision :: value
        t_sim = value
        get_t_sim = 0
      end function get_t_sim
      
      integer function set_sigma_s(value)
        double precision :: value
        sigma_s = value
        set_sigma_s = 0
      end function set_sigma_s
      
      integer function get_sigma_s(value)
        double precision :: value
        value = sigma_s
        get_sigma_s = 0
      end function get_sigma_s
      
      !LEPINE MODEL
      
      integer function set_spiral_model(value)
        integer :: value
        spiral_model = value
        set_spiral_model = 0
      end function set_spiral_model
      
      integer function get_spiral_model(value)
        integer :: value
        value = spiral_model
        get_spiral_model = 0
      end function get_spiral_model
      
      integer function set_omega_spiral2(value)
        double precision :: value
        ombra2 = value
        set_omega_spiral2 = 0
      end function set_omega_spiral2
      
      integer function get_omega_spiral2(value)
        double precision :: value
        value = ombra2
        get_omega_spiral2 = 0
      end function get_omega_spiral2
      
      
      integer function set_amplitude2(value)
        double precision :: value
        asp2 = value
        set_amplitude2 = 0
      end function set_amplitude2
      
      integer function get_amplitude2(value)
        double precision :: value
        value = asp2
        get_amplitude2 = 0
      end function get_amplitude2
      
      integer function set_tan_pitch_angle2(value)
        double precision :: value
        tani2 = value
        set_tan_pitch_angle2 = 0
      end function set_tan_pitch_angle2
      
      integer function get_tan_pitch_angle2(value)
        double precision :: value
        value = tani2
        get_tan_pitch_angle2 = 0
      end function get_tan_pitch_angle2
      
      integer function set_m2(value)
        double precision :: value
        m2 = value
        set_m2 = 0
      end function set_m2
      
      integer function get_m2(value)
        double precision :: value
        value = m2
        get_m2 = 0
      end function get_m2
      
      integer function set_phi21(value)
        double precision :: value
        phi21 = value
        set_phi21 = 0
      end function set_phi21
      
      integer function get_phi21(value)
        double precision :: value
        value = phi21
        get_phi21 = 0
      end function get_phi21
      
      
      ! AXI POTENTIAL
      
      integer function set_mass_bulge(value)
        double precision :: value
        MB = value
        set_mass_bulge = 0
      end function set_mass_bulge
      
      integer function get_mass_bulge(value)
        double precision :: value
        value = MB
        get_mass_bulge = 0
      end function get_mass_bulge
      
      integer function set_b_bulge(value)
        double precision :: value
        B1 = value
        set_b_bulge = 0
      end function set_b_bulge
      
      integer function get_b_bulge(value)
        double precision :: value
        value = B1
        get_b_bulge = 0
      end function get_b_bulge
      
      integer function set_mass_disk(value)
        double precision :: value
        MD = value
        set_mass_disk = 0
      end function set_mass_disk
      
      integer function get_mass_disk(value)
        double precision :: value
        value = MD
        get_mass_disk = 0
      end function get_mass_disk
      
      integer function set_a_disk(value)
        double precision :: value
        AD2 = value
        set_a_disk = 0
      end function set_a_disk
      
      integer function get_a_disk(value)
        double precision :: value
        value = AD2
        get_a_disk = 0
      end function get_a_disk
      
      integer function set_b_disk(value)
        double precision :: value
        BD2 = value
        set_b_disk = 0
      end function set_b_disk
      
      integer function get_b_disk(value)
        double precision :: value
        value = BD2
        get_b_disk = 0
      end function get_b_disk
      
      integer function set_mass_halo(value)
        double precision :: value
        MH = value
        set_mass_halo = 0
      end function set_mass_halo
      
      integer function get_mass_halo(value)
        double precision :: value
        value = MH
        get_mass_halo = 0
      end function get_mass_halo
      
      integer function set_a_halo(value)
        double precision :: value
        A3 = value
        set_a_halo = 0
      end function set_a_halo
      
      integer function get_a_halo(value)
        double precision :: value
        value = A3
        get_a_halo = 0
      end function get_a_halo
      
      
      ! These flags set of non axi parts will be added in the axi force and pot
      
      integer function set_spiral_contribution(flag)
        integer :: flag
        if (flag.eq.0) then
           spiral_contribution = .false.
        else
           spiral_contribution = .true.
        endif
        set_spiral_contribution = 0
      end function set_spiral_contribution
      
      integer function get_spiral_contribution(flag)
        integer :: flag
        if (spiral_contribution) then
           flag = 1
        else
           flag = 0
        endif
        get_spiral_contribution = 0
      end function get_spiral_contribution

      integer function set_bar_contribution(flag)
        integer :: flag
        if (flag.eq.0) then
           bar_contribution = .false.
        else
           bar_contribution = .true.
        endif
        set_bar_contribution = 0
      end function set_bar_contribution

      integer function get_bar_contribution(flag)
        integer :: flag
        if (bar_contribution) then
           flag = 1
        else
           flag = 0
        endif
        get_bar_contribution = 0
      end function get_bar_contribution
      
      
      integer function set_transient_spiral(flag)
        integer :: flag
        if (flag.eq.0) then
           transient_spiral = .false.
        else
           transient_spiral = .true.
        endif
        set_transient_spiral = 0
      end function set_transient_spiral
      
      integer function get_transient_spiral(flag)
        integer :: flag
        if (transient_spiral) then
           flag = 1
        else
           flag = 0
        endif
        get_transient_spiral = 0
      end function get_transient_spiral
      
      !THE FOLLOWING FUNCTIONS MAKE THE COMPUTATION OF 
      !FORCES AND POTENTIALS DUE TO BULGE,DISK, HALO,  BAR AND SARMS. 
      !ALSO THERE ARE OTHER HELP FUNCTIONS THAT COMPUTE
      ! THE CIRCULAR VELOCITY, THE EPICYCLIC FREQUENCY AND THE 
      !THE TIDAL TENSOR AT A GIVEN POINT

      integer function get_gravity_at_point(eps, x, y, z, ax, ay, az, n)
        integer :: i,n  
        real*8 :: eps(n), x(n), y(n), z(n), ax(n), ay(n), az(n)
        real*8:: fx, fy, fz
        ax=0
        ay=0
        az=0
        fx=0
        fy=0
        fz=0
        
        do i= 1,n
           call total_force(time, x(i), y(i), z(i), fx, fy, fz)
           !call forbra(x(i),y(i),z(i),fx,fy,fz)
           !write(*,*) "in grav", x(i),y(i),z(i),fx,fy,fz
           ax(i)= -fx
           ay(i)= -fy
           az(i)= -fz
        end do
        
        get_gravity_at_point = 0
      end function get_gravity_at_point
      
      
      integer function get_potential_at_point(eps, x, y, z, potential, n)
        integer :: i,n
        real*8 :: eps(n),x(n), y(n), z(n), potential(n)
        real*8 :: phi
        potential = 0
        
        do i=1,n
           phi=0
           call total_potential(time, x(i), y(i), z(i), phi)         
           
           potential(i)= phi
        enddo
        
        get_potential_at_point = 0
      end function get_potential_at_point


      integer function get_local_density(t, x,y,z, rho,n)
        integer :: i,n  
        real(8) :: t(n), x(n), y(n), z(n), rho(n)
        real(8) :: density
        rho=0
        
        do i=1,n
           call local_density(t(i), x(i),y(i),z(i), density)
           rho(i)= density
        end do
        
        get_local_density = 0
      end function get_local_density


      integer function get_spiral_density(x,y,z, density,n)
        integer :: i, n  
        real(8) :: x(n), y(n), z(n), density(n)
        real(8) :: dens
        
        do i=1,n
           call dsa_cg_anal10(0.d0, (/x(i), y(i), z(i)/), dens,rho0_scaled, eps_sp, atan(tani),r0_sp, m, ombra, 0.d0, h_sp)
           density(i)= dens
        enddo
        get_spiral_density = 0
      end function get_spiral_density
      
      integer function get_velcirc(x,y,z,vc,n)
        integer :: i,n  
        real*8 :: x(n), y(n), z(n),vc(n)
        real*8 :: velc
        vc=0
        
        do i=1,n
           call velcirc(x(i),y(i),z(i),velc)
           vc(i)= velc
        end do
        
        get_velcirc = 0
      end function get_velcirc
      
      integer function get_epifreq(x,y,z,k,n)
        integer :: i, n  
        real*8 :: x(n), y(n), z(n),k(n)
        real*8 :: epi
        
        do i=1,n
           call epic_freq(x(i),y(i),z(i), epi)
           k(i)= epi
        enddo
        get_epifreq = 0
      end function get_epifreq
      
      integer function get_tidal_tensor(t, x,y,z, fuxx, fuxy, fuxz, fuyx, fuyy, fuyz,fuzx, fuzy, fuzz, n)
        integer:: i,n
        real*8 :: x(n), y(n), z(n), t
        real*8:: fuxx(n), fuxy(n), fuxz(n)
        real*8:: fuyx(n), fuyy(n), fuyz(n), fuzx(n), fuzy(n), fuzz(n)
        real*8:: F0xx, F0yx, F0zx, F0xy,F0yy,F0zy, F0xz, F0yz, F0zz
        
        do i=1,n
           call  tidal_tensor(t, x(i),y(i),z(i),F0xx, F0yx, F0zx, F0xy,F0yy,F0zy, F0xz, F0yz, F0zz)
           fuxx(i)= F0xx
           fuyx(i)= F0yx
           fuzx(i)= F0zx 
           fuxy(i)= F0xy
           fuyy(i)= F0yy
           fuzy(i)= F0zy
           fuxz(i)= F0xz
           fuyz(i)= F0yz
           fuzz(i)= F0zz
        end do
        
        get_tidal_tensor = 0
      end function get_tidal_tensor
      
      integer function get_eigen_values(t,x,y,z, lamda1, lamda2, lamda3, n)
        integer :: i, n 
        real*8 :: x(n), y(n), z(n), lamda1(n), lamda2(n), lamda3(n), t
        real*8:: wr1, wr2, wr3
        
        do i=1,n
           call eigen_values(t,x(i),y(i),z(i), wr1, wr2, wr3)
           lamda1(i)= wr1
           lamda2(i)= wr2
           lamda3(i)= wr3
        enddo
        get_eigen_values = 0
      end function get_eigen_values
      
    
      integer function get_tidal_radius(t, x,y,z, mass_cluster, tidalr, n)
        integer:: i, n
        real*8:: t, x(n),y(n) ,z(n), rt, mass_cluster, tidalr(n)
        
        do i=1,n
           call tidal_radius(t, x(i),y(i),z(i), mass_cluster, rt)
           tidalr(i)= rt
        end do
        
        get_tidal_radius=0
      end function get_tidal_radius
          
        
end module BarAndSpiralsInterface

!_____________________ALL SOUBRUTINES___________________________

!______________ TOTAL FORCE__________________________

subroutine total_force(t, xk1, yk, zk,ffx, ffy, ffz)
  !-----------------------------------------------------------------------
  !     Computes the total force in the model
  !-----------------------------------------------------------------------
  use BarAndSpiralsInterface
  implicit none
  double precision :: xk1, yk, zk, ffx, ffy, ffz,  tmt, T, tmtbra
  double precision :: fax= 0.d0, fay= 0.d0, faz= 0.d0
  double precision :: fbx= 0.d0, fby= 0.d0, fbz= 0.d0
  double precision :: fbrax= 0.d0, fbray= 0.d0, fbraz= 0.d0
  double precision :: fx= 0.d0, fy= 0.d0, fz= 0.d0
  double precision :: theta, xsp=0.d0, ysp=0.d0, zsp=0.d0
  
  ! AXI CASE
  if ((bar_contribution.eqv. .false.).and.(spiral_contribution.eqv. .false.))then
     call forax(t,xk1, yk, zk, fax, fay, faz)
     ffx = fax
     ffy = fay
     ffz = faz
  endif
  
  !AXI+BAR
  if ((bar_contribution.eqv. .true.).and.(spiral_contribution.eqv. .false.))then
     call forax(t,xk1, yk, zk, fax, fay, faz)
     CALL TIMEF(T,TIN,TGROW,TMT)
     aux2= tmt 
     call FORBAR(xk1,yk,zk,fbx,fby,fbz)
     ffx = fax + aux2*fbx
     ffy = fay + aux2*fby
     ffz = faz + aux2*fbz
  endif

  !AXI+SPIRAL
  if ((bar_contribution.eqv. .False.).and.(spiral_contribution.eqv. .True.))then
     call forax(t,xk1, yk, zk, fax, fay, faz)
     call force_sp(t, xk1, yk, zk, fbrax, fbray, fbraz)
     ffx = fax + fbrax
     ffy = fay + fbray
     ffz = faz + fbraz
     !     write(*,*) "force:", ffx, ffy, ffz
     !     write(24,120) t*(9.78462e7/1e6), ffx, ffy, ffz, tmtbra, ombra
     !120  format(7f18.10)
  endif
 
  !AXI+BAR+SPIRAL
  if ((bar_contribution.eqv. .True.).and.(spiral_contribution.eqv. .True.))then   
    
     if (transient_spiral .eqv. .False.) then
        theta= phi0+ (ombra-omegs)*t

        call forax(t,xk1, yk, zk, fax, fay, faz)
        CALL TIMEF(T,TIN,TGROW,TMT)
        aux2=tmt
        call FORBAR(xk1,yk,zk,fbx,fby,fbz)
        ffx = fax + aux2*fbx
        ffy = fay + aux2*fby
        ffz = faz + aux2*fbz
        ! bar-> sp
        xsp= xk1*cos(theta) +yk*sin(theta)
        ysp= -xk1*sin(theta) +yk*cos(theta)
        zsp= zk
        call force_sp(t, xsp, ysp, zsp, fbrax, fbray, fbraz)
        ! sp-> bar
        fx= fbrax*cos(theta) -fbray*sin(theta)
        fy= fbrax*sin(theta) + fbray*cos(theta)
        fz= fbraz
        
        ffx= ffx+ fx 
        ffy= ffy+ fy
        ffz= ffz +fz
        !write(40, *) T, xsp, ysp, zsp, ffx, ffy, ffz
        !write(24,110) t*(9.78462e7/1e6), ffx, ffy, ffz, tmtbra, ombra
        !110  format(7f18.10)
     endif

     if (transient_spiral .eqv. .True.) then !transient structure

        call timefbra(t,tmtbra) ! should be called 1st to update phi0 and ombra
        !write(40,*) t, tmtbra

        theta= phi0+ (ombra-omegs)*t

        call forax(t,xk1, yk, zk, fax, fay, faz)
        CALL TIMEF(T,TIN,TGROW,TMT)
        !aux2=TMT*barmas/totalmass
        aux2=tmt
        call FORBAR(xk1,yk,zk,fbx,fby,fbz)
        ffx = fax + aux2*fbx
        ffy = fay + aux2*fby
        ffz = faz + aux2*fbz
        
        ! bar-> sp
        xsp= xk1*cos(theta) +yk*sin(theta)
        ysp= -xk1*sin(theta) +yk*cos(theta)
        zsp= zk
        call force_sp(t, xsp, ysp, zsp,fbrax, fbray, fbraz)
        ! sp-> bar
        fx= tmtbra*fbrax*cos(theta) - tmtbra*fbray*sin(theta)
        fy= tmtbra*fbrax*sin(theta) + tmtbra*fbray*cos(theta)
        fz= tmtbra*fbraz
        
        ffx= ffx+ fx 
        ffy= ffy+ fy
        ffz= ffz +fz

        !write(60, *) T, (ffx**2+ ffy**2 + ffz**2)**0.5
        !write(24,110) t*(9.78462e7/1e6), ffx, ffy, ffz, tmtbra, ombra
        !110  format(7f18.10)       
     endif 

  endif
   
  return
end subroutine total_force
    

!___________Axi force ______________________

subroutine forax(t,x,y,z,px,py,pz)
   !-----------------------------------------------------------------------
   ! first derivatives of the bulge and disc potentials
   !-----------------------------------------------------------------------
   use BarAndSpiralsInterface
   implicit none
   real*8 :: x,y,z,pxd,pyd,pzd,pxb,pyb,pzb,px,py,pz,xp(3),pxh,pyh,pzh
   real*8 :: aux1,t,tmt
      
   CALL TIMEF(T,TIN,TGROW,TMT)
   aux1=(totalmass-tmt*barmas)/totalmass
   !write(24,*) t, aux1
   !write(*,*) , totalmass, tmt, barmass, aux1
   !      aux2=tmt
   !      write(*,*)'forax',tmt,aux1
   
   ! HALO
   call FORHALO(x,y,z,pxh,pyh,pzh)
   ! BULGE
   call forbulge(x,y,z,pxb,pyb,pzb)
   ! DISK
   CALL FORDISK(x,y,z,pxd,pyd,pzd)
   ! BULGE + DISK
   px=aux1*pxb +pxh+pxd
   py=aux1*pyb +pyh+pyd
   pz=aux1*pzb +pzh+pzd

   !      write(*,*)'forax',aux1
   !      write(*,*)aux1*pxb,aux1*pyb,aux1*pzb
   !      write(*,*)pxh,pyh,pzh
   !      write(*,*)pxd,pyd,pzd
   !write(24,100), T, px/aux1, py/aux1, pz/aux1
   !100   format(7f18.10)
   return
 end subroutine forax

SUBROUTINE FORBULGE(X1,Y1,Z1,PX,PY,PZ)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !DIMENSION X(3)
      AM1=G*MB1 

      SUM=X1*X1+Y1*Y1+Z1*Z1+B1*B1
      RAD2=SUM**1.5D0
      FAC=AM1/RAD2
      PX=FAC*X1 
      PY=FAC*Y1
      PZ=FAC*Z1
!     write(*,*)'bulge',x1,Y1,Z1
      RETURN
    END SUBROUTINE FORBULGE


SUBROUTINE FORDISK(X1,Y1,Z1,PX,PY,PZ)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !DIMENSION X(3)
          
      AM2= G*MD

      SUM1=SQRT(Z1*Z1+BD2*BD2)
      SUM2=AD2+SUM1
      DRC=X1*X1+Y1*Y1
      SUM3=DRC+SUM2*SUM2
      RAD1=SQRT(SUM3)
      RAD2=SUM3**1.5D0
      FAC=AM2/RAD2
      PX=FAC*X1      
      PY=FAC*Y1
      PZ=FAC*Z1*SUM2/SUM1  
      RETURN
    END SUBROUTINE FORDISK

SUBROUTINE FORHALO(X1,Y1,Z1,PX,PY,PZ)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !DIMENSION X(3)

      AM3= G*MH

      DIST=SQRT(X1*X1+Y1*Y1+Z1*Z1)
      DIST3=DIST**3
      DEN1=1.D0+(RMAX/A3)**CS2
      IF(DIST.LT.RMAX) GO TO 1
      ANUM=AM3*((RMAX/A3)**CS1)
      AMHAL=ANUM/DEN1
      FAC=AMHAL/DIST3
      PX=FAC*X1
      PY=FAC*Y1
      PZ=FAC*Z1
      GO TO 2
 1    DEN2=1.D0+(DIST/A3)**CS2
      ANUM=AM3*((DIST/A3)**CS1)
      AMHALR=ANUM/DEN2
      FACT=AM3/(A3*CS2)
      IF(DIST.GT.0.0D0) GO TO 3
      PX=0.0D0
      PY=0.0D0
      PZ=0.0D0
      GO TO 2
 3    CONTINUE
      FAC=AMHALR/DIST3
      PX=FAC*X1       
      PY=FAC*Y1
      PZ=FAC*Z1
 2    CONTINUE
      RETURN
    END SUBROUTINE FORHALO

!___________________Force bar___________________

 SUBROUTINE FORBAR (X1, Y1,Z1,px,py,pz)
!*******************************************************************
!       Forces FX, FY, FZ of the Ferrers n=1 potential at X, Y, Z
!       Semi-axes A,B,C
!       Restriction : A > B > C >= 0
!       Parameters via common parsb
!*******************************************************************
   use BarAndSpiralsInterface 
   IMPLICIT NONE
   REAL*8:: X1,Y1,Z1,PX,PY,PZ
   REAL*8:: XQ(3)
   REAL*8:: X2,Y2,Z2, XL,UA3,B3,UB3,UC3
   REAL*8:: W100,W010,W001,W110,W101,W011,W200,W020,W002
   REAL*8:: XLMBD

   XQ(1)=X1*1.D0
   XQ(2)=Y1*1.D0
   XQ(3)=Z1*1.D0

   X2 = XQ(1)*XQ(1)
   Y2 = XQ(2)*XQ(2)
   Z2 = XQ(3)*XQ(3)
	
   IF (X2*UA2+Y2*UB2+Z2*UC2.LE.1.D0) THEN
      px = 2.D0*CTE1*XQ(1)*(V100-X2*V200-Y2*V110-Z2*V101)
      py = 2.D0*CTE1*XQ(2)*(V010-X2*V110-Y2*V020-Z2*V011)
      pz = 2.D0*CTE1*XQ(3)*(V001-X2*V101-Y2*V011-Z2*V002)
   ELSE
      XL = XLMBD(X2,Y2,Z2,A2,B2,C2)
      UA3 = 1.D0/(A2+XL)
      B3 = B2+XL
      UB3 = 1.D0/B3
      UC3 = 1.D0/(C2+XL)
      PHI= DASIN(DSQRT(UA3/UA2C2))
      CALL ELINT(PHI,XK,F,E)
      D2 = 2.D0*DSQRT(UA3*UB3*UC3)
           !write(24,*) c2,xl

      W100 = 2.D0*(F-E)*UA2B2*SUA2C2
      W001 = (D2*B3 - 2.D0*E*SUA2C2)*UB2C2
      W010 = D2 - W100 - W001
      
      W110 = (W010 - W100)*UA2B2
      W101 = (W001 - W100)*UA2C2
      W011 = (W001 - W010)*UB2C2
      
      W200 = (D2*UA3 - W110 - W101)*US3
      W020 = (D2*UB3 - W110 - W011)*US3
      W002 = (D2*UC3 - W011 - W101)*US3
      
      PX = 2.D0*CTE1*XQ(1)*(W100-X2*W200-Y2*W110-Z2*W101)
      PY = 2.D0*CTE1*XQ(2)*(W010-X2*W110-Y2*W020-Z2*W011)
      PZ = 2.D0*CTE1*XQ(3)*(W001-X2*W101-Y2*W011-Z2*W002)
          ! write(24,*) -X2, D2, UA3, - W110, - W101, US3 
   ENDIF
   RETURN
 END SUBROUTINE FORBAR

 !-----------------------------------------------------------------------
 ! force for the spiral arms
 !-----------------------------------------------------------------------
 ! selects model given by spiral_model   

 subroutine force_sp(t, x_sp, y_sp, z_sp, fx_sp, fy_sp, fz_sp)
   use BarAndSpiralsInterface
   implicit none
   real*8 :: t, x_sp, y_sp, z_sp, fx_sp, fy_sp, fz_sp
   real(8) :: fvec(3)
      
      if (spiral_model .eq. 0) then
         call FORBRA(x_sp, y_sp, z_sp, fx_sp, fy_sp, fx_sp)

      elseif (spiral_model .eq. 1) then
         fvec = 0.d0
         call fsa_cg(0.d0, (/x_sp,y_sp,z_sp/), fvec, &
              rho0_scaled, eps_sp, atan(tani), r0_sp, m, ombra, 0.d0, h_sp)
         ! time and phi0 set to 0 -- rotation to the coordinate sys. of the spirals
         ! happens in total_force
         ! note the minus sign for the force components
         !   -- in the get_gravity_at_point the force is multyplied by (-1)
         !   so the output of all force routines must be (-1)*force
         fx_sp = -fvec(1)
         fy_sp = -fvec(2)
         fz_sp = -fvec(3)
       
      elseif (spiral_model .eq. 2)then
         call force_lepine(t, x_sp, y_sp, z_sp,fx_sp, fy_sp, fx_sp)
        
      else
         write(*,*) ' ** error -- spiral_model', spiral_model, 'not defined! **'
         stop 'unrecognized spiral_model -- code terminated'
      endif

       
    end subroutine force_sp

!____________________FORCE LEPINE________________________________________
    !
    !FUNCTION THAT COMPUTES THE FORCE DUE TO THE SUPERPOSITION
    ! OF TWO SPIRAL PATTERNS. It can be both with m=2 or m=2 and m=4
    ! SEE MISHUROV AND ACHAROVA, MNRAS 2011 AND LEPINE ET AL. 2001
    ! Parameters primary SA:  (can be changed through the python sript)
    !   A1= asp 
    !   m1= m
    !   omega_sp1= ombra
    !   tani1= tani
    !   phi1= phi_sp
    ! Parameters secondary SA: (can be changed through the python sript)
    !   A2= asp2 
    !   m2= m2
    !   omega_sp2= ombra2
    !   tani2= tani2
    !   phi21= phi21

     SUBROUTINE force_lepine(t, x_sp, y_sp, z_sp, ftx_sp, fty_sp, ftz_sp)

      use BarAndSpiralsInterface

      double precision :: t, x_sp, y_sp, z_sp, ftx_sp, fty_sp, ftz_sp
      double precision :: fx_sp1, fy_sp1, fz_sp1, fx_sp2, fy_sp2, fz_sp2
      double precision :: fx_sp2_in1, fy_sp2_in1, fz_sp2_in1
      double precision :: x_sp2, y_sp2, z_sp2
      double precision :: theta1


      call FORBRA(x_sp, y_sp, z_sp, fx_sp1, fy_sp1, fz_sp1) 
      ftx_sp= fx_sp1
      fty_sp= fy_sp1
      ftz_sp= fz_sp1
      !sp1 -> sp2
      theta1= phi21 +(ombra2 -ombra)*t
      x_sp2= x_sp*cos(theta1) + y_sp*sin(theta1)
      y_sp2= -x_sp*sin(theta1) + y_sp*cos(theta1)
      z_sp2= z_sp
      call FORBRA2(x_sp2, y_sp2, z_sp2, fx_sp2, fy_sp2, fz_sp2)
      ! sp2 -> sp1
      fx_sp2_in1= fx_sp2*cos(theta1) - fy_sp2*sin(theta1)
      fy_sp2_in1= fx_sp2*sin(theta1) + fy_sp2*cos(theta1)
      fz_sp2_in1= fz_sp2
      
      ftx_sp= ftx_sp + fx_sp2_in1
      fty_sp= fty_sp + fy_sp2_in1
      ftz_sp= ftz_sp + fz_sp2_in1

    end SUBROUTINE FORCE_LEPINE


    SUBROUTINE FORBRA(X1,Y1,Z1,FBRX,FBRY,FBRZ)
      use BarAndSpiralsInterface
      double precision :: X1,Y1,Z1, R,PHI1, termz, COC1, COC2, SUM
      double precision :: FUNF, FUNG, DFDR, DGDR, FACS, FACC, epss
      double precision :: frbr, ffibr,fbrx,fbry,fbrz
      
      R=DSQRT(X1*X1+Y1*Y1)
      PHI1=ATAN2(Y1,X1)

      epss= 1.d0/rsigma
      COC1=(R/RSP)**(NSP-1)
      COC2=COC1*R/RSP
      SUM=1.D0+COC2
      
      FUNF= DBLE(m)*DLOG(SUM)/(DBLE(NSP)*TANI)
      FUNG=-ASP*R*DEXP(-EPSS*R)
      DFDR= DBLE(m)*COC1/(RSP*SUM*TANI)
      DGDR=-ASP*(1.D0-EPSS*R)*DEXP(-EPSS*R)
      FACS=SIN(DBLE(m)*PHI1-FUNF)
      FACC=COS(DBLE(m)*PHI1-FUNF)
      
      FRBR= FUNG*FACS*DFDR +DGDR*FACC !RADIAL FORCE
      FFIBR= -DBLE(m)*FUNG*FACS/R   !TANGENTIAL FORCE
      FBRX= FRBR*COS(PHI1)-FFIBR*SIN(PHI1)
      FBRY= FRBR*SIN(PHI1)+FFIBR*COS(PHI1)
      FBRZ= 0.D0
      
      RETURN
    END SUBROUTINE FORBRA


    SUBROUTINE FORBRA2(X1,Y1,Z1,FBRX,FBRY,FBRZ)
      use BarAndSpiralsInterface
      double precision :: X1,Y1,Z1, R,PHI1, termz, COC1, COC2, SUM
      double precision :: FUNF, FUNG, DFDR, DGDR, FACS, FACC, epss
      double precision :: frbr, ffibr,fbrx,fbry,fbrz
      
      R=DSQRT(X1*X1+Y1*Y1)
      PHI1=ATAN2(Y1,X1)

      epss= 1.d0/rsigma
      COC1=(R/RSP)**(NSP-1)
      COC2=COC1*R/RSP
      SUM=1.D0+COC2
      
      FUNF= DBLE(m2)*DLOG(SUM)/(DBLE(NSP)*TANI2)
      FUNG=-ASP2*R*DEXP(-EPSS*R)
      DFDR= DBLE(m2)*COC1/(RSP*SUM*TANI2)
      DGDR=-ASP2*(1.D0-EPSS*R)*DEXP(-EPSS*R)
      FACS=SIN(DBLE(m2)*PHI1-FUNF)
      FACC=COS(DBLE(m2)*PHI1-FUNF)
      
      FRBR= FUNG*FACS*DFDR +DGDR*FACC !RADIAL FORCE
      FFIBR= -DBLE(m2)*FUNG*FACS/R   !TANGENTIAL FORCE
      FBRX= FRBR*COS(PHI1)-FFIBR*SIN(PHI1)
      FBRY= FRBR*SIN(PHI1)+FFIBR*COS(PHI1)
      FBRZ= 0.D0
      
      RETURN
    END SUBROUTINE FORBRA2
   

!______________________log spiral force__________________________________
subroutine force_log_spiral(x1,y1,z1,fx1,fy1,fz1)
  use BarAndSpiralsInterface
  double precision:: x1,y1,z1,fx1,fy1,fz1, phi1, R
  double precision:: fm, k, expFactor, kDer
  double precision:: FR, Fphi, Fz, sigma2

  sigma2= sigma*sigma

  phi1= atan2(y1,x1)
  R= sqrt(x1*x1+ y1*y1)

  if (R .eq. 0)then
     FR=0
     Fphi=0
     Fz=0

  else
     fm= (m/tani)*log(R/rsp) + gamma
     k= m/(R*tani)
     expFactor= exp(-R*R/sigma2*(1-cos(m*phi1-fm))-epss*R-abs(k*z1))
     !write(*,*) "expF:", expFactor, "K:", k, "A0:", asp, "R", R 

     if (k .ge. 0)then
        kDer= -m/(R*R*tani)
     endif
     if (k .lt. 0)then
        kDer= m/(R*R*tani)
     endif
     FR= asp*expFactor*(1 +2*R*R*(cos(m*phi1-fm)-1)/sigma2 & 
          +R*R*m*sin(m*phi1-fm)/(sigma2*tani) -epss*R - kDer*abs(z1)*R )

     Fphi= -R*R*m*asp*expFactor*sin(m*phi1-fm)/sigma2
  
     if (z1 .gt. 0)then
        Fz= -asp*abs(k)*R*expFactor
        !write(*,*) "FZ!!!!!", Fz
     endif
     
     if(z1 .lt. 0)then
        Fz= asp*abs(k)*R*expFactor
        !write(*,*) "or FZ!!!!!", Fz
     endif

     ! Spiral force on the plane must be zero
     if (z1 .eq.0) then
        Fz=0
     endif

  endif

  fx1=  -(FR*cos(phi1) -Fphi*sin(phi1))
  fy1=  -(FR*sin(phi1)+ Fphi*cos(phi1))
  fz1=  -Fz

end subroutine force_log_spiral


!____________________________TOTAL POTENTIAL______________
!_________________________________________________________

subroutine total_potential(t,x,y,z,pot)
  !**********************************************************************
  ! returns the potential of the model
  !**********************************************************************
  use BarAndSpiralsInterface
  implicit none
  double precision :: pot,POTA,POTB,POTSP,tmt,t,x,y,z, potential, tmtbra
  double precision :: theta, xsp, ysp, zsp
  ! x the complete case, bar and spiral should be initially aligned =>
  double precision:: angle, pbulge, aux1 

  ! AXI CASE
  if ((bar_contribution.eqv. .false.).and.(spiral_contribution.eqv. .false.))then
     call potax(t,x, y, z,pota)
     pot=pota
  endif
  
  !AXI+BAR
  if ((bar_contribution.eqv. .true.).and.(spiral_contribution.eqv. .false.))then
     call potax(t,x, y, z,pota)
     CALL TIMEF(T,TIN,TGROW,TMT)
     aux2=tmt
     call potbar (x,y,z, potb)
     pot=pota+aux2*potb
     !these 4 lines are for plotting time-dep of potentials
     !call POTBULGE(X,Y,Z,pbulge)
     !aux1=(totalmass-tmt*barmas)/totalmass
     !write(60, 100) t, (-1)*(aux1*pbulge), (-1)*(aux2*potb)
     !write(50, 100) t, aux1, aux2
!100   format(7f18.10)     
  endif

  
  !AXI+SPIRAL
  if ((bar_contribution.eqv. .False.).and.(spiral_contribution.eqv. .True.))then
     call potax(t,x, y, z,pota)
     call pot_sp(t, x,y,z, potsp)
     pot= pota+ potsp
     !pot=  tmtbra*potsp
!          write(24,100) t, pota, potsp
!100  format(7f18.10)    
  endif
  
  !AXI+BAR+SPIRAL
  
   if ((bar_contribution.eqv. .True.).and.(spiral_contribution.eqv. .True.))then

      if (transient_spiral .eqv. .false.) then
         theta= phi0- (ombra-omegs)*t
         call potax(t,x, y, z,pota)
         CALL TIMEF(T,TIN,TGROW,TMT)
         aux2= tmt
         call potbar (x,y,z, potb)
         pot=pota +aux2*potb
         !bar -> sp
         xsp= x*cos(theta) + y*sin(theta)
         ysp= -x*sin(theta) + y*cos(theta)
         call pot_sp(t, xsp,ysp,z, potsp)
         pot= pot+ potsp
      endif
      
      if (transient_spiral .eqv. .True.)then !transient spirals
         call TIMEFBRA (T, TMTBRA) !this should be called 1st to update ombra and phi0 
         theta= phi0- (ombra-omegs)*t
         call potax(t,x, y, z,pota)
         CALL TIMEF(T,TIN,TGROW,TMT)
         aux2=tmt
         call potbar (x,y,z, potb)
         pot=pota +aux2*potb
         !bar -> sp
         xsp= x*cos(theta) + y*sin(theta)
         ysp= -x*sin(theta) + y*cos(theta)
         call pot_sp(t, xsp,ysp,z, potsp)
         pot= pot+ tmtbra*potsp
      endif

   endif
  
  return
end subroutine total_potential


 !__________________________________Axi potential__________________________
subroutine potax (t,x, y, z, pot)
   use BarAndSpiralsInterface
   implicit none
   real*8 :: x,y,z,pot,POTBU,POTDI,POTHA
   real*8 :: aux1,t, tmt
   
   CALL TIMEF(T,TIN,TGROW,TMT)
   aux1=(totalmass-tmt*barmas)/totalmass
   !write(*,*) totalmass, tmt, barmass
   
!     bulge
   call potbulge(x,y,z,potbu)
!     halo
   call POTHALO(x,y,z,potha)
!     disk
   call potdisk(x,y,z,potdi)
!     total
   pot = aux1*potbu + potha + potdi
   !pot= potbu
   !write(*,*) "In potax", 'aux1', aux1, 'phibulge', potbu, 'phihalo', potha , 'phidisk', potdi

   return
 end subroutine potax


SUBROUTINE POTBULGE(X1,Y1,Z1,POT)
  use BarAndSpiralsInterface
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  !DIMENSION X(3)    
  AM1=G*MB1
  
  SUM=X1*X1+Y1*Y1+Z1*Z1+B1*B1
  RAD1=SQRT(SUM)
  POT=-AM1/RAD1
  RETURN
END SUBROUTINE POTBULGE


SUBROUTINE POTDISK(X1,Y1,Z1,POT)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !DIMENSION X(3)
      
      AM2= G*MD
      
      SUM1=SQRT(Z1*Z1+BD2*BD2)
      SUM2=AD2+SUM1
      DRC=X1*X1+Y1*Y1
      SUM3=DRC+SUM2*SUM2
      RAD1=SQRT(SUM3)
      POT=-AM2/RAD1
      RETURN
    END SUBROUTINE POTDISK


SUBROUTINE POTHALO(X1,Y1,Z1,POT)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !DIMENSION X(3)
      
      AM3= G*MH

      DIST=SQRT(X1*X1+Y1*Y1+Z1*Z1)
      DIST3=DIST**3
      DEN1=1.D0+(RMAX/A3)**CS2
      IF(DIST.LT.RMAX) GO TO 1
      ANUM=AM3*((RMAX/A3)**CS1)
      AMHAL=ANUM/DEN1
      POT=-AMHAL/DIST
      GO TO 2
 1    DEN2=1.D0+(DIST/A3)**CS2
      ANUM=AM3*((DIST/A3)**CS1)
      AMHALR=ANUM/DEN2
      FACT=-AM3/(A3*CS2)
      IF(DIST.GT.0.0D0) GO TO 3
      POT=FACT*(CS2/DEN2-CS2/DEN1+LOG(DEN1/DEN2))
      GO TO 2
 3    CONTINUE
      POT=-AMHALR/DIST+FACT*(CS2/DEN2-CS2/DEN1+LOG(DEN1/DEN2))
 2    CONTINUE
      RETURN
    END SUBROUTINE POTHALO

         
!_________________________ Bar potential___________________

     SUBROUTINE POTBAR(X1,Y1,Z1,POT)
!**********************************************************************
! Given X(x,y,z,xd,yd,zd) coord of a point, computes the potential at
! this point, where the potential is Ferrers with n=1, and semiaxes
! a>b>c>=0
! NOTE: parameters of the potential must be given via common PARSB
!**********************************************************************
       use BarAndSpiralsInterface
       IMPLICIT NONE
       REAL*8:: X1,Y1,Z1,POT
       REAL*8:: XQ(3)
       REAL*8:: X2,Y2,Z2,P1,P2,P3,XL,UA3,B3,UB3,UC3
       REAL*8:: W000,W100,W010,W001,W110,W101,W011,W200,W020,W002
       REAL*8:: XLMBD, term1, term2, term3       
       
       XQ(1)=X1*1.D0
       XQ(2)=Y1*1.D0
       XQ(3)=Z1*1.D0
       
       X2=XQ(1)*XQ(1)
       Y2=XQ(2)*XQ(2)
       Z2=XQ(3)*XQ(3)
       
       IF(X2*UA2+Y2*UB2+Z2*UC2.LE.1.D0)THEN
          P1=X2*V100+Y2*V010+Z2*V001
          P2=-(x2*y2*v110+x2*z2*v101+y2*z2*v011)
          P3=-0.5Q0*(v000+x2*x2*v200+y2*y2*v020+z2*z2*v002)
          POT=CTE1*(P1+P2+P3)
          
       ELSE
          XL = XLMBD(X2,Y2,Z2,A2,B2,C2)
          UA3 = 1.D0/(A2+XL)
          B3 = B2+XL
          UB3 = 1.D0/B3
          UC3 = 1.D0/(C2+XL)
          PHI= DASIN(DSQRT(UA3/UA2C2))
          CALL ELINT(PHI,XK,F,E)
          D2 = 2.D0*DSQRT(UA3*UB3*UC3)
          
          W000 = 2.D0*F*SUA2C2 
          
          W100 = 2.D0*(F-E)*UA2B2*SUA2C2
          W001 = (D2*B3 - 2.D0*E*SUA2C2)*UB2C2
          W010 = D2 - W100 - W001
          
          W110 = (W010 - W100)*UA2B2
          W101 = (W001 - W100)*UA2C2
          W011 = (W001 - W010)*UB2C2
          
          W200 = (D2*UA3 - W110 - W101)*US3
          W020 = (D2*UB3 - W110 - W011)*US3
          W002 = (D2*UC3 - W011 - W101)*US3
          
          POT = CTE1*(-0.5D0*W000+X2*W100+Y2*W010+Z2*W001- &
             X2*Y2*W110-X2*Z2*W101-Y2*Z2*W011- &
             0.5D0*X2*X2*W200-0.5D0*Y2*Y2*W020-0.5D0*Z2*Z2*W002)
         
       ENDIF
  RETURN
END SUBROUTINE POTBAR


!-----------------------------------------------------------------------
! potential for the spiral arms
!-----------------------------------------------------------------------
! selects model given by spiral_model


     SUBROUTINE pot_sp(t, X1, Y1,Z1, POT)
      use BarAndSpiralsInterface
      IMPLICIT NONE 
      double precision :: t, X1, Y1,Z1, POT
      
      if (spiral_model .eq. 0)then
         call potbra(x1, y1, z1, POT)
         
      else if (spiral_model .eq. 1) then
         call potsa_cg(0.d0, (/x1,y1,z1/), POT, &
              rho0_scaled, eps_sp, atan(tani), r0_sp, m, ombra, 0.d0, h_sp)
         ! time and phi0 set to 0 -- rotation to the coordinate sys. of the spirals
         ! happens in total_potential
         
      elseif (spiral_model .eq. 2)then
         call pot_lepine(t, X1, Y1,Z1, POT)
         
      else
         write(*,*) ' ** error -- spiral_model', spiral_model, 'not defined! **'
         stop 'unrecognized spiral_model -- code terminated'
      endif

    end SUBROUTINE pot_sp



    SUBROUTINE POTBRA(X1, Y1,Z1, POTBR)
      use BarAndSpiralsInterface
      double precision:: X1, Y1, Z1, R, PHI1, termz, AMP, FASE, POTBR, epss
      
      R=DSQRT(X1*X1+Y1*Y1)
      PHI1=ATAN2(Y1,X1)
      
      epss= 1.d0/rsigma
      termz= 1/cosh(Z1/H)**2

      AMP=-ASP*R*DEXP(-EPSS*R)
      FASE=(DBLE(m)/(NSP*TANI))*DLOG(1.D0+(R/RSP)**NSP)

      POTBR=AMP*DCOS(DBLE(m)*PHI1-FASE)

      RETURN
    END SUBROUTINE POTBRA

    SUBROUTINE POTBRA2(X1, Y1,Z1, POTBR)
      use BarAndSpiralsInterface
      double precision:: X1, Y1, Z1, R, PHI1, termz, AMP, FASE, POTBR, epss
      
      R=DSQRT(X1*X1+Y1*Y1)
      PHI1=ATAN2(Y1,X1)
      
      epss= 1.d0/rsigma
      termz= 1/cosh(Z1/H)**2

      AMP=-ASP2*R*DEXP(-EPSS*R)
      FASE=(DBLE(m2)/(NSP*TANI2))*DLOG(1.D0+(R/RSP)**NSP)

      POTBR=AMP*DCOS(DBLE(m2)*PHI1-FASE)

      RETURN
    END SUBROUTINE POTBRA2


!____________________Lepine potential ________________________________________
    !
    !FUNCTION THAT COMPUTES THE POTENTIAL DUE TO THE SUPERPOSITION
    ! OF TWO SPIRAL PATTERNS WITH M=2 AND M=4
    ! SEE MISHUROV AND ACHAROVA, MNRAS 2011 AND LEPINE ET AL. 2001
    ! Parameters m=2:
    !   A= asp (can be changed in the script)
    !   m=2
    !   i= -7 grades
    !   omega= ombra (can be changed in the script)
    ! Parameters m=4:
    !   A= 0.8 asp
    !   m=4
    !   i= -14 grades
    !   varphi= -200 grades (w.r.to m=2) 
    !   omega= ombra_m4 (can be changed in the script)

    SUBROUTINE pot_lepine(t, X1, Y1,Z1, POTBR)
      use BarAndSpiralsInterface
      implicit none
      double precision:: t, X1, Y1, Z1, R, POTBR, potsp1
      double precision :: theta1
      double precision :: xsp2, ysp2, zsp2, potsp2
      

      call potbra(x1, y1, z1, potsp1)
      !sp1 -> sp2
      theta1= phi21 -(ombra2 -ombra)*t
      xsp2= x1*cos(theta1) + y1*sin(theta1)
      ysp2= -x1*sin(theta1) + y1*cos(theta1)
      zsp2= z1
      call potbra2(xsp2, ysp2, zsp2, potsp2 )
      potbr= potsp1+ potsp2

    end SUBROUTINE pot_lepine



   
!______________________logaritmic spiral potential_________________
subroutine pot_log_spiral(x1,y1,z1, potential)

  use BarAndSpiralsInterface
  double precision:: x1,y1,z1, phi1, R, potential
  double precision:: fm, expFactor,k, sigma2
  sigma2= sigma*sigma

  phi1= atan2(y1,x1)
  R= sqrt(x1*x1+ y1*y1)
  
  if (R .eq. 0) then
     potential= 0

  else
     fm= (m/tani)*log(R/rsp) + gamma
     k= m/(R*tani)
     expFactor= exp(- R*R/sigma2*(1-cos(m*phi1-fm))-epss*R -abs(k*z1) )
     potential= -asp*R*expFactor
   
  endif
  
end subroutine pot_log_spiral



!____________________________TOTAL POTENTIAL of the stellar component______________
!________________________________________________

subroutine total_stellar_potential(t,x,y,z,pot)
  !**********************************************************************
  ! returns the potential of the stellar component. This is to compute
  ! Local density.
  !**********************************************************************
  use BarAndSpiralsInterface
  implicit none
  real(8) :: pot,POTA,POTB,POTSA,tmt,t,x,y,z, tmtbra
  real(8) :: theta, xsp, ysp, zsp
  ! x the complete case, bar and spiral should be initially aligned =>
  real(8):: angle 

  ! AXI CASE
  if ((bar_contribution.eqv. .false.).and.(spiral_contribution.eqv. .false.))then
     call potax_stellar(t,x, y, z,pota)
     pot=pota
     !write(*,*) "pot axi", pota
  endif
  
  !AXI+BAR
  if ((bar_contribution.eqv. .true.).and.(spiral_contribution.eqv. .false.))then
     call potax_stellar(t,x, y, z,pota)
     CALL TIMEF(T,TIN,TGROW,TMT)
     aux2=tmt
     call potbar (x,y,z, potb)
     pot=pota+aux2*potb
      write(*,*) "pot axi+bar", pota, aux2, potb, pot
  endif
  
  !AXI+SPIRAL
  if ((bar_contribution.eqv. .False.).and.(spiral_contribution.eqv. .True.))then
     call potax_stellar(t,x, y, z,pota)
     call pot_sp(t, x,y,z, potsa)
     pot= pota+ potsa
     write(*,*) 'potaxi+ pot SA', pota, potsa
!          write(24,100) t, pota, potsp
!100  format(7f18.10)    
  endif
  
  !AXI+BAR+SPIRAL  
  if ((bar_contribution.eqv. .True.).and.(spiral_contribution.eqv. .True.))then
     if (transient_spiral .eqv. .false.) then
        theta= phi0- (ombra-omegs)*t
        call potax_stellar(t,x, y, z,pota)
        CALL TIMEF(T,TIN,TGROW,TMT)
        aux2= tmt
        call potbar(x,y,z, potb)
        pot=pota +aux2*potb
        !bar -> sp
        xsp= x*cos(theta) + y*sin(theta)
        ysp= -x*sin(theta) + y*cos(theta)
        call pot_sp(t, xsp,ysp,z, potsa)
        pot= pot+ potsa
        write(*,*)'all_pot', pota, potb, potsa, pot
     endif
      
     if (transient_spiral .eqv. .True.)then !transient spirals
        call TIMEFBRA (T, TMTBRA) !this should be called 1st to update ombra and phi0 
        theta= phi0- (ombra-omegs)*t
        call potax_stellar(t,x, y, z,pota)
        CALL TIMEF(T,TIN,TGROW,TMT)
        aux2=tmt
        call potbar (x,y,z, potb)
        pot=pota +aux2*potb
        !bar -> sp
        xsp= x*cos(theta) + y*sin(theta)
        ysp= -x*sin(theta) + y*cos(theta)
        call pot_sp(t, xsp,ysp,z, potsa)
        pot= pot+ tmtbra*potsa
     endif  
  endif
  return
end subroutine total_stellar_potential

subroutine potax_stellar(t,x, y, z, pot)
   use BarAndSpiralsInterface
   implicit none
   real(8) :: x,y,z,pot,POTBU,POTDI,POTHA
   real(8) :: aux1,t, tmt
   
   CALL TIMEF(T,TIN,TGROW,TMT)
   aux1=(totalmass-tmt*barmas)/totalmass
   call potbulge(x,y,z,potbu)
   call potdisk(x,y,z,potdi)
   pot = aux1*potbu + potdi
   return
 end subroutine potax_stellar

!_______________________ Local density ______________________________________


 subroutine local_density(t, x1,y1,z1, density)
   use BarAndSpiralsInterface
   implicit none
   real*8 :: t, x1,y1,z1, dx, dy, dz, density, stellar_phi
   real*8 :: pot, pot1xx, pot2xx, pot1yy, pot2yy, pot1zz, pot2zz
   real*8 :: Fxx, Fyy, Fzz, Fxx_r, Fyy_r, Fzz_r
   real*8 :: angle, xr, yr, zr
   dx= 0.0001
   dy= dx
   dz= dx
   !write(*,*) dx, dy, dz
   
   if (xflag .EQ. 1) then ! coordinates in rotating frame
      call total_stellar_potential(t,x1,y1,z1, pot)
      call total_stellar_potential(t, (x1+dx), y1, z1, pot1xx)
      call total_stellar_potential(t, (x1-dx), y1, z1, pot2xx)
      Fxx= (pot1xx-2*pot+ pot2xx)/dx**2
      call total_stellar_potential(T, x1, (y1+dy), z1, pot1yy)
      call total_stellar_potential(T, x1, (y1-dy), z1, pot2yy)
      Fyy= (pot1yy-2*pot+ pot2yy)/dy**2
      call total_stellar_potential(T, x1, y1, (z1+dz), pot1zz)
      call total_stellar_potential(T, x1, y1, (z1-dz), pot2zz)
      Fzz= (pot1zz-2*pot+ pot2zz)/dz**2
      stellar_phi= Fxx+ Fyy+ Fzz  
      density= stellar_phi/(4.d0*pi*G)
      write(*,*)'density in 1:', Fxx, Fyy, Fzz, stellar_phi, density
   endif
   
   if (xflag .EQ. 2) then
      ! inertial -> rotating
      angle= initial_phase + omega_sys*t
      xr= x1*cos(angle) + y1*sin(angle)
      yr= -x1*sin(angle) + y1*cos(angle)
      zr= z1
      
      call total_stellar_potential(t,xr,yr,zr, pot)
      call total_stellar_potential(t, (xr+dx), yr, zr, pot1xx)
      call total_stellar_potential(t, (xr-dx), yr, zr, pot2xx)
      Fxx_r= (pot1xx-2*pot+ pot2xx)/dx**2
      call total_stellar_potential(T, xr, (yr+dy), zr, pot1yy)
      call total_stellar_potential(T, xr, (yr-dy), zr, pot2yy)
      Fyy_r= (pot1yy-2*pot+ pot2yy)/dy**2
      call total_stellar_potential(T, xr, yr, (zr+dz), pot1zz)
      call total_stellar_potential(T, xr, yr, (zr-dz), pot2zz)
      Fzz_r= (pot1zz-2*pot+ pot2zz)/dz**2
      
      stellar_phi= Fxx_r+ Fyy_r+Fzz_r !checked that gradient is the same in inertial or rotating systems  
      density= stellar_phi/(4.d0*pi*G)
      write(*,*)'density in 2:', Fxx_r, Fyy_r, Fzz_r, stellar_phi, density
   end if
    
 end subroutine local_density

 



!__________________________ DERIVATIVE OF THE AXI FORCE___________________
subroutine der_forax (t,x, y, z, der_force)
   use BarAndSpiralsInterface
   implicit none
   real*8 :: x,y,z,der_force,der_fbu,der_fdisk,der_fhalo
   real*8 :: aux1,t, tmt
   
   CALL TIMEF(T,TIN,TGROW,TMT)
   aux1=(totalmass-tmt*barmas)/totalmass
   call der_forbulge(x,y,z,der_fbu)
   call der_fordisk(x,y,z,der_fdisk)
   call der_forhalo(x,y,z,der_fhalo)
   der_force = aux1*der_fbu + der_fdisk + der_fhalo
   !write(*,*) 'derivative of the axi force is:'
   !write(*,*) der_force

   return
 end subroutine der_forax

! second derivative of the bulge potential 
SUBROUTINE DER_FORBULGE(x1,y1,z1,der_force)
  use BarAndSpiralsInterface
  IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
 
  AM1=G*MB1
  R2=x1*x1+y1*y1
  factor= AM1/(R2+z1*z1+b1*b1)**1.5D0 
  der_pot= factor*(1-3*R2/(R2+z1*z1+b1*b1))
  RETURN
END SUBROUTINE DER_FORBULGE

!Second derivative of the disk potential
SUBROUTINE DER_FORDISK(X1,Y1,Z1,der_force)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      AM2= G*MD     
      R2=x1*x1+y1*y1
      SUM1=SQRT(Z1*Z1+BD2*BD2)
      SUM2=AD2+SUM1
      SUM3=R2+SUM2*SUM2
      RAD1=(SUM3)**1.5D0
      factor= AM2/RAD1
      der_pot= factor*(1- 3*R2/SUM3)
      
      RETURN
    END SUBROUTINE DER_FORDISK

!Second derivative of the halo  
SUBROUTINE DER_FORHALO(X1,Y1,Z1,der_force)
      use BarAndSpiralsInterface
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      AM3= G*MH
      r=SQRT(X1*X1+Y1*Y1+Z1*Z1)
      cte1= r/A3
      cte2= 1+d**CS2
      factor=AM3/(A3*A3*cte2) 
      der_pot= factor*(0.02D0/cte1**0.98D0 -(CS2*d**0.04D0)/cte2 )
     
      RETURN
    END SUBROUTINE DER_FORHALO



!________________________OTHER FUNCTIONS_______________________
!______________________________________________________________

    subroutine potbaref(t,x1,y1,z1,pot)
       !**********************************************************************
       ! returns the effective potential of the model
       !**********************************************************************
       use BarAndSpiralsInterface
       implicit none
       real*8 :: x1,y1,z1,pot,POTA,POTB, t,tmt
       
       CALL TIMEF(T,TIN,TGROW,TMT)
       !   aux1=(totalmass-tmt*barmass)/totalmass
       !aux2=tmt*barmas/totalmass
       aux2=tmt
       !      write(*,*)t,tmt,aux2
       call potax(t,x1,y1,z1,pota)
       if (barmas .gt. 1.D-9) THEN
          call potbar(x1,y1,z1, potb) 
          pot=pota+aux2*potb-0.5D0*omegs2*(x1*x1+ y1*y1)
          !      write(*,*)'aki',x(1),x(2),pot,aux2
       else
          pot=pota-0.5D0*omegs2*(x1*x1+ y1*y1)
       end if
       !      write(*,*)'r',x(1)   
       !      write(*,*)'axi',pota
       !      write(*,*)'bar',potb
       !      write(*,*)'pot',pota+potb
       !      write(*,*)'poteff',pot
       return
     end subroutine potbaref
     

     FUNCTION XLMBD (X,Y,Z,A,B,C)
!*********************************************************************
!       Resolution of the equation for L:
!
!       X/(A+L) + Y/(B+L) + Z/(C+L) = 1
!
!       Restriction :X/A + Y/B + Z/C > 1.
!***********************************************************************
! The input parameters are assumed to be positive
!***********************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (US3 = 1.D0/3.D0)
       DATA AA/-1.D3/,BB/-1.D3/,CC/-1.D3/
       SAVE XA,YA,ZA,AA,BB,CC,XLA
       IF (DABS(A-AA)+DABS(B-BB)+DABS(C-CC).LT.1.D-14) THEN
          IF (DABS(X-XA)+DABS(Y-YA)+DABS(Z-ZA).LT.1.D-14) THEN
             XLMBD=XLA
             RETURN
          ENDIF
       ENDIF
	
       C2 = (A+B+C-X-Y-Z)*US3
       C1 = A*(B-Y-Z) + B*(C-X-Z) + C*(A-X-Y)
       C0 = A*B*(C-Z) - (A*Y+X*B)*C
       P = C1*US3-C2*C2
       Q = C2**3 + (C0-C1*C2)*.5D0
       DE = P**3+Q*Q
	
       IF (DE .LT. 0.D0) THEN
          R = DSQRT(-P)
          XLMBD = 2.D0*R*DCOS( DACOS(-Q/R**3)*US3 ) - C2
       ELSE
          R = -Q+DSQRT(DE)
          IF (R .LT. 0.D0) THEN
             U = -(-R)**US3
          ELSE
             U = R**US3
          END IF
          R = -Q-DSQRT(DE)
          IF (R .LT. 0.D0) THEN
             V = -(-R)**US3
          ELSE
             V = R**US3
          END IF
          XLMBD = U + V - C2
       END IF
       XA=X
       YA=Y
       ZA=Z
       AA=A
       BB=B
       CC=C
       XLA=XLMBD
       
       RETURN
     END FUNCTION XLMBD


 SUBROUTINE TIMEF (T,TIN,TGROW,TMT)
      !----------------------------------------------------------------------
      !----------------------------------------------------------------------
      !     COMPUTES THE TIME FUNCTION: from Dehnen 2000
      !----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION T,TMT,TFAC, TIN, TGROW
      
      IF(T.LE.TIN)TMT=0.D0
      IF(TIN.LT.T.AND.T.LT.TGROW)THEN
         TFAC=2.D0*(T-TIN)/(TGROW-TIN)-1.D0
         TMT=3.D0*(TFAC**5)/16.D0-5.D0*(TFAC**3)/8.D0+ 15.D0*TFAC/16.D0+0.5D0
      ENDIF
      IF(T.GE.TGROW)TMT=1.D0
      !write(*,*)'evolution time',t, "tmt", tmt
      RETURN
    END SUBROUTINE TIMEF



 SUBROUTINE TIMEFBRA (T, TMTBRA)
   !----------------------------------------------------------------------
   !----------------------------------------------------------------------
   !     COMPUTES THE TIME FUNCTION TO MAKE TRANSIENT SPIRALS: 
   !     from De Simone, Wu and TREMAINE 2004
   !     TS defines time in which the spiral is a maximum
      !     SIGMA_S defines the duration of the spiral arm
   !----------------------------------------------------------------------
   use BarAndSpiralsInterface
   implicit none
   double precision:: T,TMTBRA
   integer:: i

   TMTBRA= 0D0

   if( (t .LT. ts-4*sigma_s) .or. (t .GE. ntot*ts +4*sigma_s) )then
      ombra=0
      phi0=0
      tmtbra=0
   else

      do i=1,ntot
         if( (t .LT. i*ts+ 4*sigma_s) .and. (t .GE. i*ts- 4*sigma_s) )then 
            ombra= omegas_sp(i)
            phi0= phi_sp(i) - omegs*(i*ts- 4*sigma_s)
            TMTBRA= aspi*exp(-0.5D0*((t-i*ts)/sigma_s)**2)
            !write(20,*) t, tmtbra
         endif    
      enddo
   endif
   
  
   
   !write(*,*) "in tmtbra:", tmtbra, ombra
    !write(24,110) t, tmtbra
!110 format(7f18.10)
      !write(*,*) "omegs_sp:", omegas_sp

   RETURN
 END SUBROUTINE TIMEFBRA
   
      

 SUBROUTINE ELINT(PHI,XK,F,E)
      !******************************************************************
       ! Incomplete elliptic integrals F & E
       !
       ! Method : descending Landen transformation
       ! Ref.: Abramowitz & Stegun 17.5 p. 597
       !
       ! Input parameters : PHI, XK : amplitude angle and module
       ! 0 = XK & 1
       ! Output values : F, E : incomplete elliptic integrals of the
       ! first and second kind
       !*****************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (PI=3.141592653589793D0, PI2=2.D0*PI)
       
       A = 1.D0
       C = XK
       S = C*C
       B = DSQRT(1.D0-S)
       P = PHI
       P2 = P
       T = 0.D0
       XI = 1.D0
       
10     XI = XI + XI
       P = P + DINT(P/PI+.5D0)*PI + DATAN(DTAN(P2)*B/A)
       P2 = DMOD(P,PI2)
       C = .5D0*(A-B)
       A1= .5D0*(A+B)
       B = DSQRT(A*B)
       A = A1
       S = S + C*C*XI
       T = T + C*DSIN(P2)
       IF (C.GT.1D-15) GOTO 10
       
       F = P/(XI*A)
       E = T + F*(1.D0-.5D0*S)
       
       RETURN
     END SUBROUTINE ELINT
     
!-----------------------------------------------------------------------
! 3D Spiral Arms -- subroutines
!-----------------------------------------------------------------------
!
! Cox & Gomez (2002) [CG02] -- potentail and mass density
!  rho(R,z) ~ rho_0 * exp(-(R-R_0)/R_s) * sech**2(z/H)
!
! subroutines from code mwf9 (v18), LJ, 25 July 2014
!
! parameters in interface notation:
!   am_ssa  :: rho_sp density amplitude [M_o/kpc**3]
!                     am_saa = rho_scaled = rho_0*exp(R_0/R_s) 
!                            = rho_sp*exp(eps_sa*r_sa)
!   eps_sa  :: 1/rsigma  1/(R_s) inverse scale lenght [kpc^-1]
!   al_sa   :: atan(tani)   pitch angle [rad]
!   r_sa    :: r0_sp  fiducial radius at which the rho_0 is measured [kpc]
!   m_sa    :: m      number of arms
!   om_sa   :: ombra  pattern velocity [km/s^2/kpc]
!   phi0_sa :: phi0   initial phase measured at r_0 [rad]
!   z_sa    :: h_sp   scale height     [kpc]
!
!   phi0_sa=0 and time=0 -- rotation to the coordinate sys. of the spirals
!     happens in total_potential and total_force
!
!-----------------------------------------------------------------------
 subroutine potsa_cg(t,rvec,epoo, &
                     am_ssa,eps_sa, al_sa, r_sa,m_sa,om_sa,phi0_sa,z_sa)
!-----------------------------------------------------------------------
! potential for spiral arms -- eqn. (8) CG02

  use BarAndSpiralsInterface
  implicit none

  real(8), intent(in)  :: t,rvec(3)
  real(8), intent(in)  :: am_ssa,eps_sa,al_sa,r_sa,om_sa, &
                           phi0_sa,z_sa
  integer, intent(in)   :: m_sa
  real(8), intent(out) :: epoo
  
  real(8), external :: kn,bn,dn,radi
  
  real(8) :: rad,thet
  real(8) :: kk,bb,dd,sin_a,itan_a,pom,arg,gzet
  
! transformation polar coordinates
  rad = radi(rvec(1),rvec(2))
  thet = atan2(rvec(2),rvec(1))

! auxiliar variables
  sin_a = sin(al_sa)
  itan_a = 1.d+0/tan(al_sa)

  kk = kn(1,m_sa,sin_a,rad)
  bb = bn(1,m_sa,sin_a,z_sa,rad)
  dd = dn(1,m_sa,sin_a,z_sa,rad)
  
  pom = -4.d+0*pi*G*z_sa*am_ssa
  arg = 1*m_sa*(thet - phi0_sa - om_sa*t - log(rad/r_sa)*itan_a)
  gzet = (1.d+0/cosh(kk*rvec(3)/bb))**bb

! eqn. (8)
  epoo = pom * exp(-eps_sa*rad) * (1.d+0/(kk*dd)) * cos(arg) * gzet

 end subroutine potsa_cg

!-----------------------------------------------------------------------
 subroutine fsa_cg(t,rvec,fvec, &
                   am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)
!-----------------------------------------------------------------------
! force from spiral arms model from CG02, transformation to cylindrical coordinates
! rvec is x,y,z vector of position in a frame corotating with the spiral arms
! t=phi0_sa= 0 because of the frame
! am_ssa= rho0_scaled
! eps_sa= eps_sp (1/rsigma)
! al_sa= atan(tani)
! r_sa= r0_sp = fiducial_radius
! m_sa= m   
! om_sa= ombra
! z_sa= h_sp
   

  use BarAndSpiralsInterface
  implicit none

  real(8), intent(in)  :: t,rvec(3)
  real(8), intent(in)  :: am_ssa,eps_sa,al_sa,r_sa,om_sa, &
                           phi0_sa,z_sa
  integer, intent(in)   :: m_sa
  real(8), intent(out) :: fvec(3)
  
  real(8), external :: kn,bn,dn,kndr,bndr,dndr,radi

  real(8) :: rad,thet,fr,ft,phip,sint,cost
  real(8) :: kk,kkdr,bb,bbdr,dd,dddr,sin_a,itan_a,pom,arg,za,secz,tanz

! transformation polar coordinates
  rad = radi(rvec(1),rvec(2))
  thet = atan2(rvec(2),rvec(1))

  sint = sin(thet)
  cost = cos(thet)

! auxiliar variables
  sin_a = sin(al_sa)
  itan_a = 1.d+0/tan(al_sa)

  kk = kn(1,m_sa,sin_a,rad)
  bb = bn(1,m_sa,sin_a,z_sa,rad)
  dd = dn(1,m_sa,sin_a,z_sa,rad)

  kkdr = kndr(1,m_sa,sin_a,rad)
  bbdr = bndr(1,m_sa,sin_a,z_sa,rad)
  dddr = dndr(1,m_sa,sin_a,z_sa,rad)
  
  pom = -4.d+0*pi*G*z_sa*am_ssa
  arg = m_sa*(thet - phi0_sa - om_sa*t - log(rad/r_sa)*itan_a)

  za = kk*rvec(3)/bb
  secz = 1.d+0/cosh(za)
  tanz = tanh(za)

! potential
  call potsa_cg(t,rvec,phip, &
                am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

! tangential force
  ft = (1.d+0/rad) * m_sa * &
       pom * exp(-eps_sa*rad) * (1.d+0/(kk*dd)) * sin(arg) * secz**bb

! radial force
  fr = (-phip) * (-eps_sa - kkdr/kk - kkdr*rvec(3)*tanz - dddr/dd + &
                 (log(secz)+tanz*za)*bbdr + m_sa*itan_a/rad*tan(arg))

! z-force
  fvec(3) = phip * tanz * kk

! transformation to F_x, F_y
  fvec(1) = fr*cost - ft*sint
  fvec(2) = fr*sint + ft*cost

 end subroutine fsa_cg

!-----------------------------------------------------------------------
 subroutine dsa_cg_num(t,rvec,densa, &
                     am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)
!-----------------------------------------------------------------------
! mass density of SA from CG02 model
!  calculate numerically from Poisson's equation
  
  use BarAndSpiralsInterface
  implicit none

  real(8), intent(in)  :: t,rvec(3)
  real(8), intent(in)  :: am_ssa,eps_sa,al_sa,r_sa,om_sa, &
                           phi0_sa,z_sa
  integer, intent(in)   :: m_sa
  real(8), intent(out) :: densa

  real(8) :: f1(3),f2(3),dpx,dpy,dpz,hh,rh(3)

! hh = rvec(1)*(epsilon(sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2))**0.25)

! d(F_x)/d(x)
  hh = rvec(1)*(epsilon(rvec(1))**0.25)

  rh = (/rvec(1)+hh,rvec(2),rvec(3)/)
  call fsa_cg(t,rh,f1, &
              am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

  rh = (/rvec(1)-hh,rvec(2),rvec(3)/)
  call fsa_cg(t,rh,f2, &
              am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

  dpx = -(f1(1)-f2(1))/(2.d+0*hh)

! d(F_y)/d(y)
  hh = rvec(2)*(epsilon(rvec(2))**0.25)

  rh = (/rvec(1),rvec(2)+hh,rvec(3)/)
  call fsa_cg(t,rh,f1, &
              am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

  rh = (/rvec(1),rvec(2)-hh,rvec(3)/)
  call fsa_cg(t,rh,f2, &
              am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

  dpy = -(f1(2)-f2(2))/(2.d+0*hh)

! d(F_z)/d(z)
  hh = rvec(1)*(epsilon(rvec(1))**0.25)

  rh = (/rvec(1),rvec(2),rvec(3)+hh/)
  call fsa_cg(t,rh,f1, &
              am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

  rh = (/rvec(1),rvec(2),rvec(3)-hh/)
  call fsa_cg(t,rh,f2, &
              am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)

  dpz = -(f1(3)-f2(3))/(2.d+0*hh)

! Poisson
  densa = (dpx + dpy + dpz)/(4.d+0*pi*G)

 end subroutine dsa_cg_num

!-----------------------------------------------------------------------
 subroutine dsa_cg_analit(t,rvec,densa, &
                          am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)
!-----------------------------------------------------------------------
! mass density of SA from CG02 model
!  eqs. (1)--(3) from CG02
!  rho_A(r,z) = rho_0 * exp(-(r-r0)/Rs) * sech**2(z/h)
!  rho(r,phi,z) = rho_A(r,z) cos(gamma)
!  gamma = N*[phi - phi_p(r0) - ln(r/r0)/tan(alpha)]
  
  use BarAndSpiralsInterface
  implicit none

  real(8), intent(in)  :: t,rvec(3)
  real(8), intent(in)  :: am_ssa,eps_sa,al_sa,r_sa,om_sa, &
                           phi0_sa,z_sa
  integer, intent(in)  :: m_sa
  real(8), intent(out) :: densa
  
  real(8), external :: radi

  real(8) :: f1(3),f2(3),dpx,dpy,dpz,hh,rh(3)
  real(8) :: rad,thet,z
  real(8) :: gammma,rho_a

! transformation polar coordinates
  rad = radi(rvec(1),rvec(2))
  thet = atan2(rvec(2),rvec(1))
  z = rvec(3)

! eq. (3) + rotation (om_sa*t)
  gammma = m_sa*(thet - phi0_sa - om_sa*t - log(rad/r_sa)/tan(al_sa))
!  gammma = m_sa*(thet - om_sa*t - log(rad/r_sa)/tan(al_sa))

! eq. (1)
  rho_a = am_ssa*exp(-rad*eps_sa)*(1.d+0/cosh(z/z_sa))**2

! eq. (2)
  densa = rho_a*cos(gammma)

 end subroutine dsa_cg_analit

!-----------------------------------------------------------------------
 subroutine dsa_cg_anal10(t,rvec,densa, &
                        am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)
!-----------------------------------------------------------------------
! mass density of SA from CG02 model
!  eq. (10) from CG02
  
  use BarAndSpiralsInterface
  implicit none

  real(8), intent(in)  :: t,rvec(3)
  real(8), intent(in)  :: am_ssa,eps_sa,al_sa,r_sa,om_sa, &
                           phi0_sa,z_sa
  integer, intent(in)   :: m_sa
  real(8), intent(out) :: densa
  
  real(8), external :: kn,bn,dn,radi

  real(8) :: f1(3),f2(3),dpx,dpy,dpz,hh,rh(3)
  real(8) :: rad,thet,z
  real(8) :: gammma,rho_a
  real(8) :: kk,bb,dd,sin_a,itan_a,pom,arg,gzet

! transformation polar coordinates
  rad = radi(rvec(1),rvec(2))
  thet = atan2(rvec(2),rvec(1))
  z = rvec(3)

! auxiliar variables
  sin_a = sin(al_sa)
  itan_a = 1.d+0/tan(al_sa)

  kk = kn(1,m_sa,sin_a,rad)
  bb = bn(1,m_sa,sin_a,z_sa,rad)
  dd = dn(1,m_sa,sin_a,z_sa,rad)

  gzet = (1.d+0/cosh(kk*rvec(3)/bb))**(bb+2.d+0)

! eq. (3) + rotation (om_sa*t)
  gammma = m_sa*(thet - phi0_sa - om_sa*t - log(rad/r_sa)*itan_a)
!  gammma = m_sa*(thet - om_sa*t - log(rad/r_sa)*itan_a)

  rho_a = am_ssa*exp(-rad*eps_sa)* (kk*z_sa/dd)*(bb+1.d+0)/bb *gzet

  densa = rho_a*cos(gammma)

 end subroutine dsa_cg_anal10

!-----------------------------------------------------------------------
 subroutine dsa_cg_anal_app(t,rvec,densa, &
                        am_ssa,eps_sa,al_sa,r_sa,m_sa,om_sa,phi0_sa,z_sa)
!-----------------------------------------------------------------------
! mass density of SA from CG02 model
!  eqs. from appendix

  use BarAndSpiralsInterface
  implicit none

  real(8), intent(in)  :: t,rvec(3)
  real(8), intent(in)  :: am_ssa,eps_sa,al_sa,r_sa,om_sa, &
                           phi0_sa,z_sa
  integer, intent(in)  :: m_sa
  real(8), intent(out) :: densa

  real(8), external :: kn,bn,dn,radi

  real(8) :: f1(3),f2(3),dpx,dpy,dpz,hh,rh(3)
  real(8) :: rad,thet,z
  real(8) :: gammma,rho_a
  real(8) :: kk,bb,dd,sin_a,itan_a,pom,arg,gzet,zet
  real(8) :: e_app,re_app,par1

! transformation polar coordinates
  rad = radi(rvec(1),rvec(2))
  thet = atan2(rvec(2),rvec(1))
  z = rvec(3)

! auxiliar variables
  sin_a = sin(al_sa)
  itan_a = 1.d+0/tan(al_sa)

  kk = kn(1,m_sa,sin_a,rad)
  bb = bn(1,m_sa,sin_a,z_sa,rad)
  dd = dn(1,m_sa,sin_a,z_sa,rad)

  zet = kk*z/bb
  gzet = (1.d+0/cosh(zet))
  par1 = (1.d+0+0.3*kk*z_sa)

! eq. (3) + rotation (om_sa*t)
  gammma = m_sa*(thet - phi0_sa - om_sa*t - log(rad/r_sa)*itan_a)
! gammma = m_sa*(thet - om_sa*t - log(rad/r_sa)*itan_a)

  e_app = 1.d+0 + kk*z_sa/dd*(1.d+0-0.3/par1**2) - rad*eps_sa - &
          (kk*z_sa)*(1.d+0+0.8*kk*z_sa)*log(gzet) - &
          0.4*((kk*z_sa)**2)*zet*tanh(zet)

  re_app = -kk*z_sa/dd*(1.d+0 - 0.3*(1.d+0-0.3*kk*z_sa)/par1**3) + &
            (kk*z_sa/dd*(1.d+0-0.3/par1**2))**2 - rad*eps_sa + &
            kk*z_sa*(1.d+0+1.6*kk*z_sa)*log(gzet) - ((0.4*((kk*z_sa)**2)*zet*gzet)**2)/bb + &
            1.2*(kk*z_sa)**2*zet*tanh(zet)

  densa = am_ssa*z_sa/(dd*rad)*exp(-rad*eps_sa)*(gzet**bb) * &
         ((kk*rad*(bb+1.d+0)/bb*(gzet**2) - 1/(kk*rad)*(e_app**2+re_app))*cos(gammma) - &
          2.d+0*e_app*cos(al_sa)*sin(gammma))

 end subroutine dsa_cg_anal_app

!-----------------------------------------------------------------------
 real(8) function kn(nn,mm,sin_a,rad)
!-----------------------------------------------------------------------
! function K_n -- eqn. (5) in CG02

  implicit none

  real(8), intent(in) :: sin_a,rad
  integer, intent(in)  :: nn,mm
 
  kn = nn*mm/(rad*sin_a)

 end function kn

!-----------------------------------------------------------------------
 real(8) function bn(nn,mm,sin_a,z_sa,rad)
!-----------------------------------------------------------------------
! function beta_n -- eqn. (6) in CG02

  implicit none

  real(8), intent(in) :: sin_a,z_sa,rad
  integer, intent(in)  :: nn,mm

  real(8), external :: kn

  real(8) :: kk

  kk = kn(nn,mm,sin_a,rad)

  bn = kk*z_sa*(1.d+0 + 0.4d+0*kk*z_sa)

 end function bn

!-----------------------------------------------------------------------
 real(8) function dn(nn,mm,sin_a,z_sa,rad)
!-----------------------------------------------------------------------
! function D_n -- eqn. (7) in CG02

  implicit none

  real(8), intent(in) :: sin_a,z_sa,rad
  integer, intent(in)  :: nn,mm

  real(8), external :: kn

  real(8) :: kk

  kk = kn(nn,mm,sin_a,rad)

  dn = (1.d+0 + kk*z_sa + 0.3d+0*(kk*z_sa)**2)/(1.d+0 + 0.3d+0*kk*z_sa)

 end function dn

!-----------------------------------------------------------------------
 real(8) function kndr(nn,mm,sin_a,rad)
!-----------------------------------------------------------------------
! d(K_n)/d(R)

  implicit none

  real(8), intent(in) :: sin_a,rad
  integer, intent(in)  :: nn,mm

  kndr = -nn*mm/(sin_a*rad**2)

 end function kndr

!-----------------------------------------------------------------------
 real(8) function bndr(nn,mm,sin_a,z_sa,rad)
!-----------------------------------------------------------------------
! d(Beta_n)/d(R)

  implicit none

  real(8), intent(in) :: sin_a,z_sa,rad
  integer, intent(in)  :: nn,mm

  real(8), external :: kn,kndr

  real(8) :: kk,kkdr

  kk = kn(nn,mm,sin_a,rad)
  kkdr = kndr(nn,mm,sin_a,rad)

  bndr = 0.8d+0*z_sa**2*kk*kkdr + z_sa*kkdr

 end function bndr

!-----------------------------------------------------------------------
 real(8) function dndr(nn,mm,sin_a,z_sa,rad)
!-----------------------------------------------------------------------
! d(D_n)/d(R)

  implicit none

  real(8), intent(in) :: sin_a,z_sa,rad
  integer, intent(in)  :: nn,mm

  real(8), external :: kn,kndr

  real(8) :: kk,kkdr

  kk = kn(nn,mm,sin_a,rad)
  kkdr = kndr(nn,mm,sin_a,rad)

  dndr = kkdr * (0.7d+0*z_sa + 0.6d+0*kk*z_sa**2 + 0.09*kk**2*z_sa**3) / &
                (1.d+0 + 0.3d+0*kk*z_sa)**2

 end function dndr
 
!-----------------------------------------------------------------------
 real(8) function radi(rx,ry)
!-----------------------------------------------------------------------
! radial distance in xy plane

 real(8), intent(in) :: rx,ry

 radi = sqrt(rx**2 + ry**2)

 end function radi

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

     
      
    subroutine velcirc(x1,y1,z1, vc)
       !-----------------------------------------------------------------------
       ! circular velocity corresponding to the total potential
       ! xflag = 1 (x1, y1, z1) are in a corotating system 
       ! xflag =2 (x1, y1, z1) are in an inertial system
       !-----------------------------------------------------------------------
       use BarAndSpiralsInterface
       implicit none
       real*8:: x1,y1,z1, vc, T, PX,PY,PZ,R, force
       real*8:: angle, xr,yr,zr, fx,fy,fz
       
       T=0.D0
       if (xflag .EQ. 1) then 
          call total_force(t,x1,y1,z1,px,py,pz)
          r=dsqrt(x1*x1 + Y1*Y1 + Z1*Z1)
          force= dsqrt(PX*PX+PY*PY+PZ*PZ)
          vc=dsqrt(force*r)

       elseif (xflag .EQ. 2) then
          ! inertial -> rotating
          angle= initial_phase + omega_sys*t
          xr= x1*cos(angle) + y1*sin(angle)
          yr= -x1*sin(angle) + y1*cos(angle)
          zr= z1
          call total_force(t,xr,yr,zr,px,py,pz)
          fx= px*cos(angle) -py*sin(angle)
          fy= px*sin(angle) + py*cos(angle)
          fz= pz
          r=dsqrt(x1*x1 + y1*y1 + z1*z1)
          force= dsqrt(fx*fx+ fy*fy+ fz*fz)
          vc=dsqrt(force*r)
       endif

       RETURN
     END subroutine velcirc   

    subroutine velcirc_axi(x1,y1,z1, vc)
       !-----------------------------------------------------------------------
       ! circular velocity corresponding to the A&S 
       !-----------------------------------------------------------------------
       use BarAndSpiralsInterface
       implicit none
       real*8:: x1,y1,z1, vc, T, PX,PY,PZ,R, force
       
       T=0.D0
       !write(*,*) "pos in sub", x1,y1,z1
       call forax(t,x1,y1,z1,px,py,pz)
       !write(*,*) "f in sub:", px, py, pz
       r=dsqrt(x1*x1 + Y1*Y1 + Z1*Z1)
       force= dsqrt(PX*PX+PY*PY+PZ*PZ)
       vc=dsqrt(force*r)
       
       RETURN
     END subroutine velcirc_axi

    subroutine epic_freq(x1,y1,z1, k)
       !-----------------------------------------------------------------------
       !epicyclic frecuency corresponding to the A&S 
       !-----------------------------------------------------------------------
       use BarAndSpiralsInterface
       implicit none
       real*8:: x1,y1,z1,k, T, omega, der_force
       
       T=0.D0
       call der_forax (t,x1, y1, z1, der_force)
       call velcirc_axi(x1,y1,z1, omega)
       k= dsqrt(der_force+ 3*omega*omega)
            
       RETURN
     END subroutine epic_freq


     subroutine compute_tidal_derivatives(t, x1,y1,z1, FFxx, FFyy, FFzz, FFxy, FFxz, FFyz)
       use BarAndSpiralsInterface
       implicit none
       real*8:: x1,y1,z1,dx, dy, dz, t
       real*8:: FFxx, FFyy, FFzz, FFxy, FFxz, FFyz
       real*8:: pot, pot1xx, pot2xx, pot1yy, pot2yy, pot1zz, pot2zz
       real*8:: pot1xy, pot2xy, pot1xz, pot2xz, pot1yz, pot2yz
       dx= 0.0001
       dy= dx
       dz= dx
       
       ! computation Fxx, Fyy, Fzz
       call total_potential(T, x1, y1, z1, pot)
       call total_potential(T, (x1+dx), y1, z1, pot1xx)
       call total_potential(T, (x1-dx), y1, z1, pot2xx)
       FFxx= (pot1xx-2*pot+ pot2xx)/dx**2
       call total_potential(T, x1, (y1+dy), z1, pot1yy)
       call total_potential(T, x1, (y1-dy), z1, pot2yy)
       FFyy= (pot1yy-2*pot+ pot2yy)/dy**2
       call total_potential(T, x1, y1, (z1+dz), pot1zz)
       call total_potential(T, x1, y1, (z1-dz), pot2zz)
       FFzz= (pot1zz-2*pot+ pot2zz)/dz**2
       ! Computation of the cross derivatives
       call total_potential(T, (x1+dx), (y1+dy), z1, pot1xy)
       call total_potential(T, (x1-dx), (y1-dy), z1, pot2xy)
       FFxy= (pot1xy-pot1xx-pot1yy+2*pot-pot2xx-pot2yy+pot2xy)/(2*dx*dy)
       call total_potential(T, (x1+dx), y1, (z1+dz), pot1xz)
       call total_potential(T, (x1-dx), y1, (z1-dz), pot2xz)
       FFxz= (pot1xz-pot1xx-pot1zz +2*pot-pot2xx-pot2zz +pot2xz)/(2*dx*dz)
       call total_potential(T, x1, (y1+dy), (z1+dz), pot1yz)
       call total_potential(T, x1, (y1-dy), (z1-dz), pot2yz)
       FFyz= (pot1yz-pot1yy-pot1zz +2*pot-pot2yy-pot2zz +pot2yz)/(2*dy*dz)
       
       return
     end subroutine compute_tidal_derivatives

     
     subroutine tidal_tensor(t, x1,y1,z1, Fxx, Fyx, Fzx, Fxy,Fyy,Fzy, Fxz, Fyz, Fzz)
       !-----------------------------------------------------------------------
       ! This function computes the components of the tidal tensor
       ! Tij at coordinates (t, x1, y1, z1)
       ! by calculating the second derivative of the potential
       !flag = 1 (x1, y1, z1) are in a corotating system 
       !flag =2 (x1, y1, z1) are in an inertial system
       !-----------------------------------------------------------------------
       use BarAndSpiralsInterface
       implicit none
       
       real*8:: x1,y1,z1,dx, dy, dz, t, angle
       real*8:: Fxx, Fyy, Fzz, Fxy, Fxz, Fyz
       real*8:: Fyx, Fzx, Fzy
       real*8:: Fxx_r, Fyy_r, Fzz_r, Fxy_r, Fxz_r, Fyz_r
       real*8:: xr, yr, zr, dxr, dyr, dzr
       real*8:: cos2theta, sin2theta, sincostheta
       
       if (xflag .EQ. 1) then 
          call compute_tidal_derivatives(t, x1,y1,z1, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz)
          
       elseif (xflag .EQ. 2) then
          ! inertial -> rotating
          angle= initial_phase + omega_sys*t
          xr= x1*cos(angle) + y1*sin(angle)
          yr= -x1*sin(angle) + y1*cos(angle)
          zr= z1

          call compute_tidal_derivatives(t, xr,yr,zr, Fxx_r, Fyy_r, Fzz_r, Fxy_r, Fxz_r, Fyz_r)
          
          !conversion tidal tensor to the inertial frame
          cos2theta= (1.d0+ cos(2.d0*angle))/2.d0
          sin2theta= (1.d0- cos(2.d0*angle))/2.d0
          sincostheta= cos(angle)*sin(angle)
          
          Fxx= Fxx_r*cos2theta- 2.d0*Fxy_r*sincostheta +Fyy_r*sin2theta
          Fyy= 2*Fxy_r*sincostheta +Fxx_r*sin2theta +Fyy_r*cos2theta
          Fzz= Fzz_r
          Fxy= (Fxx_r-Fyy_r)*sincostheta +Fxy_r*cos(2.d0*angle)
          Fxz= Fxz_r*cos(angle) - Fyz_r*sin(angle)
          Fyz= Fxz_r*sin(angle) + Fyz_r*cos(angle)
       else
          !write(*,*) 'Choose flag= 1 or 2'
       end if
       
       ! Build the symmetric components
       Fyx= Fxy
       Fzx= Fxz
       Fzy= Fyz
       

       RETURN
     END subroutine tidal_tensor


     subroutine eigen_values(t, x1,y1,z1, l1, l2, l3)
       !--------------------------------------------------------
       ! Compute the eigen values of the tidal tensor       !
       !flag = 1 computation in a corotating frame with the bar ->
       ! x1, y1, z1 are defined in that system
       !flag =2 computation in an inertial frame
       !--------------------------------------------------------
       use BarAndSpiralsInterface
       implicit none
       integer, parameter:: N=3, LDVL=1, LDVR=1 
       integer :: LDA, ok
       real*8:: t, x1,y1,z1
       real*8:: Tij(N,N)
       real*8:: wr(N), wi(N), VL(LDVL,N), VR(LDVR,N), WORK(3*N) !eigenvalues and vectors
       real*8:: Foxx, Foyy, Fozz, Foxy, Foyx, Foxz, Fozx, Foyz, Fozy
       real*8:: l1, l2, l3

       ! Build the tidal tensor. First columns and then arrows
       call  tidal_tensor(t, x1,y1,z1,Foxx, Foyx, Fozx, Foxy,Foyy,Fozy, Foxz, Foyz, Fozz)
       Tij(1,1)= Foxx
       Tij(1,2)= Foyx
       Tij(1,3)= Fozx
       Tij(2,1)= Foxy
       Tij(2,2)= Foyy
       Tij(2,3)= Fozy
       Tij(3,1)= Foxz
       Tij(3,2)= Foyz
       Tij(3,3)= Fozz

       
       ! Compute eigen values of Tij:
       LDA = max(1,N)
       call DGEEV('N', 'N', N, Tij, LDA, wr, wi, VL, LDVL, VR, LDVR, WORK, 3*N, ok)
       l1= wr(1)
       l2= wr(2)
       l3= wr(3)
       
       return
     end subroutine eigen_values
     
     subroutine tidal_radius(t, x1,y1,z1, Mc, rt)
       !--------------------------------------------------------
       ! Compute the tidal radius when the cluster center is at x,y,z
       ! at a time t. x,y,z are in a NON ROTATING FRAME
       !Mc is an input corresponding to
       !the initial mass of the cluster 
       !--------------------------------------------------------
       use BarAndSpiralsInterface
       implicit none
       real*8:: t, x1,y1,z1, rt, Mc, phi_bar
       real*8:: lambda1, lambda2, lambda3, lambdamax
       
       call eigen_values(t, x1,y1,z1, lambda1, lambda2, lambda3)
       lambdamax= max(lambda1, lambda2, lambda3)

       rt= (G*Mc/lambdamax)**(1.d0/3.d0)

       return
     end subroutine tidal_radius



    
