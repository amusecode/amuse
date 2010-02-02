
      subroutine initialize(z_in,
     &			    neta_in, bwind_in, hewind_in, 
     &              alpha1_in, CElambda_in,
     &              ceflag_in, tflag_in,
     &			    ifflag_in, wdflag_in, bhflag_in,
     &			    nsflag_in, mxns_in,
     &			    pts1_in, pts2_in, pts3_in,
     &              sigma_in, beta_in, xi_in, acc2_in,
     &              epsnov_in, eddfac_in, gamma_in,
     &     		    status)

      implicit none
      real*8 z_in, z, zpars(20)
      real*8 neta_in, bwind_in, hewind_in
      real*8 alpha1_in, CElambda_in
      integer ceflag_in, tflag_in
      integer ifflag_in, wdflag_in, bhflag_in
      integer nsflag_in, mxns_in
      real*8 pts1_in, pts2_in, pts3_in
      real*8 sigma_in, beta_in, xi_in, acc2_in
      real*8 epsnov_in, eddfac_in, gamma_in
      integer status
      include 'src/const_bse.h'
      common /SSE_init/ z, zpars

c     Input parameters are passed from MUSE, rather than being read here.

      z = z_in
      neta = neta_in
      bwind = bwind_in
      hewind = hewind_in
      alpha1 = alpha1_in
      lambda = CElambda_in
      ceflag = ceflag_in
      tflag = tflag_in
      ifflag = ifflag_in
      wdflag = wdflag_in
      bhflag = bhflag_in
      nsflag = nsflag_in
      mxns = mxns_in
      pts1 = pts1_in
      pts2 = pts2_in
      pts3 = pts3_in
      sigma = sigma_in
      beta = beta_in
      xi = xi_in
      acc2 = acc2_in
      epsnov = epsnov_in
      eddfac = eddfac_in
      gamma = gamma_in
      
      call zcnsts(z, zpars)
      if(idum.gt.0) idum = -idum

      status = 0
      return
      end


      subroutine evolve_binary(type1,type2,initial_mass1,initial_mass2,
     &         mass1, mass2, radius1, radius2, luminosity1, luminosity2,
     &         core_mass1, core_mass2, core_radius1, core_radius2,
     &         envelope_mass1, envelope_mass2, envelope_radius1,
     &         envelope_radius2, spin1, spin2, epoch1, epoch2,
     &         MS_lifetime1, MS_lifetime2, age,
     &         orbital_period, eccentricity, end_time)
      
      implicit none
      integer,intent(inout) :: type1, type2
      real*8,intent(inout) :: initial_mass1,initial_mass2,mass1, mass2
      real*8,intent(inout) :: radius1, radius2, luminosity1, luminosity2
      real*8,intent(inout) :: core_mass1, core_mass2, core_radius1
      real*8,intent(inout) :: core_radius2,envelope_mass1,envelope_mass2
      real*8,intent(inout) :: envelope_radius1,envelope_radius2, spin1
      real*8,intent(inout) :: spin2,epoch1,epoch2,MS_lifetime1
      real*8,intent(inout) :: MS_lifetime2, age, end_time
      real*8,intent(inout) :: orbital_period, eccentricity
      real*8 dtp, z, zpars(20)
      integer kstar(2)
      real*8 mass0(2),mass(2),rad(2),lum(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 ospin(2),epoch(2),tms(2)
      common /SSE_init/ z, zpars
      
      kstar(1) = type1
      kstar(2) = type2
      mass0(1) = initial_mass1
      mass0(2) = initial_mass2
      mass(1) = mass1
      mass(2) = mass2
      rad(1) = radius1
      rad(2) = radius2
      lum(1) = luminosity1
      lum(2) = luminosity2
      massc(1) = core_mass1
      massc(2) = core_mass2
      radc(1) = core_radius1
      radc(2) = core_radius2
      menv(1) = envelope_mass1
      menv(2) = envelope_mass2
      renv(1) = envelope_radius1
      renv(2) = envelope_radius2
      ospin(1) = spin1
      ospin(2) = spin2
      epoch(1) = epoch1
      epoch(2) = epoch2
      tms(1) = MS_lifetime1
      tms(2) = MS_lifetime2
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the bcm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = end_time+1
!      dtp = end_time
!      dtp = 0.0d0
!      dtp = age+1

c      call evolv2(type1,type2,initial_mass1,initial_mass2, 
c     &         mass1, mass2, radius1, radius2, luminosity1, luminosity2, 
c     &         core_mass1, core_mass2, core_radius1, core_radius2, 
c     &         envelope_mass1, envelope_mass2, envelope_radius1, 
c     &         envelope_radius2, spin1, spin2, epoch1, epoch2, 
c     &         MS_lifetime1, MS_lifetime2, age, end_time, 
c     &         dtp,z,zpars,
c     &         orbital_period, eccentricity)

      call evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &            menv,renv,ospin,epoch,tms,
     &            age, end_time,dtp,z,zpars,
     &            orbital_period, eccentricity)

      type1 = kstar(1)
      type2 = kstar(2)
      initial_mass1 = mass0(1)
      initial_mass2 = mass0(2)
      mass1 = mass(1)
      mass2 = mass(2)
      radius1 = rad(1)
      radius2 = rad(2)
      luminosity1 = lum(1)
      luminosity2 = lum(2)
      core_mass1 = massc(1)
      core_mass2 = massc(2)
      core_radius1 = radc(1)
      core_radius2 = radc(2)
      envelope_mass1 = menv(1)
      envelope_mass2 = menv(2)
      envelope_radius1 = renv(1)
      envelope_radius2 = renv(2)
      spin1 = ospin(1)
      spin2 = ospin(2)
      epoch1 = epoch(1)
      epoch2 = epoch(2)
      MS_lifetime1 = tms(1)
      MS_lifetime2 = tms(2)

      return
      end
      


      subroutine get_time_step(type1, type2, initial_mass1, 
     &      initial_mass2, mass1, mass2, MS_lifetime1,
     &      MS_lifetime2, epoch1, epoch2, age, time_step)
      
cf2py intent(out) dt
cf2py intent(in) kw, mass, age, mt, tm, epoch
      implicit none
      integer type1, type2
      real*8 initial_mass1, initial_mass2, mass1, mass2
      real*8 MS_lifetime1, MS_lifetime2, epoch1, epoch2
      real*8 age, time_step
      real*8 dtm1, dtr1, dtm2, dtr2
      real*8 tscls(20), lums(10), GB(10), tn
      real*8 z, zpars(20)
      common /SSE_init/ z, zpars

!     Call star fuction to get stellar parameters
      call star(type1, initial_mass1, mass1, MS_lifetime1,
     &      tn, tscls, lums, GB, zpars)
!     Call deltat function to get next timestep
      call deltat(type1, age-epoch1, MS_lifetime1,
     &      tn, tscls, dtm1, dtr1)

!     Call star fuction to get stellar parameters
      call star(type2, initial_mass2, mass2, MS_lifetime2,
     &      tn, tscls, lums, GB, zpars)
!     Call deltat function to get next timestep
      call deltat(type2, age-epoch2, MS_lifetime2,
     &      tn, tscls, dtm2, dtr2)
      
      time_step = min(dtr1, dtm1, dtr2, dtm2)

      return
      end
