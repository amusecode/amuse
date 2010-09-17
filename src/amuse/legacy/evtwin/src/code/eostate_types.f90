module eostate_types

   use real_kind
   
   implicit none
   
   ! User defined data type used to hold information about abundances
   type :: abundance
       real(double) :: xa(9)           ! Abundance by baryon fraction (usually input)
       real(double) :: na(9)           ! Abundance by particle number fraction
       real(double) :: neo             ! Number of electrons per baryon (full ion.)
       real(double) :: nio             ! Ion fraction
       real(double) :: nzz             ! Average charge^2
       real(double) :: avm             ! Mean molecular weight [amu]
       real(double) :: ne              ! Actual number of electrons per baryon
       real(double) :: xai(9, 26)      ! Ionisation states for all elements
   end type abundance

   ! User defined data type used to hold output from the equation-of-state
   ! calculation.
   ! Parts of this are filled in by estate, parts by nucrat
   type :: eostate
      real(double) :: ap         ! ln P, P in dyne/cm^2 (=1/10 Pa)
      real(double) :: p          ! Pressure, [dyne/cm^2]
      real(double) :: pr         ! Radiation pressure [dyne/cm^2]
      real(double) :: pg         ! Gas pressure (=P-Pr) [dyne/cm^2]
      real(double) :: pf         ! dlnP/dlnf at constant T [-]
      real(double) :: pt         ! dlnP/dlnT at constant F [-]

      real(double) :: rho        ! Density, [g/cm^3]
      real(double) :: arho       ! ln rho, rho in g/cm^3
      real(double) :: rf         ! dlnrho/dlnf at constant T [-]
      real(double) :: rt         ! dlnrho/dlnT at constant F [-]

      real(double) :: s          ! Specific entropy [erg/K/g]
      real(double) :: sf         ! dS/dlnf at T constant [erg/K/g]
      real(double) :: st         ! dS/dlnT at f constant [erg/K/g]

      real(double) :: scp        ! Specific heat at constant pressure [erg/K/g]
      real(double) :: grada      ! Adiabatic temperature "gradient", dlnT/dlnP
      real(double) :: delta      ! dlnrho/dlnT at constant P
      real(double) :: phi        ! dlnrho/dlnmu at constant P and T (mu is mean molecular weight)
      real(double) :: gamma1

      real(double) :: u          ! Specific internal energy, erg/g
      real(double) :: uf         ! dU/dlnf at constant T [erg/g]
      real(double) :: ut         ! dU/dlnT at constant f [erg/g]

      real(double) :: T          ! Temperature [K]
      real(double) :: wmu        ! Mean molecular weight [amu]
      real(double) :: aT         ! ln T
      real(double) :: af         ! ln f

      real(double) :: fk         ! Opacity, [cm^2/g]
      real(double) :: fkt        !> \todo dlnfk/dlnT      [FIXME: set to -3.5 for Kramers]
      real(double) :: fkr        !> \todo dlnfk/dlnrho    [FIXME: set to  1.0 for Kramers]
      real(double) :: xhi        ! Thermal diffusion coefficient [cm^2/s]
      real(double) :: nu         ! Viscosity [cm^2/s]
      real(double) :: prandtl    ! Prandtl number, nu/xhi (nu is viscosity)

      real(double) :: zt         ! ? (used for screening corrections)
      real(double) :: dv         ! Electron chemical potential, corrected for pressure ionisation and Coulomb interactions

      real(double) :: rpp        ! Rate of p(pp;nu,g,e+)D reaction per unit volume [s^-1cm^-3]
      real(double) :: r33        ! Rate of He3(He3;2p)a reaction per unit volume [s^-1cm^-3]
      real(double) :: r34        ! Rate of He3(a;g)Be7 reaction per unit volume [s^-1cm^-3]
      real(double) :: rbe        ! Rate of Be7(e-;nu)Li7(p;a)a reaction per unit volume [s^-1cm^-3]
      real(double) :: rbp        ! Rate of Be7(p;a)a reaction per unit volume [s^-1cm^-3]
      real(double) :: rpc        ! Rate of C12(p;e+,nu)C13(p;g)N14 reaction per unit volume [s^-1cm^-3]
      real(double) :: rpna       ! Rate of N14(p;e+,nu)N15(p,a)C12 reaction per unit volume [s^-1cm^-3]
      real(double) :: rpo        ! Rate of O16(p;e+,nu)O17(p,a)N14 reaction per unit volume [s^-1cm^-3]
      real(double) :: r3a        ! Rate of a(aa;g)C12 reaction per unit volume [s^-1cm^-3]
      real(double) :: rac        ! Rate of C12(a;g)O16 reaction per unit volume [s^-1cm^-3]
      real(double) :: ran        ! Rate of N14(3a/2;e+,nu)Ne20 reaction per unit volume [s^-1cm^-3]
      real(double) :: rao        ! Rate of O16(a;g)Ne20 reaction per unit volume [s^-1cm^-3]
      real(double) :: rane       ! Rate of Ne20(a;g)Mg24 reaction per unit volume [s^-1cm^-3]
      real(double) :: rcca       ! Rate of C12(C12;ag)Ne20 reaction per unit volume [s^-1cm^-3]
      real(double) :: rco        ! Rate of C12(O16;ag)Mg24 reaction per unit volume [s^-1cm^-3]
      real(double) :: roo        ! Rate of O16(O16;aa)Mg24 reaction per unit volume [s^-1cm^-3]
      real(double) :: rgne       ! Rate of Ne20(g;a)O16 reaction per unit volume [s^-1cm^-3]
      real(double) :: rgmg       ! Rate of Mg24(g;a)Ne20 reaction per unit volume [s^-1cm^-3]
      real(double) :: rccg       ! Rate of C12(C12;g)Mg24 reaction per unit volume [s^-1cm^-3]
      real(double) :: rpng       ! Rate of N14(p;e+,nu)N15(p,g)O16 reaction per unit volume [s^-1cm^-3]

      real(double) :: ex         ! Total energy generated by nuclear reactions [erg/g/s]
      real(double) :: enx        ! Total energy lost to neutrinos in reactions [erg/g/s]
      real(double) :: en         ! Total energy lost to plasma neutrinos [erg/g/s]
      real(double) :: ext        ! dlnex/dlnT
   end type eostate

end module eostate_types
