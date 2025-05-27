module reaction_rate_data
   use real_kind
   
   implicit none
   integer, parameter :: nuc_num_rates = 92

   real(double) :: crt(200,20)      ! Reaction rates, for structure code
   real(double) :: qrt(20)          ! Reaction rate Q values [MeV?]
   real(double) :: qnt(20)          ! Neutrino losses for reactions [MeV?]

   !     Nucleosynthesis reaction rates
   real(double), save :: nucsyn_crt(200,nuc_num_rates)
   real(double), save :: nucsyn_nrt(200,45)

   ! Reaction-dependent charge parameters for electron screening according
   ! to Graboske, Dewitt, Grossman & Cooper (1973),
   ! for strong (za, zb, zc) are intermediate screening (zd).
   real(double) :: cza(92)
   real(double) :: czb(92)
   real(double) :: czc(92)
   real(double) :: czd(92)
   real(double) :: vz(9)
end module reaction_rate_data

