module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew ! State vector and its update
  real(dp),allocatable,dimension(:)  ::divu,enew ! Non conservative variables
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0
  real(dp)::ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  integer::nbins
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0e+30
  real(dp)::zero = 0.0
  real(dp)::one = 1.0
  real(dp)::two = 2.0
  real(dp)::three = 3.0
  real(dp)::four = 4.0
  real(dp)::two3rd = 0.6666666666666667
  real(dp)::half = 0.5
  real(dp)::third = 0.33333333333333333
  real(dp)::forth = 0.25
  real(dp)::sixth = 0.16666666666666667
end module const

