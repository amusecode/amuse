! Subroutines for nucleosynthesis post-procesing
! originally by Richard Stancliffe (Stancliffe & al. 2005 MNRAS).

module nucleosynthesis
   use real_kind
   use mesh
   
   implicit none
   integer, parameter :: nvar_nuc = 50
   integer, parameter :: nuc_num_rates = 92
   integer, save :: kh_nuc = 0

   !     Global switch to determine whether nucleosynthesis has been switched
   !     on or not
   logical, save :: nucleosynthesis_enabled = .false.

   !     Input quantities needed to calculate the reaction rates for the
   !     nucleosynthesis code
   !     FIXME: the star index ("2") should be the last one rather than the first one
   real(double), save :: ht(2,25,NM)

   !     Nucleosynthesis data structures. These are the analogues of the H and
   !     DH arrays.
   real(double), allocatable, save :: Hnucpr(:,:,:)
   real(double), allocatable, save :: Hnuc(:,:,:)
   real(double), allocatable, save :: DHnuc(:,:,:)
   real(double), allocatable, save :: nucacc(:,:)
   real(double), allocatable, save :: gradnuc(:,:,:)
   real(double), allocatable, save :: congradnuc(:,:,:)
   real(double), allocatable, save :: diffpnuc(:,:,:)
   real(double), allocatable, save :: difftnuc(:,:,:)
   real(double), allocatable, save :: diffcnuc(:,:,:,:)

   !     Nucleosynthesis reaction rates
   real(double), save :: nucsyn_crt(200,nuc_num_rates)
   real(double), save :: nucsyn_nrt(200,45)

   !     Misc. options
   !     Short-circuit decay of meta-stable Al26. Shouldn't matter but messes
   !     up convergence some times.
   logical, parameter :: instantly_decay_Al26m = .true.
end module nucleosynthesis



! Atomic data for nucleosynthesis stuff
module nucleosynthesis_atomic_data
   use real_kind
   ! Baryon number for different isotopes
   real(double), parameter :: baryn(50) = (/                                  &
         1.d0,  1.d0,  2.d0,  3.d0,  7.d0,  7.d0, 11.d0, 13.d0, 14.d0, 15.d0, &
        17.d0, 18.d0, 19.d0, 21.d0, 22.d0, 22.d0, 23.d0, 24.d0, 25.d0, 26.d0, &
        26.d0, 26.d0, 27.d0, 28.d0, 29.d0, 30.d0, 31.d0, 32.d0, 33.d0, 34.d0, &
        56.d0, 57.d0, 58.d0, 59.d0, 60.d0, 59.d0, 58.d0, 59.d0, 60.d0, 61.d0, &
         1.d0,  4.d0, 12.d0, 14.d0, 16.d0, 20.d0, 40.d0,  1.d0,  1.d0,  1.d0/)

   ! Mass in AMU
   real(double), parameter :: aamu(50) = (/                                                                 &
       0.00000,  0.00000,  2.00912,  3.01037,  7.01035,  7.01095, 11.00602, 13.00217, 14.00210, 15.00007,   &
      16.99944, 17.99946, 18.99897, 20.99602, 21.99443, 21.99640, 22.99338, 23.99032, 24.99084, 25.98874,   &
      25.99152, 25.99152, 26.98806, 27.98507, 28.98479, 29.98303, 30.98303, 31.98193, 32.98154, 33.97921,   &
      55.95791, 56.95821, 57.95684, 58.95787, 59.95735, 58.95678, 57.95817, 58.95753, 59.95523, 60.95540,   &
       1.00506,  4.00168, 12.00000, 14.00199, 15.99671, 19.99511, 39.97580,  0.00000,  0.00000,  0.00000/)

   ! Nuclear charge
   integer, parameter :: zz(50) = (/            &
        0,  0,  1,  2,  3,  4,  5,  6,  6,  7,  &
        8,  8,  9, 10, 10, 11, 11, 12, 12, 12,  &
       13, 13, 13, 14, 14, 14, 15, 16, 16, 16,  &
       26, 26, 26, 26, 26, 27, 28, 28, 28, 28,  &
        1,  2,  6,  7,  8, 10, 20,  0,  0,  0/)

   ! Maximum charge of any isotope
   integer, parameter :: Zmax = 28

   ! Statistical weight of the ground state for a nucleus with i electrons.
   ! In principle these should be selected using Hund's rules; values here
   ! taken from STARS, but some of them look funny: the neutral state of Ca, for
   ! instance, is 1s(2) 2s(2) 2p(6) 3s(2) 3p(6) 4s(2), which is a singlet
   ! state. Yet its multiplicity is listed as 21...
   ! This list goes up to Co
   real(double) :: com(Zmax+1) = (/                                  &
          1.0,  2.0,  1.0, 2.0,  2.0,  6.0,  9.0,  4.0, 9.0,  6.0,   &
          1.0,  2.0,  1.0, 6.0,  9.0,  4.0,  9.0,  6.0, 1.0, 10.0,   &
         21.0, 28.0, 25.0, 6.0, 25.0, 28.0, 21.0, 10.0, 1.0/)

   ! Ionisation energies for each element, in eV.
   ! Indexed first by ionisation state (1 = neutral), then by atom number
   ! Note that this is different from the equivalent data structure in the
   ! structure code, which indexes by ionisation state and variable number
   ! (1-9 for H, He, C, N, O, Ne, Mg, Si and Fe)
   real(double) :: chi(1:28, 1:28) = reshape( (/                                                                              &
      13.5984, 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                            &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      24.5874, 79.0052, 0., 0., 0., 0., 0., 0., 0., 0.,                                                                       &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      5.39172, 81.03192, 203.48592, 0., 0., 0., 0., 0., 0., 0.,                                                               &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      9.3227, 27.5339, 181.4309, 399.1499, 0., 0., 0., 0., 0., 0.,                                                            &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      8.29803, 33.45283, 71.38343, 330.75843, 670.98443, 0., 0., 0., 0., 0.,                                                  &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      11.2603, 35.6436, 83.5314, 148.0253, 540.1123, 1030.1053, 0., 0., 0., 0.,                                               &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      14.5341, 44.1354, 91.5846, 169.0581, 266.9483, 819.0203, 1486.0663, 0., 0., 0.,                                         &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      13.6181, 48.7354, 103.6709, 181.0844, 294.9834, 433.1034, 1172.3934, 2043.8034, 0., 0.,                                 &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      17.4228, 52.3936, 115.102, 202.2418, 316.4848, 473.6498, 658.8358, 1612.7468, 2715.8668, 0.,                            &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      21.5646, 62.5279, 125.9779, 223.0979, 349.3079, 507.2379, 714.5139, 953.6129, 2149.4429, 3511.6429,                     &
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      5.13908, 52.42548, 124.04548, 222.95548, 361.35548, 533.53548, 742.03548, 1006.28548, 1306.14948, 2771.26948,           &
      4419.96948, 0., 0., 0., 0., 0., 0., 0., 0., 0.,                                                                         &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      7.64624, 22.68154, 102.82524, 212.09124, 353.36124, 540.12124, 765.14124, 1031.10124, 1359.16124, 1726.66124,           &
      3488.47124, 5451.13124, 0., 0., 0., 0., 0., 0., 0., 0.,                                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      5.98577, 24.81437, 53.26197, 173.25397, 327.07897, 517.56897, 759.32897, 1043.98897, 1374.11897, 1772.86897,            &
      2214.86897, 4300.84897, 6604.98897, 0., 0., 0., 0., 0., 0., 0.,                                                         &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      8.15169, 24.49749, 57.99049, 103.13229, 269.89929, 475.16929, 721.66929, 1025.20929, 1376.32929, 1777.69929,            &
      2254.05929, 2777.47929, 5215.10929, 7888.28929, 0., 0., 0., 0., 0., 0.,                                                 &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      10.4867, 30.2561, 60.4588, 111.9027, 176.9278, 397.3488, 660.9188, 970.5188, 1342.6488, 1767.0488,                      &
      2246.5088, 2807.3088, 3419.0488, 6235.9588, 9305.7988, 0., 0., 0., 0., 0.,                                              &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      10.36, 33.6979, 68.4879, 115.7099, 188.3044, 276.3574, 557.3054, 886.0554, 1265.6054, 1713.1054,                        &
      2217.9054, 2782.3454, 3434.5454, 4141.5554, 7365.3354, 10859.5254, 0., 0., 0., 0.,                                      &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      12.9676, 36.7816, 76.3916, 129.8568, 197.6568, 294.6868, 408.8828, 757.1628, 1157.2228, 1612.8528,                      &
      2142.1328, 2734.1228, 3390.8328, 4140.5928, 4949.9928, 8608.5128, 12554.8128, 0., 0., 0.,                               &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      15.7596, 43.3893, 84.1293, 143.9393, 218.9593, 309.9683, 434.2913, 577.7513, 1000.2013, 1478.8913,                      &
      2017.8513, 2636.1113, 3322.2113, 4077.9513, 4932.7213, 5850.7513, 9971.6413, 14397.8713, 0., 0.,                        &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      4.34066, 35.97066, 81.77666, 142.68666, 225.34666, 324.74666, 442.30666, 597.18666, 773.00366, 1276.80366,              &
      1841.50366, 2470.90366, 3185.50366, 3972.10366, 4833.20366, 5801.20366, 6834.60366, 11445.40366, 16379.45366, 0.,       &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      6.11316, 17.98486, 68.89796, 136.16796, 220.66796, 329.44796, 456.64796, 603.88796, 792.42796, 1003.70296,              &
      1595.60296, 2252.80296, 2979.40296, 3797.00296, 4691.50296, 5665.50296, 6752.50296, 7910.30296, 13039.10296, 18508.96296, &
      0., 0., 0., 0., 0., 0., 0., 0.,                                                                                         &
      6.5615, 19.3612, 44.1179, 117.6073, 209.2573, 319.9373, 457.9373, 616.0373, 796.0673, 1021.2473,                        &
      1271.0453, 1958.4053, 2715.1053, 3545.9053, 4473.4053, 5482.4053, 6576.4053, 7789.4053, 9077.3753, 14752.1753,          &
      20785.8853, 0., 0., 0., 0., 0., 0., 0.,                                                                                 &
      6.8281, 20.4036, 47.8953, 91.1625, 190.4625, 309.9925, 450.7925, 621.1925, 813.2925, 1029.2125,                         &
      1294.2825, 1585.7825, 2373.6225, 3236.7225, 4178.6225, 5222.6225, 6353.6225, 7574.6225, 8920.6225, 10346.0225,          &
      16595.0225, 23220.8425, 0., 0., 0., 0., 0., 0.,                                                                         &
      6.7462, 21.4062, 50.7172, 97.4262, 162.7079, 290.8379, 441.4379, 614.8379, 820.6379, 1051.1379,                         &
      1306.8379, 1614.9379, 1951.2149, 2847.2149, 3823.2149, 4883.2149, 6051.2149, 7311.2149, 8666.2149, 10152.2149,          &
      11721.8149, 18573.1149, 25819.2349, 0., 0., 0., 0., 0.,                                                                 &
      6.7665, 23.2522, 54.2122, 103.3722, 172.8322, 263.4671, 423.6471, 608.3471, 817.6471, 1062.0471,                        &
      1332.8471, 1630.8471, 1985.6471, 2369.8151, 3380.4151, 4477.4151, 5662.4151, 6961.4151, 8357.4151, 9853.4151,           &
      11487.4151, 13208.8151, 20690.5151, 28585.3251, 0., 0., 0., 0.,                                                         &
      7.43402, 23.07402, 56.74202, 107.94202, 180.34202, 275.94202, 395.14502, 589.64502, 811.44502, 1059.74502,              &
      1345.74502, 1660.14502, 2003.74502, 2406.74502, 2841.90802, 3976.60802, 5200.60802, 6517.60802, 7954.60802, 9493.60802, &
      11137.60802, 12925.60802, 14805.50802, 22946.10802, 31518.04802, 0., 0., 0.,                                            &
      7.9024, 24.0902, 54.7422, 109.5422, 184.5422, 283.6422, 408.6222, 559.6822, 793.2822, 1055.3822,                        &
      1345.5822, 1676.3822, 2037.3822, 2429.5822, 2886.5822, 3375.8382, 4641.8382, 5999.8382, 7455.8382, 9037.8382,           &
      10726.8382, 12525.8382, 14475.8382, 16498.8382, 25326.8382, 34604.5282, 0., 0.,                                         &
      7.881, 24.964, 58.464, 109.764, 189.264, 291.264, 420.164, 577.964, 764.094, 1039.494,                                  &
      1344.494, 1680.494, 2059.494, 2470.494, 2914.494, 3426.454, 3973.034, 5370.234, 6874.834, 8477.834,                     &
      10212.834, 12058.834, 14020.834, 16139.834, 18358.834, 27902.934, 37915.034, 0.,                                        &
      7.6398, 25.8086, 60.9986, 115.8986, 191.9586, 299.9586, 432.9586, 594.9586, 787.9586, 1012.5586,                        &
      1333.5586, 1685.5586, 2069.5586, 2499.5586, 2963.5586, 3462.5586, 4033.6386, 4640.6986, 6181.6986, 7829.6986,           &
      9585.6986, 11479.6986, 13490.6986, 15621.6986, 17916.6986, 20315.8986, 30604.6986, 41380.0986 /), &
   (/ 28, 28 /) )

end module nucleosynthesis_atomic_data



module nucleosynthesis_neutron_ratios
   use real_kind
   use mesh, only: NM
   ! Branching ratios for neutron reactions
   real(double) :: frac(50,NM)
end module nucleosynthesis_neutron_ratios



module zams_nucleosynthesis_abundances
   !     Default ZAMS abundance mixture, as a fraction of Z.
   !     Based on Anders & Grevesse (1989)
   use real_kind
   
   implicit none
   real(double), save :: cxD     = 2.406e-03
   real(double), save :: cxHe3   = 1.468e-03
   real(double), save :: cxLi7   = 4.687e-07
   real(double), save :: cxBe7   = 0.000e+00
   real(double), save :: cxB11   = 2.368e-07
   real(double), save :: cxC12   = 1.764e-01
   real(double), save :: cxC13   = 1.829e-03
   real(double), save :: cxC14   = 0.000e+00
   real(double), save :: cxN14   = 5.212e-02
   real(double), save :: cxN15   = 2.186e-04
   real(double), save :: cxO16   = 5.031e-01
   real(double), save :: cxO18   = 1.086e-03
   real(double), save :: cxO17   = 1.948e-04
   real(double), save :: cxF19   = 2.030e-05
   real(double), save :: cxNe21  = 2.068e-04
   real(double), save :: cxNe20  = 9.221e-02
   real(double), save :: cxNe22  = 6.525e-03
   real(double), save :: cxNa22  = 0.000e+00
   real(double), save :: cxNa23  = 1.673e-03
   real(double), save :: cxMg24  = 2.580e-02
   real(double), save :: cxMg25  = 3.391e-03
   real(double), save :: cxMg26  = 3.889e-03
   real(double), save :: cxAl26m = 0.000e+00
   real(double), save :: cxAl27  = 2.906e-03
   real(double), save :: cxAl26g = 0.000e+00
   real(double), save :: cxSi28  = 3.272e-02
   real(double), save :: cxSi30  = 1.179e-03
   real(double), save :: cxSi29  = 1.717e-03
   real(double), save :: cxP31   = 4.087e-03
   real(double), save :: cxS33   = 1.615e-04
   real(double), save :: cxS32   = 1.984e-02
   real(double), save :: cxS34   = 9.351e-04
   real(double), save :: cxFe57  = 1.431e-03
   real(double), save :: cxFe60  = 0.000e+00
   real(double), save :: cxFe56  = 5.858e-02
   real(double), save :: cxFe58  = 1.853e-04
   real(double), save :: cxFe59  = 0.000e+00
   real(double), save :: cxCo59  = 1.683e-04
   real(double), save :: cxNi61  = 4.307e-05
   real(double), save :: cxNi59  = 0.000e+00
   real(double), save :: cxNi58  = 2.478e-03
   real(double), save :: cxNi60  = 9.812e-04
   real(double), save :: cxCa40  = 3.735e-03
end module zams_nucleosynthesis_abundances



subroutine allocate_nucleosynthesis_data(kh_in)
   use real_kind
   use nucleosynthesis
   use mesh
   
   implicit none
   integer, intent(in) :: kh_in

   if (.not. nucleosynthesis_enabled) return

   if (nvar_nuc > NVAR) then
      write(*, *) 'stellar structure/solver code does not hold enough data to solve nucleosynthesis'
      write(*, *) 'increase the value of NVAR in mesh.f'
      stop
   end if
   if (allocated(Hnucpr)) deallocate(Hnucpr)
   if (allocated(Hnuc)) deallocate(Hnuc)
   if (allocated(DHnuc)) deallocate(DHnuc)
   if (allocated(gradnuc)) deallocate(gradnuc)
   if (allocated(nucacc)) deallocate(nucacc)
   if (allocated(congradnuc)) deallocate(congradnuc)
   if (allocated(diffpnuc)) deallocate(diffpnuc)
   if (allocated(difftnuc)) deallocate(difftnuc)
   if (allocated(diffcnuc)) deallocate(diffcnuc)

   ! FIXME: these arrays have their indices in the wrong order: for better performance the star index ("2") should be the last one.
   allocate(Hnucpr(2, nvar_nuc, kh_in))
   allocate(Hnuc(2, nvar_nuc, kh_in))
   allocate(DHnuc(2, nvar_nuc, kh_in))
   allocate(gradnuc(2, nvar_nuc, kh_in))
   allocate(nucacc(nvar_nuc, 2))
   allocate(congradnuc(nvar_nuc, kh_in, 2))
   allocate(diffpnuc(nvar_nuc,kh_in,2))
   allocate(difftnuc(nvar_nuc,kh_in,2))
   allocate(diffcnuc(nvar_nuc,nvar_nuc,kh_in,2))

   kh_nuc = kh_in
end subroutine allocate_nucleosynthesis_data



subroutine set_initl_nucleosynt_abundances(Jstar)
   use real_kind
   use zams_nucleosynthesis_abundances
   use nucleosynthesis
   use constants
   use settings
   use atomic_data
   use extra_elements
   
   implicit none
   integer, intent(in) :: Jstar
   real(double) :: xa(9), ya(9), avm
   integer :: jk, i

   if (.not. nucleosynthesis_enabled) return

   ! Beginning of minor variables, in mass order (except gallinoes)
   Hnuc(Jstar, 1 , 1:kh) = 0d0       ! gallinoes
   Hnuc(Jstar, 2 , 1:kh) = 0d0       ! neutrons
   Hnuc(Jstar, 3 , 1:kh) = cxD * czs
   Hnuc(Jstar, 4 , 1:kh) = cxHe3 * czs
   Hnuc(Jstar, 5 , 1:kh) = cxLi7 * czs
   Hnuc(Jstar, 6 , 1:kh) = cxBe7 * czs
   Hnuc(Jstar, 7 , 1:kh) = cxB11 * czs
   Hnuc(Jstar, 8 , 1:kh) = cxC13 * czs
   Hnuc(Jstar, 9 , 1:kh) = cxC14 * czs
   Hnuc(Jstar, 10, 1:kh) = cxN15 * czs
   Hnuc(Jstar, 11, 1:kh) = cxO17 * czs
   Hnuc(Jstar, 12, 1:kh) = cxO18 * czs
   Hnuc(Jstar, 13, 1:kh) = cxF19 * czs
   Hnuc(Jstar, 14, 1:kh) = cxNe21 * czs
   Hnuc(Jstar, 15, 1:kh) = cxNe22 * czs
   Hnuc(Jstar, 16, 1:kh) = cxNa22 * czs
   Hnuc(Jstar, 17, 1:kh) = cxNa23 * czs
   Hnuc(Jstar, 18, 1:kh) = cxMg24 * czs
   Hnuc(Jstar, 19, 1:kh) = cxMg25 * czs
   Hnuc(Jstar, 20, 1:kh) = cxMg26 * czs
   Hnuc(Jstar, 21, 1:kh) = cxAl26m * czs
   Hnuc(Jstar, 22, 1:kh) = cxAl26g * czs
   Hnuc(Jstar, 23, 1:kh) = cxAl27 * czs
   Hnuc(Jstar, 24, 1:kh) = cxSi28 * czs
   Hnuc(Jstar, 25, 1:kh) = cxSi29 * czs
   Hnuc(Jstar, 26, 1:kh) = cxSi30 * czs
   Hnuc(Jstar, 27, 1:kh) = cxP31 * czs
   Hnuc(Jstar, 28, 1:kh) = cxS32 * czs
   Hnuc(Jstar, 29, 1:kh) = cxS33 * czs
   Hnuc(Jstar, 30, 1:kh) = cxS34 * czs
   Hnuc(Jstar, 31, 1:kh) = cxFe56 * czs
   Hnuc(Jstar, 32, 1:kh) = cxFe57 * czs
   Hnuc(Jstar, 33, 1:kh) = cxFe58 * czs
   Hnuc(Jstar, 34, 1:kh) = cxFe59 * czs
   Hnuc(Jstar, 35, 1:kh) = cxFe60 * czs
   Hnuc(Jstar, 36, 1:kh) = cxCo59 * czs
   Hnuc(Jstar, 37, 1:kh) = cxNi58 * czs
   Hnuc(Jstar, 38, 1:kh) = cxNi59 * czs
   Hnuc(Jstar, 39, 1:kh) = cxNi60 * czs
   Hnuc(Jstar, 40, 1:kh) = cxNi61 * czs
   Hnuc(Jstar, 47, 1:kh) = cxCa40 * czs

   ! Remesh using inputs - not necessarily consistent with above
   do jk = 1, kh
      xa(1) = H(24*(Jstar-1)+5, jk)
      xa(2) = H(24*(Jstar-1)+9, jk)
      xa(3) = H(24*(Jstar-1)+10, jk)
      xa(4) = H(24*(Jstar-1)+16, jk)
      xa(5) = H(24*(Jstar-1)+3, jk)
      xa(6) = H(24*(Jstar-1)+11, jk)
      xa(7) = H(NXFSTAR*(Jstar-1)+nMg24, jk)
      xa(8) = H(NXFSTAR*(Jstar-1)+nSi28, jk)
      xa(9) = H(NXFSTAR*(Jstar-1)+nFe56, jk)
      ! Mass fraction abundances of Mg, Fe and Si were set globally by
      ! remesh, we don't touch these.
      ya(1:9) = xa(1:9)*can(1:9)/cbn(1:9)
      avm = sum(ya(1:9))
      do i=1,6
         Hnuc(Jstar, 40+i, jk) = ya(i)/avm
      end do
   end do
   ! For high core temperatures, log10 T > 6.7, set D, Li7 to 0 (should all
   ! have burned on the pre-main sequence)
   if (H(24*(Jstar-1)+2, kh) > 6.5*cln) then
      Hnuc(Jstar, 3 , 1:kh) = 0.0d0
      Hnuc(Jstar, 4 , 1:kh) = 0.0d0
      Hnuc(Jstar, 5 , 1:kh) = 0.0d0
      Hnuc(Jstar, 6 , 1:kh) = 0.0d0
      Hnuc(Jstar, 7 , 1:kh) = 0.0d0
   end if

   ! Default nucacc: current surface abundance
   nucacc(:, Jstar) = Hnuc(Jstar, :, 1)
end subroutine set_initl_nucleosynt_abundances



subroutine load_nucleosynthesis_rates ( rates_in, nrates_in )
   use real_kind
   use nucleosynthesis
   use reaction_rate_data
   use nucleosynthesis_atomic_data
   
   implicit none
   integer, intent(in) :: rates_in, nrates_in
990 format (1x, 10f7.3)

   read (rates_in, 990) nucsyn_crt
   read (nrates_in, 990) nucsyn_nrt
   rewind (rates_in)
   rewind (nrates_in)
   !     Make sure the first 20 reaction rates (structure and nucleosynthesis)
   !     agree
   nucsyn_crt(:,1:20) = crt(:, 1:20)
end subroutine load_nucleosynthesis_rates




! Prepare the nucleosynthesis surface accretion abundances
! For the structure code, we can do this in funcs1 for both stars simultaneously in TWIN mode, but for the nucleosynthesis code
! the stars are solved one after the other and we should call this function once before we begin the iterations
subroutine update_accretion_abundances2
   use real_kind
   use mesh
   use settings
   use nucleosynthesis
   implicit none
   integer :: Jstar

   if (.not. allocated(nucacc)) return

   do Jstar = 1, ktw
      nucacc(:, 3-Jstar) = Hnuc(Jstar, :, 1)
   end do
end subroutine update_accretion_abundances2




subroutine update_concentration_gradients()
   use real_kind
   use mesh
   use settings
   use nucleosynthesis
   use nucleosynthesis_atomic_data
   use constants
   implicit none
   integer :: Jstar, k
   real(double) :: C(50, kh)

   if (.not. allocated(nucacc)) return

   do Jstar = 1, ktw
      ! Find concentrations
      forall (k = 1:kh)
         C(:, k) = Hnuc(Jstar, :, k) * ht(Jstar, 1, k) / ( max(aamu(:), 1.0_dbl) * AMU * ht(Jstar,9,k) )
      end forall
      where (C > 0.0d0) C = log(C)
      ! Concentration gradients
      forall (k = 2:kh-1)
         congradnuc(:, k, Jstar) = ( C(:, k) - C(:, k-1) ) / ( ht(Jstar, 24, k) - ht(Jstar, 24, k+1) )
      end forall
      congradnuc(:, 1, Jstar) = 0.0d0
      congradnuc(:, kh, Jstar) = 0.0d0
   end do
end subroutine update_concentration_gradients




subroutine update_radative_accelerations2
   use real_kind
   use mesh
   use settings
   use constants
   use nucleosynthesis
   use nucleosynthesis_atomic_data
   use radiative_acceleration
   use ionisation
   implicit none
   integer :: Jstar, jk
   real(double) :: logT, logrho, T, rho, Frad, mub, dv
   real(double) :: na(50), nn(50)
   real(double) :: naelem(50), nnelem(50)
   real(double) :: grad(50), lk
   integer :: ze(0:50), ep(0:50), eflag(0:50), n, nae
   real(double) :: avm, avgz(50), nne
   
   integer :: ns, ne, lastz
   real(double) :: lasta
   real(double) :: xai(50*53/2)
   integer :: kzi(50*53/2), iz(50*53/2)
   real(double) :: aa(50*53/2), dz(50*53/2)

   real(double), allocatable :: diffp_iso(:), difft_iso(:), diffc_iso(:, :)

   logical, parameter :: full_ionisation = .true.
   logical, parameter :: average_diffusion_coefficients = .false.

   if (.not. allocated(gradnuc)) return
   gradnuc = 0.0
   diffpnuc = 0.0
   difftnuc = 0.0
   diffcnuc = 0.0

   ! Allocate memory for isotope-specific diffusion coefficients
   ! FIXME: this number should not be hard-coded here; we should rather use the result from the code that calculates different
   ! ionisation stages.
   if (.not. average_diffusion_coefficients) then
      allocate(diffp_iso(653))
      allocate(difft_iso(653))
      allocate(diffc_iso(653, 653))
   end if

   do Jstar = 1, ktw
      do jk = 1, kh
         ! Input parameters
         avm = ht(Jstar,8,jk)
         Frad = ht(Jstar, 23, jk)
         rho = ht(Jstar, 1, jk)
         logrho = log10(rho)
         logT = H(2 + 24*(Jstar-1), jk)/cln
         T = 10**logT
         dv = ht(Jstar, 25, jk)

         ! Compute elemental abundances
         na(1:50) = Hnuc(Jstar, 1:50, jk) / baryn(1:50)

         ! Convert particle number fractions to number densities
         nn(1:50) = na(1:50) * rho/(avm*amu)
         nn(1:2) = 0.0

         ! Find (average) charge for each isotope and the electron density
         if (full_ionisation .and. .not. average_diffusion_coefficients) then
            ! We don't have ionisation information for the isotopes from the nucleosynthesis code.
            ! Compute gravitational settling coefficients pretending that everything is fully ionised.
            ! It's wrong, but the best we can do without more detailed nuclear info
            avgz(1:50) = dble(zz(1:50))
            avgz(1:2) = 0.0
            nne = dot_product(avgz(3:50), na(3:50))
            nne = nne * rho/(avm*amu)
         else
            ! Calculate abundance fractions for different ionisation stages
            call find_ionisation_levels(50, Zmax, T, dv, na,   com, zz, aamu, chi,   ns, xai, kzi, aa, iz) 

            ! Calculate actual average charge Z for each isotope
            avgz = 0.0_dbl
            lastz = 0
            lasta = 0.0
            ne = 1
            do n = 1, ns
               if (lasta /= aa(n) .or. lastz /= kzi(n)) then
                  ne = 1
                  do while (.not. (aamu(ne) == aa(n) .and. zz(ne) == kzi(n)) )
                     ne = ne + 1
                  end do
                  lastz = zz(ne)
                  lasta = aamu(ne)
               end if
               if (na(ne) > 0) avgz(ne) = avgz(ne) + iz(n) * xai(n) / na(ne)
            end do

            ! Calculate electron density
            nne = dot_product(avgz(3:50), na(3:50))
            nne = nne * rho/(avm*amu)

            !lastz = 0
            !lasta = 0.0
            !ne = 1
            !do n = 1, ns
            !   if (lasta /= aa(n) .or. lastz /= kzi(n)) then
            !      ne = 1
            !      do while (.not. (aamu(ne) == aa(n) .and. zz(ne) == kzi(n)) )
            !         ne = ne + 1
            !      end do
            !      lastz = zz(ne)
            !      lasta = aamu(ne)
            !   end if
            !   print '(i3, i3, " Z",i3, " A",i3, 2g, 2e, i3, g)', &
            !      n, ne, kzi(n), int(baryn(ne)), aa(n), aamu(ne), xai(n), nn(ne), iz(n), avgz(ne)
            !end do
         end if

         if (.not. average_diffusion_coefficients) then
            ! Find diffusion coefficients, using Burgers' formalism, assuming all ions are "average"
            call diffusion_terms_burgers(rho,T, 50,nn,nne, avgz,aamu, &
                                         diffpnuc(:,jk,Jstar),difftnuc(:,jk,Jstar),diffcnuc(:,:,jk,Jstar))
         else
            ! Calculate average diffusion coefficients for each element, the hard way

            ! Set small abundances to 0, to speed up the calculation
            do n=1, ns
               if (xai(n) < 1.0e-50_dbl) xai(n) = 0.0_dbl
            end do

            ! Convert number fractions to mass densities
            xai(1:ns) = xai(1:ns) * rho/(avm*amu)
            dz = real(iz)

            ! Calculate diffusion coefficients separately for each ionisation stage
            call diffusion_terms_burgers(rho,T, ns,xai(1:ns),nne, dz(1:ns),aa(1:ns), diffp_iso,difft_iso,diffc_iso)

            ! Calculate average for each isotope
            diffpnuc(:,jk,Jstar) = 0
            difftnuc(:,jk,Jstar) = 0
            diffcnuc(:,:,jk,Jstar) = 0
            lastz = 0
            lasta = 0.0
            ne = 1
            do n=1, ns
               if (lasta /= aa(n) .or. lastz /= kzi(n)) then
                  ne = 1
                  do while (.not. (aamu(ne) == aa(n) .and. zz(ne) == kzi(n)) )
                     ne = ne + 1
                  end do
                  lastz = zz(ne)
                  lasta = aamu(ne)
               end if
               if (na(ne) > 0) then
                  diffpnuc(ne,jk,Jstar)   = diffpnuc(ne,jk,Jstar)   + diffp_iso(n)   * xai(n)/nn(ne)
                  difftnuc(ne,jk,Jstar)   = difftnuc(ne,jk,Jstar)   + difft_iso(n)   * xai(n)/nn(ne)
                  diffcnuc(ne,:,jk,Jstar) = diffcnuc(ne,:,jk,Jstar) + diffc_iso(n,:) * xai(n)/nn(ne)
               end if
            end do
         end if

         !print '(i4, 48e)', jk, logT, h(4, jk)/CMSN, diffpnuc(1:46, jk, Jstar)


         ! Don't bother with accelerations if this zone is convective
         if (ht(Jstar,2, jk) /= 0.0d0) cycle

         ! Nothing more to do if we're not calculating accelerations
         if (crlev == 0.0d0) cycle

         ! Abundance fraction by number
         nn(1:50) = na(1:50) / sum(na(1:50))
         mub = dot_product(nn, aamu)

         ! Mark all elements (unique Z values)
         ! The number of elements is at most the same as the number of isotopes
         ep = 0
         ze = 0
         eflag = 0
         forall (n=1:50) eflag(zz(n)) = 1
         nae = 0
         ! Construct permutation list of all elements
         ! Map Zi -> j = 1..#elements
         do n=1, 50
            if (eflag(n) /= 0) then    ! Element
               nae = nae + 1
               eflag(n) = nae
               ze(nae) = n
            end if
         end do
         ! Elemental abundances
         naelem = 0.0d0
         nnelem = 0.0d0
         do n=1, 50
            ep(n) = eflag(zz(n))
            if (ep(n) > 0) then
               naelem(ep(n)) = naelem(ep(n)) + na(n)
               nnelem(ep(n)) = nnelem(ep(n)) + nn(n)
            end if
         end do

         ! Calculate accelerations
         call get_radiative_accelerations(logT, logrho, mub, Frad, nae, nnelem(1:nae), ze(1:nae), grad(1:nae), lk)

         ! Accelerations for each isotope
         forall (n=1:50) gradnuc(Jstar, n, jk) = grad(ep(n)) / max(1.0d0, aamu(n))
      end do
   end do
   if (.not. average_diffusion_coefficients) then
      deallocate(diffp_iso)
      deallocate(difft_iso)
      deallocate(diffc_iso)
   end if
end subroutine update_radative_accelerations2



!> \todo Input "variable" k1=1 always
subroutine funcs2 ( Jstar, k, k1, var, varnuc, dvarnuc, fn1, ylin )
   use real_kind
   use mesh
   use constants
   use settings
   use control
   use nucleosynthesis
   use test_variables
   use current_model_properties
   use nucleosynthesis_neutron_ratios
   use nucleosynthesis_atomic_data
   
   implicit none
   integer, intent(in) :: k,k1,Jstar
   real(double), intent(in) :: var(NVAR)
   real(double), intent(in) :: varnuc(50), dvarnuc(50)
   real(double), intent(out) :: fn1(NFUNC), ylin(50, NFUNC)

   integer :: i
   real(double) :: at,dm,dmk,ratetot
   real(double) :: rho, zt
   real(double) :: ydot(50), y(50), sg, sgth, dmt, mu, apr, atr

   real(double) :: xa2(50), na2(50)

   ! Rates from nucrat2 - format is Rabc
   ! a=light reactant,b=heavy reactant,c=light product, blank if photon
   real(double) :: rpp, r33, r34, rBe, rBp, rpc, rpN, rpO, r3a, raC             !10
   real(double) :: raN, raO, raNe, rCC, rCO, rOO, rgNe, rgMg, rCCg, rpng        !20
   real(double) :: rpd,rpbe,rpli,rebe,rpC13,rpN15a,rpO17a                       !27
   real(double) :: raC13n, raN15g, raO17n, rlia, rpB11, rpC14, raC14            !34
   real(double) :: rpO18,rpO18a,raO18,raO18n,rpF19,rpF19a,raF19p,rpNe21         !42
   real(double) :: raNe21,raNe21n,rpNe22,raNe22,raNe22n,RnNa22p,RnNa22a,rpNa22  !50
   real(double) :: rpNa23,rpNa23n,rpNa23a,rpMg24,raMg24,rpMg25t,rpMg25m,rpMg25g !58
   real(double) :: raMg25,raMg25n,raMg25p,rpMg26,rpMg26tn,rpMg26mn,rpMg26gn     !65
   real(double) :: raMg26,raMg26n,raMg26p,rg26tp,Rn26tp,Rn26ta,rp26t,rg26mp     !73
   real(double) :: Rn26mp,Rn26ma,rp26m,rg26gp,Rn26gp,Rn26ga,rp26g,rpAl27        !81
   real(double) :: rpAl27a,raAl27n,raNa23nt,raNa23nm,raNa23ng,rpSi28,rpSi29     !88
   real(double) :: rpSi30,rpN15,rpO17,rpNe20                                    !92
   real(double) :: rrt2(92)

   real(double) :: RnFe56,RnFe57,RnFe58,RnCo59,RnNi58,RnNi59,RnNi60
   real(double) :: rnp,RnHe3,RnLi7,RnC12,RnC13,RnC14,RnN14,RnN15,RnO16,RnO18
   real(double) :: RnF19,RnNe20,RnNe21,RnNe22,RnNa23,RnMg24,RnMg25,RnMg26
   real(double) :: RnAl27,RnSi28,RnSi29,RnSi30,RnP31,RnS32,RnS33,RnS34,Rn26g
   real(double) :: RnS33a,RnN14p,RnNi59p,RnNi59a,RnO17a,Rn26gad,RnS32a,RnFe59
   real(double) :: RnFe60,RnS34s,RnNi61s
   real(double) :: nrate(45)


   real(double) :: rdAl26g,rdNa22,rd26m,rdFe59,rdFe60,rdNi59,rdn,rdC14
   real(double) :: drt(8)

   real(double) :: xa1(9), na(9), neo, nio, nzz, avm, ne, t, r2
   real(double) :: wpp, wpi(50), diffp(50),difft(50),diffc(50,50)
   !real(double) :: avgz(50), nn(50), nne, 
   integer :: j

   ! Set up composition for reaction calculation
   xa1(1) = var(5+24*(Jstar-1))  !X1
   xa1(2) = var(9+24*(Jstar-1))  !X4
   xa1(3) = var(10+24*(Jstar-1)) !X12
   xa1(4) = var(16+24*(Jstar-1)) !X14
   xa1(5) = var(3+24*(Jstar-1))  !X16
   xa1(6) = var(11+24*(Jstar-1)) !X20
   xa1(7) = cmg*czs
   xa1(8) = csi*czs
   xa1(9) = cfe*czs
   
   ! Temperature
   at = var(2+24*(Jstar-1))
   T = exp(at)
   r2 = exp(2.0d0*var(7+24*(Jstar-1))) - ct(8)
   
   ! Number abundance stuff from FUNCS1
   neo = ht(Jstar,5,k)
   nio = ht(Jstar,6,k)
   nzz = ht(Jstar,7,k)
   avm = ht(Jstar,8,k)
   ne = ht(Jstar,9,k)
   na(1:9) = ht(Jstar, 10:18, k)
   
   ! Blank any -ve abundances. There shouldn't be any anyway...
   do i = 1, 9
      if (xa1(i) < 0d0) then
         xa1(i) = 0d0
         na(i) = 0d0
      end if
   end do
   xa2(:) = varnuc(:)
   xa2(2) = 0d0 ! neutrons are annoying - we treat them differently
   
   ! Proton fudging: if structure says there are none, there are none
   if (xa1(1) < 1d-40) xa2(41) = 0d0
   
   ! Work out abundances for nucrat2 - RJS
   na2(1:50) = xa2(1:50)/baryn(1:50)
   if ( k <= k1 ) then
      ylin(1:50, 1:NFUNC) = 0.0d0
   end if
   rho = ht(Jstar, 1, k)
   sg = ht(Jstar, 2, k)
   sgth = ht(Jstar, 22, k)
   zt = ht(Jstar, 3, k)
   dm = var(4+24*(Jstar-1))
   dmt = ht(Jstar, 19, k)
   dmk = ht(Jstar, 4, k)
   apr = ht(Jstar, 20, k)
   atr = ht(Jstar, 21, k)

   ! Mean molecular weight, for thermohaline mixing
   ! Use values from previous iteration, as in the structure code
   !mu = 1.0 / dot_product(xa2(3:46)/aamu(3:46), 1 + zz(3:46))
   mu = 1.0 / dot_product(Hnuc(Jstar, 3:46, k)/aamu(3:46), 1 + zz(3:46))

   ! Advection terms: gravitational settling, radiative levitation
   if (cgrs > 0.0d0 .and. grs_burgers .and. jmod > 2) then
      ! We don't have ionisation information for the isotopes from the nucleosynthesis code.
      ! Compute gravitational settling coefficients pretending that everything is fully ionised.
      ! It's wrong, but the best we can do without more detailed nuclear info
      !avgz(1:50) = dble(zz(1:50))
      !avgz(1:2) = 0.0
      !nne = dot_product(avgz(3:50), na2(3:50))
      !! Convert particle number fractions to number densities
      !nn(1:50) = na2(1:50) * rho/(avm*amu)
      !nn(1:2) = 0.0
      !nne = nne * rho/(avm*amu)
      !call diffusion_terms_burgers(rho,T, 50,nn,nne, avgz,aamu, diffp,difft,diffc)

      ! The diffusion coefficients were pre-calculated
      diffp(:) = diffpnuc(:,k,Jstar)
      difft(:) = difftnuc(:,k,Jstar)
      diffc(:,:) = diffcnuc(:,:,k,Jstar)
      wpi(3:50) = 1.0d-11 * (diffp(3:50)*apr + difft(3:50)*atr)
      forall (j=3:50)
         wpi(j) = wpi(j) + 1.0d-22 * dot_product(diffc(j, 3:50), congradnuc(3:50, k, Jstar))
      end forall
      wpi(1:2) = 0.0d0

      ! Accelerations
      if (crlev>0.0d0) wpi(3:50) = wpi(3:50) + crlev*1.0d-11 * diffp(3:50) * gradnuc(Jstar,3:50,k) * aamu(3:50) / (CR*T)

      ! No net mass flux
      wpp = -dot_product(wpi(3:50), xa2(3:50))
      wpi(3:50) = wpi(3:50) + wpp

      ! NB: note difference with respect to the structure code: a factor XA2(1:50) was omitted.
      ! This is included in equns2. The reason for this change is that equns2 also needs to calculate the
      ! derivative. If we approximate the flux as FF(i) = XA2(i) * const(i), which should be ok to lowest order,
      ! then this derivative is just FF(i)/XA2(i).
      forall (j=1:50) fn1(110+j) = cgrs * rho*cpi4*r2 * wpi(j)
   else
      fn1(111:160) = 0.0d0
   end if

   
   call nucrat2 ( at, rho/avm, zt, na, neo, nio, na2,  rrt2, nrate, drt )
   
   ! Copy results to local variables
   !> \todo FIXME: should perhaps make a big table of what rates to include where instead...
   ! Normal nuclear reaction rates:
   rpp     = rrt2(1);  r33     = rrt2(2);  r34     = rrt2(3);  rBe     = rrt2(4)
   rBp     = rrt2(5);  rpc     = rrt2(6);  rpN     = rrt2(7);  rpO     = rrt2(8)
   r3a     = rrt2(9);  raC     = rrt2(10); raN     = rrt2(11); raO     = rrt2(12)
   raNe    = rrt2(13); rCC     = rrt2(14); rCO     = rrt2(15); rOO     = rrt2(16)
   rgNe    = rrt2(17); rgMg    = rrt2(18); rCCg    = rrt2(19); rpng    = rrt2(20)
   rpd     = rrt2(21); rpbe    = rrt2(22); rpli    = rrt2(23); rebe    = rrt2(24)
   rpC13   = rrt2(25); rpN15a  = rrt2(26); rpO17a  = rrt2(27); raC13n  = rrt2(28)
   raN15g  = rrt2(29); raO17n  = rrt2(30); rlia    = rrt2(31); rpB11   = rrt2(32)
   rpC14   = rrt2(33); raC14   = rrt2(34); rpO18   = rrt2(35); rpO18a  = rrt2(36)
   raO18   = rrt2(37); raO18n  = rrt2(38); rpF19   = rrt2(39); rpF19a  = rrt2(40)
   raF19p  = rrt2(41); rpNe21  = rrt2(42); raNe21  = rrt2(43); raNe21n = rrt2(44)
   rpNe22  = rrt2(45); raNe22  = rrt2(46); raNe22n = rrt2(47); RnNa22p = rrt2(48)
   RnNa22a = rrt2(49); rpNa22  = rrt2(50); rpNa23  = rrt2(51); rpNa23n = rrt2(52)
   rpNa23a = rrt2(53); rpMg24  = rrt2(54); raMg24  = rrt2(55); rpMg25t = rrt2(56)
   rpMg25m = rrt2(57); rpMg25g = rrt2(58); raMg25  = rrt2(59); raMg25n = rrt2(60)
   raMg25p = rrt2(61); rpMg26  = rrt2(62); rpMg26tn= rrt2(63); rpMg26mn= rrt2(64)
   rpMg26gn= rrt2(65); raMg26  = rrt2(66); raMg26n = rrt2(67); raMg26p = rrt2(68)
   rg26tp  = rrt2(69); Rn26tp  = rrt2(70); Rn26ta  = rrt2(71); rp26t   = rrt2(72)
   rg26mp  = rrt2(73); Rn26mp  = rrt2(74); Rn26ma  = rrt2(75); rp26m   = rrt2(76)
   rg26gp  = rrt2(77); Rn26gp  = rrt2(78); Rn26ga  = rrt2(79); rp26g   = rrt2(80)
   rpAl27  = rrt2(81); rpAl27a = rrt2(82); raAl27n = rrt2(83); raNa23nt= rrt2(84)
   raNa23nm= rrt2(85); raNa23ng= rrt2(86); rpSi28  = rrt2(87); rpSi29  = rrt2(88)
   rpSi30  = rrt2(89); rpN15   = rrt2(90); rpO17   = rrt2(91); rpNe20  = rrt2(92)

   ! Neutron rates
   RnFe56   = nrate(1);  RnFe57   = nrate(2);  RnFe58   = nrate(3);  RnCo59   = nrate(4);
   RnNi58   = nrate(5);  RnNi59   = nrate(6);  RnNi60   = nrate(7);  rnp      = nrate(8);
   RnHe3    = nrate(9);  RnLi7    = nrate(10); RnC12    = nrate(11); RnC13    = nrate(12);
   RnC14    = nrate(13); RnN14    = nrate(14); RnN15    = nrate(15); RnO16    = nrate(16);
   RnO18    = nrate(17); RnF19    = nrate(18); RnNe20   = nrate(19); RnNe21   = nrate(20);
   RnNe22   = nrate(21); RnNa23   = nrate(22); RnMg24   = nrate(23); RnMg25   = nrate(24);
   RnMg26   = nrate(25); RnAl27   = nrate(26); RnSi28   = nrate(27); RnSi29   = nrate(28);
   RnSi30   = nrate(29); RnP31    = nrate(30); RnS32    = nrate(31); RnS33    = nrate(32);
   RnS34    = nrate(33); Rn26g    = nrate(34); RnS33a   = nrate(35); RnN14p   = nrate(36);
   RnNi59p  = nrate(37); RnNi59a  = nrate(38); RnO17a   = nrate(39); Rn26gad  = nrate(40);
   RnS32a   = nrate(41); RnFe59   = nrate(42); RnFe60   = nrate(43); RnS34s   = nrate(44);
   RnNi61s  = nrate(45);

   ! Decay rates
   rdAl26g = drt(1); rdNa22  = drt(2); rd26m   = drt(3); rdFe59  = drt(4);
   rdFe60  = drt(5); rdNi59  = drt(6); rdn     = drt(7); rdC14   = drt(8);

   !  EVALUATE COMPOSITION CHANGE
   ! 1/9/03 Minor element reactions from Clayton 1963
   ydot(1) = 0d0
   ydot(2) =  rdn - raC13n - raO17n - raNe21n - raO18n  - raNe22n      &
        + RnNa22p + RnNa22a - rpNa23n - raMg25n - rpMg26mn - rpMg26gn  &
        - raMg26n + Rn26mp + Rn26ma + Rn26g                            &
        + Rn26gp + Rn26ga - raAl27n - raNa23nm - raNa23ng              &
        + RnFe56 + RnFe57 + RnFe58 + RnCo59 + RnNi58 + RnNi59          &
        + RnNi60 + rnp + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14         &
        + RnN14 + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21      &
        + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + RnAl27          &
        + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 + RnS33 + RnS34     &
        + RnS33a + RnN14p + RnNi59p + RnNi59a + RnO17a + RnS32a
   ydot(3)  = rpd - rpp - rnp
   ydot(4)  = - rpd + 2.0*r33 + r34 + RnHe3
   ydot(5)  =  rpli - rebe + rlia + RnLi7
   ydot(6)  = -r34 + rpbe + rebe
   ydot(7)  = rpB11 - rlia
   ydot(8)  = rpC13 - rpc + raC13n + RnC13 - RnC12
   ydot(9)  = rpC14 + raC14 + RnC14 - RnN14p - RnO17a - RnC13 + rdC14
   ydot(10) = rpN15a - rpN + raN15g + RnN15 -  RnC14 - RnN14 - rpO18a + rpN15 - rpC14
   ydot(11) = - rpO + rpO17a + raO17n - RnO16 + RnO17a + rpO17
   ydot(12) = - raN - raC14 + rpO18 + rpO18a + raO18 + raO18n + RnO18 - rpO17
   ydot(13) = - raN15g + rpF19 + rpF19a + raF19p - RnNa22a - RnO18 + RnF19 - rpO18
   ydot(14) = - raO18n + rpNe21 + raNe21 + raNe21n + RnNe21 - RnNe20 - rpNe20
   ydot(15) = - raO18 - raF19p + rpNe22 + raNe22 + raNe22n - RnNa22p - rdNa22 - RnNe21 + RnNe22
   ydot(16) = RnNa22p + RnNa22a + rpNa22 - rpNe21 + rdNa22
   ydot(17) = rpNa23 + rpNa23n + rpNa23a - rpNe22 - Rn26ga - Rn26ma + raNa23nm + raNa23ng - RnNe22 + RnNa23
   ydot(18) = rpMg24 + raMg24 - rCC - raNe - raNe21n - rpNa23 - rpAl27a - RnNa23 + RnMg24
   ydot(19) = rpMg25m + rpMg25g + raMg25 + raMg25n + raMg25p - raNe21 - raNe22n - rg26mp - rg26gp + RnMg25 - RnMg24
   ydot(20) = rpMg26 + rpMg26mn + rpMg26gn + raMg26 + raMg26n + raMg26p - raNe22 - Rn26mp &
        - Rn26gp - rdAl26g - rd26m + RnMg26 - RnMg25
   ydot(21) = rg26mp + Rn26mp + Rn26ma + rp26m - rpMg25m - rpMg26mn - raNa23nm + rd26m
   ydot(22) = rg26gp + Rn26gp + Rn26ga + rp26g - rpMg25g - rpMg26gn - raNa23ng + rdAl26g + Rn26g
   ydot(23) = rpAl27 + rpAl27a + raAl27n - rpMg26 - RnMg26 + RnAl27 - Rn26g
   
   ! assume decay of Si28(p,g)P29 -> Si29
   ydot(24) = rpSi28 - rCO - raMg24 - raMg25n - rpAl27 - RnAl27 + RnSi28
   
   ! assume decay of Si29(p,g)P30 -> Si30, Mg26(a,p)Al29 -> Si29
   ydot(25) = rpSi29 - raMg25 - raMg26n - rpSi28 - RnSi28 + RnSi29 - RnS32a - raMg26p
   ydot(26) = rpSi30 - raMg26 - rpSi29 - RnSi29 + RnSi30 - RnS33a
   ydot(27) = - rpSi30 - RnSi30 + RnP31
   ydot(28) = - RnP31 + RnS32 + RnS32a
   ydot(29) = - RnS32 + RnS33 + RnS33a
   ydot(30) = - RnS33 + RnS34 + RnS34s
   ydot(31) = RnFe56 - RnNi59a
   ydot(32) = RnFe57 - RnFe56
   ydot(33) = RnFe58 - RnFe57
   ydot(34) = - RnFe58 + RnFe59 + rdFe59
   ydot(35) = - RnFe59 + RnFe60 + rdFe60
   ydot(36) = RnCo59 - rdFe59 - RnNi59p - rdNi59
   ydot(37) = RnNi58
   ydot(38) = RnNi59 - RnNi58 + RnNi59a + RnNi59p + rdNi59
   ydot(39) = RnNi60 - RnNi59
   ydot(40) = RnNi61s - RnNi60
   
   ! Protons duplicate
   ydot(41) = 2.0*rpp + rpd + rnp + rpli + rpbe + rpB11 - 2.0*r33 &
        + rpc + rpC13 + rpC14 + rpN + rpN15a + rpO18a + rpO &
        + rpO17a - rdn + rpO18 + rpF19 + rpF19a - raF19p + rpNe21 &
        + rpNe22 - RnNa22p + rpNa22 + rpNa23 + rpNa23n + rpNa23a &
        + rpMg24 + rpAl27a + rpMg25m + rpMg25g + rpMg26 + rpN15 &
        + rpMg26mn + rpMg26gn - raMg26p - raMg25p &
        + rp26m + rp26g + rpAl27 + rpSi28 + rpSi29 + rpSi30 &
        + rpO17 + rpNe20
   
   ! Helium  duplicate - RPO17a is O17(p,a)N14
   ydot(42)= - r33 + r34 - 2.0*rpli - 2.0*rpbe + rlia + raC13n + raN &
        - rpN15a - rpO18a + raC14 + raN15g + raO17n + raO + 3.0*r3a &
        + raC + raO18 + raO18n - rpO17a - rpF19a + raF19p - rCC &
        + raNe21 + raNe21n + raNe22 + raNe22n - RnNa22a - rpNa23a &
        - Rn26ga - Rn26ma + raNa23nm + raNa23ng + raNe + raMg24 &
        - rpAl27a + raMg25 + raMg25n + raMg25p + raMg26 + raMg26n &
        + raMg26p + raAl27n - RnS32a - RnS33a - RnNi59a
   
   ! Carbon duplicate - rCC to alpha and Ne20 - should be ~50-50 with
   ! this and Na23 + p
   ydot(43)= - rpB11 + rpc - rpN15a - r3a + raC + 2.0*rCC + RnC12 + rCO
   
   ! Nitrogen duplicate
   ydot(44)= - rpC13 + RnN14p + rpN + raN - rpO17a
   
   ! Oxygen duplicate
   ydot(45)= - raC13n + rpO + raO - raC + RnO16 - rpF19a + rCO - rpN15
   
   ! Neon duplicate
   ydot(46)= - raO - rpF19 - RnF19 - rCC + RnNe20 + raNe - rpNa23a + rpNe20
   
   ! Make sure we don't run into problems when calculating the derivative for
   ! rates for which the abundance is 0 (the rate will be zero, but we
   ! calculate the derivative by rate/X, so X should be non-zero).
   where (xa2 <= 0.0d0) xa2 = -9.9d0
   
   !  EVALUATE DERIVATIVES
   ! YLIN(derived by,equation)
   ! RJS 1/9/03 - Added species for PP-I,II chain from Clayton 1963
   ! 1st eq gallinoes

   ! 2nd eq neutrons
   ylin(2, 2) = (rdn + RnNa22p + RnNa22a + Rn26mp + Rn26ma + &
        Rn26gp + Rn26ga + RnFe56 + RnFe57 + RnFe58 + &
        RnCo59 + RnNi58 + RnNi59 + RnNi60 + rnp + &
        RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14 + RnN14 + &
        RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21 + &
        RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + &
        RnAl27 + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 &
        + RnS33 + RnS34 + Rn26g+ RnS33a + RnN14p + RnNi59p &
        + RnNi59a + RnO17a + RnS32a)/xa2(2)
   ylin(4, 2) = RnHe3/xa2(4)
   ylin(5, 2) = RnLi7/xa2(5)
   ylin(8, 2) = (RnC13- raC13n)/xa2(8)
   ylin(9, 2) = RnC14/xa2(9)
   ylin(10,2) = RnN15/xa2(10)
   ylin(11,2) = (RnO17a - raO17n)/xa2(11)
   ylin(12,2) = (RnO18 - raO18n)/xa2(12)
   ylin(13,2) = RnF19/xa2(13)
   ylin(14,2) = (RnNe21 - raNe21n)/xa2(14)
   ylin(15,2) = (RnNe22 - raNe22n)/xa2(15)
   ylin(16,2) = (RnNa22p + RnNa22a)/xa2(16)
   ylin(17,2) = (RnNa23 - rpNa23n - raNa23nm - raNa23ng)/xa2(17)
   ylin(18,2) = RnMg24/xa2(18)
   ylin(19,2) = (RnMg25 - raMg25n)/xa2(19)
   ylin(20,2) = (RnMg26 - raMg26n - rpMg26mn - rpMg26gn)/xa2(20)
   ylin(21,2) = (Rn26mp + Rn26ma)/xa2(21)
   ylin(22,2) = (Rn26gp + Rn26ga + Rn26g)/xa2(22)
   ylin(23,2) = (- raAl27n + RnAl27)/xa2(23)
   ylin(24,2) = RnSi28/xa2(24)
   ylin(25,2) = RnSi29/xa2(25)
   ylin(26,2) = RnSi30/xa2(26)
   ylin(27,2) = RnP31/xa2(27)
   ylin(28,2) = (RnS32 + RnS32a)/xa2(28)
   ylin(29,2) = (RnS33 + RnS33a)/xa2(29)
   ylin(30,2) = RnS34/xa2(30)
   ylin(31,2) = RnFe56/xa2(31)
   ylin(32,2) = RnFe57/xa2(32)
   ylin(33,2) = RnFe58/xa2(33)
   ylin(36,2) = RnCo59/xa2(36)
   ylin(37,2) = RnNi58/xa2(37)
   ylin(38,2) = (RnNi59 + RnNi59p + RnNi59a)/xa2(38)
   ylin(39,2) = RnNi60/xa2(39)
   ylin(41,2) = (rnp - rpNa23n - rpMg26mn - rpMg26gn)/xa2(41)
   ylin(42,2) = - (raC13n + raO17n + raO18n + raNe21n + raNe22n + raNa23nm + raNa23ng + raMg25n + raMg26n + raAl27n)/xa2(42)
   ylin(43,2) = RnC12/xa2(43)
   ylin(45,2) = RnO16/xa2(45)
   ylin(46,2) = RnNe20/xa2(46)
   
   ! 3rd eq D
   ylin(2, 3) = - rnp/xa2(2)
   ylin(3, 3) = rpd/xa2(3)
   ylin(41,3) = (rpd - rnp - 2.0*rpp)/xa2(41)
   
   ! 4th eq He-3
   ylin(2, 4) = RnHe3/xa2(2)
   ylin(3, 4) = -rpd/xa2(3)
   ylin(4, 4) = (4.0*r33+r34+RnHe3)/xa2(4)
   ylin(41,4) = - rpd/xa2(41)
   ylin(42,4) = r34/xa2(42)
   
   ! 5th eq Li-7
   ylin(2, 5) = RnLi7/xa2(2)
   ylin(5, 5) = (rpli+rlia+RnLi7)/xa2(5)
   ylin(6, 5) = -rebe/xa2(6)
   ylin(41,5) = rpli/xa2(41)
   ylin(42,5) = rlia/xa2(42)
   
   ! 6th eq Be-7
   ylin(4, 6) = -r34/xa2(4)
   ylin(6, 6) = (rpbe + rebe)/xa2(6)
   ylin(41,6) = rpbe/xa2(41)
   ylin(42,6) = - r34/xa2(42)
   
   ! 7th eq B-11
   ylin(5, 7) = - rlia/xa2(5)
   ylin(7, 7) = rpB11/xa2(7)
   ylin(41,7) = rpB11/xa2(41)
   ylin(42,7) = - rlia/xa2(42)
   
   ! CNO elements. Will assume beta-decay is instantaneous (which is ok on all but the very shortest timescales)
   ! 8th eq C-13
   ylin(2, 8) = (RnC13 - RnC12)/xa2(2)
   ylin(8, 8) = (rpC13+raC13n+RnC13)/xa2(8)
   ylin(41,8) = (rpC13 - rpc)/xa2(41)
   ylin(42,8) = raC13n/xa2(42)
   ylin(43,8) = (- rpc - RnC12)/xa2(43)
   ! 9th eq C-14
   ylin(2, 9) = (RnC14 - RnC13 - RnN14p - RnO17a)/xa2(2)
   ylin(8, 9) = - RnC13/xa2(8)
   ylin(9, 9) = (rpC14 + raC14 + RnC14 + rdC14)/xa2(9)
   ylin(11,9) = - RnO17a/xa2(11)
   ylin(41,9) = rpC14/xa2(41)
   ylin(42,9) = raC14/xa2(42)
   ylin(44,9) = - RnN14p/xa2(44)
   
   ! 10th eq N-15
   ! Assume instantaneous C15 -> N15 beta decay
   ylin(2 ,10) = (RnN15 - RnC14 - RnN14)/xa2(2)
   ylin(9 ,10) = -(RnC14 + rpC14)/xa2(9)
   ylin(10,10) = (rpN15a+raN15g+RnN15+rpN15)/xa2(10)
   ylin(12,10) = - rpO18a/xa2(12)
   ylin(41,10) = (rpN15a - rpN - rpC14 + rpO18a+rpN15)/xa2(41)
   ylin(42,10) = raN15g/xa2(42)
   ylin(44,10) = (- RnN14 - rpN)/xa2(44)
   
   ! 11th eq O-17
   ylin(2 ,11) = (RnO17a - RnO16)/xa2(2)
   ylin(11,11) = (rpO17a+raO17n+RnO17a+rpO17)/xa2(11)
   ylin(41,11) = (rpO17a+rpO17)/xa2(41)
   ylin(42,11) = raO17n/xa2(42)
   ylin(45,11) = - RnO16/xa2(45)
   
   ! 12th eq O18
   ylin(2, 12) = RnO18/xa2(2)
   ylin(9 ,12) = - raC14/xa2(9)
   ylin(11,12) = - rpO17/xa2(11)
   ylin(12,12) = (rpO18 + rpO18a + raO18 + raO18n + RnO18)/xa2(12)
   ylin(41,12) = (rpO18 + rpO18a - rpO17)/xa2(41)
   ylin(42,12) = (raO18 + raO18n - raN)/xa2(42)
   ylin(44,12) = - raN/xa2(44)
   
   ! 13th eq F19
   ylin(2 ,13) = (RnF19 - RnNa22a - RnO18)/xa2(2)
   ylin(10,13) = - raN15g/xa2(10)
   ylin(12,13) = - (RnO18 + rpO18)/xa2(12)
   ylin(13,13) = (rpF19 + rpF19a + raF19p + RnF19)/xa2(13)
   ylin(16,13) = - RnNa22a/xa2(16)
   ylin(41,13) = (rpF19 + rpF19a - rpO18)/xa2(41)
   ylin(42,13) = (raF19p - raN15g)/xa2(42)
   
   ! 14th eq Ne21
   ylin(2, 14) = (RnNe21-RnNe20)/xa2(2)
   ylin(12,14) = - raO18n/xa2(12)
   ylin(14,14) = (rpNe21 + raNe21 + raNe21n + RnNe21)/xa2(14)
   ylin(41,14) = (rpNe21 - rpNe20)/xa2(41)
   ylin(42,14) = (raNe21 + raNe21n - raO18n)/xa2(42)
   ylin(46,14) = -rpNe20/xa2(46)
   
   ! 15th eq Ne22
   ylin(2 ,15) = (RnNe22 - RnNe21 - RnNa22p)/xa2(2)
   ylin(12,15) = - raO18/xa2(12)
   ylin(13,15) = - raF19p/xa2(13)
   ylin(14,15) = -RnNe21/xa2(14)
   ylin(15,15) = (rpNe22 + raNe22 + raNe22n + RnNe22)/xa2(15)
   ylin(16,15) = - (RnNa22p + rdNa22)/xa2(16)
   ylin(41,15) = rpNe22/xa2(41)
   ylin(42,15) = (raNe22 + raNe22n - raO18 - raF19p)/xa2(42)
   
   ! 16th eq Na22
   ylin(2 ,16) = (RnNa22p + RnNa22a)/xa2(2)
   ylin(14,16) = - rpNe21/xa2(14)
   ylin(16,16) = (RnNa22p + RnNa22a + rpNa22 + rdNa22)/xa2(16)
   ylin(41,16) = (rpNa22 - rpNe21)/xa2(41)
   
   ! 17th eq Na23
   ylin(2 ,17) = (RnNa23 - Rn26ma - Rn26ga - RnNe22)/xa2(2)
   ylin(15,17) = - (rpNe22+RnNe22)/xa2(15)
   ylin(17,17) = (rpNa23 + rpNa23n + rpNa23a + raNa23nm + raNa23ng + RnNa23)/xa2(17)
   ylin(21,17) = - Rn26ma/xa2(21)
   ylin(22,17) = - Rn26ga/xa2(22)
   ylin(41,17) = (rpNa23 + rpNa23n + rpNa23a - rpNe22)/xa2(41)
   ylin(42,17) = (raNa23nm + raNa23ng)/xa2(42)
   
   ! 18th eq Mg24
   ylin(2 ,18) = (RnMg24 - RnNa23)/xa2(2)
   ylin(14,18) = - raNe21n/xa2(14)
   ylin(17,18) = (- RnNa23 - rpNa23)/xa2(17)
   ylin(18,18) = (rpMg24 + raMg24 + RnMg24)/xa2(18)
   ylin(23,18) = - rpAl27a/xa2(23)
   ylin(41,18) = (rpMg24 - rpAl27a - rpNa23)/xa2(41)
   ylin(42,18) = (raMg24 - raNe - raNe21n)/xa2(42)
   ylin(46,18) = - raNe/xa2(46)
   
   ! 19th eq Mg25
   ylin(2 ,19) = (RnMg25 - RnMg24)/xa2(2)
   ylin(14,19) = - raNe21/xa2(14)
   ylin(15,19) = - raNe22n/xa2(15)
   ylin(18,19) = - RnMg24/xa2(18)
   ylin(19,19) = (rpMg25m + rpMg25g + raMg25 + raMg25n + raMg25p + RnMg25)/xa2(19)
   ylin(21,19) = - rg26mp/xa2(21)
   ylin(22,19) = - rg26gp/xa2(22)
   ylin(41,19) = (rpMg25m + rpMg25g)/xa2(41)
   ylin(42,19) = (raMg25 + raMg25n + raMg25p)/xa2(42)
   
   ! 20th eq Mg26
   ylin(2 ,20) = (RnMg26 - RnMg25 - Rn26mp - Rn26gp)/xa2(2)
   ylin(15,20) = - raNe22/xa2(15)
   ylin(19,20) = - RnMg25/xa2(19)
   ylin(20,20) = (rpMg26 + rpMg26mn + rpMg26gn + raMg26 + raMg26n + raMg26p + RnMg26)/xa2(20)
   ylin(21,20) = - (Rn26mp + rd26m)/xa2(21)
   ylin(22,20) = - (Rn26gp + rdAl26g)/xa2(22)
   ylin(41,20) = (rpMg26 + rpMg26mn + rpMg26gn)/xa2(41)
   ylin(42,20) = (raMg26 + raMg26n + raMg26p)/xa2(42)
   
   ! 21st eq Al26M
   ylin(2 ,21) = (Rn26mp + Rn26ma)/xa2(2)
   ylin(17,21) = - raNa23nm/xa2(17)
   ylin(19,21) = - rpMg25m/xa2(19)
   ylin(20,21) = - rpMg26mn/xa2(20)
   ylin(21,21) = (rg26mp + Rn26mp + Rn26ma + rp26m + rd26m)/xa2(21)
   ylin(41,21) = (rp26m - rpMg25m - rpMg26mn)/xa2(41)
   ylin(42,21) = - raNa23nm/xa2(42)
   
   ! 22nd eq Al26G
   ylin(2 ,22) = (Rn26gp + Rn26ga + Rn26g)/xa2(2)
   ylin(17,22) = - raNa23ng/xa2(17)
   ylin(19,22) = - rpMg25g/xa2(19)
   ylin(20,22) = - rpMg26gn/xa2(20)
   ylin(22,22) = (rg26gp + Rn26gp + Rn26ga + rp26g + rdAl26g + Rn26g)/xa2(22)
   ylin(41,22) = (rp26g - rpMg25g - rpMg26gn)/xa2(41)
   ylin(42,22) = - raNa23ng/xa2(42)
   
   ! 23rd eq Al27
   ylin(2 ,23) = (RnAl27 - RnMg26 - Rn26g)/xa2(2)
   ylin(20,23) = (-RnMg26 - rpMg26)/xa2(20)
   ylin(22,23) = -Rn26g/xa2(22)
   ylin(23,23) = (rpAl27 + rpAl27a + raAl27n + RnAl27)/xa2(23)
   ylin(41,23) = (rpAl27 + rpAl27a - rpMg26 - rp26g - rp26m)/xa2(41)
   ylin(42,23) = raAl27n/xa2(42)
   
   ! 24th eq Si28
   ylin(2, 24) = (RnSi28 - RnAl27)/xa2(2)
   ylin(18,24) = - raMg24/xa2(18)
   ylin(19,24) = - raMg25n/xa2(19)
   ylin(23,24) = (-RnAl27 - rpAl27)/xa2(23)
   ylin(24,24) = (RnSi28 + rpSi28)/xa2(24)
   ylin(41,24) = (rpSi28 - rpAl27)/xa2(41)
   ylin(42,24) = (- raMg24 - raMg25n)/xa2(42)
   ylin(43,24) = - rCO/xa2(43)
   
   ! 25th eq Si29
   ylin(2 ,25) = (RnSi29 - RnSi28 - RnS32a)/xa2(2)
   ylin(19,25) = - raMg25/xa2(19)
   ylin(20,25) = - (raMg26n + raMg26p)/xa2(20)
   ylin(24,25) = - (RnSi28 + rpSi28)/xa2(24)
   ylin(25,25) = (rpSi29 + RnSi29)/xa2(25)
   ylin(28,25) = - RnS32a/xa2(28)
   ylin(41,25) = (rpSi29 - rpSi28)/xa2(41)
   ylin(42,25) = - (raMg25 + raMg26n + raMg26p)/xa2(42)
   
   ! 26th eq Si30
   ylin(2 ,26) = (RnSi30 - RnSi29 - RnS33a)/xa2(2)
   ylin(20,26) = - raMg26/xa2(20)
   ylin(25,26) = - (rpSi29 + RnSi29)/xa2(25)
   ylin(26,26) = (rpSi30+RnSi30)/xa2(26)
   ylin(29,26) = - RnS33a/xa2(29)
   ylin(41,26) = (rpSi30 - rpSi29)/xa2(41)
   ylin(42,26) = - raMg26/xa2(42)
   
   ! 27th eq P31
   ylin(2 ,27) = (RnP31 - RnSi30)/xa2(2)
   ylin(26,27) = - (rpSi30 + RnSi30)/xa2(26)
   ylin(27,27) = RnP31/xa2(27)
   ylin(41,27) = - rpSi30/xa2(41)
   
   ! 28th eq S32
   ylin(2 ,28) = (RnS32 - RnP31 + RnS32a)/xa2(2)
   ylin(27,28) = - RnP31/xa2(27)
   ylin(28,28) = (RnS32 + RnS32a)/xa2(28)
   
   ! 29th eq S33
   ylin(2 ,29) = (RnS33 - RnS32 + RnS33a)/xa2(2)
   ylin(28,29) = - RnS32/xa2(28)
   ylin(29,29) = (RnS33 + RnS33a)/xa2(29)
   
   ! 30th eq S34
   ylin(2 ,30) = (RnS34 - RnS33 + RnS34s)/xa2(2)
   ylin(29,30) = - RnS33/xa2(29)
   ylin(30,30) = (RnS34 + RnS34s)/xa2(30)
   
   ! 31st eq Fe56
   ylin(2 ,31) = (RnFe56 - RnNi59a)/xa2(2)
   ylin(31,31) = RnFe56/xa2(31)
   ylin(38,31) = -RnNi59a/xa2(38)
   
   ! 32nd eq Fe57
   ylin(2 ,32) = (RnFe57 - RnFe56)/xa2(2)
   ylin(31,32) = - RnFe56/xa2(31)
   ylin(32,32) = RnFe57/xa2(32)
   
   ! 33rd eq Fe58
   ylin(2 ,33) = (RnFe58-RnFe57)/xa2(2)
   ylin(32,33) = - RnFe57/xa2(32)
   ylin(33,33) = RnFe58/xa2(33)
   
   ! 34th eq Fe59
   ylin(2 ,34) = (RnFe59 - RnFe58)/xa2(2)
   ylin(33,34) = - RnFe58/xa2(33)
   ylin(34,34) = (RnFe59 + rdFe59)/xa2(34)
   
   ! 35th eq Fe60
   ylin(2 ,35) = (RnFe60 - RnFe59)/xa2(2)
   ylin(34,35) = - RnFe59/xa2(34)
   ylin(35,35) = (RnFe60 + rdFe60)/xa2(35)
   
   ! 36th eq Co59
   ylin(2 ,36) = (RnCo59 - RnNi59p)/xa2(2)
   ylin(34,36) = - rdFe59/xa2(34)
   ylin(36,36) = RnCo59/xa2(36)
   ylin(38,36) = (-RnNi59p - rdNi59)/xa2(38)
   
   ! 37th eq Ni58
   ylin(2 ,37) = RnNi58/xa2(2)
   ylin(37,37) = RnNi58/xa2(37)
   
   ! 38th eq Ni59
   ylin(2 ,38) = (RnNi59 - RnNi58 + RnNi59p + RnNi59a)/xa2(2)
   ylin(37,38) = - RnNi58/xa2(37)
   ylin(38,38) = (RnNi59+ RnNi59p + RnNi59a + rdNi59)/xa2(38)
   
   ! 39th eq Ni60
   ylin(2 ,39) = (RnNi60 - RnNi59)/xa2(2)
   ylin(38,39) = - RnNi59/xa2(38)
   ylin(39,39) = RnNi60/xa2(39)
   
   ! 40th eq Ni61
   ylin(2 ,40) = (RnNi61s - RnNi60)/xa2(2)
   ylin(39,40) = - RnNi60/xa2(39)
   ylin(40,40) = RnNi61s/xa2(40)
   
   ! 41st eq H - duplicate
   ylin(2 ,41) = (rnp - RnN14p - rdn - Rn26mp - Rn26gp - RnNi59p)/xa2(2)
   ylin(3 ,41) = rpd/xa2(3)
   ylin(4 ,41) = - 4.0*r33/xa2(4)
   ylin(5 ,41) = rpli/xa2(5)
   ylin(6 ,41) = rpbe/xa2(6)
   ylin(7 ,41) = rpB11/xa2(7)
   ylin(8 ,41) = rpC13/xa2(8)
   ylin(9 ,41) = rpC14/xa2(9)
   ylin(10,41) = (rpN15a+rpN15)/xa2(10)
   ylin(11,41) = (rpO17a+rpO17)/xa2(11)
   ylin(12,41) = rpO18a/xa2(12)
   ylin(13,41) = (rpF19 + rpF19a  - raF19p)/xa2(13)
   ylin(14,41) = rpNe21/xa2(14)
   ylin(15,41) = rpNe22/xa2(15)
   ylin(16,41) = rpNa22/xa2(16)
   ylin(17,41) = (rpNa23 + rpNa23n + rpNa23a)/xa2(17)
   ylin(18,41) = rpMg24/xa2(18)
   ylin(19,41) = (rpMg25m + rpMg25g - raMg25p)/xa2(19)
   ylin(20,41) = (rpMg26 + rpMg26mn + rpMg26gn - raMg26p)/xa2(20)
   ylin(21,41) = (rp26m - Rn26mp - rg26mp)/xa2(21)
   ylin(22,41) = (rp26g - Rn26gp - rg26gp)/xa2(22)
   ylin(23,41) = (rpAl27 + rpAl27a)/xa2(23)
   ylin(24,41) = rpSi28/xa2(24)
   ylin(25,41) = rpSi29/xa2(25)
   ylin(26,41) = rpSi30/xa2(26)
   ylin(38,41) = - RnNi59p/xa2(38)
   ylin(41,41) = (4.0*rpp + rpd + rnp + rpli + rpbe + rpB11 - raF19p &
        + rpC13 + rpC14 + rpN + rpN15a + rpO18a + rpO + rpO17 + rpF19 &
        + rpF19a + rpNe21 + rpNe22 + rpNa22 + rpNa23 + rpNa23n &
        + rpNa23a + rpMg24 + rpAl27a + rpMg25m + rpMg25g + rpMg26 &
        + rpMg26mn + rpMg26gn + rp26m + rp26g + rpAl27 + rpSi28 &
        + rpSi29 + rpSi30 + rpN15 + rpO17 + rpNe20)/xa2(41)
   ylin(42,41) = - (raF19p + raMg25p + raMg26p)/xa2(42)
   ylin(43,41) = rpc/xa2(43)
   ylin(44,41) = (rpN - RnN14p)/xa2(44)
   ylin(45,41) = rpO/xa2(45)
   ylin(46,41) = rpNe20/xa2(46)
   
   ! 42nd eq He4 - duplicate
   ylin(2 ,42) = (- Rn26ga - Rn26ma - RnNa22a - RnS32a - RnS33a - RnNi59a)/xa2(2)
   ylin(4 ,42) = (- 4.0*r33 + r34)/xa2(4)
   ylin(5 ,42) = (rlia - 2.0*rpli)/xa2(5)
   ylin(6 ,42) = -2.0*rpbe/xa2(6)
   ylin(8 ,42) = raC13n/xa2(8)
   ylin(9 ,42) = raC14/xa2(9)
   ylin(10,42) = (raN15g - rpN15a)/xa2(10)
   ylin(11,42) = (raO17n - rpO17a)/xa2(11)
   ylin(12,42) = - rpO18a/xa2(12)
   ylin(13,42) = (raF19p - rpF19a)/xa2(13)
   ylin(14,42) = (raNe21 + raNe21n)/xa2(14)
   ylin(15,42) = (raNe22 + raNe22n)/xa2(15)
   ylin(16,42) = - RnNa22a/xa2(16)
   ylin(17,42) = (raNa23nm + raNa23ng)/xa2(17)
   ylin(18,42) = raMg24/xa2(18)
   ylin(19,42) = (raMg25 + raMg25n + raMg25p)/xa2(19)
   ylin(20,42) = (raMg26 + raMg26n + raMg26p)/xa2(20)
   ylin(21,42) = - Rn26ma/xa2(21)
   ylin(22,42) = - Rn26ga/xa2(22)
   ylin(23,42) = (raAl27n - rpAl27a)/xa2(23)
   ylin(28,42) = - RnS32a/xa2(28)
   ylin(29,42) = - RnS33a/xa2(29)
   ylin(38,42) = - RnNi59a/xa2(38)
   ylin(41,42) = (-2.0*rpli - 2.0*rpbe - rpN15a - rpO18a - rpO17a - rpF19a - rpNa23a - rpAl27a)/xa2(41)
   ylin(42,42) = (r34+rlia+raC13n+raN + raC14 + raN15g + raO17n      &
        + raO + 9.0*r3a + raC + raF19p + raNe21 + raNe21n + raNe22   &
        + raNe22n + raNa23nm + raNa23ng + raNe + raMg24 + raMg25     &
        + raMg25n + raMg25p + raMg26 + raMg26n + raMg26p + raAl27n)  &
        /xa2(42)
   ylin(43,42) = (raC - 2.0*rCC)/xa2(43)
   ylin(44,42) = raN/xa2(44)
   ylin(45,42) = raO/xa2(45)
   ylin(46,42) = raNe/xa2(46)
   
   ! 43rd eq C12 - duplicate
   ylin(7 ,43) = - rpB11/xa2(7)
   ylin(10,43) = - rpN15a/xa2(10)
   ylin(41,43) = (rpc - rpB11 - rpN15a)/xa2(41)
   ylin(42,43) = (raC - 9.0*r3a)/xa2(42)
   ylin(43,43) = (rpc + raC + 4.0*rCC + rCO)/xa2(43)
   ylin(45,43) = rCO/xa2(45)
   
   ! 44th eq N14 - duplicate
   ylin(2 ,44) = RnN14p/xa2(2)
   ylin(8 ,44) = - rpC13/xa2(8)
   ylin(11,44) = - rpO17a/xa2(11)
   ylin(41,44) = (rpN - rpC13 - rpO17a)/xa2(41)
   ylin(44,44) = (rpN + raN + RnN14p)/xa2(44)
   
   ! 45th eq O16 - duplicate
   ylin(2, 45) = RnO16/xa2(2)
   ylin(8 ,45) = - raC13n/xa2(8)
   ylin(10,45) = - rpN15/xa2(10)
   ylin(13,45) = - rpF19a/xa2(13)
   ylin(41,45) = (rpO - rpF19a - rpN15)/xa2(41)
   ylin(42,45) = (raO - raC13n - raC)/xa2(42)
   ylin(43,45) = (rCO - raC)/xa2(43)
   ylin(45,45) = (rpO + raO + RnO16 + rCO)/xa2(45)
   
   ! 46th eq Ne20 - duplicate
   ylin(2 ,46) = (RnNe20 - RnF19)/xa2(2)
   ylin(13,46) = (-rpF19 - RnF19)/xa2(13)
   ylin(17,46) = - rpNa23a/xa2(17)
   ylin(41,46) = (rpNe20 - rpF19 - rpNa23a)/xa2(41)
   ylin(42,46) = (- raO + raNe)/xa2(42)
   ylin(43,46) = - 2.0*rCC/xa2(43)
   ylin(45,46) = - raO/xa2(45)
   ylin(46,46) = (RnNe20 + raNe + rpNe20)/xa2(46)

   !     Shortcut the cycling of Al26M -> Mg26
   if (instantly_decay_Al26m) then
      ydot(20) = ydot(20) - ydot(21)
      ylin(:, 20) = ylin(:,20) - ylin(:, 21)
      ydot(21) = 0.0d0
      ylin(:, 21) = 0.0d0
   end if
   
   ! Remove any NaN issues, by blanking rates if abundance = 0
   do i=1, 50
      if (xa2(i) <= -9.9d0) then
         ylin(i,1:50) = 0d0
      end if
   end do
   y(1:50) = varnuc(1:50)
   forall (i = 1:50)
      ylin(1:50, i) = ylin(1:50, i)*baryn(i)*dmk
      ylin(i, i) = ylin(i, i) + dmk/dt
      ydot(i) = (ydot(i)*baryn(i)+dvarnuc(i)/dt)*dmk
   end forall

   ! Copy results to output variables
   fn1(1:50) = ydot(1:50)
   fn1(51:100) = y(1:50)
   fn1(101) = sg
   fn1(102) = dmt
   fn1(103) = sgth / mu    ! Thermohaline mixing coefficient (apart from grad mu term)
   fn1(104) = mu
   ! fn1(111:160) - advection flux, set above

   ! Calculate neutron out rates for use later
   !> \todo FIXME: we only use these when the nucleosynthesis has actually converged,
   !! so we might as well only bother to calculate these at that point. Should
   !! speed things up a bit.
   !<
   xa2(2) = varnuc(2)
   frac(:,k) = 0.0d0
   
   if (xa2(2) == 0.0d0) return      ! Don't bother if there are no neutrons
   
   !do i = 1, 50
   !   na2(i) = xa2(i)/baryn(i)
   !end do
   
   call nucrat2 ( at, rho/avm, zt, na, neo, nio, na2,  rrt2, nrate, drt )
   
   ! Copy results to local variables
   ! At this point, we only care about the neutron capture rates.
   RnFe56   = nrate(1);  RnFe57   = nrate(2);  RnFe58   = nrate(3);  RnCo59   = nrate(4);
   RnNi58   = nrate(5);  RnNi59   = nrate(6);  RnNi60   = nrate(7);  rnp      = nrate(8);
   RnHe3    = nrate(9);  RnLi7    = nrate(10); RnC12    = nrate(11); RnC13    = nrate(12);
   RnC14    = nrate(13); RnN14    = nrate(14); RnN15    = nrate(15); RnO16    = nrate(16);
   RnO18    = nrate(17); RnF19    = nrate(18); RnNe20   = nrate(19); RnNe21   = nrate(20);
   RnNe22   = nrate(21); RnNa23   = nrate(22); RnMg24   = nrate(23); RnMg25   = nrate(24);
   RnMg26   = nrate(25); RnAl27   = nrate(26); RnSi28   = nrate(27); RnSi29   = nrate(28);
   RnSi30   = nrate(29); RnP31    = nrate(30); RnS32    = nrate(31); RnS33    = nrate(32);
   RnS34    = nrate(33); Rn26g    = nrate(34); RnS33a   = nrate(35); RnN14p   = nrate(36);
   RnNi59p  = nrate(37); RnNi59a  = nrate(38); RnO17a   = nrate(39); Rn26gad  = nrate(40);
   RnS32a   = nrate(41); RnFe59   = nrate(42); RnFe60   = nrate(43); RnS34s   = nrate(44);
   RnNi61s  = nrate(45);

   ratetot = rdn + RnNa22p + RnNa22a + Rn26mp + Rn26ma + Rn26g&
        + Rn26gp + Rn26ga&
        + RnFe56 + RnFe57 + RnFe58 + RnCo59 + RnNi58 + RnNi59&
        + RnNi60 + rnp + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14&
        + RnN14 + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21&
        + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + RnAl27&
        + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 + RnS33 + RnS34&
        + RnS33a + RnN14p + RnNi59p + RnNi59a + RnO17a + RnS32a&
        + RnS34s + RnFe59 + RnFe60 + RnNi61s
   if (ratetot == 0d0) ratetot = 1d0
   
   frac(1 ,k) = - RnNi61s
   frac(2 ,k) = 0d0
   frac(3 ,k) = - rnp
   frac(4 ,k) = RnHe3
   frac(5 ,k) = RnLi7
   frac(6 ,k) = 0d0
   frac(7 ,k) = 0d0
   frac(8 ,k) = RnC13 - RnC12
   frac(9 ,k) = RnC14 - RnN14p - RnO17a - RnC13
   frac(10,k) = RnN15 -  RnC14 - RnN14
   frac(11,k) = - RnO16 + RnO17a
   frac(12,k) = RnO18
   frac(13,k) = - RnNa22a - RnO18 + RnF19
   frac(14,k) = RnNe21 - RnNe20
   frac(15,k) = - RnNe21 + RnNe22 - RnNa22p
   frac(16,k) = RnNa22p + RnNa22a
   frac(17,k) = - Rn26ga - Rn26ma - RnNe22 + RnNa23
   frac(18,k) = - RnNa23 + RnMg24
   frac(19,k) = RnMg25 - RnMg24
   frac(20,k) = - Rn26mp - Rn26gp + RnMg26 - RnMg25
   frac(21,k) = Rn26mp + Rn26ma
   frac(22,k) = Rn26gp + Rn26ga + Rn26g
   frac(23,k) = - RnMg26 + RnAl27 - Rn26g
   frac(24,k) = - RnAl27 + RnSi28
   frac(25,k) = - RnSi28 + RnSi29 - RnS32a
   frac(26,k) = - RnSi29 + RnSi30 - RnS33a
   frac(27,k) = - RnSi30 + RnP31
   frac(28,k) = - RnP31 + RnS32 + RnS32a
   frac(29,k) = - RnS32 + RnS33 + RnS33a
   frac(30,k) = - RnS33 + RnS34 + RnS34s
   frac(31,k) = RnFe56 - RnNi59a
   frac(32,k) = RnFe57 - RnFe56
   frac(33,k) = RnFe58 - RnFe57
   frac(34,k) = RnFe59 - RnFe58
   frac(35,k) = - RnFe59 + RnFe60
   frac(36,k) = RnCo59 - RnNi59p
   frac(37,k) = RnNi58
   frac(38,k) = RnNi59 - RnNi58 + RnNi59a + RnNi59p
   frac(39,k) = RnNi60 - RnNi59
   frac(40,k) = RnNi61s - RnNi60
   frac(41,k) = rnp - rdn - Rn26mp - Rn26gp - RnNi59p - RnN14p
   frac(42,k) = - Rn26ga - Rn26ma - RnS32a - RnNa22a - RnS33a - RnNi59a  ! NEEDS RnHe4
   frac(43,k) = RnC12
   frac(44,k) = RnN14p
   frac(45,k) = RnO16
   frac(46,k) = - RnF19
   
   frac(1:48,k) = frac(1:48,k)/ratetot
   
end subroutine funcs2



subroutine equns2 ( jk, kl, kq, keq, Jstar, fn2, dxt, equ, dequ )
   use real_kind
   use mesh
   use settings
   use nucleosynthesis
   
   implicit none
   integer, intent(in) :: jk, kl, kq, keq, Jstar
   real(double), intent(in) :: fn2(3, NFUNC), dxt(3, NVAR, NFUNC)
   real(double), intent(out) :: equ(NEQ), dequ(NVAR, 3, NEQ)

   integer :: i1
   real(double) :: xt(3, 50), x(3, 50), sg(3), dmt(3), sgth(3), mu(3), ff(3, 50)
   real(double) :: sg12,sg23, dmu12, dmu23
   real(double), external :: pstv
   
   ! Copy input data to local variables
   xt(1:3, 1:50) = fn2(1:3, 1:50)
   x(1:3, 1:50) = fn2(1:3, 51:100)
   sg(1:3) = fn2(1:3, 101)
   dmt(1:3) = fn2(1:3, 102)
   sgth(1:3) = fn2(1:3, 103)
   mu(1:3) = fn2(1:3, 104)
   ff(1:3, 1:50) = fn2(1:3, 111:160)

   dequ(1:keq,1:3,1:keq) = 0.0d0
   equ(1:keq) = 0.0d0

   if ( jk + kl == 2 ) then
      ! Next-to-surface boundary conditions
      !do i1 = 1, keq
      !   dequ(i1, 2, i1) = dequ(i1, 2, i1) - 1.0
      !   dequ(i1, 3, i1) = 1.0
      !   equ(i1) = x(3,i1) - x(2,i1)
      !end do

      if (kl == 0) then    ! (3) is nearest to centre, so we have the grid points | 3 | 2 | AC |
         sg12 = 0.5d0*(sg(1) + sg(2)) - pstv(kq*dmt(2), 0.0d0)
         sg23 = 0.5d0*(sg(2) + sg(3)) - pstv(-kq*dmt(3), 0.0d0)
         equ(1:keq) = sg23*(x(3,1:keq) - x(2,1:keq))  - sg12*(x(2,1:keq) - nucacc(1:keq, Jstar)) - xt(2, 1:keq)
         dequ(1:keq, 2, 1:keq) = -dxt(2, 1:keq, 1:keq)
         forall (i1 = 1:keq)
            dequ(i1, 1, i1) = 0.0d0
            dequ(i1, 2, i1) = dequ(i1, 2, i1) - sg12 - sg23
            dequ(i1, 3, i1) = sg23
         end forall
      else                 ! (3) is nearest to surface, we have the grid points | 1 | 2 | 3 | AC |
         sg12 = 0.5d0*(sg(2) + sg(3)) - pstv(kq*dmt(3), 0.0d0)&
             + 0.5d0*(sgth(2)+sgth(3))*pstv(kq * (mu(2)-mu(3)), 0.0d0)
         sg23 = -pstv(dmt(3), 0.0d0)

         equ(1:keq) = sg23*(nucacc(1:keq, Jstar) - x(3, 1:keq)) - sg12*(x(3, 1:keq) - x(2, 1:keq)) - xt(3, 1:keq)
         dequ(1:keq, 3, 1:keq) = -dxt(3, 1:keq, 1:keq)
         forall (i1 = 1:keq)
            dequ(i1, 1, i1) = 0.0d0
            dequ(i1, 2, i1) = sg12
            dequ(i1, 3, i1) = dequ(i1, 3, i1) - sg12 - sg23
         end forall
      end if

      ! Add gravitational settling terms
      ! The surface value of the flux should vanish
      if (cgrs > 0.0d0) then
         equ(1:keq) = equ(1:keq) + kq*(ff(2, 1:keq)*x(2, 1:keq))
         ! Derivatives; the flux is proportional to the abundance (to first order anyway)
         forall (i1=1:keq)
            dequ(i1, 2, i1) = dequ(i1, 2, i1) + kq*ff(2, i1)
         end forall
      end if

   else if ( jk + kl >= kh + 1 ) then
      ! Next-to-central boundary conditions, consistent with equns1
      dmu23 = kq * (mu(2)-mu(3))
      sg23 = 0.5d0*(sg(2) + sg(3)) + 0.5d0*(sgth(2)+sgth(3))*pstv(dmu23, 0.0d0)
      forall (i1 = 1:keq)
         dequ(1:keq, 3, i1) = dxt(3 - kl, 1:keq, i1)
         dequ(i1, 2, i1) = - sg23
         dequ(i1, 3, i1) = dequ(i1, 3, i1) + sg23
         equ(i1) = sg23*(x(3, i1) - x(2, i1)) + xt(3 - kl, i1)
      end forall
      ! Add gravitational settling terms
      ! The central value of the flux should vanish
      if (cgrs > 0.0d0) then
         equ(1:keq) = equ(1:keq) - kq*ff(3, 1:keq)*x(3, 1:keq)
         ! Derivatives; the flux is proportional to the abundance (to first order anyway)
         forall (i1=1:keq)
            dequ(i1, 3, i1) = dequ(i1, 3, i1) - kq*ff(3, i1)
         end forall
      end if
   else
      ! Interior points
      dmu12 = kq * (mu(1)-mu(2))
      dmu23 = kq * (mu(2)-mu(3))
      sg12 = 0.5d0*(sg(1) + sg(2)) + 0.5d0*(sgth(1)+sgth(2))*pstv(dmu12, 0.0d0) - pstv(kq*dmt(2), 0.0d0)
      sg23 = 0.5d0*(sg(2) + sg(3)) + 0.5d0*(sgth(2)+sgth(3))*pstv(dmu23, 0.0d0) - pstv(-kq*dmt(3), 0.0d0)
      dequ(1:keq, 2, 1:keq) = -dxt(2, 1:keq, 1:keq)
      forall (i1 = 1:keq)
         dequ(i1, 1, i1) = sg12
         dequ(i1, 2, i1) = dequ(i1, 2, i1) - sg12 - sg23
         dequ(i1, 3, i1) = sg23
         equ(i1) = sg23*(x(3, i1) - x(2, i1)) - sg12*(x(2, i1) - x(1, i1)) - xt(2, i1)
      end forall
      ! Add gravitational settling terms
      if (cgrs > 0.0d0) then
         equ(1:keq) = equ(1:keq) - kq*(ff(3, 1:keq)*x(3, 1:keq) - ff(2, 1:keq)*x(2, 1:keq))
         ! Derivatives; the flux is proportional to the abundance (to first order anyway)
         forall (i1=1:keq)
            dequ(i1, 3, i1) = dequ(i1, 3, i1) - kq*ff(3, i1)
            dequ(i1, 2, i1) = dequ(i1, 2, i1) + kq*ff(2, i1)
         end forall
      end if
   end if
end subroutine equns2



! ------------------------------------------------------------------------------
!  NEUTRON
!   Perform neutron captures
!   These are not done during the main part of the nucleosynthesis
!   calculation because the reactions are so fast compared to the other
!   reactions occurring at the same time that they are almost instantaneous
!   anyway. Trying to put them in the same matrix would cause
!   non-convergence.
! ------------------------------------------------------------------------------
subroutine neutron(Jstar)
   use real_kind
   use mesh
   use nucleosynthesis
   use nucleosynthesis_neutron_ratios
   use nucleosynthesis_atomic_data
   ! RJS - deal with neutrons, basically assuming they're in equilibrium

   implicit none
   integer, intent(in) :: Jstar

   integer :: k,i

   forall (k=1:kh, i=1:50) DHnuc(Jstar,i,k) = DHnuc(Jstar,i,k) - frac(i,k)*baryn(i)*DHnuc(Jstar,2,k)
   ! reset neutron abundance to zero
   DHnuc(Jstar,2,1:kh) = 0.d0
end subroutine neutron



! ------------------------------------------------------------------------------
!  NUCRAT2
!   Compute thermonuclear reaction rates for minor isotopes, for
!   nucleosynthesis (funcs2).
! ------------------------------------------------------------------------------
!  Input:
!     TL          - log temperature (in Kelvin)
!     RHB         - Baryon mass density: 1 amu * number of baryons per cm3
!     ZT          - Term in the expression for strong screening
!     N           - Abundances of major isotopes
!     Ne          - Electron density (per baryon)
!     Ni          - Ion density (per baryon)
!     NA(:)       - Abundances of minor isotopes
!  Output:
!     rrt2(:)     - Reaction rates
!     nrate(:)    - Neutron reaction rates
!     drt(:)      - Decay rates
! ------------------------------------------------------------------------------
subroutine nucrat2(tl, rhb, zt, n, ne, ni, na, rrt2, nrate, drt)
   use real_kind
   use nucleosynthesis
   use constants
   use reaction_rate_data
   
   implicit none
   real(double), intent(in) :: tl, rhb, zt
   real(double), intent(in) :: n(9), ne, ni
   real(double), intent(in) :: na(50)
   real(double), intent(out) :: rrt2(92), nrate(45), drt(8)

   integer :: j,it
   real(double) :: wc,wb,wa,xb,vl,za,zb,zc,zd,tf,rn,tt,tu,rr,scrn,strn,dstr
   real(double) :: fpng,rpna,fccg,rcca,cln2

   real(double) :: rpp, r33, r34, rBe, rBp
   real(double) :: rpc, rpN, rpO, r3a, raC, raN, raO, raNe, rCC, rCO, rOO, rgNe
   real(double) :: rgMg, rCCg, rpng

   real(double) :: nng, nn, nn2,nn3,nnl7,nnb7,nn11,nn13,nN14
   real(double) :: nN15,nn17,nn18,nn19,nNe21,nNe22,nNa22,nNa23,nMg24,nMg25,nMg26
   real(double) :: n26m,n26g,nAl27,nSi28,nSi29,nSi30,nP31,nS32,nS33,nS34
   real(double) :: nFe56,nFe57,nFe58,nFe59,nFe60,nCo59,nNi58,nNi59,nNi60,nNi61,nn1
   real(double) :: nn4,nn12,nnN14,nn16,nn20

   real(double) :: rdAl26g,rdNa22,rd26m,rdFe59,rdFe60,rdNi59,rdn,rdC14

   real(double) :: csa, csb, csc, csd, cxD
   data csa, csb, csc, csd, cxD /0.d624, 0.d316, 0.d460, 0.d38, 0.d86/

   real(double) :: cbrt
   
   ! Copy abundances to local variables
   nng   = na(1);  nn    = na(2);  nn2   = na(3);  nn3   = na(4);
   nnl7  = na(5);  nnb7  = na(6);  nn11  = na(7);  nn13  = na(8);
   nN14  = na(9);  nN15  = na(10); nn17  = na(11); nn18  = na(12);
   nn19  = na(13); nNe21 = na(14); nNe22 = na(15); nNa22 = na(16);
   nNa23 = na(17); nMg24 = na(18); nMg25 = na(19); nMg26 = na(20);
   n26m  = na(21); n26g  = na(22); nAl27 = na(23); nSi28 = na(24);
   nSi29 = na(25); nSi30 = na(26); nP31  = na(27); nS32  = na(28);
   nS33  = na(29); nS34  = na(30); nFe56 = na(31); nFe57 = na(32);
   nFe58 = na(33); nFe59 = na(34); nFe60 = na(35); nCo59 = na(36);
   nNi58 = na(37); nNi59 = na(38); nNi60 = na(39); nNi61 = na(40);
   nn1   = na(41); nn4   = na(42); nn12  = na(43); nnN14 = na(44);
   nn16  = na(45); nn20  = na(46);

   ! Electron screening theory from Graboske, DeWitt, Grossman & Cooper (1973),
   ! for strong (ZA, ZB, ZC) are intermediate screening (ZD). The reaction
   ! dependent charge parameters are stored in CZA ... CZD.
   wc = dot_product(n(1:9), vz(1:9))
   wc = wc/ni
   wb = ne/ni
   wa = zt*zt/(wb*wb)
   xb = cbrt(dabs(wb))
   vl = cpl*dsqrt(dabs(ni)*rhb*dexp(-3.0d0*tl))
   za = csa*xb*vl**(2.0/3.0)
   zb = csb*xb*za
   zc = csc/(xb*xb)
   zd = csd*wc*wa*(vl/(wa*zt))**cxD
   
   ! Reaction rates interpolated in T, mostly from Caughlan & Fowler (1988)
   tf = tl/cln
   do j = 1, nuc_num_rates
      rn = 0.0d0
      tt = 50.0*(tf - 6.0) + 1.0
      if ( tt  >=  1.0d0 ) then
         it = max(1, min(199, int(tt)))
         tt = tt - it
         tu = 1.0d0 - tt
         rr = tu*nucsyn_crt(it, j) + tt*nucsyn_crt(it+1, j)
         if ( rr  >=  -50.0d0 ) then
            scrn = zd*czd(j)
            strn = za*cza(j) + zb*czb(j)
            dstr = zc*czc(j)
            if (dstr  <  0.29*strn) scrn = min(scrn, strn - dstr)
            rn = exp(cln*rr + scrn)
         end if
      end if
      rrt2(j) = rn
   end do
   
   ! Sort out rates
   rpp = rrt2(1)
   r33 = rrt2(2)
   r34 = rrt2(3)
   rBe = rrt2(4)
   rBp = rrt2(5)
   rpc = rrt2(6)
   rpN = rrt2(7)
   rpO = rrt2(8)
   r3a = rrt2(9)
   raC = rrt2(10)
   raN = rrt2(11)
   raO = rrt2(12)
   raNe = rrt2(13)
   rCC = rrt2(14)
   rCO = rrt2(15)
   rOO = rrt2(16)
   rgNe = rrt2(17)
   rgMg = rrt2(18)
   rCCg = rrt2(19)
   rpng = rrt2(20)
   
   ! Multiply with density and abundances to get rates per baryon per second,
   ! note that abundances of He3 and Be7 are not needed in equilibrium
   rpp = rhb*nn1*nn1*rpp/2.0
   r33 = rhb*nn3*nn3*r33/2.0
   r34 = rhb*nn3*nn4*r34
   rBe = rhb*ne*rBe
   rBp = rhb*nn1*rBp
   rpc = rhb*nn1*nn12*rpc
   rpN = rhb*nn1*nnN14*rpN
   rpO = rhb*nn1*nn16*rpO
   r3a = rhb*rhb*nn4*nn4*nn4*r3a/6.0
   raC = rhb*nn4*nn12*raC
   raN = rhb*nn4*nnN14*raN
   raO = rhb*nn4*nn16*raO
   raNe = rhb*nn4*nn20*raNe
   rCC = rhb*nn12*nn12*rCC/2.0
   rCO = rhb*nn12*nn16*rCO
   rOO = rhb*nn16*nn16*rOO/2.0
   rgNe = nn20*rgNe
   rgMg = nMg24*rgMg
   
   ! Branching of pN and CC reactions
   fpng = 8.0d-4
   rpna = (1.0 - fpng)*rpN
   rpng = fpng*rpN
   rpN = rpna
   fccg = rCCg
   rcca = (1.0 - fccg)*rCC
   rCCg = fccg*rCC
   rCC = rcca
   
   ! Put rates back to RRT
   rrt2(1) = rpp
   rrt2(2) = r33
   rrt2(3) = r34
   rrt2(4) = rBe
   rrt2(5) = rBp
   rrt2(6) = rpc
   rrt2(7) = rpN
   rrt2(8) = rpO
   rrt2(9) = r3a
   rrt2(10) = raC
   rrt2(11) = raN
   rrt2(12) = raO
   rrt2(13) = raNe
   rrt2(14) = rCC
   rrt2(15) = rCO
   rrt2(16) = rOO
   rrt2(17) = rgNe
   rrt2(18) = rgMg
   rrt2(19) = rCCg
   rrt2(20) = rpng
   
   ! Minor variable reaction rates - 3/9/03 RJS
   rrt2(21) = rhb*nn1*nn2*rrt2(21)
   rrt2(22) = rhb*nn1*nnb7*rrt2(22)
   rrt2(23) = rhb*nn1*nnl7*rrt2(23)
   rrt2(24) = rhb*ne*nnb7*rrt2(24)
   rrt2(25) = rhb*nn1*nn13*rrt2(25)
   rrt2(26) = rhb*nn1*nN15*rrt2(26)
   rrt2(27) = rhb*nn1*nn17*rrt2(27)
   rrt2(28) = rhb*nn4*nn13*rrt2(28)
   rrt2(29) = rhb*nn4*nN15*rrt2(29)
   rrt2(30) = rhb*nn4*nn17*rrt2(30)
   rrt2(31) = rhb*nn4*nnl7*rrt2(31)
   rrt2(32) = rhb*nn1*nn11*rrt2(32)
   rrt2(91) = rhb*nn1*nn17*rrt2(91)
   rrt2(92) = rhb*nn1*nn20*rrt2(92)
   
   ! C14 reactions
   rrt2(33) = rhb*nn1*nN14*rrt2(33)
   rrt2(34) = rhb*nn4*nN14*rrt2(34)
   
   ! O18 reactions
   rrt2(35) = rhb*nn1*nn18*rrt2(35)
   rrt2(36) = rhb*nn1*nn18*rrt2(36)
   rrt2(37) = rhb*nn4*nn18*rrt2(37)
   rrt2(38) = rhb*nn4*nn18*rrt2(38)
   
   ! F19 reactions
   rrt2(39) = rhb*nn1*nn19*rrt2(39)
   rrt2(40) = rhb*nn1*nn19*rrt2(40)
   rrt2(41) = rhb*nn4*nn19*rrt2(41)
   
   ! Ne21
   rrt2(42) = rhb*nn1*nNe21*rrt2(42)
   rrt2(43) = rhb*nn4*nNe21*rrt2(43)
   rrt2(44) = rhb*nn4*nNe21*rrt2(44)
   
   ! Ne22
   rrt2(45) = rhb*nn1*nNe22*rrt2(45)
   rrt2(46) = rhb*nn4*nNe22*rrt2(46)
   rrt2(47) = rhb*nn4*nNe22*rrt2(47)
   
   ! Na22
   rrt2(48) = rhb*nn*nNa22*rrt2(48)
   rrt2(49) = rhb*nn*nNa22*rrt2(49)
   rrt2(50) = rhb*nn1*nNa22*rrt2(50)
   
   ! Na23
   rrt2(51) = rhb*nn1*nNa23*rrt2(51)
   rrt2(52) = rhb*nn1*nNa23*rrt2(52)
   rrt2(53) = rhb*nn1*nNa23*rrt2(53)
   
   ! Mg24
   rrt2(54) = rhb*nn1*nMg24*rrt2(54)
   rrt2(55) = rhb*nn4*nMg24*rrt2(55)
   
   ! Mg25
   rrt2(56) = rhb*nn1*nMg25*rrt2(56)
   rrt2(57) = rhb*nn1*nMg25*rrt2(57)
   rrt2(58) = rhb*nn1*nMg25*rrt2(58)
   rrt2(59) = rhb*nn4*nMg25*rrt2(59)
   rrt2(60) = rhb*nn4*nMg25*rrt2(60)
   rrt2(61) = rhb*nn4*nMg25*rrt2(61)
   
   ! Mg26
   rrt2(62) = rhb*nn1*nMg26*rrt2(62)
   rrt2(63) = rhb*nn1*nMg26*rrt2(63)
   rrt2(64) = rhb*nn1*nMg26*rrt2(64)
   rrt2(65) = rhb*nn1*nMg26*rrt2(65)
   rrt2(66) = rhb*nn4*nMg26*rrt2(66)
   rrt2(67) = rhb*nn4*nMg26*rrt2(67)
   rrt2(68) = rhb*nn4*nMg26*rrt2(68)
   
   ! Blank beause I put in a reaction that didn't exist!
   rrt2(69) = 0d0
   rrt2(70) = 0d0
   rrt2(71) = 0d0
   rrt2(72) = 0d0
   
   ! Al26M
   rrt2(73) = rhb*n26m*rrt2(73)
   rrt2(74) = rhb*nn*n26m*rrt2(74)
   rrt2(75) = rhb*nn*n26m*rrt2(75)
   rrt2(76) = rhb*nn1*n26m*rrt2(76)
   
   ! Al26G
   rrt2(77) = rhb*n26g*rrt2(77)
   rrt2(78) = rhb*nn*n26g*rrt2(78)
   rrt2(79) = rhb*nn*n26g*rrt2(79)
   rrt2(80) = rhb*nn1*n26g*rrt2(80)
   
   ! Al27
   rrt2(81) = rhb*nn1*nAl27*rrt2(81)
   rrt2(82) = rhb*nn1*nAl27*rrt2(82)
   rrt2(83) = rhb*nn4*nAl27*rrt2(83)
   
   ! Na23(a,n)Al26TGM
   rrt2(84) = rhb*nn4*nNa23*rrt2(84)
   rrt2(85) = rhb*nn4*nNa23*rrt2(85)
   rrt2(86) = rhb*nn4*nNa23*rrt2(86)
   
   ! Si reactions
   rrt2(87) = rhb*nn1*nSi28*rrt2(87)
   rrt2(88) = rhb*nn1*nSi29*rrt2(88)
   rrt2(89) = rhb*nn1*nSi30*rrt2(89)
   
   ! N15
   rrt2(90) = rhb*nn1*nN15*rrt2(90)

   ! Unstable particle decays
   ! Al26G t_0.5 = 0.72 Myr
   cln2 = 0.69314718
   rdAl26g = (rhb/26.0)*n26g*cln2/2.27d13
   
   ! C14 t_0.5 = 5730 yr
   rdC14 = (rhb/14.0)*nN14*cln2/1.81d11
   
   ! Na22 t_0.5 = 2.6 yr
   rdNa22 = (rhb/22.0)*nNa22*cln2/8.199d7
   
   ! Al26M t_0.5 = 6 s Should just shortcut the network...
   rd26m = (rhb/26.0)*n26m*cln2/6.0
   
   ! Fe59 t_0.5 = 44.6 d
   rdFe59 = (rhb/59.0)*nFe59*cln2/3.85d6
   
   ! Fe60 t_0.5 = 1.5 Myr
   rdFe60 = (rhb/60.0)*nFe60*cln2/4.73d13
   
   ! Ni59 t_0.5 = 0.075 Myr
   rdNi59 = (rhb/59.0)*nNi59*cln2/2.365d12
   
   ! Free n t_0.5 = 10.3 min
   rdn = (rhb/1.0)*nn*cln2/6.18d2
   
   drt(1:8) = (/ rdAl26g, rdNa22, rd26m, rdFe59, rdFe60, rdNi59, rdn, rdC14 /)

   ! (n,g) reactions
   if (nn == 0.0d0) then
      nrate(:) = 0.0d0
      return
   end if
   do j = 1, 45
      rn = 0.0d0
      tt = 50.0*(tf - 6.0) + 1.0
      if ( tt  >=  1.0d0 ) then
         it = max(1, min(199, int(tt)))
         tt = tt - it
         tu = 1.0d0 - tt
         rr = tu*nucsyn_nrt(it, j) + tt*nucsyn_nrt(it+1, j)
         if (rr >= -50.0d0) rn = exp(cln*rr)
      end if
      nrate(j) = rn
   end do
   
   nrate(1) = rhb*nn*nFe56*nrate(1)
   nrate(2) = rhb*nn*nFe57*nrate(2)
   nrate(3) = rhb*nn*nFe58*nrate(3)
   nrate(4) = rhb*nn*nCo59*nrate(4)
   nrate(5) = rhb*nn*nNi58*nrate(5)
   nrate(6) = rhb*nn*nNi59*nrate(6)
   nrate(7) = rhb*nn*nNi60*nrate(7)
   nrate(8) = rhb*nn*nn1*nrate(8)
   nrate(9) = rhb*nn*nn3*nrate(9)
   nrate(10) = rhb*nn*nnl7*nrate(10)
   nrate(11) = rhb*nn*nn12*nrate(11)
   nrate(12) = rhb*nn*nn13*nrate(12)
   nrate(13) = rhb*nn*nN14*nrate(13)
   nrate(14) = rhb*nn*nnN14*nrate(14)
   nrate(15) = rhb*nn*nN15*nrate(15)
   nrate(16) = rhb*nn*nn16*nrate(16)
   nrate(17) = rhb*nn*nn18*nrate(17)
   nrate(18) = rhb*nn*nn19*nrate(18)
   nrate(19) = rhb*nn*nn20*nrate(19)
   nrate(20) = rhb*nn*nNe21*nrate(20)
   nrate(21) = rhb*nn*nNe22*nrate(21)
   nrate(22) = rhb*nn*nNa23*nrate(22)
   nrate(23) = rhb*nn*nMg24*nrate(23)
   nrate(24) = rhb*nn*nMg25*nrate(24)
   nrate(25) = rhb*nn*nMg26*nrate(25)
   nrate(26) = rhb*nn*nAl27*nrate(26)
   nrate(27) = rhb*nn*nSi28*nrate(27)
   nrate(28) = rhb*nn*nSi29*nrate(28)
   nrate(29) = rhb*nn*nSi30*nrate(29)
   nrate(30) = rhb*nn*nP31*nrate(30)
   nrate(31) = rhb*nn*nS32*nrate(31)
   nrate(32) = rhb*nn*nS33*nrate(32)
   nrate(33) = rhb*nn*nS34*nrate(33)
   nrate(34) = rhb*nn*n26g*nrate(34)
   nrate(35) = rhb*nn*nS33*nrate(35)
   nrate(36) = rhb*nn*nnN14*nrate(36)
   nrate(37) = rhb*nn*nNi59*nrate(37)
   nrate(38) = rhb*nn*nNi59*nrate(38)
   nrate(39) = rhb*nn*nn17*nrate(39)
   nrate(40) = rhb*nn*n26g*nrate(40)
   nrate(41) = rhb*nn*nS32*nrate(41)
   nrate(42) = rhb*nn*nFe59*nrate(42)
   nrate(43) = rhb*nn*nFe60*nrate(43)
   nrate(44) = rhb*nn*nS34*nrate(44)
   nrate(45) = rhb*nn*nNi61*nrate(45)
   
end subroutine nucrat2



! ------------------------------------------------------------------------------
!  MIX_GRIDPOINT
!   Mix (homogenise) the specified grid point for the specified star.
!   Modifies Hnuc for *1 and *2
! ------------------------------------------------------------------------------
subroutine mix_gridpoint(Jstar, k)
   use real_kind
   use nucleosynthesis
   use settings
   implicit none
   integer, intent(in) :: Jstar, k
   real(double) :: avgx(50), dm
   integer :: i

   ! Total mass in shells
   dm = sum(ht(Jstar, 4, k-1:k+1))

   ! Compute average abundence around this mesh point
   forall (i=1:50)
      avgx(i) = dot_product(Hnuc(Jstar, i, k-1:k+1), ht(Jstar, 4, k-1:k+1)) / dm
   end forall

   ! Mix meshpoints
   forall (i=1:50)
      Hnuc(Jstar, i, k-1:k+1) = avgx(i)
   end forall
end subroutine mix_gridpoint



! ------------------------------------------------------------------------------
!  SMOOTH_SPIKES
!   Mix (homogenise) spikes in the composition profile, for low temperatures.
!   Modifies Hnuc for *1 and *2
! ------------------------------------------------------------------------------
subroutine smooth_spikes()
   use real_kind
   use nucleosynthesis
   use settings
   use constants
   use mesh, only: H
   implicit none
   real(double) :: dx12, dx23
   integer :: Jstar, k, i

   do Jstar = 1, ktw
      do k = 2, kh_nuc-1
         if (H(2 + 24*(Jstar - 1), k) > 6.0*CLN ) cycle
         do i = 1, 50
            dx12 = hnuc(jstar, i, k-1) - hnuc(jstar, i, k)
            dx23 = hnuc(jstar, i, k) - hnuc(jstar, i, k+1)
            if ( dx12*dx23 < 0.0 ) then
               call mix_gridpoint(Jstar, k)
               exit
            end if
         end do
      end do
   end do
end subroutine smooth_spikes




! ------------------------------------------------------------------------------
!  CHECKS2
!   Check nucleosynthesis variables for small an negative values
!   Modifies DHnuc for *1 and *2
! ------------------------------------------------------------------------------
subroutine checks2
   use real_kind
   use nucleosynthesis
   use settings
   implicit none

   ! Composition variables - should not be negative!
   where (Hnuc + DHnuc < 0.0d0)
      DHnuc = -Hnuc
   end where

   ! Abundances cannot be > 1
   where (Hnuc + DHnuc > 1.0d0)
      DHnuc = 1.0d0 - Hnuc
   end where
end subroutine checks2




! ------------------------------------------------------------------------------
!  CHECKS3
!   Check nucleosynthesis variables for conservation of all isotopes.
!   Modifies DHnuc for *1 and *2
! ------------------------------------------------------------------------------
subroutine checks3
   use real_kind
   use nucleosynthesis
   use settings
   implicit none
   integer :: ik, ij, Jstar
   real(double) :: avgdh

   ! Make sure abundances are conserved
   do Jstar = 1 , ktw
      do ik = 1, kh_nuc
         ! Find average change (aught to be 0)
         avgdh = sum(DHnuc(Jstar, 1:50, ik))

         ! Spread out error equally over all isotopes that are to be updated
         do ij = 1, nvar_nuc
            DHnuc(Jstar, ij, ik) = DHnuc(Jstar, ij, ik) - avgdh * (Hnuc(Jstar, ij, ik) + DHnuc(Jstar, ij, ik))
         end do
      end do
   end do
end subroutine checks3




! ------------------------------------------------------------------------------
!  UPDATE2
!   Update nucleosynthesis variables
! ------------------------------------------------------------------------------
subroutine update2
   use real_kind
   use nucleosynthesis
   
   implicit none
   if (allocated(Hnuc)) then
      Hnucpr(:,:,:) = Hnuc(:,:,:)
      Hnuc(:,:,:) = Hnuc(:,:,:) + DHnuc(:,:,:)
   end if
end subroutine update2




! ------------------------------------------------------------------------------
!  BACKUP2
!   Backup nucleosynthesis variables
! ------------------------------------------------------------------------------
subroutine backup2
   use real_kind
   use nucleosynthesis
   
   implicit none
   if (allocated(Hnuc)) then
      Hnuc(:,:,:) = Hnucpr(:,:,:)
      DHnuc(:,:,:) = 0.0d0
   end if
end subroutine backup2

