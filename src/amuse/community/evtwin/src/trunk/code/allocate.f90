module allocate_arrays

contains

subroutine allocate_global_arrays(max_mesh)
   use mesh
   use mesh_enc
   use nucleosynthesis
   use nucleosynthesis_neutron_ratios
   use explicit_functions
   use structure_variables
   use semi_implicit_variables
   implicit none
   integer, intent(in) :: max_mesh

   nm = max_mesh

   allocate(  h(nvar, nm))
   allocate( dh(nvar, nm))
   allocate(hpr(nvar, nm))

   allocate(mutant_h(nvar, nm))

   allocate(ht(2, ht_nvar, nm))
   allocate(frac(nvar_nuc, nm))

   allocate(expl_var(NM, num_explv, 2))
   allocate(radacc(9, nm, 2))

   allocate(TH(nvar, nm))
   allocate(THb(nvar, nm))
   allocate(THc(nvar, nm))
   allocate(THd(nvar, nm))
   allocate(menc(2, nm))
   menc = 0.0d0

   allocate(sx(npx,NM+1))

   allocate(a_rot2(nm,2))
   allocate(a_tide(nm, 2:4, 2))
   allocate(r0_over_rv(nm, 2))
   allocate(fp_star(nm, 2))
   allocate(ft_star(nm, 2))
   allocate(Vmc2_star(nm, 2))
   allocate(diff_omega(nm, 2))

end subroutine allocate_global_arrays

end module allocate_arrays
