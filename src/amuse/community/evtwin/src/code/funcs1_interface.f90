!> Because FUNCS1 has optional arguments (for output) we require an explicit
!! interface for it.
!! A better way would be to place funcs1 itself in a module, but the way the
!! code is currently setup, this leads to a circular dependence of
!! output_properties.f90 and stars_structure.f90, so we place the interface
!! specification in a separate source file.
!! 
!! \todo FIXME: correct this circular dependence!
!<
module funcs1_interface
   interface funcs1
      subroutine funcs1 ( jk, ji, var, dvar, fn1, eosout, abundout )
         use real_kind
         use mesh
         use eostate_types
         implicit none
         integer, intent(in) :: jk, ji
         real(double), intent(in) :: var(NVAR), dvar(NVAR)
         real(double), intent(out) :: fn1(NFUNC)
         type(eostate), optional, intent(out) :: eosout
         type(abundance), optional, intent(out) :: abundout
      end subroutine
   end interface funcs1
end module

