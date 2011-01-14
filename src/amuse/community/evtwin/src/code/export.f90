! Export a stellar model.
! This module is disgned to work along with the MUSE library and exports
! the data for a stellar model, to be used alongside a tool like
! MakeMeAStar.
! Order of variables stored at each meshpoint:
!  Mass coordinate [Msun], Radius [Rsun], log density [cgs],
!  log pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
! Mespoints will be passed out from *surface* to *centre*

function get_number_of_meshpoints()
  use real_kind
  use mesh
  
  implicit none
  integer :: get_number_of_meshpoints

  get_number_of_meshpoints = KH
end function get_number_of_meshpoints



subroutine export_stellar_model(jstar, model)
  use real_kind
  use mesh
  use atomic_data
  use extra_elements
  use structure_variables
  use binary_history, only: hpr
  
  implicit none
  integer, intent(in) :: jstar     ! Which component of a binary
  real(double), intent(out) :: model(13,kh)
  
  real(double) :: xa(9), na(9)
  real(double) :: avm
  
  real(double) :: hh(nvar,NM)

  integer :: ik, ikk, n
  
  
  ! Should we export the *current* or *previous* model?
  !> \todo FIXME: make this an option! At the moment we output the previous model.
  hh(:,:) = h(:,:)
  h(:,:) = hpr(:,:)

  call compute_output_quantities ( Jstar )

  do ik=1, kh
     ikk = 1 + (kh+1-ik)
     model(1, ik) = sx(9, ikk)
     model(2, ik) = sx(17, ikk)
     model(3, ik) = log(sx(3, ikk))
     model(4, ik) = log(sx(2, ikk))
     ! Convert *all* abundances to mass fractions
     xa(1) = h(5, ik)
     xa(2) = h(9, ik)
     xa(3) = h(10, ik)
     xa(4) = h(16, ik)
     xa(5) = h(3, ik)
     xa(6) = h(11, ik)
     xa(8) = h(NSi28, ik)
     xa(9) = h(NFe56, ik)
     xa(7) = 1.0d0 - sum(xa(1:6)) - sum(xa(8:9))
     do n=1, 9
        na(n) = xa(n) * can(n)/cbn(n)
     end do
     avm = sum(na(1:9))
     na(1:9) = na(1:9) / avm
     model(5:13, ik) = na(1:9)
  end do
  h(:,:) = hh(:,:)

end subroutine export_stellar_model

