!> \file plotvariables.f90  Store plot variables for plt or mdl output


!> \brief Store plot variables for plt or mdl output
!!
!! This module can be approached by FUNCS1, the subroutine were all
!! physics is calculated, and by the various subroutines in PRINTB,
!! were a selected set of physcial variables are written to output
!! files such as file.plt1 etcetera
!! 
!! If you (temporarily) need any additional output in the plt or mdl
!! files, which is available in FUNCS1, this is the module you're
!! looking for. 
!! 
!! Hist:
!! SdM 30jan09
module plotvariables
  use real_kind
  
  implicit none
  
  ! Variables previously stored in the common block PLOT:
  real(double) :: wml, mtr, dhdt, dhsp(2), dhmb(2), dhso(2), dhml, dhgw, dhmt(2), qcnv, Teff 
  
end module plotvariables


