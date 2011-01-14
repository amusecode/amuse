module ionic

  ! This module contains the interface between the doric module and the
  ! mpi_hydro program. It contains two subroutines:
  ! - init_ionic: initializes all the radiation variables, including
  !         the doric initialization
  ! - update_ionic: updates ionic quantities

  ! Version: dummy routines, no atomic physics applied

  use precision, only: dp

  implicit none

contains

  subroutine init_ionic(restart,r_interface)

    ! This subroutine initializes the atomic/ionic fractions
    ! Garrelt Mellema

    ! Version: dummy

    logical,intent(in) :: restart
    real(kind=dp),intent(in) :: r_interface

  end subroutine init_ionic

  subroutine rad_evolve3d(dt)

    ! updates the radiation connected quantities one timestep
    ! garrelt mellema
    
    ! Version: dummy for runs without atomic physics

    real(kind=dp),intent(in) :: dt
    
  end subroutine rad_evolve3d
    
end module ionic


