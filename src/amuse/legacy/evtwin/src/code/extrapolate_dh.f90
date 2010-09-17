! Estimate DH(:) by quadratic extrapolation (as opposed to linear
! extrapolation)
module extrapolate_dh
   use real_kind
   use mesh

   ! Number of models to store; after this number of models are stored the
   ! storage space will wrap around. We need 4 instread of 3 because the code
   ! may step back one timestep.
   integer, parameter :: num_stored_models = 4
   integer :: last_stored_model = 0

   ! We need to store the total number of times a model has been placed in the
   ! buffer, because if this is < 3 we don't have enough models to compute the
   ! quadratic extarpolation
   integer :: buffer_calls

   ! We can go back one model, but not more than that!
   logical :: can_backup_one_model

   ! Local copies of some variables, for convenience
   integer :: number_of_meshpoints        ! KH
   integer :: number_of_variables         ! KVB

   ! Previous values of H and the corresponding ages
   real(double), allocatable :: previous_h(:, :, :)
   real(double), allocatable :: previous_age(:)
   real(double), allocatable :: previous_timestep(:)

contains

   ! initlse_parabola_storage_space:
   ! Initialise the storage space for the various models. This function only
   ! allocates enough space for the number of variables and the number of
   ! meshpoints that are passed to it.
   ! Input:
   !  KH - The number of meshpoints to be stored
   !  KVB - The total number of variables that need to be stored
   ! The variable permutation array needs to be set up prior to calling this
   ! function
   subroutine initlse_parabola_storage_space(kh, kvb)
      use real_kind
      use settings
      
      implicit none
      integer, intent(in) :: kh, kvb

      if (allocated(previous_h)) deallocate(previous_h);
      if (allocated(previous_age)) deallocate(previous_age);
      if (allocated(previous_timestep)) deallocate(previous_timestep);

      number_of_meshpoints = kh
      number_of_variables = kvb

      ! Allocate storage space for the various models
      allocate (previous_h(0:num_stored_models, number_of_meshpoints, number_of_variables))
      allocate (previous_age(0:num_stored_models))
      allocate (previous_timestep(0:num_stored_models))
      can_backup_one_model = .false.
      buffer_calls = 0
   end subroutine initlse_parabola_storage_space

   ! store_current_model:
   ! Store current model (H(,), AGE, DTY) in the buffer, advance buffer
   ! pointer
   subroutine store_current_model
      use real_kind
      use mesh
      use constants
      use test_variables, only: dt, age
      
      implicit none
      integer :: n, nv

      ! Unnamed COMMON block (we only need H and KP_VAR)
      integer :: ke1, ke2, kev, kp_var(40)

      ke1            = id(1)
      ke2            = id(2)
      kev            = id(5)
      kp_var(1:40)   = id(11:50)

      if (allocated(previous_h) .and. allocated(previous_age)) then
         last_stored_model = mod(last_stored_model+1, num_stored_models)
         previous_age(last_stored_model) = age
         previous_timestep(last_stored_model) = dt/csy
         do n=1, ke1+ke2+kev
            nv = kp_var(n)
            previous_h(last_stored_model, 1:kh, n) = h(nv, 1:kh)
         end do
      end if
      buffer_calls = buffer_calls+1
      can_backup_one_model = .true.
   end subroutine store_current_model

   ! restore_previous_model:
   ! Mark the previous model in the buffer as the current model, but only if
   ! a new model was stored in the buffer after the last call to this
   ! function. Needed when the solver backtracks to a previous model.
   subroutine restore_previous_model
      use real_kind
      
      implicit none

      if (.not. can_backup_one_model) return
      can_backup_one_model = .false.

      if (allocated(previous_h) .and. allocated(previous_age)) then
         last_stored_model = mod(last_stored_model+num_stored_models-1, num_stored_models)
         buffer_calls = buffer_calls-1
      end if
      buffer_calls = 0
   end subroutine restore_previous_model



   ! predict_dh
   ! Predict change DH(:,:) from the last model for a timestep DT
   ! This only works if there are at least 3 previous models in the buffer
   ! (otherwise there is not enough data to do the extrapolation).
   ! Input:
   !  DT - Timestep in years
   ! Output:
   !  DH() - Corrections to independent variables from quadratic
   !  extrapolation. Unchanged if there is not enough information in the
   !  buffer to compute this. We could have set this directly in the COMMON
   !  block instead, but this allows us some finer control.
   subroutine predict_dh(dt, new_dh)
      use real_kind
      use mesh
      
      implicit none
      real(double), intent(in) :: dt
      real(double), intent(inout) :: new_dh(nvar, nm)
      ! The unnamed COMMON block, needed to ket KD
      integer :: kp_var(40)
      ! The various times and time intervals from previous models
      real(double) :: tn, tn2, t, tp, tpp, dtp, dtpp, dttpp
      ! Determinant of Vandermonde matrix, matrix elements of inverse matrix
      real(double) :: detm
      real(double) :: invm(3,3) !, m(3,3)
      ! Model numbers of current, previous and previous previous models
      integer :: n, np, npp
      integer :: i, k, kv
      real(double) :: x, dx, rhs(3), lhs(3)

      ! Abort early if things are not yet defined
      if (.not. (allocated(previous_h) .or. allocated(previous_age))) then
         return
      end if
      ! Return if not enough models in buffer
      if (buffer_calls < 3) return;
      ! Determine current, previous and previous-previous model numbers in the
      ! array
      n = last_stored_model
      np = mod(n+num_stored_models-1, num_stored_models)
      npp = mod(n+num_stored_models-2, num_stored_models)

      ! Copy times to local variables
      !t = previous_age(n)
      !tp = previous_age(np)
      !tpp = previous_age(npp)
      !tn = t + dt

      ! Better: pretend current time is t=0, which elminates problems related to
      ! numerical accuracy when dt/t becomes of the order of the machine
      ! precision
      t = 0.0d0
      tp = -previous_timestep(n)
      tpp = -previous_timestep(n)-previous_timestep(np)
      tn = dt

      dtp = previous_timestep(n)
      dtpp = previous_timestep(np)
      dttpp = dtp+dtpp

      tn2 = tn**2

      ! Compute inverse Vandermondematrix
      !M(1, 1:3) = (/ 1.0d0, t, t**2 /)
      !M(2, 1:3) = (/ 1.0d0, tp, tp**2 /)
      !M(3, 1:3) = (/ 1.0d0, tpp, tpp**2 /)

      detm = -dtp*dtpp*dttpp
      invm(1, 1:3) = (/-dtpp *  tp*tpp,   dttpp *  tpp*t,  -dtp *  tp*t/)
      invm(2, 1:3) = (/ dtpp * (tp+tpp), -dttpp * (tpp+t),  dtp * (tp+t)/)
      invm(3, 1:3) = (/-dtpp,             dttpp,           -dtp/)
      invm(:,:) = invm(:,:)/detm
      !print '(3(1P,"[",3E16.8,"]",/))', matmul(M(:,:), invM(:,:))

      kp_var(1:40)   = id(11:50)

      ! Find predicted changes in variables
      do k=1, number_of_meshpoints
         do i=1, number_of_variables
            ! Compute extrapolation coefficients
            rhs(1:3) = (/ previous_h(n, k, i), previous_h(np, k, i), previous_h(npp, k, i) /)
            lhs(1:3) = matmul(invm(:,:), rhs(:))

            ! Compute extrapolated value of Hki for t+dt
            x = lhs(3)*tn2 + lhs(2)*tn + lhs(1)
            ! And dh
            dx = x - previous_h(n, k, i)
            !print '(X,1P,3I3,3E16.8)', k, i, KP_VAR(i), previous_h(np,k,i)-previous_h(npp,k,i),dx, DH(KP_VAR(i), k)
            kv = kp_var(i)
            ! Set initial guess for compositrion variables to 0.0, unless the change is
            ! large, indicating an active species. This reduces numerical noise but keeps
            ! the advantage of using quadratic extrapolations.
            if ( ((kv>8 .and. kv<12) .or. kv==16 .or. kv==5 .or. kv==3 .or. (kv >= 41 .and. kv <= 43)) .and. dabs(dx)<1.0d-4 )&
               dx = 0.0d0
            new_dh(kv, k) = dx
         end do
      end do


   end subroutine predict_dh

end module extrapolate_dh

