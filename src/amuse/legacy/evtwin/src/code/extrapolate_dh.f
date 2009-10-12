! Estimate DH(:) by quadratic extrapolation (as opposed to linear
! extrapolation)
      module extrapolate_dh
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
      double precision, allocatable :: previous_h(:, :, :)
      double precision, allocatable :: previous_age(:)
      double precision, allocatable :: previous_timestep(:)

      contains

! initialise_parabola_storage_space:
! Initialise the storage space for the various models. This function only
! allocates enough space for the number of variables and the number of
! meshpoints that are passed to it.
! Input:
!  KH - The number of meshpoints to be stored
!  KVB - The total number of variables that need to be stored
! The variable permutation array needs to be set up prior to calling this
! function
      subroutine initialise_parabola_storage_space(KH, KVB)
      use settings
      use mesh
      implicit none
      integer, intent(in) :: KH, KVB

      if (allocated(previous_h)) deallocate(previous_h);
      if (allocated(previous_age)) deallocate(previous_age);
      if (allocated(previous_timestep)) deallocate(previous_timestep);

      number_of_meshpoints = KH
      number_of_variables = KVB

! Allocate storage space for the various models
      allocate (previous_h(0:num_stored_models, number_of_meshpoints, number_of_variables))
      allocate (previous_age(0:num_stored_models))
      allocate (previous_timestep(0:num_stored_models))
      can_backup_one_model = .false.
      buffer_calls = 0
      end subroutine

! store_current_model:
! Store current model (H(,), AGE, DTY) in the buffer, advance buffer
! pointer
      subroutine store_current_model
      use mesh
      use constants
      implicit none
      integer :: n, nv

! Unnamed COMMON block (we only need H and KP_VAR)
      DOUBLE PRECISION :: H(NVAR,NM), COMMON_DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3),KP_VAR(40), KP_EQN(40), KP_BC(40), 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3),KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)
      COMMON H, COMMON_DH, EPS, DEL, DH0, KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2

! COMMON block /TVBLES/, we need DT (actually DTY...)
      double precision :: DT, ZQ(9), AGE, BM(71), PR(81), PPR(81)
      integer :: JHOLD, JM2, JM1
      COMMON /TVBLES/ DT, ZQ, AGE, BM, PR, PPR, JHOLD, JM2, JM1

      if (allocated(previous_h) .and. allocated(previous_age)) then
         last_stored_model = mod(last_stored_model+1, num_stored_models)
         previous_age(last_stored_model) = age
         previous_timestep(last_stored_model) = dt/CSY
         do n=1, KE1+KE2+KEV
            nv = KP_VAR(n)
            previous_h(last_stored_model, 1:KH, n) = H(nv, 1:KH)
         end do
      end if
      buffer_calls = buffer_calls+1
      can_backup_one_model = .true.
      end subroutine

! restore_previous_model:
! Mark the previous model in the buffer as the current model, but only if
! a new model was stored in the buffer after the last call to this
! function. Needed when the solver backtracks to a previous model.
      subroutine restore_previous_model
      implicit none

      if (.not. can_backup_one_model) return
      can_backup_one_model = .false.

      if (allocated(previous_h) .and. allocated(previous_age)) then
         last_stored_model = mod(last_stored_model+num_stored_models-1, num_stored_models)
         buffer_calls = buffer_calls-1
      end if
      buffer_calls = 0
      end subroutine

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
      subroutine predict_dh(dt, DH)
      use mesh
      implicit none
      double precision, intent(in) :: dt
      double precision, intent(inout) :: DH(NVAR, NM)
! The unnamed COMMON block, needed to ket KD
      DOUBLE PRECISION :: H(NVAR,NM), COMMON_DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3),KP_VAR(40), KP_EQN(40), KP_BC(40), 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3),KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)
      COMMON H, COMMON_DH, EPS, DEL, DH0, KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2
! The various times and time intervals from previous models
      double precision :: tn, tn2, t, tp, tpp, dtp, dtpp, dttpp
! Determinant of Vandermonde matrix, matrix elements of inverse matrix
      double precision :: detM
      double precision :: invM(3,3), M(3,3)
! Model numbers of current, previous and previous previous models
      integer :: n, np, npp
      integer :: i, k, kv
      double precision :: x, dx, rhs(3), lhs(3)

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

      detM = -dtp*dtpp*dttpp
      invM(1, 1:3) = (/-dtpp *  tp*tpp,   dttpp *  tpp*t,  -dtp *  tp*t/)
      invM(2, 1:3) = (/ dtpp * (tp+tpp), -dttpp * (tpp+t),  dtp * (tp+t)/)
      invM(3, 1:3) = (/-dtpp,             dttpp,           -dtp/)
      invM(:,:) = invM(:,:)/detM
      !print '(3(1P,"[",3E16.8,"]",/))', matmul(M(:,:), invM(:,:))

! Find predicted changes in variables
      do k=1, number_of_meshpoints
         do i=1, number_of_variables
! Compute extrapolation coefficients
            rhs(1:3) = (/ previous_h(n, k, i), previous_h(np, k, i), previous_h(npp, k, i) /)
            lhs(1:3) = matmul(invM(:,:), rhs(:))

! Compute extrapolated value of Hki for t+dt
            x = lhs(3)*tn2 + lhs(2)*tn + lhs(1)
! And dh
            dx = x - previous_h(n, k, i)
            !print '(X,1P,3I3,3E16.8)', k, i, KP_VAR(i), previous_h(np,k,i)-previous_h(npp,k,i),dx, DH(KP_VAR(i), k)
            kv = KP_VAR(i)
! Set initial guess for compositrion variables to 0.0, unless the change is 
! large, indicating an active species. This reduces numerical noise but keeps
! the advantage of using quadratic extrapolations.
            IF ( ((kv>8 .AND. kv<12) .or. kv==16 .OR. kv==5 .OR. kv==3) .and. dabs(dx)<1.0d-4 ) dx = 0.0d0
            DH(kv, k) = dx
         end do
      end do

      
      end subroutine

      end module
