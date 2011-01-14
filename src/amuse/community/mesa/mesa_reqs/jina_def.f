! ***********************************************************************
!
!   Copyright (C) 2010  Ed Brown, Bill Paxton
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
 
      module jina_def

      use screen_def, only: Screen_Info
      use chem_def, only: iso_name_length, nuclide_data, npart
      use utils_def, only: integer_dict
      
      implicit none




      integer, parameter :: max_nreaclib=80000
      integer, parameter :: max_species_per_reaction=6
      integer, parameter :: ncoefficients=7
      integer, parameter :: nchapters=11
      integer, parameter :: ninverse_coeff = 2
      integer, parameter :: max_terms_per_rate = 20
      integer, parameter :: max_id_length = 36


      ! i/o parameters
      integer, parameter :: internal_format = 0
      integer, parameter :: pretty_print_format = 1


      ! flags for reaction types -- values are chapter numbers
      integer, parameter :: &
         r_one_one = 1, &
         r_one_two = 2, &
         r_one_three = 3, &
         r_two_one = 4, &
         r_two_two = 5, &
         r_two_three = 6, &
         r_two_four = 7, &
         r_three_one = 8, &
         r_three_two = 9, &
         r_four_two = 10, &
         r_one_four = 11


      integer, dimension(nchapters) :: Nin = (/1,1,1,2,2,2,2,3,3,4,1/)
      integer, dimension(nchapters) :: Nout = (/1,2,3,1,2,3,4,1,2,2,4/)


      type reaclib_data
      	integer,dimension(:),pointer :: chapter
      	character(len=iso_name_length),dimension(:,:),pointer :: species
      	character(len=iso_name_length),dimension(:),pointer :: label
      	character,dimension(:),pointer :: reaction_flag
      	character,dimension(:),pointer :: reverse_flag
      	double precision,dimension(:),pointer :: Qvalue
      	double precision,dimension(:,:),pointer :: coefficients
      end type reaclib_data


      type reaction_data
      	integer :: nreactions
      	integer :: nchapters_included
      	integer, dimension(nchapters) :: chapters_present
      	
      	integer :: num_from_weaklib
      	integer, dimension(:), pointer :: weaklib_ids ! (num_from_weaklib)
      	
      	integer, dimension(2,nchapters) :: bookmarks
      	character(len=max_id_length), dimension(:), pointer :: reaction_handle
      	integer, dimension(:), pointer :: category
      	integer, dimension(:), pointer :: chapter
      	integer, dimension(:,:), pointer :: pspecies
      	character, dimension(:), pointer :: reaction_flag
      	double precision, dimension(:), pointer :: weight
      	double precision, dimension(:,:), pointer :: coefficients
      	double precision, dimension(:), pointer :: weak_mask
      	double precision, dimension(:,:), pointer :: inverse_coefficients
      	integer, dimension(:), pointer :: inverse_exp
      	double precision, dimension(:,:), pointer :: inverse_part
         double precision, dimension(:), pointer :: Q
         double precision, dimension(:), pointer :: Qneu
         type (integer_dict), pointer :: reaction_dict ! keys are reaction_handles
         ! low T rates
         integer :: ipp, ipep, idpg, ihep, ihe3he3, ihe3he4, ibe7pg, ibe7em, ili7pag, ib8ep
         
         ! optional extra factors for various rates
         integer :: num_forward_rate_factors
         integer :: num_reverse_rate_factors
         integer, pointer :: forward_rate_num(:) ! (num_forward_rate_factors)
         integer, pointer :: reverse_rate_num(:) ! (num_reverse_rate_factors)
            ! value in range 1 to nreactions
         double precision, pointer :: forward_rate_factor(:) ! (num_forward_rate_factors)
         double precision, pointer :: reverse_rate_factor(:) ! (num_reverse_rate_factors)
         
      end type reaction_data


      ! container to hold locations of all terms for a given rate
      type rate_location
      	character(len=max_id_length) :: reaction_handle
      	integer :: nterms
      	integer, dimension(max_terms_per_rate) :: indices
      end type rate_location
      
      
      ! matrix solver options
      integer, parameter :: lapack = 1
      integer, parameter :: sparskit = 2
      

      type Jina_Info ! this can be shared by multiple threads
      
         ! information about the net      
            type(nuclide_data) :: nuclides
            type (integer_dict), pointer :: nuclides_dict ! maps nuclides% names(i) to i
            type(reaction_data) :: rates       
            integer, pointer :: chem_id(:)   
               ! chem_id indexed by net ios number, from 1 to g% nuclides% nnuclides
               ! value from 1 to num_chem_isos
            integer, pointer :: net_iso(:)   
               ! indexed by chem_id, from 1 to num_chem_isos
               ! value from 1 to g% nuclides% nnuclides if in current net
               ! value <= 0 if not in current net
            ! g% net_iso(g% chem_id(i)) == i, forall (i=1:g% nuclides% nnuclides)
         
         ! info for 1-zone burner            
            ! error tolerances
            double precision :: rtol, atol        
            ! which matrix solver
            integer :: which_decsol ! e.g., lapack or sparskit
            ! information about the jacobian matrix
            integer :: ijac, nzmax, isparse, mljac, mujac
            ! information about the "mass" matrix
            integer :: imas, mlmas, mumas

         ! matrix solver parameter arrays
            integer :: lrd, lid
            integer :: im ! for sparskit

         ! work arrays.
            integer :: lwork, liwork

         ! local parameter arrays.
            integer :: lrpar, lipar

         ! screening
            logical :: use_graboske_et_al_screening
            double precision, dimension(:), pointer :: z ! (nnuc)
            double precision, dimension(:), pointer :: zg1, zg2, zg3, zg4 ! (num_reactions)
            double precision, dimension(:), pointer :: zs13, zhat, zhat2, lzav, aznut, zs13inv ! (num_reactions)

         ! bookkeeping        
            integer :: handle
            logical :: in_use
         
      end type Jina_Info
      

      type Jina_Private_Info ! this is private to a single thread
      
         type (Jina_Info), pointer :: g
         
         ! rates
         	double precision, dimension(:), pointer :: ln_lambda  ! (nr) log of forward rates
         	double precision, dimension(:), pointer :: lambda  ! (nr) forward rates
         	double precision, dimension(:), pointer :: dlambda_dlnT  ! (nr)
         	double precision, dimension(:), pointer :: dlambda_dlnRho  ! (nr)
         	double precision, dimension(:), pointer :: rlambda ! (nr) reverse rates
         	double precision, dimension(:), pointer :: drlambda_dlnT ! (nr)
         	double precision, dimension(:), pointer :: drlambda_dlnRho! (nr)
         	
         ! Q values
         	double precision, dimension(:), pointer :: Q ! (nr)
         	double precision, dimension(:), pointer :: dQ_dlnT ! (nr)
         	double precision, dimension(:), pointer :: dQ_dlnRho ! (nr)
         	double precision, dimension(:), pointer :: Qneu ! (nr)
         	double precision, dimension(:), pointer :: dQneu_dlnT ! (nr)
         	double precision, dimension(:), pointer :: dQneu_dlnRho ! (nr)

         ! matrix solver parameter arrays
            double precision, pointer :: rpar_decsol(:) ! (lrd)
            integer, pointer :: ipar_decsol(:) ! (lid)

         ! work arrays.
            double precision, pointer :: work(:) ! (lwork)
            integer, pointer :: iwork(:) ! (liwork)

         ! local parameter arrays.
            double precision, pointer :: rpar(:)
            integer, pointer :: ipar(:)
         
         ! screening data
            type (Screen_Info), pointer :: sc
            double precision, pointer :: graboske_cache(:,:,:) ! (3,max_z_to_cache,max_z_to_cache)

      end type Jina_Private_Info
      

      integer, parameter :: i_caller_id = 1
      integer, parameter :: i_handle = 2
      integer, parameter :: i_sparse_format = 3
      integer, parameter :: burn_lipar = 3

      integer, parameter :: burn_lrpar = 0


      ! private to the implementation
      integer, parameter :: max_jina_burn_handles = 1000 ! make this as large as you want
      type (Jina_Info), dimension(:), target :: jina_burn_handles(max_jina_burn_handles)
      
      character (len=256) :: jina_dir
      
   	type(reaclib_data) :: reaclib
   	integer :: nreaclib
   	

      contains
      
      
      subroutine jina_def_init(data_dir)
         character (*), intent(in) :: data_dir
         integer :: i, ierr
         jina_dir = trim(data_dir) // '/jina_data'
         do i=1, max_jina_burn_handles
            jina_burn_handles(i)% handle = i
            jina_burn_handles(i)% in_use = .false.
         end do
      end subroutine jina_def_init
      
      
      integer function do_alloc_jina_burn(ierr)
         use alert_lib, only:alert
         integer, intent(out) :: ierr
         integer :: i, new_num
         ierr = 0
         do_alloc_jina_burn = -1
!$omp critical (jina_handle)
         do i = 1, max_jina_burn_handles
            if (.not. jina_burn_handles(i)% in_use) then
               jina_burn_handles(i)% in_use = .true.
               do_alloc_jina_burn = i
               exit
            end if
         end do
!$omp end critical (jina_handle)
         if (do_alloc_jina_burn == -1) then
            ierr = -1
            call alert(ierr, 'no available jina burn handle')
            return
         end if
         if (jina_burn_handles(do_alloc_jina_burn)% handle /= do_alloc_jina_burn) then
            ierr = -1
            call alert(ierr, 'broken handle for jina burn')
            return
         end if
      end function do_alloc_jina_burn
      
      

      subroutine get_jina_burner_ptr(handle, g, ierr)
         integer, intent(in) :: handle
         type (Jina_Info), pointer :: g
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_jina_burn_handles) then
            ierr = -1
            write(*,*) 'invalid handle -- did you call alloc_jina_net_handle?'
            return
         end if
         g => jina_burn_handles(handle)
         if (.not. g% in_use) then
            ierr = -1
            write(*,*) 'invalid handle -- did you call alloc_jina_net_handle?'
            return
         end if
         ierr = 0
      end subroutine get_jina_burner_ptr
      
      
      subroutine do_free_jina_burner(handle) ! contents have already been deallocated
         use chem_def, only: free_nuclide_data
         integer, intent(in) :: handle
         type (Jina_Info), pointer :: g
         integer :: ierr
         ierr = 0
         if (handle >= 1 .and. handle <= max_jina_burn_handles) then
            g => jina_burn_handles(handle)
            g% in_use = .false.
         end if
      end subroutine do_free_jina_burner


      subroutine allocate_reaclib_data(r, n, ierr)
      	type(reaclib_data), intent(inout) :: r
      	integer, intent(in) :: n
      	integer, intent(out) :: ierr
      	ierr = 0
      	allocate( &
      	   r% chapter(n),r% species(max_species_per_reaction,n),r% label(n),r% reaction_flag(n), &
      		r% reverse_flag(n),r% Qvalue(n),r% coefficients(ncoefficients,n), stat=ierr)
      end subroutine allocate_reaclib_data
      

      subroutine free_reaclib_data(reaclib, ierr)
      	use alert_lib
      	type(reaclib_data), intent(inout) :: reaclib
      	integer, intent(out) :: ierr
      	ierr = 0
      	if (associated(reaclib% chapter)) & 
      		deallocate( &
      		   reaclib% chapter,reaclib% species,reaclib% label,reaclib% reaction_flag, &
      		   reaclib% reverse_flag,reaclib% Qvalue, reaclib% coefficients, stat=ierr)
      	if (ierr /= 0) call alert(ierr,"free_reaclib: Deallocation request denied")
      end subroutine free_reaclib_data
      

      subroutine allocate_reaction_data(r,n,nweak,ierr)
      	type(reaction_data), intent(out) :: r
      	integer, intent(in) :: n ! number of rates
      	integer, intent(in) :: nweak ! number of weaklib rates
      	integer, intent(out) :: ierr
      	allocate( &
      	   r% reaction_handle(n),r% category(n),r% chapter(n),&
      	   r% weaklib_ids(nweak),r% pspecies(max_species_per_reaction,n),r% reaction_flag(n), &
      		r% weight(n),r% coefficients(ncoefficients,n), &
      		r% weak_mask(n),r% inverse_coefficients(ninverse_coeff,n), &
      		r% inverse_exp(n),r% inverse_part(npart,n),r% Q(n),r% Qneu(n),stat=ierr)
      	nullify(r% reaction_dict, r% forward_rate_num, r% reverse_rate_num, &
      	   r% forward_rate_factor, r% reverse_rate_factor)
      	r% num_forward_rate_factors = 0
      	r% num_reverse_rate_factors = 0
      end subroutine allocate_reaction_data
      

      subroutine free_reaction_data(r)
         use utils_lib, only: integer_dict_free
      	type(reaction_data), intent(inout) :: r
      	if (associated(r% chapter)) then
      		deallocate(r% reaction_handle,r% category,r% chapter,&
      		   r% weaklib_ids,r% pspecies, r% reaction_flag, r% weight, r% coefficients, r% weak_mask, &
      		   r% inverse_coefficients, r% inverse_exp, r% inverse_part, r% Q, r% Qneu)
      		if (associated(r% reaction_dict)) call integer_dict_free(r% reaction_dict)
      		if (associated(r% forward_rate_num)) deallocate(r% forward_rate_num)
      		if (associated(r% reverse_rate_num)) deallocate(r% reverse_rate_num)
      		if (associated(r% forward_rate_factor)) deallocate(r% forward_rate_factor)
      		if (associated(r% reverse_rate_factor)) deallocate(r% reverse_rate_factor)
      	end if
      end subroutine free_reaction_data
      

      end module jina_def
