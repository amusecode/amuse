! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton
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

      module kap_def
      implicit none
      
      
      integer, parameter :: kap_table_fixed_metal_form = 1
      integer, parameter :: kap_table_co_enhanced_form = 2
      
      
      ! for fixed metal tables (no enhancements)
      type Kap_Z_Table ! holds pointers to all the X tables for a particular Z
         real :: Z ! the Z for this table
         integer :: num_Xs ! number of X's for this Z
         type (Kap_X_Table), dimension(:), pointer :: x_tables ! stored in order of increasing X
      end type Kap_Z_Table

      type (Kap_Z_Table), dimension(:), pointer :: kap_z_tables ! stored in order of increasing Z
      
      type Kap_X_Table
         logical :: not_loaded_yet
         real :: X
         real :: Z
         real :: logR_min
         real :: logR_max
         integer :: num_logRs
         integer :: ili_logRs ! =1 if logRs are evenly spaced
         real, dimension(:), pointer :: logRs ! indexed from 1 to num_logRs
         real :: logT_min
         real :: logT_max
         integer :: num_logTs
         integer :: ili_logTs ! =1 if logTs are evenly spaced
         real, dimension(:), pointer :: logTs ! indexed from 1 to num_logT
         real, dimension(:,:,:), pointer :: kap
      end type Kap_X_Table
      
      integer, parameter :: num_kap_Xs = 10
      ! 0,  .1, .2, .35, .5, .7, .8, .9, .95, 1-Z
      double precision, parameter :: kap_Xs(1:num_kap_Xs) =  &
            (/ 0.00d0, 0.10d0, 0.20d0, 0.35d0, 0.50d0, 0.70d0, 0.80d0, 0.90d0, 0.95d0, -1d0 /)
         ! -1 means X = 1-Z for whatever Z
      
      integer, parameter :: num_kap_Zs = 13
      double precision, parameter :: kap_Zs(1:num_kap_Zs) =  &
            (/ 0.000d0, 0.0001d0, 0.0003d0, 0.001d0, 0.002d0, 0.004d0, &
            0.01d0, 0.02d0, 0.03d0, 0.04d0, 0.06d0, 0.08d0, 0.100d0 /)
         
      integer, parameter :: num_kap_Xs_for_this_Z(1:num_kap_Zs) =  &
            (/ num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, &
            num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, num_kap_Xs, &
            num_kap_Xs-2, num_kap_Xs-2, num_kap_Xs-2 /)      
      
      ! for C/O enhanced tables
      type Kap_CO_Z_Table
         real :: Z ! the Z for this table
         real :: log10_Z ! log10(Z)
         real :: XC_base ! dXC = C excess = XC - XC_base, where XC is c12 mass fraction
         real :: XO_base ! dXO = O excess = XO - XO_base, where XO is O16 mass fraction
         type (Kap_CO_X_Table), dimension(:), pointer :: x_tables ! stored in order of increasing X
            ! the X tables need not be equally spaced
      end type Kap_CO_Z_Table

      type (Kap_CO_Z_Table), dimension(:), pointer :: kap_co_z_tables ! stored in order of increasing Z

      integer, parameter :: num_kap_CO_Zs = 7
      double precision, parameter :: kap_CO_Zs(1:num_kap_CO_Zs) =  &
         (/ 0.000, 0.001, 0.004, 0.010, 0.020, 0.030, 0.100 /)
      
      integer, parameter :: num_kap_CO_Xs = 4
      double precision, parameter :: kap_CO_Xs(1:num_kap_CO_Xs) =  &
         (/ 0.00d0, 0.10d0, 0.35d0, 0.70d0 /)
      
      integer, parameter :: num_kap_CO_dXs = 8
      double precision, parameter :: kap_CO_dXs(num_kap_CO_dXs) =  &
         (/ 0.00d0, 0.01d0, 0.03d0, 0.10d0, 0.20d0, 0.40d0, 0.60d0, 1.0d0 /)
      
      type Kap_CO_Table
         integer :: table_num ! the table number from the data file
         real :: X
         real :: Z
         real :: dXC
         real :: dXO
         real :: dXC_lookup
         real :: dXO_lookup
         real, dimension(:,:,:), pointer :: kap
      end type Kap_CO_Table
      
      
      integer, parameter :: max_num_CO_tables = 70
      
      
      
      type Kap_CO_X_Table
      
         logical :: not_loaded_yet
         real :: X
         real :: Z
         real :: logR_min
         real :: logR_max
         integer :: num_logRs
         integer :: ili_logRs ! =1 if logRs are evenly spaced
         real, dimension(:), pointer :: logRs ! indexed from 1 to num_logRs
         real :: logT_min
         real :: logT_max
         integer :: num_logTs
         integer :: ili_logTs ! =1 if logTs are evenly spaced
         real, dimension(:), pointer :: logTs ! indexed from 1 to num_logT
         
         integer :: num_CO_tables
         ! the tables are in 3 groups
         ! 1) tables with dXC > dXO, ordered by increasing dXO, and by increasing dXC within same dXO.
         ! 2) tables with dXC = dXO, ordered by increasing value.
         ! 3) tables with dXC < dXO, ordered by increasing dXC, and by increasing dXO within same dXC.
         ! the spacing of dXC's is the same as dXO's, so there are as many tables in 3) as in 1).
         integer :: num_dXC_gt_dXO ! the number of tables with dXC > dXO
         integer :: CO_table_numbers(num_kap_CO_dXs,num_kap_CO_dXs) 
            ! entry (i,j) is the co_index for table with dXC=Xs(i) and dXO=Xs(j), or -1 if no such table.
         integer :: next_dXO_table(max_num_CO_tables) 
            ! entry (i) is the co_index for the table with same dXC and next larger dXO, or -1 if none such.
         integer :: next_dXC_table(max_num_CO_tables) 
            ! entry (i) is the co_index for the table with same dXO and next larger dXC, or -1 if none such.
         type (Kap_CO_Table), dimension(:), pointer :: co_tables
          
      end type Kap_CO_X_Table



      

      type Kap_General_Info
      
         logical :: cubic_interpolation_in_X
         logical :: cubic_interpolation_in_Z

         ! for logR > 1, we extrapolate the radiative opacities.
         ! for low T, this can run into problems, so we need to clip logT when logR > 1.
         real :: min_logT_for_logR_gt_1 ! 3.3 is the default for this
         
         ! bookkeeping
         integer :: handle
         logical :: in_use
         
      end type Kap_General_Info



      ! NOTE: in the following, "log" means base 10, "ln" means natural log, and units are cgs.

      integer, parameter :: sz_per_kap_point = 4
      !
      ! function f(x,y) with samples f(i,j) has bicubic spline fit s(x,y).
      ! compact representation of spline fit uses 4 entries as follows:
      !
      ! d(1,i,j) = s(i,j)
      ! d(2,i,j) = d2s_dx2(i,j)
      ! d(3,i,j) = d2s_dy2(i,j)
      ! d(4,i,j) = d4s_dx2_dy2(i,j)
      ! 
      ! given f(i,j), the spline fitting code can compute the other entries
      !
      ! given d(1:4,i,j), spline interpolation code can compute s(x,y)
      ! and also the partials ds_dx(x,y) and ds_dy(x,y)
      !

      
      logical :: kap_is_initialized = .false.
      
      
      ! info for on-demand loading
      character (len=256) :: kap_dir
      character (len=256) :: kap_prefix = 'gn93'
      logical :: kap_use_cache

      
      logical :: clip_to_kap_table_boundaries ! typically, this should be set true.
         ! if this is set true, then temperature and density args are
         ! clipped to the boundaries of the table.
         ! if this is false and an arg is off the table, an alert is raised.
      
      integer, parameter :: max_kap_handles = 1000
      type (Kap_General_Info), target :: kap_handles(max_kap_handles)
      

      contains
      
      
      subroutine kap_def_init
         use chem_def
         integer :: i
         clip_to_kap_table_boundaries = .true.
         do i=1,max_kap_handles
            kap_handles(i)% cubic_interpolation_in_X = .true.
            kap_handles(i)% cubic_interpolation_in_Z = .false.            
            kap_handles(i)% min_logT_for_logR_gt_1 = 3.3
            kap_handles(i)% handle = i
            kap_handles(i)% in_use = .false.
         end do
      end subroutine kap_def_init

      
      integer function do_alloc_kap(ierr)
         use alert_lib,only:alert
         integer, intent(out) :: ierr
         integer :: i
         type (Kap_General_Info), pointer :: rq
         ierr = 0
         do_alloc_kap = -1
!$omp critical (kap_handle)
         do i = 1, max_kap_handles
            if (.not. kap_handles(i)% in_use) then
               kap_handles(i)% in_use = .true.
               do_alloc_kap = i
               exit
            end if
         end do
!$omp end critical (kap_handle)
         if (do_alloc_kap == -1) then
            ierr = -1
            call alert(ierr, 'no available kap handle')
            return
         end if
         if (kap_handles(do_alloc_kap)% handle /= do_alloc_kap) then
            ierr = -1
            call alert(ierr, 'broken handle for kap')
            return
         end if
         rq => kap_handles(do_alloc_kap)
      end function do_alloc_kap
      
      
      subroutine do_free_kap(handle)
         integer, intent(in) :: handle
         if (handle >= 1 .and. handle <= max_kap_handles) then
            kap_handles(handle)% in_use = .false.
         end if
      end subroutine do_free_kap
      

      subroutine get_kap_ptr(handle,rq,ierr)
         use alert_lib,only:alert
         integer, intent(in) :: handle
         type (Kap_General_Info), pointer :: rq
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_kap_handles) then
            ierr = -1
            call alert(ierr,'invalid kap handle')
            return
         end if
         rq => kap_handles(handle)
         ierr = 0
      end subroutine get_kap_ptr
      

      end module kap_def
      
