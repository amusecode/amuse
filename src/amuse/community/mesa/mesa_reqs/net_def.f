! ***********************************************************************
!
!   Copyright (C) 2006, 2007, 2008, 2009  Bill Paxton, Frank Timmes
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

      module net_def
      
      use rates_def, only: rates_reaction_id_max
      use utils_def, only: integer_dict
      
      implicit none
      
      
      ! nets
      
         ! the 'basic net' is for standard hydrogen and helium burning.
         
            ! it includes 3 PP chains, 4 CNO cycles, triple alpha, and alpha captures up to mg24.      
            ! the 8 isotopes included are h1, he3, he4, c12, n14, o16, ne20, and mg24.
            ! see basic_net.dek for the details
               
         ! to go beyond the basics, use any combination of the following extensions.
         
            ! to add C13, see add_c13.dek
            ! to add Ne22, see add_ne22.dek
            ! to add O18 and Ne22, see add_o18_and_ne22.dek
            ! for more PP reactions, see add_pp_extras.dek
            ! for more CNO reactions including hot CNO, see add_cno_extras.dek
            ! for carbon/oxygen burning, see add_co_burn.dek
            ! for links of alpha chain to ni56, see add_alpha_s32.dek, etc.
            
         ! there are also several special purpose nets you can use.
            
            ! 'net_18_to_mg24' includes 18 isotopes from h1 to mg24.
            ! it gives you all the energy generation you need
            ! for hydrogen and helium burning, so for stars with less 
            ! than 6-8 Msun this is good to get you to the WD stage.
            ! includes basic PP chains, CNO cycles, and helium burning.
            
            ! 'net_16_to_ni56' uses only the basic 8 isotopes plus
            ! 8 more "multi-alpha" isotopes from si28 to ni56
            ! to cover evolution of massive stars from zams 
            ! through advance burning by alpha captures up to ni56.
            
            ! 'net_18_to_ni56' adds o18 and ne22 to net_16_to_ni56.
      


   ! reactions
   
         integer, parameter :: maxlen_reaction_Info = 72
         character (len=maxlen_reaction_Info) :: reaction_Info(rates_reaction_id_max)
         
         double precision, pointer :: std_reaction_Qs(:) ! (rates_reaction_id_max) 
            ! set at initialization; read-only afterwards.
            ! avg energy including neutrinos
         
         double precision, pointer :: std_reaction_neuQs(:) ! (rates_reaction_id_max) 
            ! set at initialization; read-only afterwards.
            ! avg neutrino loss
            
         integer, target :: reaction_screening_info(2,rates_reaction_id_max)
            ! reaction_screen_info(1:2,i) = [chem_id1, chem_id2] for screening.  0's if no screening.

         integer, target :: reaction_ye_rho_exponents(2,rates_reaction_id_max)
            ! multiply T dependent rate by Ye^a(i) * Rho^b(i)
            ! reaction_ye_rho_coeffs(1,i) is a(i)
            ! reaction_ye_rho_coeffs(2,i) is b(i)
            ! (0,0) for photodisintegrations and decays
            ! (0,1) for standard 2 body reactions
            ! (0,2) for 3 body reactions such as ir3a
            ! (1,1) for 2 body, including an electron such as irbe7ec
            ! (1,2) for 3 body, including an electron such as irpep
      
      
         integer, parameter :: max_num_reaction_inputs = 3
         integer, target :: reaction_inputs(2*max_num_reaction_inputs,rates_reaction_id_max)
            ! up to 5 pairs of coefficients and chem id's, terminated by 0's.
            ! e.g.,  o16(p,g)f17 would be (/ 1, io16, 1, ih1, 0 /)
            ! triple alpha would be (/ 3, ihe4, 0 /)
            ! he3(he4, g)be7(e-,nu)li7(p,a)he4 would be (/ 1, ihe3, 1, ihe4, i, ih1, 0 /)

         integer, parameter :: max_num_reaction_outputs = 4
         integer, target :: reaction_outputs(2*max_num_reaction_outputs,rates_reaction_id_max)
            ! up to 5 pairs of coefficients and chem id's, terminated by 0's.
            ! e.g.,  o16(p,g)f17 would be (/ 1, if17, 0 /)
            ! c12(a, p)n15 would be (/ 1, in15, 1, ih1, 0 /)
         


      
      ! predefined nets and extensions      
         
         integer, parameter :: bit_for_Basic = 0  ! net number bits start at 0
         integer, parameter :: num_Basic = 2**bit_for_Basic
         integer, parameter :: num_isos_for_Basic = 8
            
         ! bits 1 to 29 are for extensions to the basic net -- only some numbers are used.
         
         integer, parameter :: bit_for_PP_extras = 1  
         integer, parameter :: bit_for_CNO_extras = 2  
         integer, parameter :: bit_for_C13 = 3  
         integer, parameter :: bit_for_O18_and_Ne22 = 4  
         integer, parameter :: bit_for_CO_burn = 5  

         integer, parameter :: bit_for_alpha_to_S32 = 10
         integer, parameter :: bit_for_alpha_to_Ar36 = 11
         integer, parameter :: bit_for_alpha_to_Ca40 = 12
         integer, parameter :: bit_for_alpha_to_Ti44 = 13 
         integer, parameter :: bit_for_alpha_to_Cr48 = 14  
         integer, parameter :: bit_for_alpha_to_Fe52 = 15
         integer, parameter :: bit_for_alpha_to_Ni56 = 16

         integer, parameter :: bit_for_add_Na23_to_alpha_Mg24 = 17
         integer, parameter :: bit_for_add_Al27_to_alpha_Si28 = 18
         integer, parameter :: bit_for_add_P31_to_alpha_S32 = 19
         integer, parameter :: bit_for_add_Cl35_to_alpha_Ar36 = 20
         integer, parameter :: bit_for_add_K39_to_alpha_Ca40 = 21
         integer, parameter :: bit_for_add_Sc43_to_alpha_Ti44 = 22
         integer, parameter :: bit_for_add_V47_to_alpha_Cr48 = 23
         integer, parameter :: bit_for_add_Mn51_to_alpha_Fe52 = 24
         ! note: Co55 is automatically included in alpha_Ni56

         integer, parameter :: bit_for_specific_nets = 30 
         
         ! "special nets" are nets that are not combos of the basic net plus extensions
         ! NOTE: all nets are required to include the basic net isos.
         ! beyond that, the special nets can do anything they want.
         ! 2**bit_for_specific_nets = 2**30 = 1073741824
         
         integer, parameter :: num_19_to_Ni56 = 2**bit_for_specific_nets ! 1073741824
         integer, parameter :: num_18_to_Mg24 = 2**bit_for_specific_nets+1 ! 1073741825


      
   ! internal parameters for the implementation
      
      ! for tabular evaluation of the raw reaction rates
         double precision, parameter :: rattab_thi = 10.301029995664d0 ! log10(highest temp = 2e10)
         double precision, parameter :: rattab_tlo = 6d0 ! log10(lowest temp = 1e6)
         ! all reaction rates are set to zero for temperatures lower than 10**rattab_tlo
         
         integer, parameter :: nrattab = 8603 ! number of reaction rate table temperatures
            ! nrattab = 2000*(rattab_thi - rattab_tlo) + 1
            ! approx 2000 pts per decade of T
                  
         double precision, parameter :: rattab_tstp = (rattab_thi-rattab_tlo)/(nrattab-1)! step size
         
         
         
      type Net_General_Info ! things that are constant for the particular net
      ! it is okay to have multiple threads using the same instance of this simultaneously.

         integer :: num_isos ! total number in current net            
         integer :: num_reactions ! total number of reactions for current net

         ! isotopes
         integer, pointer :: net_iso(:) ! maps chem id to net iso number
         ! index from 1 to num_chem_isos
         ! value is 0 if the iso is not in the current net
         ! else is value between 1 and num_isos in current net
         integer, pointer :: chem_id(:) ! maps net iso number to chem id
         ! index from 1 to num_isos in current net
         ! value is between 1 and num_chem_isos         

         ! reactions

         integer :: which_rate_3a
         integer :: which_rate_c12ag
         integer :: which_rate_n14pg
         
         logical :: use_rates_fxt
                  
         integer, pointer :: net_reaction(:) ! maps reaction id to net reaction number
         ! index from 1 to num_net_reactions
         ! value is 0 if the reaction is not in the current net
         ! else is value between 1 and num_reactions in current net
         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         ! index from 1 to num_reactions in current net
         ! value is between 1 and num_net_reactions     


         ! mode switches for alpha chain reactions
         
         integer :: alpha_ap_mode ! for combination reactions like si28(a,p)p31(p,g)s32
            ! =0, off -- no (a,p)+(p,g) alpha capture links.
            ! =1, on without pa vs. pg factor.
            ! =2, on with pa vs. pg factor.
            
         integer :: alpha_gp_mode ! for combination reactions like s32(g,p)p31(p,a)si28
            ! =0, off -- no (g,p)+(p,a) alpha emission links.
            ! =1, on without pa vs. pg factor.
            ! =2, on with pa vs. pg factor.
         
         ! extra info
   
         ! the following is private info for the implementation
         
         ! tables for graboske screening
         double precision, pointer :: zg1(:) ! (num_reactions)
         double precision, pointer :: zg2(:) ! (num_reactions)
         double precision, pointer :: zg3(:) ! (num_reactions)
         double precision, pointer :: zg4(:) ! (num_reactions)
         
         ! tables for screen5
         double precision, pointer :: zs13(:) ! (num_reactions) ! zs13 = (z1+z2)**(1./3.)
         double precision, pointer :: zhat(:) ! (num_reactions)
         double precision, pointer :: zhat2(:) ! (num_reactions)
         double precision, pointer :: lzav(:) ! (num_reactions)
         double precision, pointer :: aznut(:) ! (num_reactions)
         double precision, pointer :: zs13inv(:) ! (num_reactions) ! zs13inv = 1 / zs13
   
         ! info for evaluation of the raw reaction rates
         double precision, pointer :: rattab(:,:,:) ! (num_rvs, num_reactions, nrattab)
         double precision, pointer :: ttab(:) ! (nrattab)
         double precision, pointer :: logttab(:) ! (nrattab)
         double precision, pointer :: rattab_f(:,:,:) ! (4, nrattab, num_reactions) ! for interpolation

         ! bookkeeping
         integer :: handle
         logical :: net_has_been_defined
         logical :: in_use

      end type Net_General_Info
               
               
      type Net_Info
         ! this is private working storage for the nuclear reaction calculations

         type (Net_General_Info), pointer  :: g

         double precision, pointer :: reaction_Qs(:) ! if null, use standard values         
         double precision, pointer :: reaction_neuQs(:) ! if null, use standard values

         ! molar fractions and their rates of change
         double precision, pointer :: y(:) ! units [moles/gram]     (num_isos)
         double precision, pointer :: dydt(:, :) ! units [moles/(gram-second)]  (num_rvs, num_isos)
         double precision, pointer :: d_dydt_dy(:, :) ! units [1/second] (num_isos, num_isos)

         logical :: use_graboske_et_al_screening
         double precision :: theta_e_for_graboske_et_al
         double precision, pointer :: graboske_cache(:, :, :)

         double precision, pointer :: rate_screened(:, :) ! (num_rvs, num_rates)
         ! the units here depend on the number of reactants.
         ! in all cases, the rate_screened times as many molar fractions as there are reactants
            ! gives a number with the same units as dy/dt.
         ! so for a 2-body reaction, there are 2 Y factors, each with units [moles/gram]
            ! and the rate_screened units for such a reaction are [grams/(mole-sec)], 
            ! which when multiplied by [moles/gram]^2 gives the same units as dydt.
         ! for a 1-body reaction (e.g., a decay), there is only 1 Y factor, so the units are [1/second].
         ! similarly, a 3 body reaction will have rate_screened with units of [gram^2/(mole^2-sec)].

         double precision, pointer :: rate_raw(:, :) ! (num_rvs, num_rates)
         ! raw rates are unscreened
         ! units are the same as rate_screened

         double precision, pointer :: reaction_eps_nuc(:, :) ! (num_rvs, num_rates)
         
         double precision, pointer :: d_eps_nuc_dy(:) ! (num_isos)

         double precision, pointer :: eps_nuc_categories(:, :) ! (num_rvs, num_categories)
         ! eps_nuc subtotals for each reaction category

         double precision :: eps_neu_total

      end type Net_Info
      

   ! private to the implementation
      integer, parameter :: max_net_handles = 1000
      type (Net_General_Info), target :: net_handles(max_net_handles)
      
      character (len=256) :: net_dir

      
      contains


      subroutine net_def_init(data_dir)
         use utils_lib, only: integer_dict_define, integer_dict_create_hash
         character (*), intent(in) :: data_dir
         integer :: i, ierr
         net_dir = trim(data_dir) // '/net_data'
         do i=1, max_net_handles
            net_handles(i)% handle = i
            net_handles(i)% in_use = .false.
            net_handles(i)% net_has_been_defined = .false.
            net_handles(i)% num_isos = 0
            net_handles(i)% num_reactions = 0
         end do
      end subroutine net_def_init


      integer function do_alloc_net(ierr)
         use rates_def
         use alert_lib, only:alert
         integer, intent(out) :: ierr
         integer :: i
         type (Net_General_Info), pointer :: g
         ierr = 0
         do_alloc_net = -1
!$omp critical (net_handle)
         do i = 1, max_net_handles
            if (.not. net_handles(i)% in_use) then
               net_handles(i)% in_use = .true.
               do_alloc_net = i
               exit
            end if
         end do
!$omp end critical (net_handle)
         if (do_alloc_net == -1) then
            ierr = -1
            call alert(ierr, 'no available net handle')
            return
         end if
         if (net_handles(do_alloc_net)% handle /= do_alloc_net) then
            ierr = -1
            call alert(ierr, 'broken handle for net')
            return
         end if
         g => net_handles(do_alloc_net)
         nullify(g% net_iso)
         nullify(g% chem_id)
         nullify(g% net_reaction)
         nullify(g% reaction_id)
         g% net_has_been_defined = .false.
         g% num_isos = 0
         g% num_reactions = 0
         g% which_rate_3a = use_rate_3a_NACRE
         g% which_rate_c12ag = use_rate_c12ag_NACRE
         g% which_rate_n14pg = use_rate_n14pg_NACRE
         g% use_rates_fxt = .false.
         g% alpha_ap_mode = 1 ! on without pa vs. pg factor
         g% alpha_gp_mode = 0 ! off -- no (g,p)+(p,a) alpha emission links
      end function do_alloc_net
      
      
      subroutine do_free_net(handle)
         use rates_def
         integer, intent(in) :: handle
         type (Net_General_Info), pointer :: g
         if (handle >= 1 .and. handle <= max_net_handles) then
            g => net_handles(handle)
            deallocate(g% net_iso)
            deallocate(g% chem_id)
            deallocate(g% net_reaction)
            deallocate(g% reaction_id)
            deallocate(g% zg1)
            deallocate(g% zg2)
            deallocate(g% zg3)
            deallocate(g% zg4)
            deallocate(g% zs13)
            deallocate(g% zhat)
            deallocate(g% zhat2)
            deallocate(g% lzav)
            deallocate(g% aznut)
            deallocate(g% zs13inv)
            deallocate(g% rattab)
            deallocate(g% ttab)
            deallocate(g% logttab)
            deallocate(g% rattab_f)
            g% in_use = .false.
            g% net_has_been_defined = .false.
            g% num_isos = 0
            g% num_reactions = 0
            g% which_rate_3a = use_rate_3a_NACRE ! default
            g% which_rate_c12ag = use_rate_c12ag_NACRE ! default
            g% which_rate_n14pg = use_rate_n14pg_NACRE ! default
         end if
      end subroutine do_free_net
      

      subroutine get_net_ptr(handle, g, ierr)
         use alert_lib, only:alert
         integer, intent(in) :: handle
         type (Net_General_Info), pointer :: g
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_net_handles) then
            ierr = -1
            call alert(ierr, 'invalid net handle')
            return
         end if
         g => net_handles(handle)
         ierr = 0
      end subroutine get_net_ptr
      

      end module net_def

