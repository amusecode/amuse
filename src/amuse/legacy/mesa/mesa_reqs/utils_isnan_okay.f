! ***********************************************************************
!
!   Copyright (C) 2007  Bill Paxton
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

      module utils_isnan
      implicit none
      

      contains
      
      
      logical function is_inf(x)
         double precision, intent(in) :: x
         is_inf = (2*x==x .and. x /= 0)
      end function is_inf
      
      
      logical function is_real_inf(x)
         real, intent(in) :: x
         is_real_inf = (2*x==x .and. x /= 0)
      end function is_real_inf

      
      logical function check_for_bad_num(x)
         double precision, intent(in) :: x
!         check_for_bad_num = isnan(x) .or. is_inf(x)
         check_for_bad_num = is_inf(x)
      end function check_for_bad_num

      
      logical function check_for_bad_real(x)
         real, intent(in) :: x
!         check_for_bad_real = isnan(x) .or. is_real_inf(x)
         check_for_bad_real = is_real_inf(x)
      end function check_for_bad_real


      end module utils_isnan

