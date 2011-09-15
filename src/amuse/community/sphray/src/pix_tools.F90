!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------
module pix_tools
use healpix_types
implicit none

INTEGER(KIND=i4b), private, PARAMETER :: ns_max=8192 ! 2^13 : largest nside available


interface fatal_error
   module procedure fatal_error_womsg, fatal_error_msg
end interface


contains


  subroutine fatal_error_msg (msg)
    character(len=*), intent(in) :: msg
    print *,'Fatal error: ', trim(msg)
    stop
  end subroutine fatal_error_msg
  
  subroutine fatal_error_womsg
    print *,'Fatal error'
    stop
  end subroutine fatal_error_womsg
  


  !=======================================================================
  subroutine pix2ang_ring(nside, ipix, theta, phi)
    !=======================================================================
    !     renders theta and phi coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)
    !     given the map resolution parameter nside
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN) :: ipix, nside
    REAL(KIND=DP), INTENT(OUT) ::  theta, phi

    INTEGER(KIND=I4B) ::  nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=DP) ::  fodd, hip, fihip
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error ("nside out of range")
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error ("ipix out of range")

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1

    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1*0.5_dp
       fihip = AINT ( hip , kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       theta = ACOS( 1.0_dp - iring**2 / (3.0_dp*nside**2) )
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       nl4   = 4*nside
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       theta = ACOS( (nl2 - iring) / (1.5_dp*nside) )
       phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip*0.5_dp
       fihip = AINT ( hip , kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       theta = ACOS( -1.0_dp + iring**2 / (3.0_dp*nside**2) )
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

    endif

    return
  end subroutine pix2ang_ring



  subroutine pix2vec_ring(nside, ipix, vector, vertex)
    !=======================================================================
    !     renders vector (x,y,z) coordinates of the nominal pixel center
    !     for the pixel number ipix (RING scheme)
    !     given the map resolution parameter nside
    !     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    !     in the order N,W,S,E
    !=======================================================================
    INTEGER(KIND=I4B), INTENT(IN)                             :: ipix, nside
    REAL(KIND=DP),     INTENT(OUT),dimension(1:)              :: vector
    REAL(KIND=DP),     INTENT(OUT),dimension(1:,1:), optional :: vertex

    INTEGER(KIND=I4B) :: nl2, nl4, npix, ncap, iring, iphi, ip, ipix1
    REAL(KIND=DP) ::  fact1, fact2, fodd, hip, fihip, z, sth, phi

    real(kind=DP) :: phi_nv, phi_wv, phi_sv, phi_ev
    real(kind=DP) :: z_nv, z_sv, sth_nv, sth_sv
    real(kind=DP) :: hdelta_phi
    integer(kind=I4B) :: iphi_mod, iphi_rat
    logical(kind=LGT) :: do_vertex
    !-----------------------------------------------------------------------
    if (nside<1 .or. nside>ns_max) call fatal_error("nside out of range")
    npix = 12*nside**2       ! total number of points
    if (ipix <0 .or. ipix>npix-1) call fatal_error("ipix out of range")

    ipix1 = ipix + 1 ! in {1, npix}
    nl2 = 2*nside
    nl4 = 4*nside
    ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
    fact1 = 1.5_dp*nside
    fact2 = 3.0_dp*nside**2

    do_vertex = .false.
    if (present(vertex)) then
       if (size(vertex,dim=1) >= 3 .and. size(vertex,dim=2) >= 4) then
          do_vertex = .true.
       else
          call fatal_error(" pix2vec_ring : vertex array has wrong size ")
       endif
    endif

    phi_nv = 0.0_dp
    phi_sv = 0.0_dp
    if (ipix1 <= ncap) then ! North Polar cap -------------

       hip   = ipix1/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
       iphi  = ipix1 - 2*iring*(iring - 1)

       z =  1.0_dp - iring**2 / fact2
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = 1.0_dp - (iring-1)**2 / fact2
          z_sv = 1.0_dp - (iring+1)**2 / fact2
          iphi_mod = MODULO(iphi-1, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi-1) / iring      ! in {0,1,2,3}
          if (iring > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
          phi_sv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
       endif


    elseif (ipix1 <= nl2*(5*nside+1)) then ! Equatorial region ------

       ip    = ipix1 - ncap - 1
       iring = INT( ip / nl4 ) + nside ! counted from North pole
       iphi  = MODULO(ip,nl4) + 1

       fodd  = 0.5_dp * (1 + MODULO(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
       z = (nl2 - iring) / fact1
       phi   = (real(iphi,kind=dp) - fodd) * PI /(2.0_dp*nside)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*nside)   ! half pixel width
          phi_nv = phi
          phi_sv = phi
          z_nv = (nl2 - iring +1) / fact1
          z_sv = (nl2 - iring -1) / fact1
          if (iring == nside) then ! northern transition
             z_nv = 1.0_dp - (nside-1)**2 / fact2
             iphi_mod = MODULO(iphi-1, nside) ! in {0,1,... nside-1}
             iphi_rat = (iphi-1) / nside      ! in {0,1,2,3}
             if (nside > 1) phi_nv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
          elseif (iring == 3*nside) then ! southern transition
             z_sv = -1.0_dp + (nside-1)**2 / fact2
             iphi_mod = MODULO(iphi-1, nside) ! in {0,1,... iring-1}
             iphi_rat = (iphi-1) / nside      ! in {0,1,2,3}
             if (nside > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(nside-1,kind=dp))
          endif
       endif

    else ! South Polar cap -----------------------------------

       ip    = npix - ipix1 + 1
       hip   = ip/2.0_dp
       fihip = AINT ( hip ,kind=DP)
       iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
       iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

       z = -1.0_dp + iring**2 / fact2
       phi   = (real(iphi,kind=dp) - 0.5_dp) * PI/(2.0_dp*iring)

       if (do_vertex) then
          hdelta_phi = PI/(4.0_dp*iring)   ! half pixel width
          z_nv = -1.0_dp + (iring+1)**2 / fact2
          z_sv = -1.0_dp + (iring-1)**2 / fact2
          iphi_mod = MODULO(iphi-1, iring) ! in {0,1,... iring-1}
          iphi_rat = (iphi-1) / iring      ! in {0,1,2,3}
          phi_nv                = HALFPI * (iphi_rat + (iphi_mod+1)/real(iring+1,kind=dp))
          if (iring > 1) phi_sv = HALFPI * (iphi_rat +  iphi_mod   /real(iring-1,kind=dp))
       endif

    endif

    ! pixel center
    sth = SQRT((1.0_dp-z)*(1.0_dp+z))
    vector(1) = sth * COS(phi)
    vector(2) = sth * SIN(phi)
    vector(3) = z

    if (do_vertex) then
       ! west vertex
       phi_wv      = phi - hdelta_phi
       vertex(1,2) = sth * COS(phi_wv)
       vertex(2,2) = sth * SIN(phi_wv)
       vertex(3,2) = z

       ! east vertex
       phi_ev      = phi + hdelta_phi
       vertex(1,4) = sth * COS(phi_ev)
       vertex(2,4) = sth * SIN(phi_ev)
       vertex(3,4) = z

       ! north vertex
       sth_nv = SQRT((1.0_dp-z_nv)*(1.0_dp+z_nv))
       vertex(1,1) = sth_nv * COS(phi_nv)
       vertex(2,1) = sth_nv * SIN(phi_nv)
       vertex(3,1) = z_nv

       ! south vertex
       sth_sv = SQRT((1.0_dp-z_sv)*(1.0_dp+z_sv))
       vertex(1,3) = sth_sv * COS(phi_sv)
       vertex(2,3) = sth_sv * SIN(phi_sv)
       vertex(3,3) = z_sv
    endif


    return
  end subroutine pix2vec_ring



end module pix_tools
