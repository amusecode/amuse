!************************************************************************
PROGRAM generate_coolrates

! Generate files containing tables of cooling and heating rates in
! the implicit cooling-module of RAMSES, to verify they are correct.
!------------------------------------------------------------------------
  use rt_icooling_module
  implicit none
  integer,parameter::nT=500,nXi=500
  real(dp),dimension(nT) :: a_T  ! Temperature array
  real(dp),dimension(nXi) :: a_X  ! Ionization array
  real(dp)::T0, T1, X0, X1, dlogT, dX
  real(dp),dimension(nT) :: a_Aa, a_Ab, a_B, a_dAadT, a_dAbdT, a_dBdT
  real(dp),dimension(nT,nXi) :: a_L, a_dLdT, a_dLdX 
  integer::i, j, ilun_aa, ilun_ab, ilun_b, ilun_daadt, ilun_dabdt,      &
           ilun_dbdt,ncols 
!------------------------------------------------------------------------

  T0=100 ; T1=1.d8
  dlogT = (log10(T1)-log10(T0))/(nT-1)
  do i = 0, nT-1
     a_T(i+1) = 10**(log10(T0) + dlogT * i)
  end do
  X0=0 ;  X1=1
  dX = (X1-X0)/(nXi-1)
  do i = 0, nX-1
     a_X(i+1) = X0 + dX * i
  end do

  do i = 1, nT
     a_Aa(i)    = comp_AlphaA(a_T(i))
     a_Ab(i)    = comp_AlphaB(a_T(i))
     a_B(i)     = comp_Beta(a_T(i))
     a_dAadT(i) = comp_dAlphaA_dT(a_T(i))
     a_dAbdT(i) = comp_dAlphaB_dT(a_T(i))
     a_dBdT(i)  = comp_dBeta_dT(a_T(i))
     do j = 1,nXi
        a_L(i,j)  = comp_coolrate(a_T(i), a_X(j), 1.d0, 1.d0, a_dLdT(i,j), a_dLdX(i,j))
     end do
  end do
  ilun_aa=21; ilun_ab=22; ilun_b=23; ilun_daadt=24; ilun_dabdt=24; ilun_dbdt=25
  open(ilun_aa,   file='alpha_a.list',    status='unknown')
  open(ilun_ab,   file='alpha_b.list',    status='unknown')
  open(ilun_b,    file='beta.list',       status='unknown')
  open(ilun_daadt,file='dalpha_a_dt.list',status='unknown')
  open(ilun_dabdt,file='dalpha_b_dt.list',status='unknown')
  open(ilun_dbdt, file='dbeta_dt.list',   status='unknown')

  write(ilun_aa,*) ''
  write(ilun_ab,*) ''
  write(ilun_b,*) ''
  write(ilun_daadt,*) ''
  write(ilun_dabdt,*) ''
  write(ilun_dbdt,*) ''

  write(ilun_aa,*) ''
  write(ilun_ab,*) ''
  write(ilun_b,*) ''
  write(ilun_daadt,*) ''
  write(ilun_dabdt,*) ''
  write(ilun_dbdt,*) ''

  ncols=2
  write(ilun_aa,*) ncols,nT
  write(ilun_ab,*) ncols,nT
  write(ilun_b,*) ncols,nT
  write(ilun_daadt,*) ncols,nT
  write(ilun_dabdt,*) ncols,nT
  write(ilun_dbdt,*) ncols,nT

  write(ilun_aa,*) ''
  write(ilun_ab,*) ''
  write(ilun_b,*) ''
  write(ilun_daadt,*) ''
  write(ilun_dabdt,*) ''
  write(ilun_dbdt,*) ''

  do i = 1,nT
     write(ilun_aa,   '(f21.6, d21.10)') a_T(i), a_Aa(i)
     write(ilun_ab,   '(f21.6, d21.10)') a_T(i), a_Ab(i)
     write(ilun_b,    '(f21.6, d21.10)') a_T(i), a_B(i)
     write(ilun_daadt,'(f21.6, d21.10)') a_T(i), a_dAadT(i)
     write(ilun_dabdt,'(f21.6, d21.10)') a_T(i), a_dAbdT(i)
     write(ilun_dbdt, '(f21.6, d21.10)') a_T(i), a_dBdt(i)
  end do

  close(ilun_aa)
  close(ilun_ab)
  close(ilun_b)
  close(ilun_daadt)
  close(ilun_dabdt)
  close(ilun_dbdt)


END PROGRAM generate_coolrates
