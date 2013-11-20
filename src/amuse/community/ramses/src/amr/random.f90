!=======================================================================
!
!   P A R A L L E L   R A N D O M   N U M B E R   G E N E R A T O R
!
!=======================================================================
! Here's how to use these functions:
!       Initialization: call Rans( N, StartVal, Seeds )
! This returns an array of N seeds;
! the returned seeds have values that partition the basic
! ranf cycle of 2**46 pseudorandom reals in [0,1] into independent sub-
! sequences. The second argument, StartVal, can be a cycle-starting seed
! for the full-period generator; if it is zero, a special seed will be
! used that produces statistically desirable behavior.
!       Use: call Ranf( Seed, RandNum )
! This returns a pseudorandom real and a seed to be passed back in the
! next invocation. The returned seed carries the state normally hidden
! in imperative generators.
!=======================================================================
module random
  integer,parameter :: IRandNumSize = 4, IBinarySize = 48
  integer,parameter :: Mod4096DigitSize = 12
  integer,parameter :: IZero = 0 
  integer,parameter :: NPoissonLimit = 10

  integer,dimension( IRandNumSize ) :: &
       & Multiplier = (/373, 3707, 1442, 647/)
  integer,dimension( IRandNumSize ) :: &
       & DefaultSeed = (/ 3281, 4041, 595, 2376 /)
  real(kind=8),dimension( IRandNumSize ) :: &
       & Divisor = (/281474976710656.0,68719476736.0,16777216.0,4096.0/)

  integer :: IGauss = 0
  real(kind=8) :: GaussBak = 0.0d0
    
contains
  !=======================================================================
  subroutine ranf( Seed, RandNum )
    implicit none
    integer,dimension( IRandNumSize ) :: Seed
    real(kind=8) :: RandNum
    
    integer,dimension( IRandNumSize ) :: OutSeed

    RandNum = float( Seed( 4 ) ) / Divisor( 4 ) + &
         &    float( Seed( 3 ) ) / Divisor( 3 ) + &
         &    float( Seed( 2 ) ) / Divisor( 2 ) + &
         &    float( Seed( 1 ) ) / Divisor( 1 )
    
    call ranfmodmult( Multiplier, Seed, OutSeed )
    
    Seed = OutSeed

    return
  end subroutine ranf
  !=======================================================================

  !=======================================================================
  subroutine poissdev( Seed, AverNum , PoissNum)
    implicit none
    integer,dimension( IRandNumSize ) :: Seed
    real(kind=8) :: AverNum
    integer :: PoissNum

    real(kind=8) :: Norm, Repar, Proba
    real(kind=8) :: RandNum, GaussNum

    if(AverNum <= DBLE(NPoissonLimit)) then
       Norm=exp(-AverNum) 
       Repar=1.0d0
       PoissNum=0
       Proba=1.0d0
       call ranf(Seed,RandNum)
       do while(Repar*Norm <= RandNum .and. PoissNum <= 10*NPoissonLimit )
          PoissNum=PoissNum+1
          Proba=Proba*AverNum/PoissNum
          Repar=Repar+Proba
       end do
    else
       call gaussdev(Seed,GaussNum)
       GaussNum=GaussNum*sqrt(AverNum)-0.5+AverNum
       if(GaussNum<=0.0d0)GaussNum=0.0d0
       PoissNum=nint(GaussNum)
    endif

    return
  end subroutine poissdev
  !=======================================================================

  !=======================================================================
  subroutine gaussdev( Seed, GaussNum )
    implicit none
    integer,dimension( IRandNumSize ) :: Seed
    real(kind=8) :: GaussNum
    
    real(kind=8) :: fac,rsq,v1,v2
    
    if (IGauss.eq.IZero) then
       rsq=0.0d0
       do while (rsq.ge.1.0d0.or.rsq.le.0.0d0)
          call ranf(Seed,v1)
          call ranf(Seed,v2)
          v1=2.0d0*v1-1.0d0
          v2=2.0d0*v2-1.0d0
          rsq=v1**2+v2**2
       end do
       fac=sqrt(-2.0d0*log(rsq)/rsq)
       GaussBak=v1*fac
       GaussNum=v2*fac
       IGauss=1
    else
       GaussNum=GaussBak
       IGauss=0
    endif
    return
  END subroutine gaussdev
  !=======================================================================

  !=======================================================================
  subroutine rans( N, StartVal, SeedArray )
    use amr_commons,only:ncpu
    implicit none
    integer :: N
    integer :: StartVal
    integer,dimension( ncpu, IRandNumSize ) :: SeedArray

    integer,dimension( IRandNumSize ) :: atothek
    integer,dimension( IRandNumSize ) :: K
    integer,dimension( IRandNumSize ) :: InSeed
    integer,dimension( IRandNumSize ) :: OutSeed
    integer,dimension( IBinarySize )  :: KBinary

    integer :: I

    if( StartVal .eq. IZero ) then
       SeedArray( 1, 1 ) = DefaultSeed( 1 )
       SeedArray( 1, 2 ) = DefaultSeed( 2 )
       SeedArray( 1, 3 ) = DefaultSeed( 3 )
       SeedArray( 1, 4 ) = DefaultSeed( 4 )
    else
       SeedArray( 1, 1 ) = abs( StartVal )
       SeedArray( 1, 2 ) = IZero
       SeedArray( 1, 3 ) = IZero
       SeedArray( 1, 4 ) = IZero
    endif
    
    if( N .eq. 1 ) then
       atothek( 1 ) = Multiplier( 1 )
       atothek( 2 ) = Multiplier( 2 )
       atothek( 3 ) = Multiplier( 3 )
       atothek( 4 ) = Multiplier( 4 )
    else
       call ranfk( N, K )
       call ranfkbinary( K, KBinary )
       call ranfatok( Multiplier, KBinary, atothek )
       do I = 2, N
          InSeed( 1 ) = SeedArray( I-1, 1 )
          InSeed( 2 ) = SeedArray( I-1, 2 )
          InSeed( 3 ) = SeedArray( I-1, 3 )
          InSeed( 4 ) = SeedArray( I-1, 4 )
          call ranfmodmult( InSeed, atothek, OutSeed )
          SeedArray( I, 1 ) = OutSeed( 1 )
          SeedArray( I, 2 ) = OutSeed( 2 )
          SeedArray( I, 3 ) = OutSeed( 3 )
          SeedArray( I, 4 ) = OutSeed( 4 )
       end do
    endif
    
    return
  end subroutine rans
  !=======================================================================
  
  !=======================================================================
  subroutine ranfatok( a, Kbinary, atothek )
    implicit none
    integer,dimension( IRandNumSize ) :: a
    integer,dimension( IBinarySize )  :: KBinary
    integer,dimension( IRandNumSize ) :: atothek

    integer,dimension( IRandNumSize ) :: asubi
    
    integer :: I
    asubi( 1 ) = a( 1 )
    asubi( 2 ) = a( 2 )
    asubi( 3 ) = a( 3 )
    asubi( 4 ) = a( 4 )
    
    atothek( 1 ) = 1
    atothek( 2 ) = IZero
    atothek( 3 ) = IZero
    atothek( 4 ) = IZero
    
    do I = 1, 45
       if( KBinary( I ) .ne. IZero ) then
          call ranfmodmult( atothek, asubi, atothek )
       endif
       call ranfmodmult( asubi, asubi, asubi )
    end do
    
    return
  end subroutine ranfatok
  !=======================================================================
  
  !=======================================================================
  subroutine ranfk( N, K )
    implicit none
    integer :: N
    integer, dimension( IRandNumSize ) :: K

    integer :: nn, r4, r3, r2, q4, q3, q2, q1
    
    nn = N + iranfeven( N )
    
    q4 = 1024 / nn
    r4 = 1024 - (nn * q4)
    q3 = (r4 * 4096) / nn
    r3 = (r4 * 4096) - (nn * q3)
    q2 = (r3 * 4096) / nn
    r2 = (r3 * 4096) - (nn * q2)
    q1 = (r2 * 4096) / nn
    
    K( 1 ) = q1
    K( 2 ) = q2
    K( 3 ) = q3
    K( 4 ) = q4
    
    return
  end subroutine ranfk
  !=======================================================================
  
  !=======================================================================
  subroutine ranfkbinary( K, KBinary )
    implicit none
    integer, dimension( IRandNumSize ) :: K
    integer, dimension( IBinarySize ) ::  KBinary

    integer, dimension( Mod4096DigitSize ) :: Bits
    integer X, I, J
    
    do I = 1, 4
       X = K( I ) / 2
       Bits( 1 ) = iranfodd( K( I ) )     
       do J = 2, Mod4096DigitSize 
          Bits( J ) = iranfodd( X )
          X = X / 2
       end do
       do J = 1, Mod4096DigitSize
          KBinary( (I-1)*Mod4096DigitSize + J ) = Bits( J )
       end do
    end do
    
    return
  end subroutine ranfkbinary
  !=======================================================================
  
  !=======================================================================
  subroutine ranfmodmult( A, B, C )
    implicit none
    integer, dimension( IRandNumSize ) :: A, B, C

    integer :: j1, j2, j3, j4, k1, k2, k3, k4
        
    j1 = A( 1 ) * B( 1 )
    j2 = A( 1 ) * B( 2 ) + A( 2 ) * B( 1 )
    j3 = A( 1 ) * B( 3 ) + A( 2 ) * B( 2 ) + A( 3 ) * B( 1 )
    j4 = A( 1 ) * B( 4 ) + A( 2 ) * B( 3 ) + A( 3 ) * B( 2 ) + A( 4 ) * B( 1 )
    
    k1 = j1
    k2 = j2 + k1 / 4096
    k3 = j3 + k2 / 4096
    k4 = j4 + k3 / 4096
    
    C( 1 ) = mod( k1, 4096 )
    C( 2 ) = mod( k2, 4096 )
    C( 3 ) = mod( k3, 4096 )
    C( 4 ) = mod( k4, 4096 )
    
    return
  end subroutine ranfmodmult
  !=======================================================================
  
  !=======================================================================
  function iranfodd( N )
    integer N
    if( mod( N, 2 ) .eq. 0 ) then
       iranfodd = 0
    else
       iranfodd = 1
    endif
    
    return
  end function iranfodd
  !=======================================================================
  !=======================================================================
  function iranfeven( N )
    integer N
    if( mod( N, 2 ) .eq. 0 ) then
       iranfeven = 1
    else
       iranfeven = 0
    endif
    return
  end function iranfeven
  !=======================================================================

end module random
