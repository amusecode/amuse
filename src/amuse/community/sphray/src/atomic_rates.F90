!> \file atomic_rates.F90

!> \brief the module that handles atomic rates interpolation tables

module atomic_rates_mod
use myf03_mod
use hui_gnedin_atomic_rates_mod
use cen_atomic_rates_mod
implicit none

  character(14), parameter :: ratenames(25) = &
       (/ "LogT          ",                                                        &
          "HIci          ", "HeIci         ", "HeIIci        ",                    &
          "HIIrcA        ", "HeIIrcA       ", "HeIIIrcA      ",                    &
          "HIIrcB        ", "HeIIrcB       ", "HeIIIrcB      ", "HeDrc         ",  &
          "HIIrccA       ", "HeIIrccA      ", "HeIIIrccA     ",                    &
          "HIIrccB       ", "HeIIrccB      ", "HeIIIrccB     ", "HeDrcc        ",  &
          "HIcic         ", "HeIcic        ", "HeIIcic       ", "He23cic       ",  &
          "HIcec         ", "HeIcec        ", "HeIIcec       "                     /) !< rate names


  integer(i8b), private, parameter :: RateHeaderLines = 0 !< lines in header


!> atomic rates for all processes at a single temperature
!----------------------------------------------------------
type atomic_rates_type
   real(r4b) :: logT   !< log10 temperature

   real(r4b) :: HIci      !< HI   collisional ionization rate 
   real(r4b) :: HeIci     !< HeI  collisional ionization rate 
   real(r4b) :: HeIIci    !< HeII collisional ionization rate 

   real(r4b) :: HIIrcA    !< HII   recombination rate (case A)
   real(r4b) :: HeIIrcA   !< HeII  recombination rate (case A)
   real(r4b) :: HeIIIrcA  !< HeIII recombination rate (case A)

   real(r4b) :: HIIrcB    !< HII   recombination rate (case B)
   real(r4b) :: HeIIrcB   !< HeII  recombination rate (case B)
   real(r4b) :: HeIIIrcB  !< HeIII recombination rate (case B)
   real(r4b) :: HeDrc     !< dielectronic He recombination rate   

   real(r4b) :: HIIrccA   !< HII   recombination cooling rate (case A)
   real(r4b) :: HeIIrccA  !< HeII  recombination cooling rate (case A)
   real(r4b) :: HeIIIrccA !< HeIII recombination cooling rate (case A)

   real(r4b) :: HIIrccB   !< HII   recombination cooling rate (case B)
   real(r4b) :: HeIIrccB  !< HeII  recombination cooling rate (case B)
   real(r4b) :: HeIIIrccB !< HeIII recombination cooling rate (case B)
   real(r4b) :: HeDrcc    !< dielectronic He recombination cooling rate
   
   real(r4b) :: HIcic     !< HI   collisional ionization cooling rate
   real(r4b) :: HeIcic    !< HeI  collisional ionization cooling rate
   real(r4b) :: HeIIcic   !< HeII collisional ionization cooling rate
   real(r4b) :: He23cic   !< He23 collisional ionization cooling rate
   
   real(r4b) :: HIcec     !< HI   collisional excitation cooling rate
   real(r4b) :: HeIcec    !< HeI  collisional excitation cooling rate
   real(r4b) :: HeIIcec   !< HeII collisional excitation cooling rate
end type atomic_rates_type


!> atomic rates for a single process
!------------------------------------
type atomic_rate_type
   character(14) :: source             !< reference for the rate fit
   real(r4b), allocatable :: rate(:)   !< allocate to number of temperature bins
end type atomic_rate_type


!> atomic rates table
!---------------------------
type atomic_rates_table_type

   real(r4b) :: logT1
   real(r4b) :: logT2
   integer   :: Tbins
   real(r4b) :: dlogT

   real(r4b), allocatable :: logT(:)   !< log10 temperature

   type(atomic_rate_type) :: HIci      !< HI   collisional ionization rate 
   type(atomic_rate_type) :: HeIci     !< HeI  collisional ionization rate 
   type(atomic_rate_type) :: HeIIci    !< HeII collisional ionization rate 

   type(atomic_rate_type) :: HIIrcA    !< HII   recombination rate (case A)
   type(atomic_rate_type) :: HeIIrcA   !< HeII  recombination rate (case A)
   type(atomic_rate_type) :: HeIIIrcA  !< HeIII recombination rate (case A)

   type(atomic_rate_type) :: HIIrcB    !< HII   recombination rate (case B)
   type(atomic_rate_type) :: HeIIrcB   !< HeII  recombination rate (case B)
   type(atomic_rate_type) :: HeIIIrcB  !< HeIII recombination rate (case B)
   type(atomic_rate_type) :: HeDrc     !< dielectronic He recombination rate   

   type(atomic_rate_type) :: HIIrccA   !< HII   recombination cooling rate (case A)
   type(atomic_rate_type) :: HeIIrccA  !< HeII  recombination cooling rate (case A)
   type(atomic_rate_type) :: HeIIIrccA !< HeIII recombination cooling rate (case A)

   type(atomic_rate_type) :: HIIrccB   !< HII   recombination cooling rate (case B)
   type(atomic_rate_type) :: HeIIrccB  !< HeII  recombination cooling rate (case B)
   type(atomic_rate_type) :: HeIIIrccB !< HeIII recombination cooling rate (case B)
   type(atomic_rate_type) :: HeDrcc    !< dielectronic He recombination cooling rate
   
   type(atomic_rate_type) :: HIcic     !< HI   collisional ionization cooling rate
   type(atomic_rate_type) :: HeIcic    !< HeI  collisional ionization cooling rate
   type(atomic_rate_type) :: HeIIcic   !< HeII collisional ionization cooling rate
   type(atomic_rate_type) :: He23cic   !< He23 collisional ionization cooling rate
   
   type(atomic_rate_type) :: HIcec     !< HI   collisional excitation cooling rate
   type(atomic_rate_type) :: HeIcec    !< HeI  collisional excitation cooling rate
   type(atomic_rate_type) :: HeIIcec   !< HeII collisional excitation cooling rate
      
end type atomic_rates_table_type




  
contains



  !> allocates space for a table
  !---------------------------------------------
  subroutine allocate_atomic_rates_table( table, n_bins )
    type(atomic_rates_table_type) :: table
    integer :: n_bins

    allocate( table%logT(n_bins) )

    allocate( table%HIci%rate(n_bins) )
    allocate( table%HeIci%rate(n_bins) )
    allocate( table%HeIIci%rate(n_bins) )
    
    allocate( table%HIIrcA%rate(n_bins) )
    allocate( table%HeIIrcA%rate(n_bins) )
    allocate( table%HeIIIrcA%rate(n_bins) )
    
    allocate( table%HIIrcB%rate(n_bins) )
    allocate( table%HeIIrcB%rate(n_bins) )
    allocate( table%HeIIIrcB%rate(n_bins) )
    allocate( table%HeDrc%rate(n_bins) )
    
    allocate( table%HIIrccA%rate(n_bins) )
    allocate( table%HeIIrccA%rate(n_bins) )
    allocate( table%HeIIIrccA%rate(n_bins) )
    
    allocate( table%HIIrccB%rate(n_bins) )
    allocate( table%HeIIrccB%rate(n_bins) )
    allocate( table%HeIIIrccB%rate(n_bins) )
    allocate( table%HeDrcc%rate(n_bins) )
    
    allocate( table%HIcic%rate(n_bins) )
    allocate( table%HeIcic%rate(n_bins) )
    allocate( table%HeIIcic%rate(n_bins) )
    allocate( table%He23cic%rate(n_bins) )
    
    allocate( table%HIcec%rate(n_bins) )
    allocate( table%HeIcec%rate(n_bins) )
    allocate( table%HeIIcec%rate(n_bins) )
    
  end subroutine allocate_atomic_rates_table


  !> open and read in the atomic rate file
  !---------------------------------------------
  subroutine read_atomic_rates_file(table, AtomicRatesFile)
    type(atomic_rates_table_type) :: table
    character(clen) :: AtomicRatesFile

    character(clen), parameter :: myname = "read_atomic_rates_file"
    integer, parameter :: verb=1
    logical, parameter :: crash = .true.

    integer(i4b) :: lun,err,i
    logical :: fthere
    character(14) :: tags(25)
    character(14) :: srcs(24)
    character(clen) :: str



    write(str,'(A,T27,A)') 'using rates file: ', trim(AtomicRatesFile)
    call mywrite(str, verb) 
    call mywrite('', verb)

    inquire(file=AtomicRatesFile, exist=fthere)
    if (.not. fthere) call myerr('cant find atomic rates file',myname,crash)

    call open_formatted_file_r(AtomicRatesFile,lun)
    do i = 1,RateHeaderLines
       read(lun,*)
    end do

    read(lun,*) table%logT1, table%logT2, table%Tbins, table%dlogT
    read(lun,*) tags
    read(lun,*) srcs

    table%HIci%source = srcs(1)       
    table%HeIci%source = srcs(2)      
    table%HeIIci%source = srcs(3)     
    
    table%HIIrcA%source = srcs(4)     
    table%HeIIrcA%source = srcs(5)    
    table%HeIIIrcA%source = srcs(6)   
    
    table%HIIrcB%source = srcs(7)     
    table%HeIIrcB%source = srcs(8)    
    table%HeIIIrcB%source = srcs(9)     
    table%HeDrc%source = srcs(10)      
    
    table%HIIrccA%source = srcs(11)    
    table%HeIIrccA%source = srcs(12)    
    table%HeIIIrccA%source = srcs(13)  
    
    table%HIIrccB%source = srcs(14)    
    table%HeIIrccB%source = srcs(15)   
    table%HeIIIrccB%source = srcs(16)      
    table%HeDrcc%source = srcs(17) 
    
    table%HIcic%source = srcs(18)      
    table%HeIcic%source = srcs(19)     
    table%HeIIcic%source = srcs(20)    
    table%He23cic%source = srcs(21)    
    
    table%HIcec%source = srcs(22)      
    table%HeIcec%source = srcs(23)     
    table%HeIIcec%source = srcs(24)   
    

    call allocate_atomic_rates_table( table, table%Tbins )

    do i = 1, table%Tbins

     read(lun,*) &

          table%logT(i), & 
          
          table%HIci%rate(i), &       
          table%HeIci%rate(i), &      
          table%HeIIci%rate(i), &     
          
          table%HIIrcA%rate(i), &     
          table%HeIIrcA%rate(i), &    
          table%HeIIIrcA%rate(i), &   
          
          table%HIIrcB%rate(i), &     
          table%HeIIrcB%rate(i), &    
          table%HeIIIrcB%rate(i), &     
          table%HeDrc%rate(i), &      
          
          table%HIIrccA%rate(i), &    
          table%HeIIrccA%rate(i), &    
          table%HeIIIrccA%rate(i), &  
          
          table%HIIrccB%rate(i), &    
          table%HeIIrccB%rate(i), &   
          table%HeIIIrccB%rate(i), &      
          table%HeDrcc%rate(i), & 
          
          table%HIcic%rate(i), &      
          table%HeIcic%rate(i), &     
          table%HeIIcic%rate(i), &    
          table%He23cic%rate(i), &    
          
          table%HIcec%rate(i), &      
          table%HeIcec%rate(i), &     
          table%HeIIcec%rate(i)   

  end do

 

     
  end subroutine read_atomic_rates_file











  !-------------------------------------------------------------------------------
  !> Calculates collisional ionization equilibrium for all species 
  !! at given temperature using the rates table
  subroutine calc_colion_eq_table(table, Tin, caseA, xvec)
    type(atomic_rates_table_type), intent(in) :: table
    real(r8b), intent(in) :: Tin
    logical, intent(in) :: caseA(2)
    real(r4b), intent(out) :: xvec(5)
    
    type(atomic_rates_type) :: k
    real(r8b) :: GGHI, GGHeI, GGHeII
    real(r8b) :: RRHII, RRHeII, RRHeIII
    real(r8b) :: den
    
    call get_atomic_rates(Tin, table,  k)
    
    GGHI   = k%HIci 
    GGHeI  = k%HeIci
    GGHeII = k%HeIIci
    
    if (caseA(1)) then
       RRHII   = k%HIIrcA 
    else
       RRHII   = k%HIIrcB 
    end if
    
    if (caseA(2)) then
       RRHeII  = k%HeIIrcA 
       RRHeIII = k%HeIIIrcA 
    else
       RRHeII  = k%HeIIrcB 
       RRHeIII = k%HeIIIrcB 
    end if
    
    den = (RRHII + GGHI)
    xvec(1) = RRHII / den
    xvec(2) = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xvec(3) = RRHeII * RRHeIII / den
    xvec(4) = RRHeIII * GGHeI  / den 
    xvec(5) = GGHeI * GGHeII   / den
    
  end subroutine calc_colion_eq_table
  


  !> Calculates collisional ionization equilibrium for all species 
  !! at given temperature using the fits directly
  !-------------------------------------------------------------------------------
  subroutine calc_colion_eq_fits(fit, Tin, caseA, xvec)
    character(*), intent(in) :: fit
    real(r8b), intent(in) :: Tin
    logical, intent(in) :: caseA(2)
    real(r8b), intent(out) :: xvec(5)
    
    real(r8b) :: T
    real(r8b) :: GGHI, GGHeI, GGHeII
    real(r8b) :: RRHII, RRHeII, RRHeIII
    real(r8b) :: den


    T = Tin
    
    select case(trim(fit))

       case("cen")
          GGHI   = Cen_HI_col_ion(T)
          GGHeI  = Cen_HeI_col_ion(T)
          GGHeII = Cen_HeII_col_ion(T)

          if (caseA(1)) then
             RRHII   = Cen_HII_recombA(T)
          else
             write(*,*) "no case B rates from Cen"
             stop
          end if

          if (caseA(2)) then
             RRHeII  = Cen_HeII_recombA(T)
             RRHeIII = Cen_HeIII_recombA(T)
          else
             write(*,*) "no case B rates from Cen"
             stop
          end if


       case("hui")        
          GGHI   = Hui_HI_col_ion(T)
          GGHeI  = Hui_HeI_col_ion(T)
          GGHeII = Hui_HeII_col_ion(T)
    
          if (caseA(1)) then
             RRHII   = Hui_HII_recombA(T)
          else
             RRHII   = Hui_HII_recombB(T)
          end if
    
          if (caseA(2)) then
             RRHeII  = Hui_HeII_recombA(T)
             RRHeIII = Hui_HeIII_recombA(T)
          else
             RRHeII  = Hui_HeII_recombB(T)
             RRHeIII = Hui_HeIII_recombB(T)
          end if
    
       case default

          write(*,*) "fit type ", trim(fit), " not recognized"
          stop
    
    end select


    den = (RRHII + GGHI)
    xvec(1) = RRHII / den
    xvec(2) = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xvec(3) = RRHeII * RRHeIII / den
    xvec(4) = RRHeIII * GGHeI  / den 
    xvec(5) = GGHeI * GGHeII   / den
    
  end subroutine calc_colion_eq_fits



  !--------------------------------------------------------------------
  !> takes a temperature and calculates the atomic rate table index and
  !! how far the true temperature is past the table temperature
  subroutine get_Tindx_and_Tremain(table, T, Tindx, Tremain)
    type(atomic_rates_table_type), intent(in) :: table  
    real(r8b), intent(in) :: T         !< input temperature
    integer(i8b), intent(out) :: Tindx !< atomic rate table index
    real(r8b), intent(out) :: Tremain  !< rates%(T) + Tremain = T
    real(r8b) :: Tnum
    
    if (T < 0.0) stop "T < 0.0 in get_Tindx_and_Tremain"

    if (log10(T) < table%logT1) stop "log(T) < table%logT1"

    Tnum = ( log10(T) - table%logT1 ) / table%dlogT
    Tindx = ceiling(Tnum)
    
    if (Tindx == 0) then
       Tindx = 1
       Tremain = 0.0d0
    end if
    Tremain = Tnum - (Tindx-1)
    
    if(Tindx < 1)             stop "Tindx < 1 "
    if(Tindx > table%Tbins-1) stop "Tindx > Tbins - 1 "
    
  end subroutine get_Tindx_and_Tremain




  

  !=================================
  !> get atomic rates from table 
  subroutine get_atomic_rates(T, table, k)

    real(r8b), intent(in) :: T                             !< input temperature
    type(atomic_rates_table_type), intent(in) :: table     !< interpolation table
    type(atomic_rates_type), intent(out) :: k              !< returns rates
  
    integer(i8b) :: Tindx
    real(r8b) :: Tremain, rdif

    call get_Tindx_and_Tremain(table, T, Tindx, Tremain)

    ! Collisional Ionization

    rdif = ( table%HIci%rate(Tindx+1) - table%HIci%rate(Tindx) ) * Tremain
    k%HIci = table%HIci%rate(Tindx) + rdif 
    
    rdif = ( table%HeIci%rate(Tindx+1) - table%HeIci%rate(Tindx) ) * Tremain
    k%HeIci = table%HeIci%rate(Tindx) + rdif 

    rdif = ( table%HeIIci%rate(Tindx+1) - table%HeIIci%rate(Tindx) ) * Tremain
    k%HeIIci = table%HeIIci%rate(Tindx) + rdif 

        
    ! Recombination
    
    rdif = ( table%HIIrcA%rate(Tindx+1) - table%HIIrcA%rate(Tindx) ) * Tremain
    k%HIIrcA = table%HIIrcA%rate(Tindx) + rdif 

    rdif = ( table%HIIrcB%rate(Tindx+1) - table%HIIrcB%rate(Tindx) ) * Tremain
    k%HIIrcB = table%HIIrcB%rate(Tindx) + rdif 

    rdif = ( table%HeIIrcA%rate(Tindx+1) - table%HeIIrcA%rate(Tindx) ) * Tremain
    k%HeIIrcA = table%HeIIrcA%rate(Tindx) + rdif 

    rdif = ( table%HeIIrcB%rate(Tindx+1) - table%HeIIrcB%rate(Tindx) ) * Tremain
    k%HeIIrcB = table%HeIIrcB%rate(Tindx) + rdif 

    rdif = ( table%HeIIIrcA%rate(Tindx+1) - table%HeIIIrcA%rate(Tindx) ) * Tremain
    k%HeIIIrcA = table%HeIIIrcA%rate(Tindx) + rdif 

    rdif = ( table%HeIIIrcB%rate(Tindx+1) - table%HeIIIrcB%rate(Tindx) ) * Tremain
    k%HeIIIrcB = table%HeIIIrcB%rate(Tindx) + rdif 

    rdif = (table%HeDrc%rate(Tindx+1) - table%HeDrc%rate(Tindx)) * Tremain
    k%HeDrc = table%HeDrc%rate(Tindx) + rdif 


    
    ! Recombination Cooling
    
    rdif = ( table%HIIrccA%rate(Tindx+1) - table%HIIrccA%rate(Tindx) ) * Tremain
    k%HIIrccA = table%HIIrccA%rate(Tindx) + rdif 

    rdif = ( table%HIIrccB%rate(Tindx+1) - table%HIIrccB%rate(Tindx) ) * Tremain
    k%HIIrccB = table%HIIrccB%rate(Tindx) + rdif 

    rdif = ( table%HeIIrccA%rate(Tindx+1) - table%HeIIrccA%rate(Tindx) ) * Tremain
    k%HeIIrccA = table%HeIIrccA%rate(Tindx) + rdif 

    rdif = ( table%HeIIrccB%rate(Tindx+1) - table%HeIIrccB%rate(Tindx) ) * Tremain
    k%HeIIrccB = table%HeIIrccB%rate(Tindx) + rdif 

    rdif = ( table%HeIIIrccA%rate(Tindx+1) - table%HeIIIrccA%rate(Tindx) ) * Tremain
    k%HeIIIrccA = table%HeIIIrccA%rate(Tindx) + rdif 

    rdif = ( table%HeIIIrccB%rate(Tindx+1) - table%HeIIIrccB%rate(Tindx) ) * Tremain
    k%HeIIIrccB = table%HeIIIrccB%rate(Tindx) + rdif 

    rdif = ( table%HeDrcc%rate(Tindx+1) - table%HeDrcc%rate(Tindx) ) * Tremain
    k%HeDrcc = table%HeDrcc%rate(Tindx) + rdif 

    
    ! Collisional Ionization Cooling

    rdif = ( table%HIcic%rate(Tindx+1) - table%HIcic%rate(Tindx) ) * Tremain
    k%HIcic = table%HIcic%rate(Tindx) + rdif 
    
    rdif = ( table%HeIcic%rate(Tindx+1) - table%HeIcic%rate(Tindx) ) * Tremain
    k%HeIcic = table%HeIcic%rate(Tindx) + rdif 

    rdif = ( table%HeIIcic%rate(Tindx+1) - table%HeIIcic%rate(Tindx) ) * Tremain
    k%HeIIcic = table%HeIIcic%rate(Tindx) + rdif 
    
    rdif = ( table%He23cic%rate(Tindx+1) - table%He23cic%rate(Tindx) ) * Tremain
    k%He23cic = table%He23cic%rate(Tindx) + rdif 
    
    ! Collisional Excitation Cooling
    
    rdif = ( table%HIcec%rate(Tindx+1) - table%HIcec%rate(Tindx) ) * Tremain
    k%HIcec = table%HIcec%rate(Tindx) + rdif 
    
    rdif = ( table%HeIcec%rate(Tindx+1) - table%HeIcec%rate(Tindx) ) * Tremain
    k%HeIcec = table%HeIcec%rate(Tindx) + rdif 

    rdif = ( table%HeIIcec%rate(Tindx+1) - table%HeIIcec%rate(Tindx) ) * Tremain
    k%HeIIcec = table%HeIIcec%rate(Tindx) + rdif 
    
  end subroutine get_atomic_rates


  !> writes an atomic rates table to file
  !---------------------------------------------
  subroutine write_atomic_rates_table_to_file( table, file )
    type(atomic_rates_table_type) :: table
    character(*) :: file

    character(clen) :: fmt
    integer :: j

  
    open(unit=10,file=file)

    fmt = "(2ES12.5,I8,ES12.5)"
    write(10,fmt) table%logT1, table%logT2, table%Tbins, table%dlogT
    
    fmt = "(1X,25(A))"
    write(10,fmt) ratenames
    
    write(10,fmt) &
         "              ", & 
         
         table%HIci%source, &       
         table%HeIci%source, &      
         table%HeIIci%source, &     
         
         table%HIIrcA%source, &     
         table%HeIIrcA%source, &    
         table%HeIIIrcA%source, &   
         
         table%HIIrcB%source, &     
         table%HeIIrcB%source, &    
         table%HeIIIrcB%source, &     
         table%HeDrc%source, &      
         
         table%HIIrccA%source, &    
         table%HeIIrccA%source, &    
         table%HeIIIrccA%source, &  
         
         table%HIIrccB%source, &    
         table%HeIIrccB%source, &   
         table%HeIIIrccB%source, &      
         table%HeDrcc%source, & 
         
         table%HIcic%source, &      
         table%HeIcic%source, &     
         table%HeIIcic%source, &    
         table%He23cic%source, &    
         
         table%HIcec%source, &      
         table%HeIcec%source, &     
         table%HeIIcec%source   
    
    
    fmt = "(25(ES12.5,2X))"
    
    do j = 1, table%Tbins
       write(10,fmt) &
            
            table%logT(j), & 
            
            table%HIci%rate(j), &       
            table%HeIci%rate(j), &      
            table%HeIIci%rate(j), &     
            
            table%HIIrcA%rate(j), &     
            table%HeIIrcA%rate(j), &    
            table%HeIIIrcA%rate(j), &   
            
            table%HIIrcB%rate(j), &     
            table%HeIIrcB%rate(j), &    
            table%HeIIIrcB%rate(j), &     
            table%HeDrc%rate(j), &      
            
            table%HIIrccA%rate(j), &    
            table%HeIIrccA%rate(j), &    
            table%HeIIIrccA%rate(j), &  
            
            table%HIIrccB%rate(j), &    
            table%HeIIrccB%rate(j), &   
            table%HeIIIrccB%rate(j), &      
            table%HeDrcc%rate(j), & 
            
            table%HIcic%rate(j), &      
            table%HeIcic%rate(j), &     
            table%HeIIcic%rate(j), &    
            table%He23cic%rate(j), &    
            
            table%HIcec%rate(j), &      
            table%HeIcec%rate(j), &     
            table%HeIIcec%rate(j)   
       
    end do
    
  end subroutine write_atomic_rates_table_to_file



  !> writes some table values to log file
  !---------------------------------------------
  subroutine write_atomic_rates_to_log_file( table, OutputDir )
    type(atomic_rates_table_type) :: table
    character(*) :: OutputDir
    character(clen), parameter :: myname = "write_atomic_rates_to_log"
    integer, parameter :: verb=1
    logical, parameter :: crash = .true.

    character(clen) :: logfile
    integer(i4b) :: loglun
    type(atomic_rates_type) :: k

    real(r8b) :: Tdum
    real(r8b) :: GGHI, RRHII
    real(r8b) :: GGHeI, RRHeII
    real(r8b) :: GGHeII, RRHeIII
    real(r8b) :: xHI,xHII,xHeI,xHeII,xHeIII,den

    logfile = trim(OutputDir) // '/' // 'atomic_rates.log'
    call open_formatted_file_w(logfile,loglun)

    write(loglun,'(A)') 'atomic rates header'
    write(loglun,'(A,ES12.4)') "log10 of temperature floor   = ", table%logT1
    write(loglun,'(A,ES12.4)') "log10 of temperature ceiling = ", table%logT2
    write(loglun,'(A,I5)')     "rate file entries            = ", table%Tbins
    write(loglun,'(A,ES12.4)') "log10 spacing of temperature = ", table%dlogT
    write(loglun,*)    


    ! write out characteristic rates @ T = 10,000 K
    Tdum = 1.0E4
    call get_atomic_rates(Tdum, table, k)

    100 format(A,ES12.5,A)
    101 format(A,":",T32,ES12.5,A)
    102 format(A,T15,ES12.3)

    write(loglun,*) 
    write(loglun,*) "-----------------------------------------------------------"
    write(loglun,100) "Atomic Rates @ ", Tdum," K from rate tables"
    write(loglun,*) 
    write(loglun,101) "HI Col Ion", k%HIci,   " cm^3/s"

    write(loglun,101) "HeI Col Ion", k%HeIci,  " cm^3/s"
    write(loglun,101) "HeII Col Ion", k%HeIIci, " cm^3/s"

    write(loglun,*) 
    write(loglun,101) "HII Rec (A)", k%HIIrcA,   " cm^3/s"
    write(loglun,101) "HII Rec (B)", k%HIIrcB,   " cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeII Rec (A)", k%HeIIrcA,  " cm^3/s"
    write(loglun,101) "HeII Rec (B)", k%HeIIrcB,  " cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeIII Rec (A)", k%HeIIIrcA, " cm^3/s"
    write(loglun,101) "HeIII Rec (B)", k%HeIIIrcB, " cm^3/s"
    write(loglun,*) 
    write(loglun,101) "He D Rec", k%HeDrc, "cm^3/s"
    write(loglun,*) 

    write(loglun,101) "HII Rec Cool (A)", k%HIIrccA, " ergs cm^3/s"
    write(loglun,101) "HII Rec Cool (B)", k%HIIrccB, " ergs cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeII Rec Cool (A)", k%HeIIrccA,  " ergs cm^3/s"
    write(loglun,101) "HeII Rec Cool (B)", k%HeIIrccB,  " ergs cm^3/s"
    write(loglun,*) 
    write(loglun,101) "HeIII Rec Cool (A)", k%HeIIIrccA, " ergs cm^3/s"
    write(loglun,101) "HeIII Rec Cool (B)", k%HeIIIrccB, " ergs cm^3/s"
    write(loglun,*) 
    write(loglun,101) "He D Rec Cool", k%HeDrc, "cm^3/s"
    write(loglun,*) 

    write(loglun,101) "HI Col Ion Cool", k%HIcic, " ergs cm^3/s"
    write(loglun,101) "HeI Col Ion Cool", k%HeIcic, " ergs cm^3/s"
    write(loglun,101) "HeII Col Ion Cool", k%HeIIcic, " ergs cm^3/s"
    write(loglun,101) "He 23 Col Ion Cool", k%He23cic, " ergs cm^3/s"
    write(loglun,*) 


    write(loglun,101) "HI Col Ext Cool", k%HIcec, " ergs cm^3/s"
    write(loglun,101) "HeI Col Ext Cool", k%HeIcec, " ergs cm^6/s"
    write(loglun,101) "HeII Col Ext Cool", k%HeIIcec, " ergs cm^3/s"

    write(loglun,*) 

      
    write(loglun,'(A)') "Collisional Equilibrium"
    write(loglun,*) 
      
    GGHI   = k%HIci 
    GGHeI  = k%HeIci
    GGHeII = k%HeIIci
    
    RRHII   = k%HIIrcA
    RRHeII  = k%HeIIrcA
    RRHeIII = k%HeIIIrcA
    
    den = (RRHII + GGHI)
    xHI = RRHII / den
    xHII = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xHeI   = RRHeII * RRHeIII / den
    xHeII  = RRHeIII * GGHeI  / den 
    xHeIII = GGHeI * GGHeII   / den
    
    write(loglun,*)   " CASE A"
    write(loglun,102) "xHIeq    = ", xHI
    write(loglun,102) "xHIIeq   = ", xHII
    write(loglun,102) "xHeIeq   = ", xHeI
    write(loglun,102) "xHeIIeq  = ", xHeII
    write(loglun,102) "xHeIIIeq = ", xHeIII
    write(loglun,*) 
    
    RRHII   = k%HIIrcB
    RRHeII  = k%HeIIrcB
    RRHeIII = k%HeIIIrcB
    
    den = (RRHII + GGHI)
    xHI = RRHII / den
    xHII = GGHI / den
    
    den = (RRHeII * RRHeIII + RRHeIII * GGHeI + GGHeI * GGHeII)
    xHeI   = RRHeII * RRHeIII / den
    xHeII  = RRHeIII * GGHeI  / den 
    xHeIII = GGHeI * GGHeII   / den
    
    write(loglun,*)   " CASE B"
    write(loglun,102) "xHIeq    = ", xHI
    write(loglun,102) "xHIIeq   = ", xHII
    write(loglun,102) "xHeIeq   = ", xHeI
    write(loglun,102) "xHeIIeq  = ", xHeII
    write(loglun,102) "xHeIIIeq = ", xHeIII
    write(loglun,*) 



  end subroutine write_atomic_rates_to_log_file



end module atomic_rates_mod
