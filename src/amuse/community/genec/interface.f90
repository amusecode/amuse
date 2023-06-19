module AmuseInterface
    use storage, only:&
            GenecStar,genec_star
    use helpers, only:&
            copy_to_genec_star,&
            copy_namelists_from_genec_star,&
            copy_from_genec_star,&
            copy_structure_from_genec_star,&
            copy2_from_genec_star
    use evol, only: kindreal,ldi,npondcouche

    type(genec_star) :: BackupBackupGenecStar
    type(genec_star) :: BackupGenecStar
    integer :: number_of_stars = 0
    integer :: evolve_steps
    public:: BackupGenecStar,number_of_stars
contains

subroutine init_or_restore_star(Star)
  ! This function replaces read_parameters and initialise_star in genec.f90
  use genec
  implicit none
  type(genec_star), intent(in) :: Star

  call copy_namelists_from_genec_star(Star)

  if (idebug > 0) then
    verbose = .true.
  endif

  if (idebug > 1) then
    write(*,*) 'initialisations...'
  endif

  tzero=999999999.d0

  idern=0
  id1=0
  id2=5
  dm_lost=0.d0

  agdp = agdr    ! )
  agds = agdr    ! ) bounds on the corrections in henyey
  agdt = agdr    ! )

  if (.not. Star%initialised) then
  dgrp = dgrp*um ! maximum allowed variation in Ln P
  dgrl = dgrl*um ! maximum allowed variation in Ln S
  endif

  !if (nwseq == 1) then
  !  restart = 0
  !else
  !  restart = nwseq
  !endif

  !if (modanf == 0) then
  !  if (idebug > 1) then
  !    write(*,*) 'initial values check and corrections'
  !  endif
  !  if (faktor /= 1.d0) then
  !    faktor = 1.d0
  !    write(*,*) 'First model: faktor set to 1'
  !  endif
  !  if (phase /= 1) then
  !    phase = 1
  !    write(*,*) 'First model: phase set to 1'
  !  endif
  !  if (irot == 1 .and. isol /= 1) then
  !    isol = 1
  !    write(*,*) 'First model: isol set to 1'
  !  endif
  !endif

  if (idebug > 1) then
    write(*,*) 'restoring network'
  endif
  call restore_network(z)

  if (ialflu == 1) then
    xnetalu = Star%xnetalu
    zabelx = zabelx-xnetalu(1)-xnetalu(2)-xnetalu(3)-xnetalu(4)
  endif

  !if (isugi >= 1 .and. nwseq  ==  1) then
  !  nsugi=mmax
  !endif

  !inum = 0
  !if (nzmod > 1) then
  !  modell = mod(nwseq,nzmod)
  !else
  !  modell = 1
  !endif
  !nzmodini = nzmod
  !if (.not. libgenec) then
  !write(*,*) "\n\n\n\n ******FFFFFFFFF", veryFirst, "\n\n", nwmd, nwseq, "\n\n\n"
  if (nwmd == 0) then
    nwmd = nwseq
  endif
  !nwseqini = nwseq

  ! Initial model
  if (modanf == 0) then
! security if initial file is missing the iprezams parameter
    if (vwant>epsilon(vwant) .and. iprezams==0) then
      write(*,*) 'VWANT/=0 --> IPREZAMS set to 1'
      iprezams=1
    endif
    if (idebug > 1) then
      write(*,*) 'Reading of initial structure'
    endif
    call copy_structure_from_genec_star(Star)
    !read(*,nml=IniStruc)
    xmini=summas
    zams_radius = 0.d0
    if (bintide) then
      period = periodini*day
    endif
    if (irot == 1 .and. isol>=1 .and. omega /= omegi(1)) then
      omegi(:) = omega
    endif
    if (alter == 0.d0) then
      firstmods = .true.
    else
      firstmods = .false.
    endif

    if (fitm /= exphi(q(1))) then
      q(1) = log10(1.d0 - fitm)
    endif

    if (ialflu == 1) then
      xf19(:)=xnetalu(1)
      xne21(:)=xnetalu(2)
      xna23(:)=xnetalu(3)
      xal26(:)=0.d0
      xal27(:)=xnetalu(4)
      xsi28(:)=xnetalu(5)
      xneut(:)=0.d0
      xprot(:)=0.d0
      xc14(:)=0.d0
      xf18(:)=0.d0
      xbid1(:)=0.d0

      xsi28(:)=0.d0
      bibib = &
              1.d0-x(1)-y(1)-y3(1)-xc12(1)-xc13(1)-xn14(1)-xn15(1)-xo16(1)-xo17(1)-xo18(1)-xne22(1)&
              -xmg24(1)-xmg25(1)-xmg26(1)-xne20(1)-xf19(1)-xne21(1)-xal27(1)-xsi28(1)-xna23(1)

      do ii=1,nbelx
        bibib=bibib-abels(ii)
      enddo

      xbid(1:m)=bibib
    else
      xf19(:)=0.d0
      xne21(:)=0.d0
      xna23(:)=0.d0
      xal26(:)=0.d0
      xal27(:)=0.d0
      xsi28(:)=0.d0
      xneut(:)=0.d0
      xprot(:)=0.d0
      xc14(:)=0.d0
      xf18(:)=0.d0
      xbid1(:)=0.d0
    endif

! for each shell give same value
    zabelx=z
    do ii=1,nbelx
      abelx(ii,:)=abels(ii)
      zabelx=zabelx-abels(ii)
    enddo
    if (ialflu == 1) then
      zabelx=zabelx-xf19(1)-xne21(1)-xna23(1)-xal27(1)
    endif

    if (isugi >= 1) then
      nsugi = m
    endif

    NPcoucheEff = npondcouche

    !call write4

    ab=ab*um
    q(:)=q(:)*um
    p(:)=p(:)*um
    t(:)=t(:)*um
    r(:)=r(:)*um
    s(:)=s(:)*um
    vp(:)=vp(:)*um
    vt(:)=vt(:)*um
    vr(:)=vr(:)*um
    vs(:)=vs(:)*um
    veryFirst = .true.

  else ! modanf > 0

    call copy_structure_from_genec_star(Star)
    !write(*,*) 's1/vs1 (helpers): ', s(1), vs(1)

    !vvsuminenv = vsuminenv

    ! if (phase > 1 .and. CorrOmega(npondcouche) < -100.d0) then
    !   NPcoucheEff = npondcoucheAdv
    ! else
    !   NPcoucheEff = npondcouche
    ! endif

    ! if (irot /= 0) then
    !   if (idebug > 1) then
    !     write(*,*) 'call momevo'
    !   endif
    !   call momevo(r,vomegi,xltod,CorrOmega,.true.)
    ! endif

    ! if (x(1) > 7.d-1) then
    !   xini = x(1)
    !   if (x(m) > (xini - 5.d-3)) then
    !     firstmods = .true.
    !   else
    !     firstmods = .false.
    !   endif
    ! else
    !   firstmods = .false.
    ! endif

    ! if (irot==1 .and. isol>=1 .and. abs(vwant)>1.0d-5) then
    !   omegi(1:m)=sqrt(xfom)*omegi(1:m)
    ! endif

    !call write4

    !write(*,*) 'call fitmshift'
    !call fitmshift

  endif ! modanf

  ! ftfp initialisation
  call initgeo

  if (ialflu==0 .and. xmini<=9.d0) then
    ichem = 1
  endif

end subroutine init_or_restore_star


integer function restore_star()
    implicit none
    restore_star = 0
end function

integer function initialize_code()
    !use WriteSaveClose, only: quitafterclosing
    use genec, only: initialise_genec
    use evol, only: input_dir
    use inputparam, only: libgenec
    use io_definitions
    implicit none

    !io_runfile = 6
    !io_logs = 6
    libgenec = .true.
    !io_logs = 6
    input_dir = "./src/GENEC/code"
    call initialise_genec()

    initialize_code = 0
end function

integer function read_genec_model(index_of_the_star, cardfilename)
    ! This should only be called if no star has been initialised yet!
    implicit none
    integer:: index_of_the_star
    character(256):: cardfilename
    read_genec_model = 0
end function

integer function cleanup_code()
    implicit none
    cleanup_code = 0
end function

! **** Parameters

integer function get_model_number(model_number)
    implicit none
    integer:: model_number
    model_number = GenecStar%nwmd
    get_model_number = 0
end function

integer function set_model_number(model_number)
    implicit none
    integer:: model_number
    GenecStar%nwmd = model_number
    set_model_number = 0
end function

integer function get_min_timestep_stop_condition(min_timestep_stop_condition)
    implicit none
    real(kindreal):: min_timestep_stop_condition
    get_min_timestep_stop_condition = 0
end function

function finalize_stellar_model()
    implicit none
    integer:: finalize_stellar_model
    !write(*,*) "copy to GenecStar"
    !call copy_to_genec_star(GenecStar)
    !write(*,*) GenecStar
    finalize_stellar_model = 0
end function

function commit_parameters()
    implicit none
    integer:: commit_parameters
    commit_parameters = 0
end function

subroutine restore_network(z)
!----------------------------------------------------------------------
  use evol, only: input_dir
  use inputparam,only: idebug,libgenec
  use abundmod,only: mbelx,abels,xlostneu,nbzel,nbael,zabelx,nbelx
  use storage, only: GenecStar

  implicit none

  integer:: i,ii,ierror
  real(kindreal),intent(in):: z
  character (256):: namenet,namereac
!----------------------------------------------------------------------
! Reading network information (elements, ...)
! first add elements to the program
  ierror = 0
  i = 1
  xlostneu = GenecStar%xlostneu
  nbzel = GenecStar%nbzel
  nbael = GenecStar%nbael
  abels = GenecStar%abels
  do while (nbael(i) > 0)
   if (GenecStar%verbose) then
     write(*,*) nbzel(i), nbael(i), abels(i)
   endif
   i = i+1
  enddo

  nbelx = i-1

  zabelx = z
  do ii = 1, nbelx
   zabelx = zabelx - abels(ii)
  enddo

  if (nbelx > mbelx) then
    write(*,*) 'nbelx= ',nbelx,' > mbelx= ',mbelx
    stop 'stop in restore_network'
  endif

  if (.false.) then  ! skip this for now - fixme later --SR
! then decide which element are followed in netnewr.f
  if (GenecStar%phase < 4) then
    namenet=trim(input_dir)//'inputs/netinit.inCNE'
    namereac=trim(input_dir)//'inputs/vit.datCNE'
  else
    namenet=trim(input_dir)//'inputs/netinit.inCNEO'
    namereac=trim(input_dir)//'inputs/vit.datCNEO'
  endif

  if (idebug > 0) then
    write(*,*) 'call readnetZA'
  endif
  call  readnetZA
  if (idebug > 0) then
    write(*,*) 'call readreac'
  endif
  call readreac
  endif ! .false.

  return

end subroutine restore_network

function commit_particles()
    use genec, only: initialise_star
    use genec, only: evolve, modell, finalise, veryFirst
    implicit none
    integer:: commit_particles
    ! makeini will actually override some things from set_defaults now! FIXME

    call copy_from_genec_star(GenecStar)
    call init_or_restore_star(GenecStar)
    GenecStar%initialised = .true.
    call evolve()
    call finalise()
    call copy_to_genec_star(GenecStar)
    veryFirst = .false.
    commit_particles = 0
end function

function delete_star(index_of_the_star)
    implicit none
    integer:: index_of_the_star
    integer:: delete_star
    delete_star = -1  ! not supported
end function

function evolve_model(end_time)
    use timestep, only: alter
    implicit none
    real(kindreal):: end_time
    integer:: evolve_model
    !stopping_condition = ""

    do while (alter < end_time)
        !if (stopping_condition == "") then
        write(*,*) "Current time: ", alter, ", evolving to: ", end_time
        evolve_model = evolve_one_step(0)
    end do
    !call copy_to_genec_star(GenecStar)
    evolve_model = 0
end function

function evolve_for(index_of_the_star, time)
    ! get current time
    ! set max time to current time plus argument
    ! evolve
    use timestep, only: alter
    implicit none
    integer:: index_of_the_star
    real(kindreal):: time, end_time
    integer:: evolve_for
    end_time = alter+time
    do while (alter < end_time)
        write(*,*) "Current time: ", alter, ", evolving to: ", end_time
        evolve_for = evolve_one_step(index_of_the_star)
    end do
    !call copy_to_genec_star(GenecStar)

    evolve_for = 0
end function

function evolve_one_step(index_of_the_star)
    use timestep, only: alter
    use WriteSaveClose, only: OpenAll
    use genec, only: evolve, modell, finalise, veryFirst
    use inputparam,only: modanf,nwseq,nzmod,end_at_phase,end_at_model
    use genec, only: n_snap
    use genec, only: m
    use abundmod, only: x
    implicit none
    integer:: index_of_the_star
    integer:: evolve_one_step
    integer:: original_nzmod

    nzmod = 1
    modell = 1
    n_snap = 0
    call evolve()
    call finalise()
    veryFirst = .false.
    call copy_to_genec_star(GenecStar)
    evolve_one_step = 0
end function

function write_genec_model()
    !use inputparam, only: modanf
    !use WriteSaveClose, only: OpenAll
    !use genec, only: finalise
    !use helpers, only: initialise_star
    implicit none
    integer:: write_genec_model
    !call finalise()
    ! modanf = 0
    !call OpenAll()
    !call initialise_star()
    !call finalise()
    !call OpenAll()
    !call initialise_star()
    write_genec_model = -1
end function

function get_age(index_of_the_star, age)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: age
    integer:: get_age
    age = GenecStar%alter
    get_age = 0
end function

function set_age(index_of_the_star, age)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: age
    integer:: set_age
    GenecStar%alter = age
    set_age = 0
end function

function get_surface_velocity(index_of_the_star, surface_velocity)
    use genec, only: vequat
    implicit none
    integer:: index_of_the_star
    real(kindreal):: surface_velocity
    integer:: get_surface_velocity
    surface_velocity = vequat
    get_surface_velocity = 0
end function

function get_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: rho_i
    integer:: get_density_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        rho_i = exp(rho(i))
    end if
    get_density_at_zone = 0
end function

function set_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: rho_i
    integer:: set_density_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        rho(i) = log(rho_i)
    end if
    set_density_at_zone = 0
end function

function get_luminosity(index_of_the_star, luminosity)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: luminosity
    integer:: get_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    luminosity = GenecStar%gls
    get_luminosity = 0
end function

function set_luminosity(index_of_the_star, luminosity)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: luminosity
    integer:: set_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    GenecStar%gls = luminosity
    set_luminosity = 0
end function

function get_radius_at_zone(index_of_the_star, zone, R_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: R_i
    integer:: get_radius_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        R_i = exp(GenecStar%r(i))  ! in cm
    end if
    get_radius_at_zone = 0
end function

function set_radius_at_zone(index_of_the_star, zone, R_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: R_i
    integer:: set_radius_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        GenecStar%r(i) = log(R_i)
    end if
    set_radius_at_zone = 0
end function

function get_luminosity_at_zone(index_of_the_star, zone, lum_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: lum_i
    integer:: get_luminosity_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        !lum_i = exp(s(zone+1))
        lum_i = exp(GenecStar%s(i)) - 1
    end if
    get_luminosity_at_zone = 0
end function

function set_luminosity_at_zone(index_of_the_star, zone, lum_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: lum_i
    integer:: set_luminosity_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        GenecStar%s(i) = log(lum_i + 1)
    end if
    set_luminosity_at_zone = 0
end function

function get_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: dq_i
    integer:: get_mass_fraction_at_zone
    i = GenecStar%m - zone
    if (i == 1) then
        dq_i = exp(GenecStar%q(i))
    else if (i <= GenecStar%m) then
        !dq_i = 1-exp(q(zone+1))
        dq_i = (exp(GenecStar%q(i)) - exp(GenecStar%q(i-1)))
    end if
    get_mass_fraction_at_zone = 0
end function

function set_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: dq_i
    integer:: set_mass_fraction_at_zone
    i = GenecStar%m - zone
    if (i == 1) then
        GenecStar%q(i) = log(dq_i)
    else if (i <= GenecStar%m) then
        !dq_i = 1-exp(q(zone+1))
        GenecStar%q(i) = log(exp(GenecStar%q(i-1)) + dq_i)  ! this won't do
    end if
    set_mass_fraction_at_zone = -1 ! This function is incomplete!
end function

function get_mass(index_of_the_star, mass)
    implicit none
    real(kindreal):: mass
    integer:: get_mass, index_of_the_star
    mass = GenecStar%gms
    get_mass = 0
end function

function set_mass(index_of_the_star, mass)
    implicit none
    integer:: set_mass, index_of_the_star
    real(kindreal):: mass
    if (.not.GenecStar%initialised) then
        GenecStar%initial_mass = mass
        set_mass = 0
    else
        write(*,*) "This function should not be called when the star is already initialised"
        set_mass = -2
    endif
end function

function get_mass_of_species(index_of_the_star, species, species_mass)
    implicit none
    integer:: index_of_the_star
    integer:: species
    real(kindreal):: species_mass
    integer:: get_mass_of_species
    get_mass_of_species = 0
    select case(species)
    case(1)
        species_mass = 1.
    case(2)
        species_mass = 3.
    case(3)
        species_mass = 4.
    case(4)
        species_mass = 12.
    case(5)
        species_mass = 13.
    case(6)
        species_mass = 14.
    case(7)
        species_mass = 15.
    case(8)
        species_mass = 16.
    case(9)
        species_mass = 17.
    case(10)
        species_mass = 18.
    case(11)
        species_mass = 20.
    case(12)
        species_mass = 22.
    case(13)
        species_mass = 24.
    case(14)
        species_mass = 25.
    case(15)
        species_mass = 26.
    case(16)
        species_mass = 14.
    case(17)
        species_mass = 18.
    case(18)
        species_mass = 19.
    case(19)
        species_mass = 21.
    case(20)
        species_mass = 23.
    case(21)
        species_mass = 26.
    case(22)
        species_mass = 27.
    case(23)
        species_mass = 28.
    case(24)
        species_mass = 1. !neut(i)
    case(25)
        species_mass = 1. !prot(i)
!    case(26)
!        species_mass = bid(i)
!    case(27)
!        species_mass = bid1(i)
    case default
        species_mass = 0.
    end select
end function

function get_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_species_at_zone
    i = GenecStar%m-zone
    if (zone <= GenecStar%m) then
      select case(species)
      case(1)
          Xj_i = GenecStar%x(i)
      case(2)
          Xj_i = GenecStar%y3(i)
      case(3)
          Xj_i = GenecStar%y(i)
      case(4)
          Xj_i = GenecStar%xc12(i)
      case(5)
          Xj_i = GenecStar%xc13(i)
      case(6)
          Xj_i = GenecStar%xn14(i)
      case(7)
          Xj_i = GenecStar%xn15(i)
      case(8)
          Xj_i = GenecStar%xo16(i)
      case(9)
          Xj_i = GenecStar%xo17(i)
      case(10)
          Xj_i = GenecStar%xo18(i)
      case(11)
          Xj_i = GenecStar%xne20(i)
      case(12)
          Xj_i = GenecStar%xne22(i)
      case(13)
          Xj_i = GenecStar%xmg24(i)
      case(14)
          Xj_i = GenecStar%xmg25(i)
      case(15)
          Xj_i = GenecStar%xmg26(i)
      case(16)
          Xj_i = GenecStar%xc14(i)
      case(17)
          Xj_i = GenecStar%xf18(i)
      case(18)
          Xj_i = GenecStar%xf19(i)
      case(19)
          Xj_i = GenecStar%xne21(i)
      case(20)
          Xj_i = GenecStar%xna23(i)
      case(21)
          Xj_i = GenecStar%xal26(i)
      case(22)
          Xj_i = GenecStar%xal27(i)
      case(23)
          Xj_i = GenecStar%xsi28(i)
      case(24)
          Xj_i = GenecStar%xneut(i)
      case(25)
          Xj_i = GenecStar%xprot(i)
      case(26)
          Xj_i = GenecStar%xbid(i)
      case(27)
          Xj_i = GenecStar%xbid1(i)
      case default
          Xj_i = 0
      end select
    end if

    get_mass_fraction_of_species_at_zone = 0
end function

function set_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_species_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        select case(species)
      case(1)
          GenecStar%x(i) = Xj_i
      case(2)
          GenecStar%y3(i) = Xj_i
      case(3)
          GenecStar%y(i) = Xj_i
      case(4)
          GenecStar%xc12(i) = Xj_i
      case(5)
          GenecStar%xc13(i) = Xj_i
      case(6)
          GenecStar%xn14(i) = Xj_i
      case(7)
          GenecStar%xn14(i) = Xj_i
      case(8)
          GenecStar%xo16(i) = Xj_i
      case(9)
          GenecStar%xo17(i) = Xj_i
      case(10)
          GenecStar%xo18(i) = Xj_i
      case(11)
          GenecStar%xne20(i) = Xj_i
      case(12)
          GenecStar%xne22(i) = Xj_i
      case(13)
          GenecStar%xmg24(i) = Xj_i
      case(14)
          GenecStar%xmg25(i) = Xj_i
      case(15)
          GenecStar%xmg26(i) = Xj_i
      case(16)
          GenecStar%xc14(i) = Xj_i
      case(17)
          GenecStar%xf18(i) = Xj_i
      case(18)
          GenecStar%xf19(i) = Xj_i
      case(19)
          GenecStar%xne21(i) = Xj_i
      case(20)
          GenecStar%xna23(i) = Xj_i
      case(21)
          GenecStar%xal26(i) = Xj_i
      case(22)
          GenecStar%xal27(i) = Xj_i
      case(23)
          GenecStar%xsi28(i) = Xj_i
      case(24)
          GenecStar%xneut(i) = Xj_i
      case(25)
          GenecStar%xprot(i) = Xj_i
      case(26)
          GenecStar%xbid(i) = Xj_i
      case(27)
          GenecStar%xbid1(i) = Xj_i
      !case default
      end select
    end if

    set_mass_fraction_of_species_at_zone = 0
end function

function get_metallicity(metallicity)
    implicit none
    real(kindreal):: metallicity
    integer:: get_metallicity
    metallicity = GenecStar%initial_metallicity
    get_metallicity = 0
end function

function set_metallicity(metallicity)
    implicit none
    real(kindreal):: metallicity
    integer:: set_metallicity
    GenecStar%initial_metallicity = metallicity
    set_metallicity = 0
end function

function get_mu_at_zone(index_of_the_star, zone, mu_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: mu_i, X, Y3, Y
    integer:: err
    integer:: get_mu_at_zone
    i = (GenecStar%m - zone)
    if (zone <= GenecStar%m) then
        mu_i = 0.0d0
        err = get_mass_fraction_of_species_at_zone(index_of_the_star, 1, zone, X)
        err = get_mass_fraction_of_species_at_zone(index_of_the_star, 2, zone, Y3)
        err = get_mass_fraction_of_species_at_zone(index_of_the_star, 3, zone, Y)

        mu_i = (4.0 / (2.0 + 6 * X - (Y+Y3)))
    end if
    get_mu_at_zone = 0
end function

function get_name_of_species(index_of_the_star, species, species_name)
    implicit none
    integer:: index_of_the_star
    integer:: species
    character(len = 6):: species_name
    integer:: get_name_of_species

    character(len = 6), dimension(27):: species_names

    species_names(1)  = 'h'
    species_names(2)  = 'he3'
    species_names(3)  = 'he'
    species_names(4)  = 'c12'
    species_names(5)  = 'c13'
    species_names(6)  = 'n14'
    species_names(7)  = 'n15'
    species_names(8)  = 'o16'
    species_names(9)  = 'o17'
    species_names(10) = 'o18'
    species_names(11) = 'ne20'
    species_names(12) = 'ne22'
    species_names(13) = 'mg24'
    species_names(14) = 'mg25'
    species_names(15) = 'mg26'
    species_names(16) = 'c14'
    species_names(17) = 'f18'
    species_names(18) = 'f19'
    species_names(19) = 'ne21'
    species_names(20) = 'na23'
    species_names(21) = 'al26'
    species_names(22) = 'al27'
    species_names(23) = 'si28'
    species_names(24) = 'neut'
    species_names(25) = 'prot'
    species_names(26) = 'bid'
    species_names(27) = 'bid1'
    species_name = species_names(species)
    !x: H1
    !y: He4
    !y3: He3
    !xneut: neutron
    !xprot: proton
    !xc12: C12
    !xXYY: XYY - X element YY mass number
    !xXXYY: as above
    get_name_of_species = 0
end function

function get_number_of_particles(n)
    implicit none
    integer:: n
    integer:: get_number_of_particles
    n = 1
    get_number_of_particles = 0
end function

function get_number_of_species(index_of_the_star, n_species)
    implicit none
    integer:: index_of_the_star
    integer:: n_species
    integer:: get_number_of_species
    if (GenecStar%ialflu==1) then
       n_species = 27
    else
        n_species = 15
    end if
    get_number_of_species = 0
end function

function get_firstlast_species_number(first, last)
    implicit none
    !integer:: index_of_the_star
    integer:: first, last
    integer:: get_firstlast_species_number
    first = 1
    if (GenecStar%ialflu==1) then
        last = 27
    else
        last = 15
    end if
    get_firstlast_species_number = 0
end function

function get_number_of_zones(index_of_the_star, n_zones)
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: get_number_of_zones
    n_zones = GenecStar%m
    get_number_of_zones = 0
end function

function set_number_of_zones(index_of_the_star, n_zones)
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: set_number_of_zones
    GenecStar%m = n_zones
    set_number_of_zones = 0
end function

function get_firstlast_zone(index_of_the_star, first, last)
    implicit none
    integer:: index_of_the_star
    integer:: first, last
    integer:: get_firstlast_zone
    first = 0
    last = GenecStar%m-1
    get_firstlast_zone = 0
end function

function get_pressure_at_zone(index_of_the_star, zone, P_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: P_i
    integer:: get_pressure_at_zone
    if (zone <= GenecStar%m) then
        i = GenecStar%m - zone
        P_i = exp(GenecStar%p(i))
    end if
    get_pressure_at_zone = 0
end function

!function set_pressure_at_zone(index_of_the_star, zone, P_i)
!    use strucmod, only: p, m
!    implicit none
!    integer:: index_of_the_star
!    integer:: zone, i
!    real(kindreal):: P_i
!    integer:: set_pressure_at_zone
!    if (zone <= m) then
!        i = m - zone
!        p(i) = log(P_i)
!    end if
!    set_pressure_at_zone = 0
!end function

function get_radius(index_of_the_star, am_radius)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: am_radius
    integer:: get_radius
    am_radius = 10**GenecStar%radius
    !radius = exp(r(1))  ! in cm
    get_radius = 0
end function

function set_radius(index_of_the_star, am_radius)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: am_radius
    integer:: set_radius
    GenecStar%radius = log10(am_radius)
    set_radius = 0
end function

function get_stellar_type(index_of_the_star, stellar_type)
    implicit none
    integer:: index_of_the_star
    integer:: stellar_type
    integer:: get_stellar_type
    get_stellar_type = -1
end function

function get_temperature(index_of_the_star, temperature)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: temperature
    integer:: get_temperature
    temperature = GenecStar%teff
    get_temperature = 0
end function

function get_temperature_at_zone(index_of_the_star, zone, T_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: T_i
    integer:: get_temperature_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        T_i = exp(GenecStar%t(i))
    else
        get_temperature_at_zone = -2
        return
    end if
    get_temperature_at_zone = 0
end function

function get_time_step(index_of_the_star, time_step)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: time_step
    integer:: get_time_step
    get_time_step = get_dzeit(index_of_the_star, time_step)
    !time_step = GenecStar%dzeitj
    !time_step = GenecStar%dzeit
    !get_time_step = 0
end function

function set_time_step(index_of_the_star, time_step)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: time_step
    integer:: set_time_step
    set_time_step = set_dzeit(index_of_the_star, time_step)
end function

function get_time(time)
    implicit none
    real(kindreal):: time
    integer:: get_time
    time = GenecStar%alter
    get_time = 0
end function

function new_particle(index_of_the_star, initial_mass, initial_metallicity, zams_velocity, star_name)
    use makeini, only: make_initial_star
    implicit none
    integer:: index_of_the_star, key
    real(kindreal):: initial_mass, initial_metallicity, zams_velocity
    integer:: new_particle
    character(len=12):: star_name
    number_of_stars = number_of_stars + 1
    GenecStar%star_name = star_name
    GenecStar%index_of_the_star = number_of_stars
    GenecStar%initial_mass = initial_mass
    GenecStar%initial_metallicity = initial_metallicity
    GenecStar%zams_velocity = zams_velocity
    GenecStar%idefaut = 1
    index_of_the_star = GenecStar%index_of_the_star

    call make_initial_star()
    GenecStar%nzmod = 1
    GenecStar%modell = 1
    GenecStar%n_snap = 0
    call copy_from_genec_star(GenecStar)
    call copy_to_genec_star(GenecStar)
    new_particle = 0
end function

!function new_stellar_model(&
!      integer_of_the_star,&
!      initialised, starname, nwseq, modanf, nzmod, end_at_phase, end_at_model, irot, isol, imagn, ialflu, ianiso, ipop3,&
!      ibasnet, phase, var_rates, bintide, binm2, periodini, const_per, iprezams, zinit, zsol, z, iopac, ikappa, idiff, iadvec,&
!      istati, icoeff, fenerg, richac, igamma, frein, K_Kawaler, Omega_saturation, rapcrilim, vwant, xfom, omega, xdial, idialo,&
!      idialu, Add_Flux, diff_only, B_initial, add_diff, n_mag, alpha_F, nsmooth, qminsmooth, imloss, fmlos, ifitm, fitm, fitmi,&
!      deltal, deltat, nndr, RSG_Mdot, SupraEddMdot, Be_mdotfrac, start_mdot, iledou, idifcon, iover, elph, my, dovhp, iunder,&
!      dunder, gkorm, alph, agdr, faktor, dgrp, dgrl, dgry, dgrc, dgro, dgr20, nbchx, nrband, xcn, islow, icncst, tauH_fit,&
!      display_plot, iauto, iprn, iout, itmin, xyfiles, idebug, itests, verbose, stop_deg, n_snap,&
!      m,gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,dm_lost,&
!      mbelx,xtefflast,xllast,xrholast,xclast,xtclast,inum,nsugi,period,r_core,vna,vnr,&
!      q,p,t,r,s,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,&
!      xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
!      vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,vxal26,vxneut,vxprot,vomegi,&
!      vxbid,vxbid1,&
!      abelx,vabelx&
!      )
!    implicit none
!    integer:: integer_of_the_star
!    logical, intent(in):: &
!            initialised,var_rates,bintide,const_per,Add_Flux,diff_only,qminsmooth,SupraEddMdot,display_plot,xyfiles,verbose,&
!            stop_deg
!    integer, intent(in):: &
!            nwseq,modanf,nzmod,end_at_phase,end_at_model,irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,iprezams,iopac,ikappa,&
!            idiff,iadvec,istati,icoeff,igamma,idialo,idialu,n_mag,nsmooth,imloss,ifitm,nndr,RSG_Mdot,iledou,idifcon,iover,my,&
!            iunder,nbchx,nrband,islow,icncst,tauH_fit,iauto,iprn,iout,itmin,idebug,itests,n_snap
!    real(kindreal), intent(in):: &
!            binm2,periodini,zinit,zsol,z,fenerg,richac,frein,K_Kawaler,Omega_saturation,rapcrilim,vwant,xfom,omega,xdial,&
!            B_initial,add_diff,alpha_F,fmlos,fitm,fitmi,deltal,deltat,Be_mdotfrac,start_mdot,elph,dovhp,dunder,gkorm,alph,&
!            agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,xcn
!    character(len=200), intent(in):: &
!            starname
!    integer, intent(in) :: m
!    real(kindreal), intent(in) :: &
!            gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,&
!            dm_lost
!    !real(kindreal), dimension(ldi) :: &
!    integer, intent(in):: mbelx
!    integer, intent(in):: inum,nsugi
!    real(kindreal), intent(in):: &
!            xtefflast,xllast,xrholast,xclast,xtclast,period,r_core,vna,vnr
!    real(kindreal), dimension(m), intent(in):: &
!            q,p,t,r,s,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,&
!            xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
!            vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,vxal26,vxneut,vxprot,vomegi,&
!            vxbid,vxbid1
!    real(kindreal), dimension(mbelx,m), intent(in):: &
!            abelx,vabelx
!    integer:: new_stellar_model
!
!    GenecStar%initialised      = initialised
!    GenecStar%starname         = starname
!    !GenecStar%nwmd             = nwmd
!    GenecStar%nwseq            = nwseq            
!    GenecStar%modanf           = modanf
!    GenecStar%nzmod            = nzmod
!    GenecStar%end_at_phase     = end_at_phase
!    GenecStar%end_at_model     = end_at_model
!    GenecStar%irot             = irot
!    GenecStar%isol             = isol
!    GenecStar%imagn            = imagn
!    GenecStar%ialflu           = ialflu
!    GenecStar%ianiso           = ianiso
!    GenecStar%ipop3            = ipop3
!    GenecStar%ibasnet          = ibasnet
!    GenecStar%phase            = phase
!    GenecStar%iprezams         = iprezams
!    GenecStar%var_rates        = var_rates
!    GenecStar%bintide          = bintide
!    GenecStar%binm2            = binm2
!    GenecStar%periodini        = periodini
!    GenecStar%const_per        = const_per
!    GenecStar%iopac            = iopac
!    GenecStar%ikappa           = ikappa
!    GenecStar%zinit            = zinit
!    GenecStar%zsol             = zsol
!    GenecStar%z                = z
!    GenecStar%idiff            = idiff
!    GenecStar%iadvec           = iadvec
!    GenecStar%istati           = istati
!    GenecStar%icoeff           = icoeff
!    GenecStar%igamma           = igamma
!    GenecStar%idialo           = idialo
!    GenecStar%idialu           = idialu
!    GenecStar%n_mag            = n_mag
!    GenecStar%nsmooth          = nsmooth
!    GenecStar%fenerg           = fenerg
!    GenecStar%richac           = richac
!    GenecStar%frein            = frein
!    GenecStar%K_Kawaler        = K_Kawaler
!    GenecStar%Omega_saturation = Omega_saturation
!    GenecStar%rapcrilim        = rapcrilim
!    GenecStar%vwant            = vwant
!    GenecStar%xfom             = xfom
!    GenecStar%omega            = omega
!    GenecStar%xdial            = xdial
!    GenecStar%B_initial        = B_initial
!    GenecStar%add_diff         = add_diff
!    GenecStar%alpha_F          = alpha_F
!    GenecStar%Add_Flux         = Add_Flux
!    GenecStar%diff_only        = diff_only
!    GenecStar%qminsmooth       = qminsmooth
!    GenecStar%imloss           = imloss
!    GenecStar%ifitm            = ifitm
!    GenecStar%nndr             = nndr
!    GenecStar%RSG_Mdot         = RSG_Mdot
!    GenecStar%fmlos            = fmlos
!    GenecStar%fitm             = fitm
!    GenecStar%fitmi            = fitmi
!    GenecStar%deltal           = deltal
!    GenecStar%deltat           = deltat
!    GenecStar%Be_mdotfrac      = Be_mdotfrac
!    GenecStar%start_mdot       = start_mdot
!    GenecStar%SupraEddMdot     = SupraEddMdot
!    GenecStar%iledou           = iledou
!    GenecStar%idifcon          = idifcon
!    GenecStar%my               = my
!    GenecStar%iover            = iover
!    GenecStar%iunder           = iunder
!    GenecStar%elph             = elph
!    GenecStar%dovhp            = dovhp
!    GenecStar%dunder           = dunder
!    GenecStar%nbchx            = nbchx
!    GenecStar%nrband           = nrband
!    GenecStar%gkorm            = gkorm
!    GenecStar%alph             = alph
!    GenecStar%agdr             = agdr
!    GenecStar%faktor           = faktor
!    GenecStar%dgrp             = dgrp
!    GenecStar%dgrl             = dgrl
!    GenecStar%dgry             = dgry
!    GenecStar%dgrc             = dgrc
!    GenecStar%dgro             = dgro
!    GenecStar%dgr20            = dgr20
!    GenecStar%islow            = islow
!    GenecStar%icncst           = icncst
!    GenecStar%tauH_fit         = tauH_fit
!    GenecStar%xcn              = xcn
!    GenecStar%iauto            = iauto
!    GenecStar%iprn             = iprn
!    GenecStar%iout             = iout
!    GenecStar%itmin            = itmin
!    GenecStar%idebug           = idebug
!    GenecStar%itests           = itests
!    GenecStar%n_snap           = n_snap
!    GenecStar%display_plot     = display_plot
!    GenecStar%xyfiles          = xyfiles
!    GenecStar%verbose          = verbose
!    GenecStar%stop_deg         = stop_deg
!
!    GenecStar%m                = m
!    GenecStar%gms              = gms
!    GenecStar%alter            = alter
!    GenecStar%gls              = gls
!    GenecStar%teff             = teff
!    GenecStar%glsv             = glsv
!    GenecStar%teffv            = teffv
!    GenecStar%dzeitj           = dzeitj
!    GenecStar%dzeit            = dzeit
!    GenecStar%dzeitv           = dzeitv
!    GenecStar%summas           = summas
!    GenecStar%xmini            = xmini
!    GenecStar%ab               = ab
!    GenecStar%dm_lost          = dm_lost
!
!    GenecStar%mbelx            = mbelx
!    GenecStar%xteffprev        = xteffprev
!    GenecStar%xlprev           = xlprev   
!    GenecStar%xrhoprev         = xrhoprev 
!    GenecStar%xcprev           = xcprev   
!    GenecStar%xtcprev          = xtcprev  
!    GenecStar%inum             = inum     
!    GenecStar%nsugi            = nsugi    
!    GenecStar%period           = period   
!    GenecStar%r_core           = r_core   
!    GenecStar%vna              = vna      
!    GenecStar%vnr              = vnr      
!
!    GenecStar%q                = q
!    GenecStar%p                = p
!    GenecStar%t                = t
!    GenecStar%r                = r
!    GenecStar%s                = s
!    GenecStar%x                = x
!    GenecStar%y                = y
!    GenecStar%xc12             = xc12
!    GenecStar%vp               = vp
!    GenecStar%vt               = vt
!    GenecStar%vr               = vr
!    GenecStar%vs               = vs
!    GenecStar%xo16             = xo16
!    GenecStar%vx               = vx
!    GenecStar%vy               = vy
!    GenecStar%vxc12            = vxc12
!    GenecStar%vxo16            = vxo16
!    GenecStar%y3               = y3
!    GenecStar%xc13             = xc13
!    GenecStar%xn14             = xn14
!    GenecStar%xn15             = xn15
!    GenecStar%xo17             = xo17
!    GenecStar%xo18             = xo18
!    GenecStar%vy3              = vy3
!    GenecStar%vxc13            = vxc13
!    GenecStar%vxn14            = vxn14
!    GenecStar%vxn15            = vxn15
!    GenecStar%vxo17            = vxo17
!    GenecStar%vxo18            = vxo18
!    GenecStar%xne20            = xne20
!    GenecStar%xne22            = xne22
!    GenecStar%xmg24            = xmg24
!    GenecStar%xmg25            = xmg25
!    GenecStar%xmg26            = xmg26
!    GenecStar%vxne20           = vxne20
!    GenecStar%vxne22           = vxne22
!    GenecStar%vxmg24           = vxmg24
!    GenecStar%vxmg25           = vxmg25
!    GenecStar%vxmg26           = vxmg26
!    GenecStar%omegi            = omegi
!    GenecStar%vomegi           = vomegi
!    GenecStar%xf19             = xf19
!    GenecStar%xne21            = xne21
!    GenecStar%xna23            = xna23
!    GenecStar%xal26            = xal26
!    GenecStar%xal27            = xal27
!    GenecStar%xsi28            = xsi28
!    GenecStar%vxf19            = vxf19
!    GenecStar%vxne21           = vxne21
!    GenecStar%vxna23           = vxna23
!    GenecStar%vxal26           = vxal26
!    GenecStar%vxal27           = vxal27
!    GenecStar%vxsi28           = vxsi28
!    GenecStar%xneut            = xneut
!    GenecStar%xprot            = xprot
!    GenecStar%xc14             = xc14
!    GenecStar%xf18             = xf18
!    GenecStar%xbid             = xbid
!    GenecStar%xbid1            = xbid1
!    GenecStar%vxneut           = vxneut
!    GenecStar%vxprot           = vxprot
!    GenecStar%vxc14            = vxc14
!    GenecStar%vxf18            = vxf18
!    GenecStar%vxbid            = vxbid
!    GenecStar%vxbid1           = vxbid1
!
!    GenecStar%abelx            = abelx
!    GenecStar%vabelx           = vabelx
!
!    new_stellar_model = 0
!end function

function recommit_parameters()
    implicit none
    integer:: recommit_parameters
    recommit_parameters = 0
end function

function recommit_particles()
    !use genec, only: initialise_star
    implicit none
    integer:: recommit_particles
    !write(*,*) "copy from GenecStar"
    !call copy_from_genec_star(GenecStar)

    ! "read" GENEC star namelists
    !write(*,*) "call copy_namelists_from_genec_star(GenecStar)"
    !call copy_namelists_from_genec_star(GenecStar)
    !! Initialise values not in GENEC star
    !! but this also resets some values that are read, so...
    !write(*,*) 'call initialise_star()'
    !call initialise_star()
    !! "read" the GENEC star structure
    !write(*,*) 'call copy_structure_from_genec_star(GenecStar)'
    !call copy_structure_from_genec_star(GenecStar)
    !! and then copy back - some things may have changed?
    !write(*,*) 'call copy_to_genec_star(GenecStar)'
    !call copy_to_genec_star(GenecStar)
    write(*,*) "recommit particles after a change"
    call copy_from_genec_star(GenecStar)
    recommit_particles = 0
end function

function set_genec_path(path)
    use evol, only: input_dir
    implicit none
    character(len=200):: path
    integer:: set_genec_path

    input_dir = path
    set_genec_path = 0
end function

function set_temperature_at_zone(index_of_the_star, zone, T_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: T_i
    integer:: set_temperature_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        GenecStar%t(i) = log(T_i)
    end if
    
    set_temperature_at_zone = 0
end function





!!!!!!!!!
integer function get_modell(index_of_the_particle, modell)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: modell
    modell = GenecStar%modell
    get_modell = 0
end function get_modell

integer function set_modell(index_of_the_particle, modell)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: modell
    GenecStar%modell = modell
    set_modell = 0
end function set_modell

integer function get_veryFirst(index_of_the_particle, veryFirst)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: veryFirst
    veryFirst = GenecStar%veryFirst
    get_veryFirst = 0
end function get_veryFirst

integer function set_veryFirst(index_of_the_particle, veryFirst)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: veryFirst
    GenecStar%veryFirst = veryFirst
    set_veryFirst = 0
end function set_veryFirst

integer function get_initialised(index_of_the_particle, initialised)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: initialised
    initialised = GenecStar%initialised
    get_initialised = 0
end function get_initialised

integer function set_initialised(index_of_the_particle, initialised)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: initialised
    GenecStar%initialised = initialised
    set_initialised = 0
end function set_initialised

integer function get_star_name(index_of_the_particle, star_name)
    implicit none
    integer, intent(in):: index_of_the_particle
    character(256), intent(out):: star_name
    star_name = GenecStar%star_name
    get_star_name = 0
end function get_star_name

integer function set_star_name(index_of_the_particle, star_name)
    implicit none
    integer, intent(in):: index_of_the_particle
    character(256), intent(in):: star_name
    GenecStar%star_name = star_name
    set_star_name = 0
end function set_star_name

integer function get_nwmd(index_of_the_particle, nwmd)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nwmd
    nwmd = GenecStar%nwmd
    get_nwmd = 0
end function get_nwmd

integer function set_nwmd(index_of_the_particle, nwmd)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nwmd
    GenecStar%nwmd = nwmd
    set_nwmd = 0
end function set_nwmd

integer function get_nwseq(index_of_the_particle, nwseq)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nwseq
    nwseq = GenecStar%nwseq
    get_nwseq = 0
end function get_nwseq

integer function set_nwseq(index_of_the_particle, nwseq)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nwseq
    GenecStar%nwseq = nwseq
    set_nwseq = 0
end function set_nwseq

integer function get_modanf(index_of_the_particle, modanf)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: modanf
    modanf = GenecStar%modanf
    get_modanf = 0
end function get_modanf

integer function set_modanf(index_of_the_particle, modanf)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: modanf
    GenecStar%modanf = modanf
    set_modanf = 0
end function set_modanf

integer function get_nzmod(index_of_the_particle, nzmod)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nzmod
    nzmod = GenecStar%nzmod
    get_nzmod = 0
end function get_nzmod

integer function set_nzmod(index_of_the_particle, nzmod)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nzmod
    GenecStar%nzmod = nzmod
    set_nzmod = 0
end function set_nzmod

integer function get_end_at_phase(index_of_the_particle, end_at_phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: end_at_phase
    end_at_phase = GenecStar%end_at_phase
    get_end_at_phase = 0
end function get_end_at_phase

integer function set_end_at_phase(index_of_the_particle, end_at_phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: end_at_phase
    GenecStar%end_at_phase = end_at_phase
    set_end_at_phase = 0
end function set_end_at_phase

integer function get_end_at_model(index_of_the_particle, end_at_model)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: end_at_model
    end_at_model = GenecStar%end_at_model
    get_end_at_model = 0
end function get_end_at_model

integer function set_end_at_model(index_of_the_particle, end_at_model)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: end_at_model
    GenecStar%end_at_model = end_at_model
    set_end_at_model = 0
end function set_end_at_model

integer function get_irot(index_of_the_particle, irot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: irot
    irot = GenecStar%irot
    get_irot = 0
end function get_irot

integer function set_irot(index_of_the_particle, irot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: irot
    GenecStar%irot = irot
    set_irot = 0
end function set_irot

integer function get_isol(index_of_the_particle, isol)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: isol
    isol = GenecStar%isol
    get_isol = 0
end function get_isol

integer function set_isol(index_of_the_particle, isol)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: isol
    GenecStar%isol = isol
    set_isol = 0
end function set_isol

integer function get_imagn(index_of_the_particle, imagn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: imagn
    imagn = GenecStar%imagn
    get_imagn = 0
end function get_imagn

integer function set_imagn(index_of_the_particle, imagn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: imagn
    GenecStar%imagn = imagn
    set_imagn = 0
end function set_imagn

integer function get_ialflu(index_of_the_particle, ialflu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ialflu
    ialflu = GenecStar%ialflu
    get_ialflu = 0
end function get_ialflu

integer function set_ialflu(index_of_the_particle, ialflu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ialflu
    GenecStar%ialflu = ialflu
    set_ialflu = 0
end function set_ialflu

integer function get_ianiso(index_of_the_particle, ianiso)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ianiso
    ianiso = GenecStar%ianiso
    get_ianiso = 0
end function get_ianiso

integer function set_ianiso(index_of_the_particle, ianiso)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ianiso
    GenecStar%ianiso = ianiso
    set_ianiso = 0
end function set_ianiso

integer function get_ipop3(index_of_the_particle, ipop3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ipop3
    ipop3 = GenecStar%ipop3
    get_ipop3 = 0
end function get_ipop3

integer function set_ipop3(index_of_the_particle, ipop3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ipop3
    GenecStar%ipop3 = ipop3
    set_ipop3 = 0
end function set_ipop3

integer function get_ibasnet(index_of_the_particle, ibasnet)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ibasnet
    ibasnet = GenecStar%ibasnet
    get_ibasnet = 0
end function get_ibasnet

integer function set_ibasnet(index_of_the_particle, ibasnet)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ibasnet
    GenecStar%ibasnet = ibasnet
    set_ibasnet = 0
end function set_ibasnet

integer function get_phase(index_of_the_particle, phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: phase
    phase = GenecStar%phase
    get_phase = 0
end function get_phase

integer function set_phase(index_of_the_particle, phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: phase
    GenecStar%phase = phase
    set_phase = 0
end function set_phase

integer function get_var_rates(index_of_the_particle, var_rates)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: var_rates
    var_rates = GenecStar%var_rates
    get_var_rates = 0
end function get_var_rates

integer function set_var_rates(index_of_the_particle, var_rates)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: var_rates
    GenecStar%var_rates = var_rates
    set_var_rates = 0
end function set_var_rates

integer function get_bintide(index_of_the_particle, bintide)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: bintide
    bintide = GenecStar%bintide
    get_bintide = 0
end function get_bintide

integer function set_bintide(index_of_the_particle, bintide)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: bintide
    GenecStar%bintide = bintide
    set_bintide = 0
end function set_bintide

integer function get_binm2(index_of_the_particle, binm2)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: binm2
    binm2 = GenecStar%binm2
    get_binm2 = 0
end function get_binm2

integer function set_binm2(index_of_the_particle, binm2)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: binm2
    GenecStar%binm2 = binm2
    set_binm2 = 0
end function set_binm2

integer function get_periodini(index_of_the_particle, periodini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: periodini
    periodini = GenecStar%periodini
    get_periodini = 0
end function get_periodini

integer function set_periodini(index_of_the_particle, periodini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: periodini
    GenecStar%periodini = periodini
    set_periodini = 0
end function set_periodini

integer function get_const_per(index_of_the_particle, const_per)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: const_per
    const_per = GenecStar%const_per
    get_const_per = 0
end function get_const_per

integer function set_const_per(index_of_the_particle, const_per)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: const_per
    GenecStar%const_per = const_per
    set_const_per = 0
end function set_const_per

integer function get_iprezams(index_of_the_particle, iprezams)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iprezams
    iprezams = GenecStar%iprezams
    get_iprezams = 0
end function get_iprezams

integer function set_iprezams(index_of_the_particle, iprezams)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iprezams
    GenecStar%iprezams = iprezams
    set_iprezams = 0
end function set_iprezams

integer function get_initial_metallicity(index_of_the_particle, initial_metallicity)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: initial_metallicity
    initial_metallicity = GenecStar%initial_metallicity
    get_initial_metallicity = 0
end function get_initial_metallicity

integer function set_initial_metallicity(index_of_the_particle, initial_metallicity)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: initial_metallicity
    GenecStar%initial_metallicity = initial_metallicity
    set_initial_metallicity = 0
end function set_initial_metallicity

integer function get_zsol(index_of_the_particle, zsol)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: zsol
    zsol = GenecStar%zsol
    get_zsol = 0
end function get_zsol

integer function set_zsol(index_of_the_particle, zsol)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: zsol
    GenecStar%zsol = zsol
    set_zsol = 0
end function set_zsol

integer function get_z(index_of_the_particle, z)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: z
    z = GenecStar%z
    get_z = 0
end function get_z

integer function set_z(index_of_the_particle, z)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: z
    GenecStar%z = z
    set_z = 0
end function set_z

integer function get_iopac(index_of_the_particle, iopac)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iopac
    iopac = GenecStar%iopac
    get_iopac = 0
end function get_iopac

integer function set_iopac(index_of_the_particle, iopac)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iopac
    GenecStar%iopac = iopac
    set_iopac = 0
end function set_iopac

integer function get_ikappa(index_of_the_particle, ikappa)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ikappa
    ikappa = GenecStar%ikappa
    get_ikappa = 0
end function get_ikappa

integer function set_ikappa(index_of_the_particle, ikappa)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ikappa
    GenecStar%ikappa = ikappa
    set_ikappa = 0
end function set_ikappa

integer function get_idiff(index_of_the_particle, idiff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idiff
    idiff = GenecStar%idiff
    get_idiff = 0
end function get_idiff

integer function set_idiff(index_of_the_particle, idiff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idiff
    GenecStar%idiff = idiff
    set_idiff = 0
end function set_idiff

integer function get_iadvec(index_of_the_particle, iadvec)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iadvec
    iadvec = GenecStar%iadvec
    get_iadvec = 0
end function get_iadvec

integer function set_iadvec(index_of_the_particle, iadvec)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iadvec
    GenecStar%iadvec = iadvec
    set_iadvec = 0
end function set_iadvec

integer function get_istati(index_of_the_particle, istati)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: istati
    istati = GenecStar%istati
    get_istati = 0
end function get_istati

integer function set_istati(index_of_the_particle, istati)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: istati
    GenecStar%istati = istati
    set_istati = 0
end function set_istati

integer function get_icoeff(index_of_the_particle, icoeff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: icoeff
    icoeff = GenecStar%icoeff
    get_icoeff = 0
end function get_icoeff

integer function set_icoeff(index_of_the_particle, icoeff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: icoeff
    GenecStar%icoeff = icoeff
    set_icoeff = 0
end function set_icoeff

integer function get_fenerg(index_of_the_particle, fenerg)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fenerg
    fenerg = GenecStar%fenerg
    get_fenerg = 0
end function get_fenerg

integer function set_fenerg(index_of_the_particle, fenerg)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fenerg
    GenecStar%fenerg = fenerg
    set_fenerg = 0
end function set_fenerg

integer function get_richac(index_of_the_particle, richac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: richac
    richac = GenecStar%richac
    get_richac = 0
end function get_richac

integer function set_richac(index_of_the_particle, richac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: richac
    GenecStar%richac = richac
    set_richac = 0
end function set_richac

integer function get_igamma(index_of_the_particle, igamma)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: igamma
    igamma = GenecStar%igamma
    get_igamma = 0
end function get_igamma

integer function set_igamma(index_of_the_particle, igamma)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: igamma
    GenecStar%igamma = igamma
    set_igamma = 0
end function set_igamma

integer function get_frein(index_of_the_particle, frein)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: frein
    frein = GenecStar%frein
    get_frein = 0
end function get_frein

integer function set_frein(index_of_the_particle, frein)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: frein
    GenecStar%frein = frein
    set_frein = 0
end function set_frein

integer function get_K_Kawaler(index_of_the_particle, K_Kawaler)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: K_Kawaler
    K_Kawaler = GenecStar%K_Kawaler
    get_K_Kawaler = 0
end function get_K_Kawaler

integer function set_K_Kawaler(index_of_the_particle, K_Kawaler)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: K_Kawaler
    GenecStar%K_Kawaler = K_Kawaler
    set_K_Kawaler = 0
end function set_K_Kawaler

integer function get_Omega_saturation(index_of_the_particle, Omega_saturation)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: Omega_saturation
    Omega_saturation = GenecStar%Omega_saturation
    get_Omega_saturation = 0
end function get_Omega_saturation

integer function set_Omega_saturation(index_of_the_particle, Omega_saturation)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: Omega_saturation
    GenecStar%Omega_saturation = Omega_saturation
    set_Omega_saturation = 0
end function set_Omega_saturation

integer function get_rapcrilim(index_of_the_particle, rapcrilim)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rapcrilim
    rapcrilim = GenecStar%rapcrilim
    get_rapcrilim = 0
end function get_rapcrilim

integer function set_rapcrilim(index_of_the_particle, rapcrilim)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rapcrilim
    GenecStar%rapcrilim = rapcrilim
    set_rapcrilim = 0
end function set_rapcrilim

integer function get_zams_velocity(index_of_the_particle, zams_velocity)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: zams_velocity
    zams_velocity = GenecStar%zams_velocity
    get_zams_velocity = 0
end function get_zams_velocity

integer function set_zams_velocity(index_of_the_particle, zams_velocity)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: zams_velocity
    if (.not.GenecStar%initialised) then
        GenecStar%zams_velocity = zams_velocity
        set_zams_velocity = 0
    else
        write(*,*) "This function should not be called when the star is already initialised"
        set_zams_velocity = -2
    endif
end function set_zams_velocity

integer function get_xfom(index_of_the_particle, xfom)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xfom
    xfom = GenecStar%xfom
    get_xfom = 0
end function get_xfom

integer function set_xfom(index_of_the_particle, xfom)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xfom
    GenecStar%xfom = xfom
    set_xfom = 0
end function set_xfom

integer function get_omega(index_of_the_particle, omega)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: omega
    omega = GenecStar%omega
    get_omega = 0
end function get_omega

integer function set_omega(index_of_the_particle, omega)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: omega
    GenecStar%omega = omega
    set_omega = 0
end function set_omega

integer function get_xdial(index_of_the_particle, xdial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xdial
    xdial = GenecStar%xdial
    get_xdial = 0
end function get_xdial

integer function set_xdial(index_of_the_particle, xdial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xdial
    GenecStar%xdial = xdial
    set_xdial = 0
end function set_xdial

integer function get_idialo(index_of_the_particle, idialo)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idialo
    idialo = GenecStar%idialo
    get_idialo = 0
end function get_idialo

integer function set_idialo(index_of_the_particle, idialo)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idialo
    GenecStar%idialo = idialo
    set_idialo = 0
end function set_idialo

integer function get_idialu(index_of_the_particle, idialu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idialu
    idialu = GenecStar%idialu
    get_idialu = 0
end function get_idialu

integer function set_idialu(index_of_the_particle, idialu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idialu
    GenecStar%idialu = idialu
    set_idialu = 0
end function set_idialu

integer function get_Add_Flux(index_of_the_particle, Add_Flux)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: Add_Flux
    Add_Flux = GenecStar%Add_Flux
    get_Add_Flux = 0
end function get_Add_Flux

integer function set_Add_Flux(index_of_the_particle, Add_Flux)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: Add_Flux
    GenecStar%Add_Flux = Add_Flux
    set_Add_Flux = 0
end function set_Add_Flux

integer function get_diff_only(index_of_the_particle, diff_only)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: diff_only
    diff_only = GenecStar%diff_only
    get_diff_only = 0
end function get_diff_only

integer function set_diff_only(index_of_the_particle, diff_only)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: diff_only
    GenecStar%diff_only = diff_only
    set_diff_only = 0
end function set_diff_only

integer function get_B_initial(index_of_the_particle, B_initial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: B_initial
    B_initial = GenecStar%B_initial
    get_B_initial = 0
end function get_B_initial

integer function set_B_initial(index_of_the_particle, B_initial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: B_initial
    GenecStar%B_initial = B_initial
    set_B_initial = 0
end function set_B_initial

integer function get_add_diff(index_of_the_particle, add_diff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: add_diff
    add_diff = GenecStar%add_diff
    get_add_diff = 0
end function get_add_diff

integer function set_add_diff(index_of_the_particle, add_diff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: add_diff
    GenecStar%add_diff = add_diff
    set_add_diff = 0
end function set_add_diff

integer function get_n_mag(index_of_the_particle, n_mag)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: n_mag
    n_mag = GenecStar%n_mag
    get_n_mag = 0
end function get_n_mag

integer function set_n_mag(index_of_the_particle, n_mag)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: n_mag
    GenecStar%n_mag = n_mag
    set_n_mag = 0
end function set_n_mag

integer function get_alpha_F(index_of_the_particle, alpha_F)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: alpha_F
    alpha_F = GenecStar%alpha_F
    get_alpha_F = 0
end function get_alpha_F

integer function set_alpha_F(index_of_the_particle, alpha_F)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: alpha_F
    GenecStar%alpha_F = alpha_F
    set_alpha_F = 0
end function set_alpha_F

integer function get_nsmooth(index_of_the_particle, nsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nsmooth
    nsmooth = GenecStar%nsmooth
    get_nsmooth = 0
end function get_nsmooth

integer function set_nsmooth(index_of_the_particle, nsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nsmooth
    GenecStar%nsmooth = nsmooth
    set_nsmooth = 0
end function set_nsmooth

integer function get_qminsmooth(index_of_the_particle, qminsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: qminsmooth
    qminsmooth = GenecStar%qminsmooth
    get_qminsmooth = 0
end function get_qminsmooth

integer function set_qminsmooth(index_of_the_particle, qminsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: qminsmooth
    GenecStar%qminsmooth = qminsmooth
    set_qminsmooth = 0
end function set_qminsmooth

integer function get_imloss(index_of_the_particle, imloss)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: imloss
    imloss = GenecStar%imloss
    get_imloss = 0
end function get_imloss

integer function set_imloss(index_of_the_particle, imloss)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: imloss
    GenecStar%imloss = imloss
    set_imloss = 0
end function set_imloss

integer function get_fmlos(index_of_the_particle, fmlos)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fmlos
    fmlos = GenecStar%fmlos
    get_fmlos = 0
end function get_fmlos

integer function set_fmlos(index_of_the_particle, fmlos)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fmlos
    GenecStar%fmlos = fmlos
    set_fmlos = 0
end function set_fmlos

integer function get_ifitm(index_of_the_particle, ifitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ifitm
    ifitm = GenecStar%ifitm
    get_ifitm = 0
end function get_ifitm

integer function set_ifitm(index_of_the_particle, ifitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ifitm
    GenecStar%ifitm = ifitm
    set_ifitm = 0
end function set_ifitm

integer function get_fitm(index_of_the_particle, fitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fitm
    fitm = GenecStar%fitm
    get_fitm = 0
end function get_fitm

integer function set_fitm(index_of_the_particle, fitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fitm
    GenecStar%fitm = fitm
    set_fitm = 0
end function set_fitm

integer function get_fitmi(index_of_the_particle, fitmi)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fitmi
    fitmi = GenecStar%fitmi
    get_fitmi = 0
end function get_fitmi

integer function set_fitmi(index_of_the_particle, fitmi)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fitmi
    GenecStar%fitmi = fitmi
    set_fitmi = 0
end function set_fitmi

integer function get_fitmi_default(index_of_the_particle, fitmi_default)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fitmi_default
    fitmi_default = GenecStar%fitmi_default
    get_fitmi_default = 0
end function get_fitmi_default

integer function set_fitmi_default(index_of_the_particle, fitmi_default)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fitmi_default
    GenecStar%fitmi_default = fitmi_default
    set_fitmi_default = 0
end function set_fitmi_default

integer function get_deltal(index_of_the_particle, deltal)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: deltal
    deltal = GenecStar%deltal
    get_deltal = 0
end function get_deltal

integer function set_deltal(index_of_the_particle, deltal)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: deltal
    GenecStar%deltal = deltal
    set_deltal = 0
end function set_deltal

integer function get_deltat(index_of_the_particle, deltat)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: deltat
    deltat = GenecStar%deltat
    get_deltat = 0
end function get_deltat

integer function set_deltat(index_of_the_particle, deltat)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: deltat
    GenecStar%deltat = deltat
    set_deltat = 0
end function set_deltat

integer function get_nndr(index_of_the_particle, nndr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nndr
    nndr = GenecStar%nndr
    get_nndr = 0
end function get_nndr

integer function set_nndr(index_of_the_particle, nndr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nndr
    GenecStar%nndr = nndr
    set_nndr = 0
end function set_nndr

integer function get_RSG_Mdot(index_of_the_particle, RSG_Mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: RSG_Mdot
    RSG_Mdot = GenecStar%RSG_Mdot
    get_RSG_Mdot = 0
end function get_RSG_Mdot

integer function set_RSG_Mdot(index_of_the_particle, RSG_Mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: RSG_Mdot
    GenecStar%RSG_Mdot = RSG_Mdot
    set_RSG_Mdot = 0
end function set_RSG_Mdot

integer function get_SupraEddMdot(index_of_the_particle, SupraEddMdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: SupraEddMdot
    SupraEddMdot = GenecStar%SupraEddMdot
    get_SupraEddMdot = 0
end function get_SupraEddMdot

integer function set_SupraEddMdot(index_of_the_particle, SupraEddMdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: SupraEddMdot
    GenecStar%SupraEddMdot = SupraEddMdot
    set_SupraEddMdot = 0
end function set_SupraEddMdot

integer function get_Be_mdotfrac(index_of_the_particle, Be_mdotfrac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: Be_mdotfrac
    Be_mdotfrac = GenecStar%Be_mdotfrac
    get_Be_mdotfrac = 0
end function get_Be_mdotfrac

integer function set_Be_mdotfrac(index_of_the_particle, Be_mdotfrac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: Be_mdotfrac
    GenecStar%Be_mdotfrac = Be_mdotfrac
    set_Be_mdotfrac = 0
end function set_Be_mdotfrac

integer function get_start_mdot(index_of_the_particle, start_mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: start_mdot
    start_mdot = GenecStar%start_mdot
    get_start_mdot = 0
end function get_start_mdot

integer function set_start_mdot(index_of_the_particle, start_mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: start_mdot
    GenecStar%start_mdot = start_mdot
    set_start_mdot = 0
end function set_start_mdot

integer function get_iledou(index_of_the_particle, iledou)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iledou
    iledou = GenecStar%iledou
    get_iledou = 0
end function get_iledou

integer function set_iledou(index_of_the_particle, iledou)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iledou
    GenecStar%iledou = iledou
    set_iledou = 0
end function set_iledou

integer function get_idifcon(index_of_the_particle, idifcon)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idifcon
    idifcon = GenecStar%idifcon
    get_idifcon = 0
end function get_idifcon

integer function set_idifcon(index_of_the_particle, idifcon)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idifcon
    GenecStar%idifcon = idifcon
    set_idifcon = 0
end function set_idifcon

integer function get_iover(index_of_the_particle, iover)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iover
    iover = GenecStar%iover
    get_iover = 0
end function get_iover

integer function set_iover(index_of_the_particle, iover)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iover
    GenecStar%iover = iover
    set_iover = 0
end function set_iover

integer function get_elph(index_of_the_particle, elph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: elph
    elph = GenecStar%elph
    get_elph = 0
end function get_elph

integer function set_elph(index_of_the_particle, elph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: elph
    GenecStar%elph = elph
    set_elph = 0
end function set_elph

integer function get_my(index_of_the_particle, my)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: my
    my = GenecStar%my
    get_my = 0
end function get_my

integer function set_my(index_of_the_particle, my)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: my
    GenecStar%my = my
    set_my = 0
end function set_my

integer function get_dovhp(index_of_the_particle, dovhp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dovhp
    dovhp = GenecStar%dovhp
    get_dovhp = 0
end function get_dovhp

integer function set_dovhp(index_of_the_particle, dovhp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dovhp
    GenecStar%dovhp = dovhp
    set_dovhp = 0
end function set_dovhp

integer function get_iunder(index_of_the_particle, iunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iunder
    iunder = GenecStar%iunder
    get_iunder = 0
end function get_iunder

integer function set_iunder(index_of_the_particle, iunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iunder
    GenecStar%iunder = iunder
    set_iunder = 0
end function set_iunder

integer function get_dunder(index_of_the_particle, dunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dunder
    dunder = GenecStar%dunder
    get_dunder = 0
end function get_dunder

integer function set_dunder(index_of_the_particle, dunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dunder
    GenecStar%dunder = dunder
    set_dunder = 0
end function set_dunder

integer function get_gkorm(index_of_the_particle, gkorm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: gkorm
    gkorm = GenecStar%gkorm
    get_gkorm = 0
end function get_gkorm

integer function set_gkorm(index_of_the_particle, gkorm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: gkorm
    GenecStar%gkorm = gkorm
    set_gkorm = 0
end function set_gkorm

integer function get_alph(index_of_the_particle, alph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: alph
    alph = GenecStar%alph
    get_alph = 0
end function get_alph

integer function set_alph(index_of_the_particle, alph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: alph
    GenecStar%alph = alph
    set_alph = 0
end function set_alph

integer function get_agdr(index_of_the_particle, agdr)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: agdr
    agdr = GenecStar%agdr
    get_agdr = 0
end function get_agdr

integer function set_agdr(index_of_the_particle, agdr)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: agdr
    GenecStar%agdr = agdr
    set_agdr = 0
end function set_agdr

integer function get_faktor(index_of_the_particle, faktor)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: faktor
    faktor = GenecStar%faktor
    get_faktor = 0
end function get_faktor

integer function set_faktor(index_of_the_particle, faktor)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: faktor
    GenecStar%faktor = faktor
    set_faktor = 0
end function set_faktor

integer function get_dgrp(index_of_the_particle, dgrp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgrp
    dgrp = GenecStar%dgrp
    get_dgrp = 0
end function get_dgrp

integer function set_dgrp(index_of_the_particle, dgrp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgrp
    GenecStar%dgrp = dgrp
    set_dgrp = 0
end function set_dgrp

integer function get_dgrl(index_of_the_particle, dgrl)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgrl
    dgrl = GenecStar%dgrl
    get_dgrl = 0
end function get_dgrl

integer function set_dgrl(index_of_the_particle, dgrl)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgrl
    GenecStar%dgrl = dgrl
    set_dgrl = 0
end function set_dgrl

integer function get_dgry(index_of_the_particle, dgry)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgry
    dgry = GenecStar%dgry
    get_dgry = 0
end function get_dgry

integer function set_dgry(index_of_the_particle, dgry)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgry
    GenecStar%dgry = dgry
    set_dgry = 0
end function set_dgry

integer function get_dgrc(index_of_the_particle, dgrc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgrc
    dgrc = GenecStar%dgrc
    get_dgrc = 0
end function get_dgrc

integer function set_dgrc(index_of_the_particle, dgrc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgrc
    GenecStar%dgrc = dgrc
    set_dgrc = 0
end function set_dgrc

integer function get_dgro(index_of_the_particle, dgro)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgro
    dgro = GenecStar%dgro
    get_dgro = 0
end function get_dgro

integer function set_dgro(index_of_the_particle, dgro)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgro
    GenecStar%dgro = dgro
    set_dgro = 0
end function set_dgro

integer function get_dgr20(index_of_the_particle, dgr20)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgr20
    dgr20 = GenecStar%dgr20
    get_dgr20 = 0
end function get_dgr20

integer function set_dgr20(index_of_the_particle, dgr20)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgr20
    GenecStar%dgr20 = dgr20
    set_dgr20 = 0
end function set_dgr20

integer function get_nbchx(index_of_the_particle, nbchx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nbchx
    nbchx = GenecStar%nbchx
    get_nbchx = 0
end function get_nbchx

integer function set_nbchx(index_of_the_particle, nbchx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nbchx
    GenecStar%nbchx = nbchx
    set_nbchx = 0
end function set_nbchx

integer function get_nrband(index_of_the_particle, nrband)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nrband
    nrband = GenecStar%nrband
    get_nrband = 0
end function get_nrband

integer function set_nrband(index_of_the_particle, nrband)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nrband
    GenecStar%nrband = nrband
    set_nrband = 0
end function set_nrband

integer function get_xcn(index_of_the_particle, xcn)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xcn
    xcn = GenecStar%xcn
    get_xcn = 0
end function get_xcn

integer function set_xcn(index_of_the_particle, xcn)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xcn
    GenecStar%xcn = xcn
    set_xcn = 0
end function set_xcn

integer function get_islow(index_of_the_particle, islow)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: islow
    islow = GenecStar%islow
    get_islow = 0
end function get_islow

integer function set_islow(index_of_the_particle, islow)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: islow
    GenecStar%islow = islow
    set_islow = 0
end function set_islow

integer function get_icncst(index_of_the_particle, icncst)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: icncst
    icncst = GenecStar%icncst
    get_icncst = 0
end function get_icncst

integer function set_icncst(index_of_the_particle, icncst)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: icncst
    GenecStar%icncst = icncst
    set_icncst = 0
end function set_icncst

integer function get_tauH_fit(index_of_the_particle, tauH_fit)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: tauH_fit
    tauH_fit = GenecStar%tauH_fit
    get_tauH_fit = 0
end function get_tauH_fit

integer function set_tauH_fit(index_of_the_particle, tauH_fit)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: tauH_fit
    GenecStar%tauH_fit = tauH_fit
    set_tauH_fit = 0
end function set_tauH_fit

integer function get_display_plot(index_of_the_particle, display_plot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: display_plot
    display_plot = GenecStar%display_plot
    get_display_plot = 0
end function get_display_plot

integer function set_display_plot(index_of_the_particle, display_plot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: display_plot
    GenecStar%display_plot = display_plot
    set_display_plot = 0
end function set_display_plot

integer function get_iauto(index_of_the_particle, iauto)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iauto
    iauto = GenecStar%iauto
    get_iauto = 0
end function get_iauto

integer function set_iauto(index_of_the_particle, iauto)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iauto
    GenecStar%iauto = iauto
    set_iauto = 0
end function set_iauto

integer function get_iprn(index_of_the_particle, iprn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iprn
    iprn = GenecStar%iprn
    get_iprn = 0
end function get_iprn

integer function set_iprn(index_of_the_particle, iprn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iprn
    GenecStar%iprn = iprn
    set_iprn = 0
end function set_iprn

integer function get_iout(index_of_the_particle, iout)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iout
    iout = GenecStar%iout
    get_iout = 0
end function get_iout

integer function set_iout(index_of_the_particle, iout)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iout
    GenecStar%iout = iout
    set_iout = 0
end function set_iout

integer function get_itmin(index_of_the_particle, itmin)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: itmin
    itmin = GenecStar%itmin
    get_itmin = 0
end function get_itmin

integer function set_itmin(index_of_the_particle, itmin)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: itmin
    GenecStar%itmin = itmin
    set_itmin = 0
end function set_itmin

integer function get_xyfiles(index_of_the_particle, xyfiles)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: xyfiles
    xyfiles = GenecStar%xyfiles
    get_xyfiles = 0
end function get_xyfiles

integer function set_xyfiles(index_of_the_particle, xyfiles)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: xyfiles
    GenecStar%xyfiles = xyfiles
    set_xyfiles = 0
end function set_xyfiles

integer function get_idebug(index_of_the_particle, idebug)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idebug
    idebug = GenecStar%idebug
    get_idebug = 0
end function get_idebug

integer function set_idebug(index_of_the_particle, idebug)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idebug
    GenecStar%idebug = idebug
    set_idebug = 0
end function set_idebug

integer function get_itests(index_of_the_particle, itests)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: itests
    itests = GenecStar%itests
    get_itests = 0
end function get_itests

integer function set_itests(index_of_the_particle, itests)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: itests
    GenecStar%itests = itests
    set_itests = 0
end function set_itests

integer function get_verbose(index_of_the_particle, verbose)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: verbose
    verbose = GenecStar%verbose
    get_verbose = 0
end function get_verbose

integer function set_verbose(index_of_the_particle, verbose)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: verbose
    GenecStar%verbose = verbose
    set_verbose = 0
end function set_verbose

integer function get_stop_deg(index_of_the_particle, stop_deg)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: stop_deg
    stop_deg = GenecStar%stop_deg
    get_stop_deg = 0
end function get_stop_deg

integer function set_stop_deg(index_of_the_particle, stop_deg)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: stop_deg
    GenecStar%stop_deg = stop_deg
    set_stop_deg = 0
end function set_stop_deg

integer function get_n_snap(index_of_the_particle, n_snap)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: n_snap
    n_snap = GenecStar%n_snap
    get_n_snap = 0
end function get_n_snap

integer function set_n_snap(index_of_the_particle, n_snap)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: n_snap
    GenecStar%n_snap = n_snap
    set_n_snap = 0
end function set_n_snap

integer function get_gms(index_of_the_particle, gms)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: gms
    gms = GenecStar%gms
    get_gms = 0
end function get_gms

integer function set_gms(index_of_the_particle, gms)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: gms
    GenecStar%gms = gms
    set_gms = 0
end function set_gms

integer function get_alter(index_of_the_particle, alter)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: alter
    alter = GenecStar%alter
    get_alter = 0
end function get_alter

integer function set_alter(index_of_the_particle, alter)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: alter
    GenecStar%alter = alter
    set_alter = 0
end function set_alter

integer function get_gls(index_of_the_particle, gls)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: gls
    gls = GenecStar%gls
    get_gls = 0
end function get_gls

integer function set_gls(index_of_the_particle, gls)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: gls
    GenecStar%gls = gls
    set_gls = 0
end function set_gls

integer function get_teff(index_of_the_particle, teff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: teff
    teff = GenecStar%teff
    get_teff = 0
end function get_teff

integer function set_teff(index_of_the_particle, teff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: teff
    GenecStar%teff = teff
    set_teff = 0
end function set_teff

integer function get_glsv(index_of_the_particle, glsv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: glsv
    glsv = GenecStar%glsv
    get_glsv = 0
end function get_glsv

integer function set_glsv(index_of_the_particle, glsv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: glsv
    GenecStar%glsv = glsv
    set_glsv = 0
end function set_glsv

integer function get_teffv(index_of_the_particle, teffv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: teffv
    teffv = GenecStar%teffv
    get_teffv = 0
end function get_teffv

integer function set_teffv(index_of_the_particle, teffv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: teffv
    GenecStar%teffv = teffv
    set_teffv = 0
end function set_teffv

integer function get_dzeitj(index_of_the_particle, dzeitj)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dzeitj
    dzeitj = GenecStar%dzeitj
    get_dzeitj = 0
end function get_dzeitj

integer function set_dzeitj(index_of_the_particle, dzeitj)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dzeitj
    GenecStar%dzeitj = dzeitj
    set_dzeitj = 0
end function set_dzeitj

integer function get_dzeit(index_of_the_particle, dzeit)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dzeit
    dzeit = GenecStar%dzeit
    get_dzeit = 0
end function get_dzeit

integer function set_dzeit(index_of_the_particle, dzeit)
    use const, only: year
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dzeit
    GenecStar%dzeit = dzeit
    GenecStar%dzeitj = dzeit / year
    set_dzeit = 0
end function set_dzeit

integer function get_dzeitv(index_of_the_particle, dzeitv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dzeitv
    dzeitv = GenecStar%dzeitv
    get_dzeitv = 0
end function get_dzeitv

integer function set_dzeitv(index_of_the_particle, dzeitv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dzeitv
    GenecStar%dzeitv = dzeitv
    set_dzeitv = 0
end function set_dzeitv

integer function get_xmini(index_of_the_particle, xmini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xmini
    xmini = GenecStar%xmini
    get_xmini = 0
end function get_xmini

integer function set_xmini(index_of_the_particle, xmini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xmini
    GenecStar%xmini = xmini
    set_xmini = 0
end function set_xmini

integer function get_ab(index_of_the_particle, ab)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: ab
    ab = GenecStar%ab
    get_ab = 0
end function get_ab

integer function set_ab(index_of_the_particle, ab)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: ab
    GenecStar%ab = ab
    set_ab = 0
end function set_ab

integer function get_dm_lost(index_of_the_particle, dm_lost)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dm_lost
    dm_lost = GenecStar%dm_lost
    get_dm_lost = 0
end function get_dm_lost

integer function set_dm_lost(index_of_the_particle, dm_lost)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dm_lost
    GenecStar%dm_lost = dm_lost
    set_dm_lost = 0
end function set_dm_lost

integer function get_m(index_of_the_particle, m)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: m
    m = GenecStar%m
    get_m = 0
end function get_m

integer function set_m(index_of_the_particle, m)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: m
    GenecStar%m = m
    set_m = 0
end function set_m

integer function get_summas(index_of_the_particle, summas)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: summas
    summas = GenecStar%summas
    get_summas = 0
end function get_summas

integer function set_summas(index_of_the_particle, summas)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: summas
    GenecStar%summas = summas
    set_summas = 0
end function set_summas

integer function get_dk(index_of_the_particle, dk)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dk
    dk = GenecStar%dk
    get_dk = 0
end function get_dk

integer function set_dk(index_of_the_particle, dk)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dk
    GenecStar%dk = dk
    set_dk = 0
end function set_dk

integer function get_rlp(index_of_the_particle, rlp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rlp
    rlp = GenecStar%rlp
    get_rlp = 0
end function get_rlp

integer function set_rlp(index_of_the_particle, rlp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rlp
    GenecStar%rlp = rlp
    set_rlp = 0
end function set_rlp

integer function get_rlt(index_of_the_particle, rlt)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rlt
    rlt = GenecStar%rlt
    get_rlt = 0
end function get_rlt

integer function set_rlt(index_of_the_particle, rlt)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rlt
    GenecStar%rlt = rlt
    set_rlt = 0
end function set_rlt

integer function get_rlc(index_of_the_particle, rlc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rlc
    rlc = GenecStar%rlc
    get_rlc = 0
end function get_rlc

integer function set_rlc(index_of_the_particle, rlc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rlc
    GenecStar%rlc = rlc
    set_rlc = 0
end function set_rlc

integer function get_rrp(index_of_the_particle, rrp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rrp
    rrp = GenecStar%rrp
    get_rrp = 0
end function get_rrp

integer function set_rrp(index_of_the_particle, rrp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rrp
    GenecStar%rrp = rrp
    set_rrp = 0
end function set_rrp

integer function get_rrt(index_of_the_particle, rrt)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rrt
    rrt = GenecStar%rrt
    get_rrt = 0
end function get_rrt

integer function set_rrt(index_of_the_particle, rrt)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rrt
    GenecStar%rrt = rrt
    set_rrt = 0
end function set_rrt

integer function get_rrc(index_of_the_particle, rrc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rrc
    rrc = GenecStar%rrc
    get_rrc = 0
end function get_rrc

integer function set_rrc(index_of_the_particle, rrc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rrc
    GenecStar%rrc = rrc
    set_rrc = 0
end function set_rrc

integer function get_rtp(index_of_the_particle, rtp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rtp
    rtp = GenecStar%rtp
    get_rtp = 0
end function get_rtp

integer function set_rtp(index_of_the_particle, rtp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rtp
    GenecStar%rtp = rtp
    set_rtp = 0
end function set_rtp

integer function get_rtt(index_of_the_particle, rtt)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rtt
    rtt = GenecStar%rtt
    get_rtt = 0
end function get_rtt

integer function set_rtt(index_of_the_particle, rtt)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rtt
    GenecStar%rtt = rtt
    set_rtt = 0
end function set_rtt

integer function get_rtc(index_of_the_particle, rtc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rtc
    rtc = GenecStar%rtc
    get_rtc = 0
end function get_rtc

integer function set_rtc(index_of_the_particle, rtc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rtc
    GenecStar%rtc = rtc
    set_rtc = 0
end function set_rtc

integer function get_tdiff(index_of_the_particle, tdiff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: tdiff
    tdiff = GenecStar%tdiff
    get_tdiff = 0
end function get_tdiff

integer function set_tdiff(index_of_the_particle, tdiff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: tdiff
    GenecStar%tdiff = tdiff
    set_tdiff = 0
end function set_tdiff

integer function get_suminenv(index_of_the_particle, suminenv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: suminenv
    suminenv = GenecStar%suminenv
    get_suminenv = 0
end function get_suminenv

integer function set_suminenv(index_of_the_particle, suminenv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: suminenv
    GenecStar%suminenv = suminenv
    set_suminenv = 0
end function set_suminenv

integer function get_xltotbeg(index_of_the_particle, xltotbeg)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xltotbeg
    xltotbeg = GenecStar%xltotbeg
    get_xltotbeg = 0
end function get_xltotbeg

integer function set_xltotbeg(index_of_the_particle, xltotbeg)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xltotbeg
    GenecStar%xltotbeg = xltotbeg
    set_xltotbeg = 0
end function set_xltotbeg

integer function get_dlelexprev(index_of_the_particle, dlelexprev)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dlelexprev
    dlelexprev = GenecStar%dlelexprev
    get_dlelexprev = 0
end function get_dlelexprev

integer function set_dlelexprev(index_of_the_particle, dlelexprev)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dlelexprev
    GenecStar%dlelexprev = dlelexprev
    set_dlelexprev = 0
end function set_dlelexprev

integer function get_zams_radius(index_of_the_particle, zams_radius)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: zams_radius
    zams_radius = GenecStar%zams_radius
    get_zams_radius = 0
end function get_zams_radius

integer function set_zams_radius(index_of_the_particle, zams_radius)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: zams_radius
    GenecStar%zams_radius = zams_radius
    set_zams_radius = 0
end function set_zams_radius

integer function get_mbelx(index_of_the_particle, mbelx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: mbelx
    mbelx = GenecStar%mbelx
    get_mbelx = 0
end function get_mbelx

integer function set_mbelx(index_of_the_particle, mbelx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: mbelx
    GenecStar%mbelx = mbelx
    set_mbelx = 0
end function set_mbelx

integer function get_xtefflast(index_of_the_particle, xtefflast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xtefflast
    xtefflast = GenecStar%xtefflast
    get_xtefflast = 0
end function get_xtefflast

integer function set_xtefflast(index_of_the_particle, xtefflast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xtefflast
    GenecStar%xtefflast = xtefflast
    set_xtefflast = 0
end function set_xtefflast

integer function get_xllast(index_of_the_particle, xllast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xllast
    xllast = GenecStar%xllast
    get_xllast = 0
end function get_xllast

integer function set_xllast(index_of_the_particle, xllast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xllast
    GenecStar%xllast = xllast
    set_xllast = 0
end function set_xllast

integer function get_xrholast(index_of_the_particle, xrholast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xrholast
    xrholast = GenecStar%xrholast
    get_xrholast = 0
end function get_xrholast

integer function set_xrholast(index_of_the_particle, xrholast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xrholast
    GenecStar%xrholast = xrholast
    set_xrholast = 0
end function set_xrholast

integer function get_xclast(index_of_the_particle, xclast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xclast
    xclast = GenecStar%xclast
    get_xclast = 0
end function get_xclast

integer function set_xclast(index_of_the_particle, xclast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xclast
    GenecStar%xclast = xclast
    set_xclast = 0
end function set_xclast

integer function get_xtclast(index_of_the_particle, xtclast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xtclast
    xtclast = GenecStar%xtclast
    get_xtclast = 0
end function get_xtclast

integer function set_xtclast(index_of_the_particle, xtclast)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xtclast
    GenecStar%xtclast = xtclast
    set_xtclast = 0
end function set_xtclast

integer function get_inum(index_of_the_particle, inum)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: inum
    inum = GenecStar%inum
    get_inum = 0
end function get_inum

integer function set_inum(index_of_the_particle, inum)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: inum
    GenecStar%inum = inum
    set_inum = 0
end function set_inum

integer function get_nsugi(index_of_the_particle, nsugi)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nsugi
    nsugi = GenecStar%nsugi
    get_nsugi = 0
end function get_nsugi

integer function set_nsugi(index_of_the_particle, nsugi)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nsugi
    GenecStar%nsugi = nsugi
    set_nsugi = 0
end function set_nsugi

integer function get_period(index_of_the_particle, period)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: period
    period = GenecStar%period
    get_period = 0
end function get_period

integer function set_period(index_of_the_particle, period)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: period
    GenecStar%period = period
    set_period = 0
end function set_period

integer function get_r_core(index_of_the_particle, r_core)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: r_core
    r_core = GenecStar%r_core
    get_r_core = 0
end function get_r_core

integer function set_r_core(index_of_the_particle, r_core)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: r_core
    GenecStar%r_core = r_core
    set_r_core = 0
end function set_r_core

integer function get_vna(index_of_the_particle, vna)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: vna
    vna = GenecStar%vna
    get_vna = 0
end function get_vna

integer function set_vna(index_of_the_particle, vna)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: vna
    GenecStar%vna = vna
    set_vna = 0
end function set_vna

integer function get_vnr(index_of_the_particle, vnr)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: vnr
    vnr = GenecStar%vnr
    get_vnr = 0
end function get_vnr

integer function set_vnr(index_of_the_particle, vnr)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: vnr
    GenecStar%vnr = vnr
    set_vnr = 0
end function set_vnr

integer function get_xlostneu(index_of_the_particle, xlostneu)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xlostneu
    xlostneu = GenecStar%xlostneu
    get_xlostneu = 0
end function get_xlostneu

integer function set_xlostneu(index_of_the_particle, xlostneu)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xlostneu
    GenecStar%xlostneu = xlostneu
    set_xlostneu = 0
end function set_xlostneu

integer function get_q_at_zone(index_of_the_particle, zone, q)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: q
    integer:: i
    get_q_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        q = GenecStar%q(i)
        get_q_at_zone = 0
    end if
end function get_q_at_zone

integer function set_q_at_zone(index_of_the_particle, zone, q)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: q
    integer:: i
    set_q_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%q(i) = q
        set_q_at_zone = 0
    end if
end function set_q_at_zone

integer function get_p_at_zone(index_of_the_particle, zone, p)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: p
    integer:: i
    get_p_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        p = GenecStar%p(i)
        get_p_at_zone = 0
    end if
end function get_p_at_zone

integer function set_p_at_zone(index_of_the_particle, zone, p)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: p
    integer:: i
    set_p_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%p(i) = p
        set_p_at_zone = 0
    end if
end function set_p_at_zone

integer function get_t_at_zone(index_of_the_particle, zone, t)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: t
    integer:: i
    get_t_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        t = GenecStar%t(i)
        get_t_at_zone = 0
    end if
end function get_t_at_zone

integer function set_t_at_zone(index_of_the_particle, zone, t)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: t
    integer:: i
    set_t_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%t(i) = t
        set_t_at_zone = 0
    end if
end function set_t_at_zone

integer function get_r_at_zone(index_of_the_particle, zone, r)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: r
    integer:: i
    get_r_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        r = GenecStar%r(i)
        get_r_at_zone = 0
    end if
end function get_r_at_zone

integer function set_r_at_zone(index_of_the_particle, zone, r)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: r
    integer:: i
    set_r_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%r(i) = r
        set_r_at_zone = 0
    end if
end function set_r_at_zone

integer function get_s_at_zone(index_of_the_particle, zone, s)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: s
    integer:: i
    get_s_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        s = GenecStar%s(i)
        get_s_at_zone = 0
    end if
end function get_s_at_zone

integer function set_s_at_zone(index_of_the_particle, zone, s)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: s
    integer:: i
    set_s_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%s(i) = s
        set_s_at_zone = 0
    end if
end function set_s_at_zone

integer function get_x_at_zone(index_of_the_particle, zone, x)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: x
    integer:: i
    get_x_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        x = GenecStar%x(i)
        get_x_at_zone = 0
    end if
end function get_x_at_zone

integer function set_x_at_zone(index_of_the_particle, zone, x)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: x
    integer:: i
    set_x_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%x(i) = x
        set_x_at_zone = 0
    end if
end function set_x_at_zone

integer function get_y3_at_zone(index_of_the_particle, zone, y3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: y3
    integer:: i
    get_y3_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        y3 = GenecStar%y3(i)
        get_y3_at_zone = 0
    end if
end function get_y3_at_zone

integer function set_y3_at_zone(index_of_the_particle, zone, y3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: y3
    integer:: i
    set_y3_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%y3(i) = y3
        set_y3_at_zone = 0
    end if
end function set_y3_at_zone

integer function get_y_at_zone(index_of_the_particle, zone, y)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: y
    integer:: i
    get_y_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        y = GenecStar%y(i)
        get_y_at_zone = 0
    end if
end function get_y_at_zone

integer function set_y_at_zone(index_of_the_particle, zone, y)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: y
    integer:: i
    set_y_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%y(i) = y
        set_y_at_zone = 0
    end if
end function set_y_at_zone

integer function get_xc12_at_zone(index_of_the_particle, zone, xc12)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xc12
    integer:: i
    get_xc12_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xc12 = GenecStar%xc12(i)
        get_xc12_at_zone = 0
    end if
end function get_xc12_at_zone

integer function set_xc12_at_zone(index_of_the_particle, zone, xc12)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xc12
    integer:: i
    set_xc12_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xc12(i) = xc12
        set_xc12_at_zone = 0
    end if
end function set_xc12_at_zone

integer function get_xc13_at_zone(index_of_the_particle, zone, xc13)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xc13
    integer:: i
    get_xc13_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xc13 = GenecStar%xc13(i)
        get_xc13_at_zone = 0
    end if
end function get_xc13_at_zone

integer function set_xc13_at_zone(index_of_the_particle, zone, xc13)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xc13
    integer:: i
    set_xc13_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xc13(i) = xc13
        set_xc13_at_zone = 0
    end if
end function set_xc13_at_zone

integer function get_xn14_at_zone(index_of_the_particle, zone, xn14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xn14
    integer:: i
    get_xn14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xn14 = GenecStar%xn14(i)
        get_xn14_at_zone = 0
    end if
end function get_xn14_at_zone

integer function set_xn14_at_zone(index_of_the_particle, zone, xn14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xn14
    integer:: i
    set_xn14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xn14(i) = xn14
        set_xn14_at_zone = 0
    end if
end function set_xn14_at_zone

integer function get_xn15_at_zone(index_of_the_particle, zone, xn15)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xn15
    integer:: i
    get_xn15_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xn15 = GenecStar%xn15(i)
        get_xn15_at_zone = 0
    end if
end function get_xn15_at_zone

integer function set_xn15_at_zone(index_of_the_particle, zone, xn15)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xn15
    integer:: i
    set_xn15_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xn15(i) = xn15
        set_xn15_at_zone = 0
    end if
end function set_xn15_at_zone

integer function get_xo16_at_zone(index_of_the_particle, zone, xo16)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xo16
    integer:: i
    get_xo16_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xo16 = GenecStar%xo16(i)
        get_xo16_at_zone = 0
    end if
end function get_xo16_at_zone

integer function set_xo16_at_zone(index_of_the_particle, zone, xo16)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xo16
    integer:: i
    set_xo16_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xo16(i) = xo16
        set_xo16_at_zone = 0
    end if
end function set_xo16_at_zone

integer function get_xo17_at_zone(index_of_the_particle, zone, xo17)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xo17
    integer:: i
    get_xo17_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xo17 = GenecStar%xo17(i)
        get_xo17_at_zone = 0
    end if
end function get_xo17_at_zone

integer function set_xo17_at_zone(index_of_the_particle, zone, xo17)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xo17
    integer:: i
    set_xo17_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xo17(i) = xo17
        set_xo17_at_zone = 0
    end if
end function set_xo17_at_zone

integer function get_xo18_at_zone(index_of_the_particle, zone, xo18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xo18
    integer:: i
    get_xo18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xo18 = GenecStar%xo18(i)
        get_xo18_at_zone = 0
    end if
end function get_xo18_at_zone

integer function set_xo18_at_zone(index_of_the_particle, zone, xo18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xo18
    integer:: i
    set_xo18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xo18(i) = xo18
        set_xo18_at_zone = 0
    end if
end function set_xo18_at_zone

integer function get_xne20_at_zone(index_of_the_particle, zone, xne20)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xne20
    integer:: i
    get_xne20_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xne20 = GenecStar%xne20(i)
        get_xne20_at_zone = 0
    end if
end function get_xne20_at_zone

integer function set_xne20_at_zone(index_of_the_particle, zone, xne20)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xne20
    integer:: i
    set_xne20_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xne20(i) = xne20
        set_xne20_at_zone = 0
    end if
end function set_xne20_at_zone

integer function get_xne22_at_zone(index_of_the_particle, zone, xne22)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xne22
    integer:: i
    get_xne22_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xne22 = GenecStar%xne22(i)
        get_xne22_at_zone = 0
    end if
end function get_xne22_at_zone

integer function set_xne22_at_zone(index_of_the_particle, zone, xne22)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xne22
    integer:: i
    set_xne22_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xne22(i) = xne22
        set_xne22_at_zone = 0
    end if
end function set_xne22_at_zone

integer function get_xmg24_at_zone(index_of_the_particle, zone, xmg24)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xmg24
    integer:: i
    get_xmg24_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xmg24 = GenecStar%xmg24(i)
        get_xmg24_at_zone = 0
    end if
end function get_xmg24_at_zone

integer function set_xmg24_at_zone(index_of_the_particle, zone, xmg24)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xmg24
    integer:: i
    set_xmg24_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xmg24(i) = xmg24
        set_xmg24_at_zone = 0
    end if
end function set_xmg24_at_zone

integer function get_xmg25_at_zone(index_of_the_particle, zone, xmg25)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xmg25
    integer:: i
    get_xmg25_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xmg25 = GenecStar%xmg25(i)
        get_xmg25_at_zone = 0
    end if
end function get_xmg25_at_zone

integer function set_xmg25_at_zone(index_of_the_particle, zone, xmg25)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xmg25
    integer:: i
    set_xmg25_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xmg25(i) = xmg25
        set_xmg25_at_zone = 0
    end if
end function set_xmg25_at_zone

integer function get_xmg26_at_zone(index_of_the_particle, zone, xmg26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xmg26
    integer:: i
    get_xmg26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xmg26 = GenecStar%xmg26(i)
        get_xmg26_at_zone = 0
    end if
end function get_xmg26_at_zone

integer function set_xmg26_at_zone(index_of_the_particle, zone, xmg26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xmg26
    integer:: i
    set_xmg26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xmg26(i) = xmg26
        set_xmg26_at_zone = 0
    end if
end function set_xmg26_at_zone

integer function get_xf19_at_zone(index_of_the_particle, zone, xf19)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xf19
    integer:: i
    get_xf19_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xf19 = GenecStar%xf19(i)
        get_xf19_at_zone = 0
    end if
end function get_xf19_at_zone

integer function set_xf19_at_zone(index_of_the_particle, zone, xf19)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xf19
    integer:: i
    set_xf19_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xf19(i) = xf19
        set_xf19_at_zone = 0
    end if
end function set_xf19_at_zone

integer function get_xne21_at_zone(index_of_the_particle, zone, xne21)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xne21
    integer:: i
    get_xne21_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xne21 = GenecStar%xne21(i)
        get_xne21_at_zone = 0
    end if
end function get_xne21_at_zone

integer function set_xne21_at_zone(index_of_the_particle, zone, xne21)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xne21
    integer:: i
    set_xne21_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xne21(i) = xne21
        set_xne21_at_zone = 0
    end if
end function set_xne21_at_zone

integer function get_xna23_at_zone(index_of_the_particle, zone, xna23)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xna23
    integer:: i
    get_xna23_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xna23 = GenecStar%xna23(i)
        get_xna23_at_zone = 0
    end if
end function get_xna23_at_zone

integer function set_xna23_at_zone(index_of_the_particle, zone, xna23)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xna23
    integer:: i
    set_xna23_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xna23(i) = xna23
        set_xna23_at_zone = 0
    end if
end function set_xna23_at_zone

integer function get_xal27_at_zone(index_of_the_particle, zone, xal27)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xal27
    integer:: i
    get_xal27_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xal27 = GenecStar%xal27(i)
        get_xal27_at_zone = 0
    end if
end function get_xal27_at_zone

integer function set_xal27_at_zone(index_of_the_particle, zone, xal27)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xal27
    integer:: i
    set_xal27_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xal27(i) = xal27
        set_xal27_at_zone = 0
    end if
end function set_xal27_at_zone

integer function get_xsi28_at_zone(index_of_the_particle, zone, xsi28)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xsi28
    integer:: i
    get_xsi28_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xsi28 = GenecStar%xsi28(i)
        get_xsi28_at_zone = 0
    end if
end function get_xsi28_at_zone

integer function set_xsi28_at_zone(index_of_the_particle, zone, xsi28)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xsi28
    integer:: i
    set_xsi28_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xsi28(i) = xsi28
        set_xsi28_at_zone = 0
    end if
end function set_xsi28_at_zone

integer function get_xc14_at_zone(index_of_the_particle, zone, xc14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xc14
    integer:: i
    get_xc14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xc14 = GenecStar%xc14(i)
        get_xc14_at_zone = 0
    end if
end function get_xc14_at_zone

integer function set_xc14_at_zone(index_of_the_particle, zone, xc14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xc14
    integer:: i
    set_xc14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xc14(i) = xc14
        set_xc14_at_zone = 0
    end if
end function set_xc14_at_zone

integer function get_xf18_at_zone(index_of_the_particle, zone, xf18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xf18
    integer:: i
    get_xf18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xf18 = GenecStar%xf18(i)
        get_xf18_at_zone = 0
    end if
end function get_xf18_at_zone

integer function set_xf18_at_zone(index_of_the_particle, zone, xf18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xf18
    integer:: i
    set_xf18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xf18(i) = xf18
        set_xf18_at_zone = 0
    end if
end function set_xf18_at_zone

integer function get_xal26_at_zone(index_of_the_particle, zone, xal26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xal26
    integer:: i
    get_xal26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xal26 = GenecStar%xal26(i)
        get_xal26_at_zone = 0
    end if
end function get_xal26_at_zone

integer function set_xal26_at_zone(index_of_the_particle, zone, xal26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xal26
    integer:: i
    set_xal26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xal26(i) = xal26
        set_xal26_at_zone = 0
    end if
end function set_xal26_at_zone

integer function get_xneut_at_zone(index_of_the_particle, zone, xneut)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xneut
    integer:: i
    get_xneut_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xneut = GenecStar%xneut(i)
        get_xneut_at_zone = 0
    end if
end function get_xneut_at_zone

integer function set_xneut_at_zone(index_of_the_particle, zone, xneut)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xneut
    integer:: i
    set_xneut_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xneut(i) = xneut
        set_xneut_at_zone = 0
    end if
end function set_xneut_at_zone

integer function get_xprot_at_zone(index_of_the_particle, zone, xprot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xprot
    integer:: i
    get_xprot_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xprot = GenecStar%xprot(i)
        get_xprot_at_zone = 0
    end if
end function get_xprot_at_zone

integer function set_xprot_at_zone(index_of_the_particle, zone, xprot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xprot
    integer:: i
    set_xprot_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xprot(i) = xprot
        set_xprot_at_zone = 0
    end if
end function set_xprot_at_zone

integer function get_omegi_at_zone(index_of_the_particle, zone, omegi)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: omegi
    integer:: i
    get_omegi_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        omegi = GenecStar%omegi(i)
        get_omegi_at_zone = 0
    end if
end function get_omegi_at_zone

integer function set_omegi_at_zone(index_of_the_particle, zone, omegi)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: omegi
    integer:: i
    set_omegi_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%omegi(i) = omegi
        set_omegi_at_zone = 0
    end if
end function set_omegi_at_zone

integer function get_xbid_at_zone(index_of_the_particle, zone, xbid)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xbid
    integer:: i
    get_xbid_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xbid = GenecStar%xbid(i)
        get_xbid_at_zone = 0
    end if
end function get_xbid_at_zone

integer function set_xbid_at_zone(index_of_the_particle, zone, xbid)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xbid
    integer:: i
    set_xbid_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xbid(i) = xbid
        set_xbid_at_zone = 0
    end if
end function set_xbid_at_zone

integer function get_xbid1_at_zone(index_of_the_particle, zone, xbid1)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: xbid1
    integer:: i
    get_xbid1_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        xbid1 = GenecStar%xbid1(i)
        get_xbid1_at_zone = 0
    end if
end function get_xbid1_at_zone

integer function set_xbid1_at_zone(index_of_the_particle, zone, xbid1)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: xbid1
    integer:: i
    set_xbid1_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%xbid1(i) = xbid1
        set_xbid1_at_zone = 0
    end if
end function set_xbid1_at_zone

integer function get_vp_at_zone(index_of_the_particle, zone, vp)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vp
    integer:: i
    get_vp_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vp = GenecStar%vp(i)
        get_vp_at_zone = 0
    end if
end function get_vp_at_zone

integer function set_vp_at_zone(index_of_the_particle, zone, vp)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vp
    integer:: i
    set_vp_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vp(i) = vp
        set_vp_at_zone = 0
    end if
end function set_vp_at_zone

integer function get_vt_at_zone(index_of_the_particle, zone, vt)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vt
    integer:: i
    get_vt_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vt = GenecStar%vt(i)
        get_vt_at_zone = 0
    end if
end function get_vt_at_zone

integer function set_vt_at_zone(index_of_the_particle, zone, vt)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vt
    integer:: i
    set_vt_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vt(i) = vt
        set_vt_at_zone = 0
    end if
end function set_vt_at_zone

integer function get_vr_at_zone(index_of_the_particle, zone, vr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vr
    integer:: i
    get_vr_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vr = GenecStar%vr(i)
        get_vr_at_zone = 0
    end if
end function get_vr_at_zone

integer function set_vr_at_zone(index_of_the_particle, zone, vr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vr
    integer:: i
    set_vr_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vr(i) = vr
        set_vr_at_zone = 0
    end if
end function set_vr_at_zone

integer function get_vs_at_zone(index_of_the_particle, zone, vs)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vs
    integer:: i
    get_vs_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vs = GenecStar%vs(i)
        get_vs_at_zone = 0
    end if
end function get_vs_at_zone

integer function set_vs_at_zone(index_of_the_particle, zone, vs)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vs
    integer:: i
    set_vs_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vs(i) = vs
        set_vs_at_zone = 0
    end if
end function set_vs_at_zone

integer function get_vx_at_zone(index_of_the_particle, zone, vx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vx
    integer:: i
    get_vx_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vx = GenecStar%vx(i)
        get_vx_at_zone = 0
    end if
end function get_vx_at_zone

integer function set_vx_at_zone(index_of_the_particle, zone, vx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vx
    integer:: i
    set_vx_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vx(i) = vx
        set_vx_at_zone = 0
    end if
end function set_vx_at_zone

integer function get_vy_at_zone(index_of_the_particle, zone, vy)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vy
    integer:: i
    get_vy_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vy = GenecStar%vy(i)
        get_vy_at_zone = 0
    end if
end function get_vy_at_zone

integer function set_vy_at_zone(index_of_the_particle, zone, vy)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vy
    integer:: i
    set_vy_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vy(i) = vy
        set_vy_at_zone = 0
    end if
end function set_vy_at_zone

integer function get_vy3_at_zone(index_of_the_particle, zone, vy3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vy3
    integer:: i
    get_vy3_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vy3 = GenecStar%vy3(i)
        get_vy3_at_zone = 0
    end if
end function get_vy3_at_zone

integer function set_vy3_at_zone(index_of_the_particle, zone, vy3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vy3
    integer:: i
    set_vy3_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vy3(i) = vy3
        set_vy3_at_zone = 0
    end if
end function set_vy3_at_zone

integer function get_vxc12_at_zone(index_of_the_particle, zone, vxc12)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxc12
    integer:: i
    get_vxc12_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxc12 = GenecStar%vxc12(i)
        get_vxc12_at_zone = 0
    end if
end function get_vxc12_at_zone

integer function set_vxc12_at_zone(index_of_the_particle, zone, vxc12)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxc12
    integer:: i
    set_vxc12_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxc12(i) = vxc12
        set_vxc12_at_zone = 0
    end if
end function set_vxc12_at_zone

integer function get_vxc13_at_zone(index_of_the_particle, zone, vxc13)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxc13
    integer:: i
    get_vxc13_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxc13 = GenecStar%vxc13(i)
        get_vxc13_at_zone = 0
    end if
end function get_vxc13_at_zone

integer function set_vxc13_at_zone(index_of_the_particle, zone, vxc13)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxc13
    integer:: i
    set_vxc13_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxc13(i) = vxc13
        set_vxc13_at_zone = 0
    end if
end function set_vxc13_at_zone

integer function get_vxn14_at_zone(index_of_the_particle, zone, vxn14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxn14
    integer:: i
    get_vxn14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxn14 = GenecStar%vxn14(i)
        get_vxn14_at_zone = 0
    end if
end function get_vxn14_at_zone

integer function set_vxn14_at_zone(index_of_the_particle, zone, vxn14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxn14
    integer:: i
    set_vxn14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxn14(i) = vxn14
        set_vxn14_at_zone = 0
    end if
end function set_vxn14_at_zone

integer function get_vxn15_at_zone(index_of_the_particle, zone, vxn15)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxn15
    integer:: i
    get_vxn15_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxn15 = GenecStar%vxn15(i)
        get_vxn15_at_zone = 0
    end if
end function get_vxn15_at_zone

integer function set_vxn15_at_zone(index_of_the_particle, zone, vxn15)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxn15
    integer:: i
    set_vxn15_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxn15(i) = vxn15
        set_vxn15_at_zone = 0
    end if
end function set_vxn15_at_zone

integer function get_vxo16_at_zone(index_of_the_particle, zone, vxo16)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxo16
    integer:: i
    get_vxo16_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxo16 = GenecStar%vxo16(i)
        get_vxo16_at_zone = 0
    end if
end function get_vxo16_at_zone

integer function set_vxo16_at_zone(index_of_the_particle, zone, vxo16)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxo16
    integer:: i
    set_vxo16_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxo16(i) = vxo16
        set_vxo16_at_zone = 0
    end if
end function set_vxo16_at_zone

integer function get_vxo17_at_zone(index_of_the_particle, zone, vxo17)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxo17
    integer:: i
    get_vxo17_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxo17 = GenecStar%vxo17(i)
        get_vxo17_at_zone = 0
    end if
end function get_vxo17_at_zone

integer function set_vxo17_at_zone(index_of_the_particle, zone, vxo17)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxo17
    integer:: i
    set_vxo17_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxo17(i) = vxo17
        set_vxo17_at_zone = 0
    end if
end function set_vxo17_at_zone

integer function get_vxo18_at_zone(index_of_the_particle, zone, vxo18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxo18
    integer:: i
    get_vxo18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxo18 = GenecStar%vxo18(i)
        get_vxo18_at_zone = 0
    end if
end function get_vxo18_at_zone

integer function set_vxo18_at_zone(index_of_the_particle, zone, vxo18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxo18
    integer:: i
    set_vxo18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxo18(i) = vxo18
        set_vxo18_at_zone = 0
    end if
end function set_vxo18_at_zone

integer function get_vxne20_at_zone(index_of_the_particle, zone, vxne20)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxne20
    integer:: i
    get_vxne20_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxne20 = GenecStar%vxne20(i)
        get_vxne20_at_zone = 0
    end if
end function get_vxne20_at_zone

integer function set_vxne20_at_zone(index_of_the_particle, zone, vxne20)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxne20
    integer:: i
    set_vxne20_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxne20(i) = vxne20
        set_vxne20_at_zone = 0
    end if
end function set_vxne20_at_zone

integer function get_vxne22_at_zone(index_of_the_particle, zone, vxne22)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxne22
    integer:: i
    get_vxne22_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxne22 = GenecStar%vxne22(i)
        get_vxne22_at_zone = 0
    end if
end function get_vxne22_at_zone

integer function set_vxne22_at_zone(index_of_the_particle, zone, vxne22)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxne22
    integer:: i
    set_vxne22_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxne22(i) = vxne22
        set_vxne22_at_zone = 0
    end if
end function set_vxne22_at_zone

integer function get_vxmg24_at_zone(index_of_the_particle, zone, vxmg24)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxmg24
    integer:: i
    get_vxmg24_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxmg24 = GenecStar%vxmg24(i)
        get_vxmg24_at_zone = 0
    end if
end function get_vxmg24_at_zone

integer function set_vxmg24_at_zone(index_of_the_particle, zone, vxmg24)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxmg24
    integer:: i
    set_vxmg24_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxmg24(i) = vxmg24
        set_vxmg24_at_zone = 0
    end if
end function set_vxmg24_at_zone

integer function get_vxmg25_at_zone(index_of_the_particle, zone, vxmg25)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxmg25
    integer:: i
    get_vxmg25_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxmg25 = GenecStar%vxmg25(i)
        get_vxmg25_at_zone = 0
    end if
end function get_vxmg25_at_zone

integer function set_vxmg25_at_zone(index_of_the_particle, zone, vxmg25)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxmg25
    integer:: i
    set_vxmg25_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxmg25(i) = vxmg25
        set_vxmg25_at_zone = 0
    end if
end function set_vxmg25_at_zone

integer function get_vxmg26_at_zone(index_of_the_particle, zone, vxmg26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxmg26
    integer:: i
    get_vxmg26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxmg26 = GenecStar%vxmg26(i)
        get_vxmg26_at_zone = 0
    end if
end function get_vxmg26_at_zone

integer function set_vxmg26_at_zone(index_of_the_particle, zone, vxmg26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxmg26
    integer:: i
    set_vxmg26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxmg26(i) = vxmg26
        set_vxmg26_at_zone = 0
    end if
end function set_vxmg26_at_zone

integer function get_vxf19_at_zone(index_of_the_particle, zone, vxf19)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxf19
    integer:: i
    get_vxf19_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxf19 = GenecStar%vxf19(i)
        get_vxf19_at_zone = 0
    end if
end function get_vxf19_at_zone

integer function set_vxf19_at_zone(index_of_the_particle, zone, vxf19)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxf19
    integer:: i
    set_vxf19_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxf19(i) = vxf19
        set_vxf19_at_zone = 0
    end if
end function set_vxf19_at_zone

integer function get_vxne21_at_zone(index_of_the_particle, zone, vxne21)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxne21
    integer:: i
    get_vxne21_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxne21 = GenecStar%vxne21(i)
        get_vxne21_at_zone = 0
    end if
end function get_vxne21_at_zone

integer function set_vxne21_at_zone(index_of_the_particle, zone, vxne21)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxne21
    integer:: i
    set_vxne21_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxne21(i) = vxne21
        set_vxne21_at_zone = 0
    end if
end function set_vxne21_at_zone

integer function get_vxna23_at_zone(index_of_the_particle, zone, vxna23)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxna23
    integer:: i
    get_vxna23_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxna23 = GenecStar%vxna23(i)
        get_vxna23_at_zone = 0
    end if
end function get_vxna23_at_zone

integer function set_vxna23_at_zone(index_of_the_particle, zone, vxna23)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxna23
    integer:: i
    set_vxna23_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxna23(i) = vxna23
        set_vxna23_at_zone = 0
    end if
end function set_vxna23_at_zone

integer function get_vxal27_at_zone(index_of_the_particle, zone, vxal27)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxal27
    integer:: i
    get_vxal27_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxal27 = GenecStar%vxal27(i)
        get_vxal27_at_zone = 0
    end if
end function get_vxal27_at_zone

integer function set_vxal27_at_zone(index_of_the_particle, zone, vxal27)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxal27
    integer:: i
    set_vxal27_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxal27(i) = vxal27
        set_vxal27_at_zone = 0
    end if
end function set_vxal27_at_zone

integer function get_vxsi28_at_zone(index_of_the_particle, zone, vxsi28)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxsi28
    integer:: i
    get_vxsi28_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxsi28 = GenecStar%vxsi28(i)
        get_vxsi28_at_zone = 0
    end if
end function get_vxsi28_at_zone

integer function set_vxsi28_at_zone(index_of_the_particle, zone, vxsi28)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxsi28
    integer:: i
    set_vxsi28_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxsi28(i) = vxsi28
        set_vxsi28_at_zone = 0
    end if
end function set_vxsi28_at_zone

integer function get_vxc14_at_zone(index_of_the_particle, zone, vxc14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxc14
    integer:: i
    get_vxc14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxc14 = GenecStar%vxc14(i)
        get_vxc14_at_zone = 0
    end if
end function get_vxc14_at_zone

integer function set_vxc14_at_zone(index_of_the_particle, zone, vxc14)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxc14
    integer:: i
    set_vxc14_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxc14(i) = vxc14
        set_vxc14_at_zone = 0
    end if
end function set_vxc14_at_zone

integer function get_vxf18_at_zone(index_of_the_particle, zone, vxf18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxf18
    integer:: i
    get_vxf18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxf18 = GenecStar%vxf18(i)
        get_vxf18_at_zone = 0
    end if
end function get_vxf18_at_zone

integer function set_vxf18_at_zone(index_of_the_particle, zone, vxf18)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxf18
    integer:: i
    set_vxf18_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxf18(i) = vxf18
        set_vxf18_at_zone = 0
    end if
end function set_vxf18_at_zone

integer function get_vxal26_at_zone(index_of_the_particle, zone, vxal26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxal26
    integer:: i
    get_vxal26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxal26 = GenecStar%vxal26(i)
        get_vxal26_at_zone = 0
    end if
end function get_vxal26_at_zone

integer function set_vxal26_at_zone(index_of_the_particle, zone, vxal26)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxal26
    integer:: i
    set_vxal26_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxal26(i) = vxal26
        set_vxal26_at_zone = 0
    end if
end function set_vxal26_at_zone

integer function get_vxneut_at_zone(index_of_the_particle, zone, vxneut)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxneut
    integer:: i
    get_vxneut_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxneut = GenecStar%vxneut(i)
        get_vxneut_at_zone = 0
    end if
end function get_vxneut_at_zone

integer function set_vxneut_at_zone(index_of_the_particle, zone, vxneut)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxneut
    integer:: i
    set_vxneut_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxneut(i) = vxneut
        set_vxneut_at_zone = 0
    end if
end function set_vxneut_at_zone

integer function get_vxprot_at_zone(index_of_the_particle, zone, vxprot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxprot
    integer:: i
    get_vxprot_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxprot = GenecStar%vxprot(i)
        get_vxprot_at_zone = 0
    end if
end function get_vxprot_at_zone

integer function set_vxprot_at_zone(index_of_the_particle, zone, vxprot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxprot
    integer:: i
    set_vxprot_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxprot(i) = vxprot
        set_vxprot_at_zone = 0
    end if
end function set_vxprot_at_zone

integer function get_vomegi_at_zone(index_of_the_particle, zone, vomegi)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vomegi
    integer:: i
    get_vomegi_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vomegi = GenecStar%vomegi(i)
        get_vomegi_at_zone = 0
    end if
end function get_vomegi_at_zone

integer function set_vomegi_at_zone(index_of_the_particle, zone, vomegi)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vomegi
    integer:: i
    set_vomegi_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vomegi(i) = vomegi
        set_vomegi_at_zone = 0
    end if
end function set_vomegi_at_zone

integer function get_vxbid_at_zone(index_of_the_particle, zone, vxbid)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxbid
    integer:: i
    get_vxbid_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxbid = GenecStar%vxbid(i)
        get_vxbid_at_zone = 0
    end if
end function get_vxbid_at_zone

integer function set_vxbid_at_zone(index_of_the_particle, zone, vxbid)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxbid
    integer:: i
    set_vxbid_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxbid(i) = vxbid
        set_vxbid_at_zone = 0
    end if
end function set_vxbid_at_zone

integer function get_vxbid1_at_zone(index_of_the_particle, zone, vxbid1)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: vxbid1
    integer:: i
    get_vxbid1_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vxbid1 = GenecStar%vxbid1(i)
        get_vxbid1_at_zone = 0
    end if
end function get_vxbid1_at_zone

integer function set_vxbid1_at_zone(index_of_the_particle, zone, vxbid1)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(in):: vxbid1
    integer:: i
    set_vxbid1_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vxbid1(i) = vxbid1
        set_vxbid1_at_zone = 0
    end if
end function set_vxbid1_at_zone

integer function get_abelx_at_zone(index_of_the_particle, zone, element, abelx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    integer, intent(in):: element
    real(kindreal), intent(out):: abelx
    integer:: i
    get_abelx_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        abelx = GenecStar%abelx(element,i)
        get_abelx_at_zone = 0
    end if
end function get_abelx_at_zone

integer function set_abelx_at_zone(index_of_the_particle, zone, element, abelx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    integer, intent(in):: element
    real(kindreal), intent(in):: abelx
    integer:: i
    set_abelx_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%abelx(element,i) = abelx
        set_abelx_at_zone = 0
    end if
end function set_abelx_at_zone

integer function get_vabelx_at_zone(index_of_the_particle, zone, element, vabelx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    integer, intent(in):: element
    real(kindreal), intent(out):: vabelx
    integer:: i
    get_vabelx_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        vabelx = GenecStar%vabelx(element,i)
        get_vabelx_at_zone = 0
    end if
end function get_vabelx_at_zone

integer function set_vabelx_at_zone(index_of_the_particle, zone, element, vabelx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    integer, intent(in):: element
    real(kindreal), intent(in):: vabelx
    integer:: i
    set_vabelx_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        GenecStar%vabelx(element,i) = vabelx
        set_vabelx_at_zone = 0
    end if
end function set_vabelx_at_zone

integer function get_drl(index_of_the_particle, number, drl)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(out):: drl
    drl = GenecStar%drl(number)
    get_drl = 0
end function get_drl

integer function set_drl(index_of_the_particle, number, drl)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(in):: drl
    GenecStar%drl(number) = drl
    set_drl = 0
end function set_drl

integer function get_drte(index_of_the_particle, number, drte)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(out):: drte
    drte = GenecStar%drte(number)
    get_drte = 0
end function get_drte

integer function set_drte(index_of_the_particle, number, drte)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(in):: drte
    GenecStar%drte(number) = drte
    set_drte = 0
end function set_drte

integer function get_drp(index_of_the_particle, number, drp)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(out):: drp
    drp = GenecStar%drp(number)
    get_drp = 0
end function get_drp

integer function set_drp(index_of_the_particle, number, drp)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(in):: drp
    GenecStar%drp(number) = drp
    set_drp = 0
end function set_drp

integer function get_drt(index_of_the_particle, number, drt)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(out):: drt
    drt = GenecStar%drt(number)
    get_drt = 0
end function get_drt

integer function set_drt(index_of_the_particle, number, drt)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(in):: drt
    GenecStar%drt(number) = drt
    set_drt = 0
end function set_drt

integer function get_drr(index_of_the_particle, number, drr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(out):: drr
    drr = GenecStar%drr(number)
    get_drr = 0
end function get_drr

integer function set_drr(index_of_the_particle, number, drr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: number
    real(kindreal), intent(in):: drr
    GenecStar%drr(number) = drr
    set_drr = 0
end function set_drr

integer function get_eps_at_zone(index_of_the_particle, zone, eps)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps
    integer:: i
    get_eps_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps = GenecStar%eps(i)
        get_eps_at_zone = 0
    end if
end function get_eps_at_zone

integer function get_epsy_at_zone(index_of_the_particle, zone, epsy)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: epsy
    integer:: i
    get_epsy_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        epsy = GenecStar%epsy(i)
        get_epsy_at_zone = 0
    end if
end function get_epsy_at_zone

integer function get_eps_c_adv_at_zone(index_of_the_particle, zone, eps_c_adv)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps_c_adv
    integer:: i
    get_eps_c_adv_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps_c_adv = GenecStar%eps_c_adv(i)
        get_eps_c_adv_at_zone = 0
    end if
end function get_eps_c_adv_at_zone

integer function get_eps_ne_adv_at_zone(index_of_the_particle, zone, eps_ne_adv)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps_ne_adv
    integer:: i
    get_eps_ne_adv_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps_ne_adv = GenecStar%eps_ne_adv(i)
        get_eps_ne_adv_at_zone = 0
    end if
end function get_eps_ne_adv_at_zone

integer function get_eps_o_adv_at_zone(index_of_the_particle, zone, eps_o_adv)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps_o_adv
    integer:: i
    get_eps_o_adv_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps_o_adv = GenecStar%eps_o_adv(i)
        get_eps_o_adv_at_zone = 0
    end if
end function get_eps_o_adv_at_zone

integer function get_eps_si_adv_at_zone(index_of_the_particle, zone, eps_si_adv)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps_si_adv
    integer:: i
    get_eps_si_adv_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps_si_adv = GenecStar%eps_si_adv(i)
        get_eps_si_adv_at_zone = 0
    end if
end function get_eps_si_adv_at_zone

integer function get_eps_grav_at_zone(index_of_the_particle, zone, eps_grav)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps_grav
    integer:: i
    get_eps_grav_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps_grav = GenecStar%eps_grav(i)
        get_eps_grav_at_zone = 0
    end if
end function get_eps_grav_at_zone

integer function get_eps_nu_at_zone(index_of_the_particle, zone, eps_nu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: eps_nu
    integer:: i
    get_eps_nu_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        eps_nu = GenecStar%eps_nu(i)
        get_eps_nu_at_zone = 0
    end if
end function get_eps_nu_at_zone

integer function get_nabla_rad_at_zone(index_of_the_particle, zone, nabla_rad)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: nabla_rad
    integer:: i
    get_nabla_rad_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        nabla_rad = GenecStar%nabla_rad(i)
        get_nabla_rad_at_zone = 0
    end if
end function get_nabla_rad_at_zone

integer function get_nabla_ad_at_zone(index_of_the_particle, zone, nabla_ad)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: nabla_ad
    integer:: i
    get_nabla_ad_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        nabla_ad = GenecStar%nabla_ad(i)
        get_nabla_ad_at_zone = 0
    end if
end function get_nabla_ad_at_zone

integer function get_nabla_mu_at_zone(index_of_the_particle, zone, nabla_mu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: zone
    real(kindreal), intent(out):: nabla_mu
    integer:: i
    get_nabla_mu_at_zone = -1
    i = GenecStar%m-zone
    if (.true.) then
        nabla_mu = GenecStar%nabla_mu(i)
        get_nabla_mu_at_zone = 0
    end if
end function get_nabla_mu_at_zone

integer function get_nbzel(index_of_the_particle, i, nbzel)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    integer, intent(out):: nbzel
    get_nbzel = -1
    if ((i > 0) .and. (i < 9)) then
        nbzel = GenecStar%nbzel(i)
        get_nbzel = 0
    end if
end function get_nbzel

integer function set_nbzel(index_of_the_particle, i, nbzel)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    integer, intent(in):: nbzel
    set_nbzel = -1
    if ((i > 0) .and. (i < 9)) then
        GenecStar%nbzel(i) = nbzel
        set_nbzel = 0
    end if
end function set_nbzel

integer function get_nbael(index_of_the_particle, i, nbael)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    integer, intent(out):: nbael
    get_nbael = -1
    if ((i > 0) .and. (i < 9)) then
        nbael = GenecStar%nbael(i)
        get_nbael = 0
    end if
end function get_nbael

integer function set_nbael(index_of_the_particle, i, nbael)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    integer, intent(in):: nbael
    set_nbael = -1
    if ((i > 0) .and. (i < 9)) then
        GenecStar%nbael(i) = nbael
        set_nbael = 0
    end if
end function set_nbael

integer function get_abels(index_of_the_particle, i, abels)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    real(kindreal), intent(out):: abels
    get_abels = -1
    if ((i > 0) .and. (i < 9)) then
        abels = GenecStar%abels(i)
        get_abels = 0
    end if
end function get_abels

integer function set_abels(index_of_the_particle, i, abels)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    real(kindreal), intent(in):: abels
    set_abels = -1
    if ((i > 0) .and. (i < 9)) then
        GenecStar%abels(i) = abels
        set_abels = 0
    end if
end function set_abels

integer function get_xnetalu(index_of_the_particle, i, xnetalu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    real(kindreal), intent(out):: xnetalu
    get_xnetalu = -1
    if ((i > 0) .and. (i < 6)) then
        xnetalu = GenecStar%xnetalu(i)
        get_xnetalu = 0
    end if
end function get_xnetalu

integer function set_xnetalu(index_of_the_particle, i, xnetalu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: i
    real(kindreal), intent(in):: xnetalu
    set_xnetalu = -1
    if ((i > 0) .and. (i < 6)) then
        GenecStar%xnetalu(i) = xnetalu
        set_xnetalu = 0
    end if
end function set_xnetalu

function new_stellar_model(&
      index_of_the_star,&
      modell,veryFirst,&
      initialised,star_name,nwmd,nwseq,modanf,nzmod,end_at_phase,end_at_model,&
      irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,var_rates,bintide,binm2,periodini,const_per,iprezams,&
      initial_metallicity,zsol,z,iopac,ikappa,&
      idiff,iadvec,istati,icoeff,fenerg,richac,igamma,frein,K_Kawaler,&
      Omega_saturation,rapcrilim,zams_velocity,xfom,omega,xdial,idialo,idialu,Add_Flux,diff_only,B_initial,&
      add_diff,n_mag,alpha_F,nsmooth,qminsmooth,&
      imloss,fmlos,ifitm,fitm,fitmi,fitmi_default,deltal,deltat,nndr,RSG_Mdot,SupraEddMdot,Be_mdotfrac,start_mdot,&
      iledou,idifcon,iover,elph,my,dovhp,iunder,dunder,&
      gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,nbchx,nrband,&
      xcn,islow,icncst,tauH_fit,&
      display_plot,iauto,iprn,iout,itmin,xyfiles,idebug,itests,verbose,stop_deg,n_snap,&
      gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,ab,dm_lost,m,summas,&
      dk,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,tdiff,suminenv,xltotbeg,dlelexprev,radius,zams_radius,& ! new!
      mbelx,xtefflast,xllast,xrholast,xclast,xtclast,inum,nsugi,period,r_core,vna,vnr,&
      q,p,t,r,s,x,y3,y,xc12,xc13,&
      xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,xc14,&
      xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,&
      vxo17,vxo18,vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,&
      vxal26,vxneut,vxprot,vomegi,vxbid,vxbid1,&
      xlostneu,&
      n&
      )
    implicit none
    integer:: index_of_the_star, n

    integer, intent(in):: &
             modell
    logical, intent(in):: &
             initialised, veryFirst
    character(256), intent(in):: &
             star_name
    integer, intent(in):: &
             nwmd,nwseq,modanf,nzmod,end_at_phase,end_at_model,irot,isol,imagn,&
             ialflu,ianiso,ipop3,ibasnet,phase
    logical, intent(in):: &
             var_rates,bintide
    real(kindreal), intent(in):: &
             binm2,periodini
    logical, intent(in):: &
             const_per
    integer, intent(in):: &
             iprezams
    real(kindreal), intent(in):: &
             initial_metallicity,&
             zsol,z
    integer, intent(in):: &
             iopac,ikappa,idiff,iadvec,istati,icoeff
    real(kindreal), intent(in):: &
             fenerg,richac
    integer, intent(in):: &
             igamma
    real(kindreal), intent(in):: &
             frein,K_Kawaler,&
             Omega_saturation,rapcrilim,zams_velocity,xfom,omega,xdial
    integer, intent(in):: &
             idialo,idialu
    logical, intent(in):: &
             Add_Flux,diff_only
    real(kindreal), intent(in):: &
             B_initial,&
             add_diff
    integer, intent(in):: &
             n_mag
    real(kindreal), intent(in):: &
             alpha_F
    integer, intent(in):: &
             nsmooth
    logical, intent(in):: &
             qminsmooth
    integer, intent(in):: &
             imloss
    real(kindreal), intent(in):: &
             fmlos
    integer, intent(in):: &
             ifitm
    real(kindreal), intent(in):: &
             fitm,fitmi,fitmi_default,deltal,deltat
    integer, intent(in):: &
             nndr,&
             RSG_Mdot
    logical, intent(in):: &
             SupraEddMdot
    real(kindreal), intent(in):: &
             Be_mdotfrac,start_mdot
    integer, intent(in):: &
             iledou,idifcon,iover
    real(kindreal), intent(in):: &
             elph
    integer, intent(in):: &
             my
    real(kindreal), intent(in):: &
             dovhp
    integer, intent(in):: &
             iunder
    real(kindreal), intent(in):: &
             dunder,&
             gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20
    integer, intent(in):: &
             nbchx,nrband
    real(kindreal), intent(in):: &
             xcn
    integer, intent(in):: &
             islow,icncst,&
             tauH_fit
    logical, intent(in):: &
             display_plot
    integer, intent(in):: &
             iauto,iprn,iout,itmin
    logical, intent(in):: &
             xyfiles
    integer, intent(in):: &
             idebug,itests
    logical, intent(in):: &
             verbose,stop_deg
    integer, intent(in):: &
             n_snap,m
    real(kindreal), intent(in):: &
             gms,&
             alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,dm_lost
    real(kindreal), intent(in):: &  ! FIXME
            dk,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,tdiff,suminenv,xltotbeg,dlelexprev,radius,zams_radius
    integer, intent(in):: &
             mbelx
    real(kindreal), intent(in):: &
             xtefflast,&
             xllast,xrholast,xclast,xtclast
    integer, intent(in):: &
             inum,nsugi
    real(kindreal), intent(in):: &
             period,r_core,vna,vnr
    real(kindreal), dimension(m), intent(in):: &
             q,p,t,r,s,x,y3,y,xc12,&
             xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,&
             xsi28,xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,&
             vxn14,vxn15,vxo16,vxo17,vxo18,vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,&
             vxal27,vxsi28,vxc14,vxf18,vxal26,vxneut,vxprot,vomegi,vxbid,vxbid1
    real(kindreal), intent(in):: &
            xlostneu
            
    integer:: new_stellar_model

    write(*,*) "*****N: ", n

    GenecStar%modell = modell
    GenecStar%veryFirst = veryFirst
    GenecStar%initialised = initialised
    GenecStar%star_name = star_name
    GenecStar%nwmd = nwmd
    GenecStar%nwseq = nwseq
    GenecStar%modanf = modanf
    GenecStar%nzmod = nzmod
    GenecStar%end_at_phase = end_at_phase
    GenecStar%end_at_model = end_at_model
    !----
    GenecStar%irot = irot
    GenecStar%isol = isol
    GenecStar%imagn = imagn
    GenecStar%ialflu = ialflu
    GenecStar%ianiso = ianiso
    GenecStar%ipop3 = ipop3
    GenecStar%ibasnet = ibasnet
    GenecStar%phase = phase
    GenecStar%var_rates = var_rates
    GenecStar%bintide = bintide
    GenecStar%binm2 = binm2
    GenecStar%periodini = periodini
    GenecStar%const_per = const_per
    GenecStar%iprezams = iprezams
    GenecStar%initial_metallicity = initial_metallicity
    GenecStar%zsol = zsol
    GenecStar%z = z
    GenecStar%iopac = iopac
    GenecStar%ikappa = ikappa
    !----
    GenecStar%idiff = idiff
    GenecStar%iadvec = iadvec
    GenecStar%istati = istati
    GenecStar%icoeff = icoeff
    GenecStar%fenerg = fenerg
    GenecStar%richac = richac
    GenecStar%igamma = igamma
    GenecStar%frein = frein
    GenecStar%K_Kawaler = K_Kawaler
    !----
    GenecStar%Omega_saturation = Omega_saturation
    GenecStar%rapcrilim = rapcrilim
    GenecStar%zams_velocity = zams_velocity
    GenecStar%xfom = xfom
    GenecStar%omega = omega
    GenecStar%xdial = xdial
    GenecStar%idialo = idialo
    GenecStar%idialu = idialu
    GenecStar%Add_Flux = Add_Flux
    GenecStar%diff_only = diff_only
    GenecStar%B_initial = B_initial
    !----
    GenecStar%add_diff = add_diff
    GenecStar%n_mag = n_mag
    GenecStar%alpha_F = alpha_F
    GenecStar%nsmooth = nsmooth
    GenecStar%qminsmooth = qminsmooth
    !----
    GenecStar%imloss = imloss
    GenecStar%fmlos = fmlos
    GenecStar%ifitm = ifitm
    GenecStar%fitm = fitm
    GenecStar%fitmi = fitmi
    GenecStar%fitmi_default = fitmi_default
    GenecStar%deltal = deltal
    GenecStar%deltat = deltat
    GenecStar%nndr = nndr
    GenecStar%RSG_Mdot = RSG_Mdot
    GenecStar%SupraEddMdot = SupraEddMdot
    GenecStar%Be_mdotfrac = Be_mdotfrac
    GenecStar%start_mdot = start_mdot
    !----
    GenecStar%iledou = iledou
    GenecStar%idifcon = idifcon
    GenecStar%iover = iover
    GenecStar%elph = elph
    GenecStar%my = my
    GenecStar%dovhp = dovhp
    GenecStar%iunder = iunder
    GenecStar%dunder = dunder
    !----
    GenecStar%gkorm = gkorm
    GenecStar%alph = alph
    GenecStar%agdr = agdr
    GenecStar%faktor = faktor
    GenecStar%dgrp = dgrp
    GenecStar%dgrl = dgrl
    GenecStar%dgry = dgry
    GenecStar%dgrc = dgrc
    GenecStar%dgro = dgro
    GenecStar%dgr20 = dgr20
    GenecStar%nbchx = nbchx
    GenecStar%nrband = nrband
    !----
    GenecStar%xcn = xcn
    GenecStar%islow = islow
    GenecStar%icncst = icncst
    GenecStar%tauH_fit = tauH_fit
    !----
    GenecStar%display_plot = display_plot
    GenecStar%iauto = iauto
    GenecStar%iprn = iprn
    GenecStar%iout = iout
    GenecStar%itmin = itmin
    GenecStar%xyfiles = xyfiles
    GenecStar%idebug = idebug
    GenecStar%itests = itests
    GenecStar%verbose = verbose
    GenecStar%stop_deg = stop_deg
    GenecStar%n_snap = n_snap
    !----
    GenecStar%gms = gms
    GenecStar%alter = alter
    GenecStar%gls = gls
    GenecStar%teff = teff
    GenecStar%glsv = glsv
    GenecStar%teffv = teffv
    GenecStar%dzeitj = dzeitj
    GenecStar%dzeit = dzeit
    GenecStar%dzeitv = dzeitv
    GenecStar%xmini = xmini
    GenecStar%ab = ab
    GenecStar%dm_lost = dm_lost
    GenecStar%m = m
    GenecStar%summas = summas
    !----
    GenecStar%dk = dk
    GenecStar%rlp = rlp
    GenecStar%rlt = rlt
    GenecStar%rlc = rlc
    GenecStar%rrp = rrp
    GenecStar%rrt = rrt
    GenecStar%rrc = rrc
    GenecStar%rtp = rtp
    GenecStar%rtt = rtt
    GenecStar%rtc = rtc
    GenecStar%tdiff = tdiff
    GenecStar%suminenv = suminenv
    GenecStar%xltotbeg = xltotbeg
    GenecStar%dlelexprev = dlelexprev
    GenecStar%radius = log10(radius)
    GenecStar%zams_radius = zams_radius
    !----

    GenecStar%mbelx = mbelx
    GenecStar%xtefflast = xtefflast
    GenecStar%xllast = xllast
    GenecStar%xrholast = xrholast
    GenecStar%xclast = xclast
    GenecStar%xtclast = xtclast
    GenecStar%inum = inum
    GenecStar%nsugi = nsugi
    GenecStar%period = period
    GenecStar%r_core = r_core
    GenecStar%vna = vna
    GenecStar%vnr = vnr
    !----
    GenecStar%q(:m) = q
    GenecStar%p(:m) = p
    GenecStar%t(:m) = t
    GenecStar%r(:m) = r
    GenecStar%s(:m) = s
    GenecStar%x(:m) = x
    GenecStar%y3(:m) = y3
    GenecStar%y(:m) = y
    GenecStar%xc12(:m) = xc12
    GenecStar%xc13(:m) = xc13
    GenecStar%xn14(:m) = xn14
    GenecStar%xn15(:m) = xn15
    GenecStar%xo16(:m) = xo16
    GenecStar%xo17(:m) = xo17
    GenecStar%xo18(:m) = xo18
    GenecStar%xne20(:m) = xne20
    GenecStar%xne22(:m) = xne22
    GenecStar%xmg24(:m) = xmg24
    GenecStar%xmg25(:m) = xmg25
    GenecStar%xmg26(:m) = xmg26
    GenecStar%xf19(:m) = xf19
    GenecStar%xne21(:m) = xne21
    GenecStar%xna23(:m) = xna23
    GenecStar%xal27(:m) = xal27
    GenecStar%xsi28(:m) = xsi28
    GenecStar%xc14(:m) = xc14
    GenecStar%xf18(:m) = xf18
    GenecStar%xal26(:m) = xal26
    GenecStar%xneut(:m) = xneut
    GenecStar%xprot(:m) = xprot
    GenecStar%omegi(:m) = omegi
    !----
    GenecStar%xbid(:m) = xbid
    GenecStar%xbid1(:m) = xbid1
    GenecStar%vp(:m) = vp
    GenecStar%vt(:m) = vt
    GenecStar%vr(:m) = vr
    GenecStar%vs(:m) = vs
    GenecStar%vx(:m) = vx
    GenecStar%vy(:m) = vy
    GenecStar%vy3(:m) = vy3
    GenecStar%vxc12(:m) = vxc12
    GenecStar%vxc13(:m) = vxc13
    GenecStar%vxn14(:m) = vxn14
    GenecStar%vxn15(:m) = vxn15
    GenecStar%vxo16(:m) = vxo16
    GenecStar%vxo17(:m) = vxo17
    GenecStar%vxo18(:m) = vxo18
    GenecStar%vxne20(:m) = vxne20
    GenecStar%vxne22(:m) = vxne22
    GenecStar%vxmg24(:m) = vxmg24
    GenecStar%vxmg25(:m) = vxmg25
    GenecStar%vxmg26(:m) = vxmg26
    GenecStar%vxf19(:m) = vxf19
    GenecStar%vxne21(:m) = vxne21
    GenecStar%vxna23(:m) = vxna23
    GenecStar%vxal27(:m) = vxal27
    GenecStar%vxsi28(:m) = vxsi28
    GenecStar%vxc14(:m) = vxc14
    GenecStar%vxf18(:m) = vxf18
    GenecStar%vxal26(:m) = vxal26
    GenecStar%vxneut(:m) = vxneut
    GenecStar%vxprot(:m) = vxprot
    GenecStar%vomegi(:m) = vomegi
    GenecStar%vxbid(:m) = vxbid
    GenecStar%vxbid1(:m) = vxbid1
    GenecStar%xlostneu = xlostneu
    !GenecStar%abelx(:,:) = 0.
    !GenecStar%vabelx(:,:) = 0.

    !call copy_from_genec_star(GenecStar)
    new_stellar_model = 0
end function new_stellar_model


!!!!!!!!!
end module AmuseInterface
