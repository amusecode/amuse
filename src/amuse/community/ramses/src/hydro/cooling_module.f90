!=============================================================================
!                   SUPERMEGAGIGAMODULECOOLINGQUIDEPOTE
!=============================================================================
! Les subroutines et fonctions d'interet general sont :
!
! ROUTINES A APPELER PAR LE CODE HYDRO
!
!    subroutine set_model(...) : 
!          pour choisir le modele de cooling et ses parametres
!
!    subroutine set_table(aexp) : 
!          pour creer la table avec les parametres par defaut 
!          Plus pratique a appeler que 
!          cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
!
!    subroutine solve_cooling(...) : 
!          pour calculer le cooling
!
! ROUTINE A MODIFIER SI NECESSAIRE
!
!    function J0simple(aexp) : 
!          donne le J0 en fonction du redshift dans les modeles Teyssier 
!          ou Theuns
!
!=============================================================================
module cooling_module
  use amr_parameters
  implicit none
  logical :: verbose_cooling=.false.

  real(kind=8),parameter ::smallnum_cooling= 1d-30
  real(kind=8),parameter ::twopi   = 6.2831853d0
  real(kind=8),parameter ::hplanck = 6.6262000d-27
  real(kind=8),parameter ::eV      = 1.6022000d-12
  real(kind=8),parameter ::kB      = 1.3806200d-16
  real(kind=8),parameter ::clight  = 2.9979250d+10
  real(kind=8),parameter ::Gyr     = 3.1536000d+16
  real(kind=8)           ::X       = 0.76
  real(kind=8)           ::Y       = 0.24 
  real(kind=8),parameter ::rhoc    = 1.8800000d-29
  real(kind=8),parameter ::mH      = 1.6600000d-24
  real(kind=8),parameter ::mu_mol  = 1.2195D0
  integer,parameter::HI      = 1
  integer,parameter::HEI     = 2
  integer,parameter::HEII    = 3

  ! Les parametres de la table par defaut
  integer,parameter     :: nbin_T_fix=101
  integer,parameter     :: nbin_n_fix=161
  real(kind=8),parameter:: nH_min_fix=1.d-10
  real(kind=8),parameter:: nH_max_fix=1.d+6
  real(kind=8),parameter:: T2_min_fix=1.d-2
  real(kind=8),parameter:: T2_max_fix=1.d+9
  
  type cooling_table
     integer::n1
     integer::n2
     real(kind=8),dimension(:)    ,pointer::nH
     real(kind=8),dimension(:)    ,pointer::T2
     real(kind=8),dimension(:,:)  ,pointer::cool
     real(kind=8),dimension(:,:)  ,pointer::heat
     real(kind=8),dimension(:,:)  ,pointer::cool_com
     real(kind=8),dimension(:,:)  ,pointer::heat_com
     real(kind=8),dimension(:,:)  ,pointer::metal
     real(kind=8),dimension(:,:)  ,pointer::cool_prime
     real(kind=8),dimension(:,:)  ,pointer::heat_prime
     real(kind=8),dimension(:,:)  ,pointer::cool_com_prime
     real(kind=8),dimension(:,:)  ,pointer::heat_com_prime
     real(kind=8),dimension(:,:)  ,pointer::metal_prime
     real(kind=8),dimension(:,:)  ,pointer::mu
     real(kind=8),dimension(:,:,:),pointer::n_spec
  end type cooling_table

  type(cooling_table)::table,table2
  ! Utilisation de table%n_spec si necessaire
  logical, parameter :: if_species_abundances=.true. 
  ! Facteur correctif de Theuns et al.
  real(kind=8),parameter :: dumfac_ion_theuns=2.d0 
  real(kind=8),parameter :: dumfac_rec_theuns=0.75D0  ! idem
  real(kind=8) :: dumfac_ion=dumfac_ion_theuns
  real(kind=8) :: dumfac_rec=dumfac_rec_theuns

  ! On DOIT AVOIR OU teyssier OU theuns OU madau 
  ! OU weinbergint OU courty avec un OU exclusif
  logical :: teyssier=.false.         
  logical :: theuns=.false.           
  logical :: madau=.false.           
  logical :: weinberg=.false. 
  logical :: weinbergint=.false. 
  logical :: courty=.true.  ! Default model

  ! Si teyssier ou theuns :
  real(kind=8) :: J0in=1.d-22  ! J0 default 
  real(kind=8) :: J0min=1.d-29 ! Valeur minimale du J0 
  logical :: force_j0_one=.false. ! Force constant UV bkg
  ! (saturation a grand redshift)
  real(kind=8) :: aexp_ref=0.0001        
  real(kind=8) :: J0min_ref=2.77168510365299962D-25 ! J0min_ref precalcule pour
  ! H0=70, omegab=0.04, omega0=0.3, omegaL=0.7
  logical :: high_z_realistic_ne=.true. ! Calcul du J0min de telle sorte 
  ! que le n_e soit realiste a grand z. J0min=J0min_ref/(aexp/aexp_ref)^2
  real(kind=8) :: alpha=1.d0   ! J(nu) \propto nu^{-alpha} 
  ! Si madau ou weinbergint :
  real(kind=8) :: normfacJ0=0.74627   ! Facteur de normalisation pour J0 
  ! pour un J(nu,z) de type haardt et Madau
  ! Ce facteur la est celui utilise par Dave et al. pour LCDM
  ! Sauvegarde des termes de cooling/heating dans les
  logical, parameter :: if_cooling_functions=.true. 
  ! variables en dessous
  real(kind=8)::cb1s,cb2s,cb3s,ci1s,ci2s,ci3s,cr1s,cr2s,cr3s,cds
  real(kind=8)::ce1s,ce3s,ch1s,ch2s,ch3s,cocs,cohs
  real(kind=8)::cool_out, heat_out

  ! Les heating et photoionization rates de Dave et al. 
  ! pour le J0 derniere version de HM (weinberg ou weinbergint si
  ! if_read_weinberg=.true. (voir plus bas) dans ce dernier cas)
  real(kind=8),allocatable, dimension(:,:)::table_weinberg 
  ! Table d'interpolation en input
  character(len=128), parameter :: table_weinberg_name='TREECOOL' 
  ! Nom du fichier avec les donnees
  integer,parameter :: luweinberg=21                        
  ! unit pour lire le fichier
  integer :: Nweinberg                                    
  ! Nombre de bins en redshift

  ! Les coefficients d'interpolation des heating rates de Dave et al.
  ! (weinbergint)
  logical,parameter :: if_read_weinberg=.false. 
  ! .true. pour lire le fichier table_weinberg_name
  ! puis interpoler par un polynome
  ! .false. pour utiliser les valeurs des coefficients
  ! precalcules listes plus bas
  integer,parameter :: Norderweinberg=7       
  ! Ordre+1 du polynome d'interpolation (NE PAS CHANGER)
  real(kind=8) :: coefweinberg(Norderweinberg,6)= reshape( &
 &                    (/ -0.31086729929951613D+002, 0.34803667059463761D+001,-0.15145716066316397D+001, &
 &                        0.54649951450632972D+000,-0.16395924120387340D+000, 0.25197466148524143D-001, &
 &                       -0.15352763785487806D-002, &
 &                       -0.31887274113252204D+002, 0.44178493140927095D+001,-0.20158132553082293D+001, &  
 &                        0.64080497292269134D+000,-0.15981267091909040D+000, 0.22056900050237707D-001, &
 &                       -0.12837570029562849D-002, &
 &                       -0.35693331167978656D+002, 0.20207245722165794D+001,-0.76856976101363744D-001, &
 &                       -0.75691470654320359D-001,-0.54502220282734729D-001, 0.20633345104660583D-001, & 
 &                       -0.18410307456285177D-002, &
 &                       -0.56967559787460921D+002, 0.38601174525546353D+001,-0.18318926655684415D+001, &
 &                        0.67360594266440688D+000,-0.18983466813215341D+000, 0.27768907786915147D-001, &
 &                       -0.16330066969315893D-002, &
 &                       -0.56977907250821026D+002, 0.38686249565302266D+001,-0.13330942368518774D+001, &  
 &                        0.33988839029092172D+000,-0.98997915675929332D-001, 0.16781612113050747D-001, &
 &                       -0.11514328893746039D-002, &
 &                       -0.59825233828609278D+002, 0.21898162706563347D+001,-0.42982055888598525D+000, &
 &                        0.50312144291614215D-001,-0.61550639239553132D-001, 0.18017109270959387D-001, & 
 &                       -0.15438891584271634D-002 /), (/Norderweinberg,6/) )

  real(kind=8) :: zreioniz=8.5d0
  integer,parameter :: Nordercourty=7       
  ! Ordre+1 du polynome d'interpolation (NE PAS CHANGER)
  real(kind=8) :: coefcourty(0:Nordercourty,6)= reshape( &
                      (/ -13.5857,  1.24475,    0.187739, &
                        -0.430409, 0.152544,  -0.0246448, &
                         0.00192622, -5.89772e-05, &
                        -14.0242,  1.99211,   -0.490766, &
                        -0.122646, 0.0776501, -0.0146310, &
                         0.00123335, -3.96066e-05, &
                        -15.6627,  0.128240,   1.65633, &
                        -1.23799,  0.372157,  -0.0561687, &
                         0.00422696, -0.000126344, &
                        -24.8422,  1.50750,   -0.0699428, &
                        -0.308682, 0.122196,  -0.0205179, &
                         0.00163695, -5.08050e-05, &
                        -25.0252,  1.79577,   -0.159054, &
                        -0.300924, 0.125343,  -0.0214598, &
                         0.00173377, -5.43576e-05, &
                        -26.4168,  0.0479454,  1.70948, &
                        -1.26395,  0.378922,  -0.0570957, &
                         0.00428897, -0.000127909 /),(/Nordercourty+1,6/) )
  real(kind=8),dimension(6)    :: coef_fit= (/ 20., 20., 20., 20., 20., 20. /) 
  integer,dimension(6) :: beta_fit= (/  6,   6,   8,   6,   6,  8  /)

contains 
!=======================================================================
subroutine set_model(Nmodel,J0in_in,J0min_in,alpha_in,normfacJ0_in,zreioniz_in, &
 &                   correct_cooling,realistic_ne, &
 &                   h,omegab,omega0,omegaL,astart_sim,T2_sim)
!=======================================================================
! Nmodel(integer) =1 : Teyssier : ancien choix de l'evolution et de la forme du J(nu,z)
!                  2 : Theuns   : pareil mais avec les fonctions interpolees de Theuns (+ rapide)
!                  3 : Madau    : J(nu,z) de Theuns et al. 1998 avec les anciennes mesures de 
!                                 Haardt et Madau (HM)
!                  4 : Weinberg : J(nu,z) de Dave et al. 1999 avec les nouvelles mesure de HM 
!                                 lues dans le fichier table_weinberg_name (inactive)
!                  5 : idem 4 mais interpole interpole de maniere polynomiale : RECOMMANDE
!                  6 : Courty
!                 -1 : defaut defini dans le module 
! J0in_in (dble) : valeur du J0 utilisee pour Teyssier et Theuns
!            Exemple : J0in_in=1.d-22
!            J0in_in <= 0 utilise le defaut defini dans le module
! J0min_in (dble) : valeur du J0min ou J0min_ref (voir option realistic_ne) 
!            utilisee dans tous les modeles a grand redshift 
!            Exemple : J0min_in=1.d-29
!            J0min_in <= 0 utilise le defaut defini dans le module
! alpha_in (dble) : valeur de l'indice spectral du J(nu) \propto nu^{-alpha}
!            Exemple : alpha=1.
!            alpha_in < 0 utilise le defaut defini dans le module
! zreioniz_in (dble) : valeur du redshift de reionisation
!            Exemple : zerion=10.
!            zreioniz_in < 0 utilise le defaut defini dans le module
! normfacJ0_in (dble) : valeur du facteur de normalisation dans le cas des
!            spectres de Haardt et Madau. C'est un nombre de l'ordre de
!            l'unite en general plus petit que 1.
!            Exemple : normfacJ0_in=0.74627
!            normfacJ0_in prend le defaut defini dans le module
! correct_cooling (integer) : 0 : pas de correction
!                             1 : correction de Theuns et al 98
!                            -1 : defaut defini dans le module
! realistic_ne (integer) : 0 : pas de n_e realiste a grand redshift :
!                              Le J0min reste le meme quel que soit le redshift
!                              (J0min=J0min_in si celui-ci est > 0)
!                          1 : n_e realiste a grand redshift : J0min proportionnel a 1/a^2 
!                              egal initialement a J0min_ref pour a=aexp_ref=0.0001
!                              (J0min_ref=J0min_in si celui-ci est > 0)
!                          2 : RECOMMANDE : pareil que 1, mais J0min_ref est calcule de 
!                              maniere iterative pour avoir le bon n_e a z=19. 
!                              Le J0min_in n'est pas relevant dans ce cas la. 
! h (dble)          : H0/100
! omegab (dble)     : omega baryons
! omega0 (dble)     : omega matiere total
! omegaL (dble)     : omega Lambda
! astart_sim (dble) : redshift auquel on veut commencer la simulation
! T2_sim     (dble) : ce sera en output, le T/mu en K a ce redshift pour des regions de contraste
!                     de densite nul. 
!
! NOTE :
! Dans les cas madau, ou weinberg ou weinbergint, le J0 a grand redshift est calcule comme 
! dans l'option theuns :
!   madau :       pour z >= 15 ou quand le taux trouve est plus petit que celui donne par 
!                 l'option theuns=.true.
!   weinberg :    quand on sort de la table des taux
!   weinbergint : pour z >= 8.5 ou quand le taux trouve est plus petit que celui donne 
!                 par l'option theuns=.true.
!   courty : 
!=======================================================================
  implicit none
  real(kind=8) :: J0in_in,zreioniz_in,J0min_in,alpha_in,normfacJ0_in,astart_sim,T2_sim
  real(kind=8) :: J0min_ref_calc,h,omegab,omega0,omegaL
  integer :: Nmodel,correct_cooling,realistic_ne
  real(kind=8) :: astart,aend,dasura,T2end,mu,ne,minus1
  if (Nmodel /= -1) then
     teyssier=.false.
     theuns=.false.
     madau=.false.
     weinberg=.false.
     weinbergint=.false.
     courty=.false.
     if (Nmodel==1) then
        teyssier=.true.
     elseif (Nmodel==2) then
        theuns=.true.
     elseif (Nmodel==3) then
        madau=.true.
     elseif (Nmodel==4) then
        weinberg=.true.
     elseif (Nmodel==5) then
        weinbergint=.true.
     elseif (Nmodel==6) then
        courty=.true.
     else
        write(*,*) 'ERROR in set_model : wrong value of Nmodel'
        write(*,*) 'Nmodel =',Nmodel
        STOP
     endif
  endif
  if (J0in_in >= 0.0) J0in=J0in_in
  if (zreioniz_in >= 0.0) zreioniz=zreioniz_in
  if (alpha_in > 0.0) alpha=alpha_in
  if (normfacJ0_in > 0.0) normfacJ0=normfacJ0_in
  if (correct_cooling == 0) then
     dumfac_ion=1.d0
     dumfac_rec=1.d0
  elseif (correct_cooling == 1) then
     dumfac_ion=dumfac_ion_theuns
     dumfac_rec=dumfac_rec_theuns
  elseif (correct_cooling /= -1) then
     write(*,*) 'ERROR in set_model : wrong value of correct_cooling'
     write(*,*) 'correct_cooling =',correct_cooling
     STOP
  endif
  if (realistic_ne == 0) then
     astart=5.d-4
     high_z_realistic_ne=.false.
     if (J0min_in > 0.d0) J0min=J0min_in
  elseif (realistic_ne == 1) then
     astart=aexp_ref
     high_z_realistic_ne=.true.
     if (J0min_in > 0.d0) J0min_ref=J0min_in
  elseif (realistic_ne == 2) then
     astart=aexp_ref
     high_z_realistic_ne=.true.
     call compute_J0min(h,omegab,omega0,omegaL,J0min_ref_calc)
     J0min_ref=J0min_ref_calc
  else
     write(*,*) 'ERROR in set_model : wrong value of realistic_ne'
     write(*,*) 'realistic_ne =',realistic_ne
     STOP
  endif
  if (astart_sim < astart) then
     write(*,*) 'ERROR in set_model : astart_sim is too small.'
     write(*,*) 'astart     =',astart
     write(*,*) 'astart_sim =',astart_sim
     STOP
  endif
  ! Calcul de la temperature initiale
  aend=astart_sim
  dasura=0.02d0
  minus1=-1.0
  call evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL,minus1,T2end,mu,ne,.false.)
  if (verbose_cooling) write(*,*) 'Starting temperature in K :',T2end*mu
  T2_sim=T2end 
end subroutine set_model
!=======================================================================
subroutine set_table(aexp)
!=======================================================================
  implicit none
  real(kind=8) :: aexp
  integer :: nbin_n,nbin_T
  real(kind=8) :: nH_min,nH_max,T2_min,T2_max
  nH_min=nH_min_fix
  nH_max=nH_max_fix
  T2_min=T2_min_fix
  T2_max=T2_max_fix
  nbin_n=nbin_n_fix
  nbin_T=nbin_T_fix
  call cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
end subroutine set_table
!=======================================================================
subroutine output_cool(filename)
!=======================================================================
  implicit none
  character(LEN=80)::filename
  open(unit=10,file=filename,form='unformatted')
  write(10)table%n1,table%n2
  write(10)table%nH
  write(10)table%T2
  write(10)table%cool
  write(10)table%heat
  write(10)table%cool_com
  write(10)table%heat_com
  write(10)table%metal
  write(10)table%cool_prime
  write(10)table%heat_prime
  write(10)table%cool_com_prime
  write(10)table%heat_com_prime
  write(10)table%metal_prime
  write(10)table%mu
  if (if_species_abundances) write(10)table%n_spec
  close(10)
end subroutine output_cool
!=======================================================================
subroutine evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL, &
 &                          J0min_in,T2end,mu,ne,if_write_result)
!=======================================================================
! astart : valeur du facteur d'expansion au debut du calcul
! aend   : valeur du facteur d'expansion a la fin du calcul
! dasura : la valeur de da/a entre 2 pas de temps
! h      : la valeur de H0/100 
! omegab : la valeur de Omega baryons
! omega0 : la valeur de Omega matiere (total)
! omegaL : la valeur de Omega Lambda
! J0min_in : la valeur du J0min a injecter :
!          Si high_z_realistic_ne alors c'est J0min a a=astart qui
!          est considere
!          Sinon, c'est le J0min habituel.
!          Si J0min_in <=0, les parametres par defaut ou predefinis
!          auparavant sont pris pour le J0min.
! T2end  : Le T/mu en output
! mu     : le poids moleculaire en output
! ne     : le ne en output
! if_write_result : .true. pour ecrire l'evolution de la temperature
!          et de n_e sur l'ecran.
!=======================================================================
  implicit none
  real(kind=8)::astart,aend,T2end,h,omegab,omega0,omegaL,J0min_in,ne,dasura
  logical :: if_write_result
  real(kind=8)::aexp,daexp,dt_cool,coeff,coeff2
  real(kind=8)::T2_com,T2_old,T2,T2_left,T2_right,err_T2
  real(kind=8)::nH_com,nH  
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) ::mu
  real(kind=8) ::cool_tot,heat_tot,cool_com,heat_com
  real(kind=8) ::diff
  integer::niter
  real(kind=8) :: n_spec(1:6)
  if (J0min_in > 0.0) then
     if (high_z_realistic_ne) then
        J0min_ref = J0min_in
        aexp_ref = astart
     else
        J0min = J0min_in
     endif
  endif
  aexp = astart
  T2_com = 2.726d0 / aexp * aexp**2 / mu_mol
  nH_com = omegab*rhoc*h**2*X/mH
  do while (aexp < aend)
     daexp = dasura*aexp
     dt_cool=daexp/(aexp*100.*h*3.2408608e-20*HsurH0(1.0/aexp-1.,omega0,omegaL,1.-omega0-omegaL))
     
     nH = nH_com/aexp**3
     T2_old = T2_com/aexp**2

     ! Compute radiative ionization and heating rates
     call set_rates(t_rad_spec,h_rad_spec,aexp)
     
     ! Iteration to find new T2
     err_T2=1.
     T2_left=1.d-2
     T2_right=1.d8
     niter=0
     coeff = 2.*nH*X/3./kB
     coeff2 = 2.*X/3./kB
     do while (err_T2 > 1.d-10.and.niter <= 100)
        T2=0.5*(T2_left+T2_right)        
        call cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu,aexp,n_spec)
        diff = coeff*(heat_tot-cool_tot) + coeff2*(heat_com-cool_com) + (T2_old-T2)/dt_cool
        if(diff>0.)then 
           T2_left =0.5*(T2_left+T2_right)
           T2_right=T2_right
        else
           T2_left =T2_left
           T2_right=0.5*(T2_left+T2_right)
        end if
        err_T2=abs(T2_right-T2_left)/T2_left
        niter=niter+1
     end do
     if (niter > 100) then
        write(*,*) 'ERROR in evol_single_cell : too many iterations'
        STOP
     endif
     T2_com=T2*aexp**2
     aexp = aexp + daexp
     if (if_write_result) write(*,'(4(1pe10.3))')aexp,nH,T2_com*mu/aexp**2,n_spec(1)/nH
  end do
  T2end=T2
  ne=n_spec(1)/nH
end subroutine evol_single_cell
!=======================================================================
subroutine compute_J0min(h,omegab,omega0,omegaL,J0min_in)
!=======================================================================
  implicit none
  real(kind=8) :: omega0,omegaL,h,omegab,ne_to_find,mu
  real(kind=8) :: h0,astart,aend,J0min_in,T2end,ne
  real(kind=8) :: J0min_left,J0min_right,err_J0min,diff,xval,dasura
  integer :: niter
  logical :: if_write_result=.false.

  xval=sqrt(omega0)/(h*omegab)
  ne_to_find=1.2d-5*xval ! From the book of Peebles p. 173
  astart=aexp_ref
  aend=MIN(0.05d0,0.5d0/(1d0+zreioniz)) ! Always end before reionization
  dasura=0.05
  err_J0min=1.
  J0min_left=1d-20
  J0min_right=1d-30
  niter=0
  do while (err_J0min > 1.d-3 .and. niter <= 100)
     J0min_in=0.5*(J0min_left+J0min_right)     
     call evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL,J0min_in,T2end,mu,ne,if_write_result)
     diff=ne-ne_to_find
     if (diff>0.d0) then
        J0min_left=0.5*(J0min_left+J0min_right)
        J0min_right=J0min_right
     else
        J0min_left=J0min_left
        J0min_right=0.5*(J0min_left+J0min_right)
     endif
     err_J0min=abs(J0min_right-J0min_left)/J0min_left
     niter=niter+1
  enddo
  if (niter > 100) then
     write(*,*) 'ERROR in compute_J0min : too many iterations'
     STOP
  endif
  if (verbose_cooling)  write(*,*) 'J0min found ',J0min_in
end subroutine compute_J0min
!=======================================================================
subroutine solve_cooling(nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none  
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost
    
  real(kind=8)::facT,dlog_nH,dlog_T2,coeff,precoeff,h,h2,h3
  real(kind=8)::metal,cool,heat,cool_com,heat_com,w1T,w2T,w11,w12,w21,w22,err,yy,yy2,yy3
  real(kind=8)::metal_prime,cool_prime,heat_prime,cool_com_prime,heat_com_prime,wcool
  real(kind=8)::lambda,lambda_prime,logT2max
  real(kind=8)::fa,fb,fprimea,fprimeb,alpha,beta,gamma
  real(kind=8),dimension(1:ncell)::rgt,lft,tau,tau_old
  real(kind=8),dimension(1:ncell)::time,time_old,facH,zzz,tau_ini
  real(kind=8),dimension(1:ncell)::w1H,w2H,wmax,time_max
  real(kind=8)::varmax=4d0
  integer::i,i_T2,iter,n,n_active
  integer,dimension(1:ncell)::ind,iii,i_nH
  logical::tau_negative
  
  ! Initializations
  logT2max=log10(T2_max_fix)
  dlog_nH=dble(table%n1-1)/(table%nH(table%n1)-table%nH(1))
  dlog_T2=dble(table%n2-1)/(table%T2(table%n2)-table%T2(1))
  h=1d0/dlog_T2
  h2=h*h
  h3=h2*h
  precoeff=2d0*X/(3d0*kB)
  do i=1,ncell
     zzz(i)=zsolar(i)
     facH(i)=MIN(MAX(log10(nH(i)/boost(i)),table%nH(1)),table%nH(table%n1))
     i_nH(i)=MIN(MAX(int((facH(i)-table%nH(1))*dlog_nH)+1,1),table%n1-1)
     w1H(i)=(table%nH(i_nH(i)+1)-facH(i))*dlog_nH
     w2H(i)=(facH(i)-table%nH(i_nH(i)  ))*dlog_nH
     tau(i)=T2(i)
     tau_ini(i)=T2(i)
     time_max(i)=dt*precoeff*nH(i)
     time(i)=0d0
     wmax(i)=1d0/time_max(i)
     ind(i)=i
  end do
  
  ! Check positivity 
  tau_negative=.false.
  do i=1,ncell
     if(tau(i)<=0.)tau_negative=.true.
  end do  
  if (tau_negative) then
     write(*,*)'ERROR in solve_cooling :'
     write(*,*)'Initial temperature is negative'
     STOP
  endif
  
  ! Loop over active cells
  iter=0
  n=ncell
  do while(n>0)
     
     iter=iter+1
     if (iter > 500) then
        write(*,*) 'Too many iterations in solve_cooling',iter,n
        do i=1,n
           write(*,*)i,tau(ind(i)),T2(ind(i)),nH(ind(i)),i_nH(ind(i))
        end do
        STOP
     endif
     
     n_active=0
     do i=1,n
        facT=log10(tau(ind(i)))

        if(facT.le.logT2max)then

           i_T2=MIN(MAX(int((facT-table%T2(1))*dlog_T2)+1,1),table%n2-1)
           yy=facT-table%T2(i_T2)
           yy2=yy*yy
           yy3=yy2*yy

           ! Cooling
           fa=table%cool(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%cool(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=table%cool(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%cool(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=table%cool_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%cool_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=table%cool_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%cool_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           cool=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           cool_prime=cool/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Heating
           fa=table%heat(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%heat(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=table%heat(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%heat(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=table%heat_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%heat_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=table%heat_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%heat_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           heat=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           heat_prime=heat/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Compton cooling
           fa=table%cool_com(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%cool_com(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=table%cool_com(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%cool_com(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=table%cool_com_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%cool_com_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=table%cool_com_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%cool_com_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           cool_com=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           cool_com_prime=cool_com/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Compton heating
           fa=table%heat_com(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%heat_com(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=table%heat_com(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%heat_com(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=table%heat_com_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%heat_com_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=table%heat_com_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%heat_com_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           heat_com=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           heat_com_prime=heat_com/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Metal cooling
           fa=table%metal(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%metal(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=table%metal(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%metal(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=table%metal_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+table%metal_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=table%metal_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+table%metal_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           metal=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           metal_prime=metal/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Total net cooling
           lambda=cool+zzz(ind(i))*metal-heat+(cool_com-heat_com)/nH(ind(i))
           lambda_prime=cool_prime+zzz(ind(i))*metal_prime-heat_prime+(cool_com_prime-heat_com_prime)/nH(ind(i))

        else

           lambda=1.42*1d-27*sqrt(tau(ind(i)))*1.1
           lambda_prime=lambda/2./tau(ind(i))

        endif

        wcool=MAX(abs(lambda)/tau(ind(i))*varmax,wmax(ind(i)),-lambda_prime*varmax)

        tau_old(ind(i))=tau(ind(i))
        tau(ind(i))=tau(ind(i))*(1d0+lambda_prime/wcool-lambda/tau(ind(i))/wcool)/(1d0+lambda_prime/wcool)
        time_old(ind(i))=time(ind(i))
        time(ind(i))=time(ind(i))+1d0/wcool

!!$        if(i==1)then
!!$           write(10,'(I5,10(1PE10.3,1X))')iter,tau_old(ind(i)),cool+zzz(ind(i))*metal,heat,lambda
!!$        endif

        if(time(ind(i))<time_max(ind(i)))then
           n_active=n_active+1
           ind(n_active)=ind(i)
        end if
        
     end do
     n=n_active
  end do
  ! End loop over active cells

  ! Compute exact time solution
  do i=1,ncell
     tau(i)=tau(i)*(time_max(i)-time_old(i))/(time(i)-time_old(i))+tau_old(i)*(time(i)-time_max(i))/(time(i)-time_old(i))
  end do

  ! Check positivity 
  tau_negative=.false.
  do i=1,ncell
     if (tau(i)<=0.)tau_negative=.true.
  end do  
  if (tau_negative) then
     write(*,*)'ERROR in solve_cooling :'
     write(*,*)'Final temperature is negative'
     STOP
  endif

  ! Compute delta T
  do i=1,ncell
     deltaT2(i)=tau(i)-tau_ini(i)
  end do
  
end subroutine solve_cooling
!=======================================================================
function J0simple(aexp)
!=======================================================================
! Le J0 dans le cas teyssier ou theuns
!=======================================================================
  real(kind=8) :: J0simple,aexp
  if (aexp .lt. 1.d0/(1d0+zreioniz)) then
     J0simple=0.d0
  elseif (aexp .lt. 1.d0/4.d0)then
     J0simple=4.d0*aexp
  elseif (aexp .lt. 1.d0/3.d0)then
     J0simple=1.d0
  else
     J0simple=1.d0/(3.*aexp)**3
  endif
  if(force_j0_one)J0simple=1.0
  J0simple=max(J0simple*J0in,J0min)
  return
end function J0simple
!=======================================================================
subroutine cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
!=======================================================================
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(kind=8)::nH_min,nH_max,T2_min,T2_max,aexp,tmp
  integer::nbin_n,nbin_T
  integer::myid,ncpu,ierr
  integer::i_n,i_T
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  integer :: i,j,n1,n2
  logical,save:: first=.true.
  
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif

  if(.not.first)then
     deallocate(table%cool)
     deallocate(table%heat)
     deallocate(table%cool_com)
     deallocate(table%heat_com)
     deallocate(table%metal)
     deallocate(table%cool_prime)
     deallocate(table%heat_prime)
     deallocate(table%cool_com_prime)
     deallocate(table%heat_com_prime)
     deallocate(table%metal_prime)
     deallocate(table%mu)
     deallocate(table%T2)
     deallocate(table%nH)
     if (if_species_abundances) deallocate(table%n_spec)
  else
     first=.false.
  endif
  
  table%n1=nbin_n
  table%n2=nbin_T
  allocate(table%cool(nbin_n,nbin_T))
  allocate(table%heat(nbin_n,nbin_T))
  allocate(table%cool_com(nbin_n,nbin_T))
  allocate(table%heat_com(nbin_n,nbin_T))
  allocate(table%metal(nbin_n,nbin_T))
  allocate(table%cool_prime(nbin_n,nbin_T))
  allocate(table%heat_prime(nbin_n,nbin_T))
  allocate(table%cool_com_prime(nbin_n,nbin_T))
  allocate(table%heat_com_prime(nbin_n,nbin_T))
  allocate(table%metal_prime(nbin_n,nbin_T))
  allocate(table%mu  (nbin_n,nbin_T))
  allocate(table%nH  (1:nbin_n))
  allocate(table%T2  (1:nbin_T))
  if (if_species_abundances) allocate(table%n_spec(nbin_n,nbin_T,1:6))
#ifndef WITHOUTMPI
  allocate(table2%cool(nbin_n,nbin_T))
  allocate(table2%heat(nbin_n,nbin_T))
  allocate(table2%cool_com(nbin_n,nbin_T))
  allocate(table2%heat_com(nbin_n,nbin_T))
  allocate(table2%metal(nbin_n,nbin_T))
  allocate(table2%cool_prime(nbin_n,nbin_T))
  allocate(table2%heat_prime(nbin_n,nbin_T))
  allocate(table2%cool_com_prime(nbin_n,nbin_T))
  allocate(table2%heat_com_prime(nbin_n,nbin_T))
  allocate(table2%metal_prime(nbin_n,nbin_T))
  allocate(table2%mu  (nbin_n,nbin_T))
  if (if_species_abundances) allocate(table2%n_spec(nbin_n,nbin_T,1:6))
#endif  
  do i_n=1,nbin_n
     tmp=log10(nH_min)+(dble(i_n)-1d0)/(dble(nbin_n)-1d0)*(log10(nH_max)-log10(nH_min))
     table%nH(i_n)=tmp
  end do
  do i_T=1,nbin_T
     tmp=log10(T2_min)+(dble(i_T)-1d0)/(dble(nbin_T)-1d0)*(log10(T2_max)-log10(T2_min))
     table%T2(i_T)=tmp
  end do

  ! Compute radiative ionization and heating rates
  call set_rates(t_rad_spec,h_rad_spec,aexp)

  ! Create the table
  table%mu=0.0
  table%cool=0.0
  table%heat=0.0
  table%cool_com=0.0
  table%heat_com=0.0
  table%metal=0.0
  table%cool_prime=0.0
  table%heat_prime=0.0
  table%cool_com_prime=0.0
  table%heat_com_prime=0.0
  table%metal_prime=0.0
  if (if_species_abundances) table%n_spec=0.0
  do i_n = myid+1,nbin_n,ncpu
     call iterate(i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(table%mu,table2%mu,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%cool,table2%cool,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%heat,table2%heat,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%cool_com,table2%cool_com,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%heat_com,table2%heat_com,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%metal,table2%metal,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%cool_prime,table2%cool_prime,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%heat_prime,table2%heat_prime,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%cool_com_prime,table2%cool_com_prime,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%heat_com_prime,table2%heat_com_prime,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%metal_prime,table2%metal_prime,nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (if_species_abundances)then
     call MPI_ALLREDUCE(table%n_spec,table2%n_spec,6*nbin_n*nbin_T,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if
  table%cool = table2%cool
  table%heat = table2%heat
  table%cool_com = table2%cool_com
  table%heat_com = table2%heat_com
  table%metal = table2%metal
  table%cool_prime = table2%cool_prime
  table%heat_prime = table2%heat_prime
  table%cool_com_prime = table2%cool_com_prime
  table%heat_com_prime = table2%heat_com_prime
  table%metal_prime = table2%metal_prime
  table%mu = table2%mu
  if (if_species_abundances)then
     table%n_spec=table2%n_spec
  endif
#endif

#ifndef WITHOUTMPI
  deallocate(table2%cool)
  deallocate(table2%heat)
  deallocate(table2%cool_com)
  deallocate(table2%heat_com)
  deallocate(table2%metal)
  deallocate(table2%cool_prime)
  deallocate(table2%heat_prime)
  deallocate(table2%cool_com_prime)
  deallocate(table2%heat_com_prime)
  deallocate(table2%metal_prime)
  deallocate(table2%mu)
  if (if_species_abundances)then
     deallocate(table2%n_spec)
  endif
#endif

end subroutine cmp_table
!=======================================================================
subroutine set_rates(t_rad_spec,h_rad_spec,aexp)
!=======================================================================
  implicit none
  real(kind=8),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(kind=8) :: J0,z,aexp
  logical :: first=.true.
  save first

  z=1.d0/aexp-1.D0
  if (high_z_realistic_ne) J0min=J0min_ref/(aexp/aexp_ref)**2
  if (teyssier) then
     J0=J0simple(aexp)
     t_rad_spec(HI  ) = taux_rad(HI  ,J0)
     t_rad_spec(HEI ) = taux_rad(HEI ,J0)
     t_rad_spec(HEII) = taux_rad(HEII,J0)
     h_rad_spec(HI  ) = heat_rad(HI  ,J0)
     h_rad_spec(HEI ) = heat_rad(HEI ,J0)
     h_rad_spec(HEII) = heat_rad(HEII,J0)
  elseif (theuns) then
     J0=J0simple(aexp)
     t_rad_spec(HI  ) = taux_rad_theuns(HI  ,J0)
     t_rad_spec(HEI ) = taux_rad_theuns(HEI ,J0)
     t_rad_spec(HEII) = taux_rad_theuns(HEII,J0)
     h_rad_spec(HI  ) = heat_rad_theuns(HI  ,J0)
     h_rad_spec(HEI ) = heat_rad_theuns(HEI ,J0)
     h_rad_spec(HEII) = heat_rad_theuns(HEII,J0)
  elseif (madau) then
     z=1.d0/aexp-1.D0
     t_rad_spec(HI  ) = taux_rad_madau(HI  ,z)
     t_rad_spec(HEI ) = taux_rad_madau(HEI ,z)
     t_rad_spec(HEII) = taux_rad_madau(HEII,z)
     h_rad_spec(HI  ) = heat_rad_madau(HI  ,z)
     h_rad_spec(HEI ) = heat_rad_madau(HEI ,z)
     h_rad_spec(HEII) = heat_rad_madau(HEII,z)     
  elseif (weinbergint) then
     t_rad_spec(HI  ) = taux_rad_weinbergint(HI  ,z)
     t_rad_spec(HEI ) = taux_rad_weinbergint(HEI ,z)
     t_rad_spec(HEII) = taux_rad_weinbergint(HEII,z)
     h_rad_spec(HI  ) = heat_rad_weinbergint(HI  ,z)
     h_rad_spec(HEI ) = heat_rad_weinbergint(HEI ,z)
     h_rad_spec(HEII) = heat_rad_weinbergint(HEII,z)     
  elseif (courty) then
     t_rad_spec(HI  ) = taux_rad_courty(HI  ,z)
     t_rad_spec(HEI ) = taux_rad_courty(HEI ,z)
     t_rad_spec(HEII) = taux_rad_courty(HEII,z)
     h_rad_spec(HI  ) = heat_rad_courty(HI  ,z)
     h_rad_spec(HEI ) = heat_rad_courty(HEI ,z)
     h_rad_spec(HEII) = heat_rad_courty(HEII,z)
  endif  
end subroutine set_rates
!=======================================================================
subroutine iterate(i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
!=======================================================================
  implicit none
  integer :: i_n
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) :: aexp
  integer::nbin_T    
  integer::i_T
  real(kind=8) ::T2,T2_eps,nH
  real(kind=8) ::mu,mu_eps
  real(kind=8) ::T2_left,T2_right,err_T2
  real(kind=8) ::cool_tot,heat_tot,cool_com,heat_com,metal_tot,metal_prime
  real(kind=8) ::cool_tot_eps,heat_tot_eps,cool_com_eps,heat_com_eps,metal_tot_eps
  real(kind=8) ::diff
  integer::niter
  real(kind=8),dimension(1:6) :: n_spec,n_spec_eps
  
  nH=10d0**table%nH(i_n)         
  do i_T = 1,nbin_T
     T2=10d0**table%T2(i_T)
     ! Compute cooling, heating and mean molecular weight
     call cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu,aexp,n_spec)
     table%cool(i_n,i_T)=log10(cool_tot)
     table%heat(i_n,i_T)=log10(heat_tot)
     table%cool_com(i_n,i_T)=log10(cool_com)
     table%heat_com(i_n,i_T)=log10(heat_com)
     table%mu(i_n,i_T)=mu
     if (if_species_abundances)then
        table%n_spec(i_n,i_T,1:6)=log10(n_spec(1:6))
     endif
     ! Compute cooling and heating derivatives
     T2_eps=10d0**(table%T2(i_T)+0.01d0)
     call cmp_cooling(T2_eps,nH,t_rad_spec,h_rad_spec,cool_tot_eps,heat_tot_eps,cool_com_eps,heat_com_eps,mu_eps,aexp,n_spec_eps)
     table%cool_prime(i_n,i_T)=(log10(cool_tot_eps)-log10(cool_tot))/0.01
     table%heat_prime(i_n,i_T)=(log10(heat_tot_eps)-log10(heat_tot))/0.01
     table%cool_com_prime(i_n,i_T)=(log10(cool_com_eps)-log10(cool_com))/0.01
     table%heat_com_prime(i_n,i_T)=(log10(heat_com_eps)-log10(heat_com))/0.01
     ! Compute metal contribution for solar metallicity
     call cmp_metals(T2,nH,mu,metal_tot,metal_prime,aexp)
     table%metal(i_n,i_T)=log10(metal_tot)
     table%metal_prime(i_n,i_T)=metal_prime
  end do
end subroutine iterate
!=======================================================================
subroutine cmp_metals(T2,nH,mu,metal_tot,metal_prime,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T2,nH,mu,metal_tot,metal_prime,aexp
  ! Compute cooling enhancement due to metals                                                
  ! Sutherland and Dopita (93) at solar metalicity                                           
  real(kind=8),dimension(1:91) :: temperature_sd93 = (/ &
       & 4.00,4.05,4.10,4.15,4.20,4.25,4.30,4.35,4.40,4.45,4.50,4.55,4.60, &
       & 4.65,4.70,4.75,4.80,4.85,4.90,4.95,5.00,5.05,5.10,5.15,5.20,5.25, &
       & 5.30,5.35,5.40,5.45,5.50,5.55,5.60,5.65,5.70,5.75,5.80,5.85,5.90, &
       & 5.95,6.00,6.05,6.10,6.15,6.20,6.25,6.30,6.35,6.40,6.45,6.50,6.55, &
       & 6.60,6.65,6.70,6.75,6.80,6.85,6.90,6.95,7.00,7.05,7.10,7.15,7.20, &
       & 7.25,7.30,7.35,7.40,7.45,7.50,7.55,7.60,7.65,7.70,7.75,7.80,7.85, &
       & 7.90,7.95,8.00,8.05,8.10,8.15,8.20,8.25,8.30,8.35,8.40,8.45,8.50  /)
  real(kind=8),dimension(1:91) :: excess_cooling_sd93 = (/ &
       & -25.8772,-24.4777,-23.6389,-22.9812,-22.5772,-22.3998,-22.3194, &
       & -22.2163,-22.0605,-21.9099,-21.7450,-21.6143,-21.4835,-21.3623, &
       & -21.2572,-21.1564,-21.0694,-20.9940,-20.9351,-20.8923,-20.8885, &
       & -20.9153,-20.9224,-20.8994,-20.8669,-20.8556,-20.8446,-20.8439, &
       & -20.8736,-21.0144,-21.2366,-21.4396,-21.5513,-21.5916,-21.6013, &
       & -21.6008,-21.6516,-21.7543,-21.8264,-21.8468,-21.8572,-21.8572, &
       & -21.8468,-21.8364,-21.8364,-21.8681,-21.9734,-22.1119,-22.2315, &
       & -22.3230,-22.3814,-22.4178,-22.4549,-22.4950,-22.5342,-22.5645, &
       & -22.5960,-22.5991,-22.5791,-22.5723,-22.5756,-22.5962,-22.6461, &
       & -22.7149,-22.7740,-22.8215,-22.8739,-22.9121,-22.9331,-22.9689, &
       & -22.9721,-23.0007,-23.0063,-22.9863,-22.9929,-22.9729,-22.9994, &
       & -22.9794,-22.9594,-22.9696,-22.9712,-22.9512,-22.9312,-22.9112, &
       & -22.9145,-22.8945,-22.8745,-22.8798,-22.8598,-22.8398,-22.8472  /)
  real(kind=8),dimension(1:91) :: excess_prime_sd93 = (/ &
       & 33.5968475,22.3829498,14.9650421,10.6169891, 5.8140259, 2.5779724, &
       & 1.8350220, 2.5890045, 3.0639954, 3.1549835, 2.9560089, 2.6150055, &
       & 2.5199890, 2.2629852, 2.0589905, 1.8779907, 1.6240082, 1.3430023, &
       & 1.0169983, 0.4660034,-0.2300110,-0.3390045, 0.1589813, 0.5549927, &
       & 0.4380035, 0.2229919, 0.1170044,-0.2899933,-1.7050018,-3.6300049, &
       & -4.2519836,-3.1469879,-1.5200043,-0.4999847,-0.0919800,-0.5030060, &
       & -1.5350037,-1.7480164,-0.9250031,-0.3079987,-0.1040039, 0.1040039, &
       & 0.2080078, 0.1040039,-0.3169861,-1.3700104,-2.4380188,-2.5809937, &
       & -2.1109924,-1.4989929,-0.9480133,-0.7350159,-0.7720032,-0.7930145, &
       & -0.6950073,-0.6180115,-0.3460083, 0.1690063, 0.2679901, 0.0350037, &
       & -0.2390137,-0.7050018,-1.1869659,-1.2790070,-1.0660248,-0.9989929, &
       & -0.9059906,-0.5919952,-0.5680084,-0.3899994,-0.3179932,-0.3419952, &
       & 0.1439972, 0.1339722, 0.1339874,-0.0649872,-0.0650024, 0.3999939, &
       & 0.0980072,-0.1180115, 0.1840057, 0.4000092, 0.4000092, 0.1670074, &
       & 0.1669769, 0.3999939, 0.1470032, 0.1470032, 0.4000244, 0.1260071, &
       & 0.0000000 /)
  ! Compute cooling enhancement due to metals
  ! Cloudy at solar metalicity
  real(kind=8),dimension(1:91) :: temperature_cc07 = (/ &
       & 3.9684,4.0187,4.0690,4.1194,4.1697,4.2200,4.2703, &
       & 4.3206,4.3709,4.4212,4.4716,4.5219,4.5722,4.6225, &
       & 4.6728,4.7231,4.7734,4.8238,4.8741,4.9244,4.9747, &
       & 5.0250,5.0753,5.1256,5.1760,5.2263,5.2766,5.3269, &
       & 5.3772,5.4275,5.4778,5.5282,5.5785,5.6288,5.6791, &
       & 5.7294,5.7797,5.8300,5.8804,5.9307,5.9810,6.0313, &
       & 6.0816,6.1319,6.1822,6.2326,6.2829,6.3332,6.3835, &
       & 6.4338,6.4841,6.5345,6.5848,6.6351,6.6854,6.7357, &
       & 6.7860,6.8363,6.8867,6.9370,6.9873,7.0376,7.0879, &
       & 7.1382,7.1885,7.2388,7.2892,7.3395,7.3898,7.4401, &
       & 7.4904,7.5407,7.5911,7.6414,7.6917,7.7420,7.7923, &
       & 7.8426,7.8929,7.9433,7.9936,8.0439,8.0942,8.1445, &
       & 8.1948,8.2451,8.2955,8.3458,8.3961,8.4464,8.4967 /)
  real(kind=8),dimension(1:91) :: excess_cooling_cc07 = (/ &
       & -24.9949,-24.7270,-24.0473,-23.0713,-22.2907,-21.8917,-21.8058, &
       & -21.8501,-21.9142,-21.9553,-21.9644,-21.9491,-21.9134,-21.8559, &
       & -21.7797,-21.6863,-21.5791,-21.4648,-21.3640,-21.2995,-21.2691, &
       & -21.2658,-21.2838,-21.2985,-21.2941,-21.2845,-21.2809,-21.2748, &
       & -21.2727,-21.3198,-21.4505,-21.5921,-21.6724,-21.6963,-21.6925, &
       & -21.6892,-21.7142,-21.7595,-21.7779,-21.7674,-21.7541,-21.7532, &
       & -21.7679,-21.7866,-21.8052,-21.8291,-21.8716,-21.9316,-22.0055, &
       & -22.0800,-22.1600,-22.2375,-22.3126,-22.3701,-22.4125,-22.4353, &
       & -22.4462,-22.4450,-22.4406,-22.4337,-22.4310,-22.4300,-22.4356, &
       & -22.4455,-22.4631,-22.4856,-22.5147,-22.5444,-22.5718,-22.5904, &
       & -22.6004,-22.5979,-22.5885,-22.5728,-22.5554,-22.5350,-22.5159, &
       & -22.4955,-22.4781,-22.4600,-22.4452,-22.4262,-22.4089,-22.3900, &
       & -22.3722,-22.3529,-22.3339,-22.3137,-22.2936,-22.2729,-22.2521 /)
  real(kind=8),dimension(1:91) :: excess_prime_cc07 = (/ & 
       &   2.0037,  4.7267, 12.2283, 13.5820,  9.8755,  4.8379,  1.8046, &
       &   1.4574,  1.8086,  2.0685,  2.2012,  2.2250,  2.2060,  2.1605, &
       &   2.1121,  2.0335,  1.9254,  1.7861,  1.5357,  1.1784,  0.7628, &
       &   0.1500, -0.1401,  0.1272,  0.3884,  0.2761,  0.1707,  0.2279, &
       &  -0.2417, -1.7802, -3.0381, -2.3511, -0.9864, -0.0989,  0.1854, &
       &  -0.1282, -0.8028, -0.7363, -0.0093,  0.3132,  0.1894, -0.1526, &
       &  -0.3663, -0.3873, -0.3993, -0.6790, -1.0615, -1.4633, -1.5687, &
       &  -1.7183, -1.7313, -1.8324, -1.5909, -1.3199, -0.8634, -0.5542, &
       &  -0.1961, -0.0552,  0.0646, -0.0109, -0.0662, -0.2539, -0.3869, &
       &  -0.6379, -0.8404, -1.1662, -1.3930, -1.6136, -1.5706, -1.4266, &
       &  -1.0460, -0.7244, -0.3006, -0.1300,  0.1491,  0.0972,  0.2463, &
       &   0.0252,  0.1079, -0.1893, -0.1033, -0.3547, -0.2393, -0.4280, &
       &  -0.2735, -0.3670, -0.2033, -0.2261, -0.0821, -0.0754,  0.0634 /)
  real(kind=8),dimension(1:50)::z_courty=(/ &
       & 0.00000,0.04912,0.10060,0.15470,0.21140,0.27090,0.33330,0.39880, &
       & 0.46750,0.53960,0.61520,0.69450,0.77780,0.86510,0.95670,1.05300, &
       & 1.15400,1.25900,1.37000,1.48700,1.60900,1.73700,1.87100,2.01300, &
       & 2.16000,2.31600,2.47900,2.64900,2.82900,3.01700,3.21400,3.42100, &
       & 3.63800,3.86600,4.10500,4.35600,4.61900,4.89500,5.18400,5.48800, &
       & 5.80700,6.14100,6.49200,6.85900,7.24600,7.65000,8.07500,8.52100, &
       & 8.98900,9.50000 /)
  real(kind=8),dimension(1:50)::phi_courty=(/ &
       & 0.0499886,0.0582622,0.0678333,0.0788739,0.0915889,0.1061913,0.1229119, &
       & 0.1419961,0.1637082,0.1883230,0.2161014,0.2473183,0.2822266,0.3210551, &
       & 0.3639784,0.4111301,0.4623273,0.5172858,0.5752659,0.6351540,0.6950232, &
       & 0.7529284,0.8063160,0.8520859,0.8920522,0.9305764,0.9682031,1.0058810, &
       & 1.0444020,1.0848160,1.1282190,1.1745120,1.2226670,1.2723200,1.3231350, &
       & 1.3743020,1.4247480,1.4730590,1.5174060,1.5552610,1.5833640,1.5976390, &
       & 1.5925270,1.5613110,1.4949610,1.3813710,1.2041510,0.9403100,0.5555344, & 
       & 0.0000000 /)
  real(kind=8)::TT,lTT,deltaT,lcool,lcool1,lcool2,lcool1_prime,lcool2_prime
  real(kind=8)::ZZ,deltaZ
  real(kind=8)::c1=0.4,c2=10.0,TT0=1d5,TTC=1d6,alpha1=0.15
  real(kind=8)::ux,g_courty,f_courty=1d0,g_courty_prime,f_courty_prime
  integer::iT,iZ

  ZZ=1d0/aexp-1d0
  TT=T2*mu
  lTT=log10(TT)

  ! This is a simple model to take into account the ionization background
  ! on metal cooling (calibrated using CLOUDY). 
  if(madau.or.weinbergint.or.courty)then
     if(ZZ.le.0.0.or.ZZ.ge.z_courty(50))then
        ux=0.0
     else
        iZ=1+int(ZZ/z_courty(50)*49.)
        iZ=min(iZ,49)
        iZ=max(iZ,1)
        deltaZ=z_courty(iZ+1)-z_courty(iZ)
        ux=1d-4*(phi_courty(iZ+1)*(ZZ-z_courty(iZ))/deltaZ & 
             & + phi_courty(iZ)*(z_courty(iZ+1)-ZZ)/deltaZ )/nH
     endif
  else ! Theuns or Teyssier
     ux=1d-4*J0simple(aexp)/1d-22/nH
  endif
  g_courty=c1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)
  g_courty_prime=(c1*alpha1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)*TTC/TT)/TT
  f_courty=1d0/(1d0+ux/g_courty)
  f_courty_prime=ux/g_courty/(1d0+ux/g_courty)**2*g_courty_prime/g_courty

!  if(lTT.ge.temperature_sd93(91))then
  if(lTT.ge.temperature_cc07(91))then
     metal_tot=1d-100
     metal_prime=0d0
  else if(lTT.ge.1.0)then
     lcool1=-100d0
     lcool1_prime=0d0
!      if(lTT.ge.temperature_sd93(1))then
      if(lTT.ge.temperature_cc07(1))then
!        iT=1+int((lTT-temperature_sd93(1))/(temperature_sd93(91)-temperature_sd93(1))*90.0)
        iT=1+int((lTT-temperature_cc07(1))/(temperature_cc07(91)-temperature_cc07(1))*90.0)
        iT=min(iT,90)
        iT=max(iT,1)
!        deltaT=temperature_sd93(iT+1)-temperature_sd93(iT)
!        lcool1 = excess_cooling_sd93(iT+1)*(lTT-temperature_sd93(iT))/deltaT  &
!             & + excess_cooling_sd93(iT)*(temperature_sd93(iT+1)-lTT)/deltaT
!        lcool1_prime  = excess_prime_sd93(iT+1)*(lTT-temperature_sd93(iT))/deltaT  &
!                    & + excess_prime_sd93(iT)*(temperature_sd93(iT+1)-lTT)/deltaT
        deltaT=temperature_cc07(iT+1)-temperature_cc07(iT)
        lcool1 = excess_cooling_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT  &
             & + excess_cooling_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT 
        lcool1_prime  = excess_prime_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT  &
                    & + excess_prime_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT 
     endif
     ! Fine structure cooling from infrared lines
     lcool2=-31.522879+2.0*lTT-20.0/TT-TT*4.342944d-5
     lcool2_prime=2d0+(20d0/TT-TT*4.342944d-5)*log(10d0)
     ! Total metal cooling and temperature derivative
     metal_tot=10d0**lcool1+10d0**lcool2
     metal_prime=(10d0**lcool1*lcool1_prime+10d0**lcool2*lcool2_prime)/metal_tot
     metal_prime=metal_prime*f_courty+metal_tot*f_courty_prime
     metal_tot=metal_tot*f_courty
  else
     metal_tot=1d-100
     metal_prime=0d0
  endif

end subroutine cmp_metals
!=======================================================================
subroutine cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu_out,aexp,n_spec)
!=======================================================================
  implicit none
  
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) ::T2,nH,cool_tot,heat_tot,cool_com,heat_com,mu_out,aexp
  real(kind=8) ::mu,mu_old,err_mu,mu_left,mu_right
  real(kind=8) ::T
  real(kind=8) ::n_E,n_HI,n_HII,n_HEI,n_HEII,n_HEIII,n_TOT
  real(kind=8),dimension(1:6)::n_spec
  real(kind=8) ::cb1,cb2,cb3,ci1,ci2,ci3,cr1,cr2,cr3,cd,ce1,ce2,ce3,ch1,ch2,ch3,coc,coh
  integer :: niter
  
  ! Iteration to find mu
  err_mu=1.
  mu_left=0.5
  mu_right=1.3
  niter=0
  do while (err_mu > 1.d-4 .and. niter <= 50)
     mu_old=0.5*(mu_left+mu_right)
     T = T2*mu_old
     call cmp_chem_eq(T,nH,t_rad_spec,n_spec,n_TOT,mu)
     err_mu = (mu-mu_old)/mu_old
     if(err_mu>0.)then 
        mu_left =0.5*(mu_left+mu_right)
        mu_right=mu_right
     else
        mu_left =mu_left
        mu_right=0.5*(mu_left+mu_right)
     end if
     err_mu=ABS(err_mu)
     niter=niter+1
  end do
  if (niter > 50) then
     write(*,*) 'ERROR in cmp_cooling : too many iterations.'
     STOP
  endif
    
  ! Get equilibrium abundances
  n_E     = n_spec(1) ! electrons
  n_HI    = n_spec(2) ! H
  n_HII   = n_spec(3) ! H+
  n_HEI   = n_spec(4) ! He
  n_HEII  = n_spec(5) ! He+
  n_HEIII = n_spec(6) ! He++
  ! Bremstrahlung
  cb1 = cool_bre(HI  ,T)*n_E*n_HII  /nH**2
  cb2 = cool_bre(HEI ,T)*n_E*n_HEII /nH**2
  cb3 = cool_bre(HEII,T)*n_E*n_HEIII/nH**2
  ! Ionization cooling
  ci1 = cool_ion(HI  ,T)*n_E*n_HI   /nH**2
  ci2 = cool_ion(HEI ,T)*n_E*n_HEI  /nH**2
  ci3 = cool_ion(HEII,T)*n_E*n_HEII /nH**2
  ! Recombination cooling
  cr1 = cool_rec(HI  ,T)*n_E*n_HII  /nH**2
  cr2 = cool_rec(HEI ,T)*n_E*n_HEII /nH**2
  cr3 = cool_rec(HEII,T)*n_E*n_HEIII/nH**2
  ! Dielectric recombination cooling
  cd  = cool_die(T     )*n_E*n_HEII /nH**2
  ! Line cooling
  ce1 = cool_exc(HI  ,T)*n_E*n_HI   /nH**2
  ce2 = cool_exc(HEI, T)*n_E*n_HEI  /nH**2
  ce3 = cool_exc(HEII,T)*n_E*n_HEII /nH**2
  ! Radiative heating
  ch1 = h_rad_spec(HI  )    *n_HI   /nH**2
  ch2 = h_rad_spec(HEI )    *n_HEI  /nH**2
  ch3 = h_rad_spec(HEII)    *n_HEII /nH**2
  ! Total cooling and heating rates
  heat_tot = ch1+ch2+ch3
  cool_tot = cb1+cb2+cb3+ci1+ci2+ci3+cr1+cr2+cr3+cd+ce1+ce2+ce3
  ! Compton cooling
  coc = cool_compton(T,aexp)*n_E/nH
  cool_com = coc
  ! Compton heating
  coh = heat_compton(T,aexp)*n_E/nH
  heat_com = coh
  ! Mean molecular weight
  mu_out = mu
  
  if (if_cooling_functions) then
     cool_out=max(cool_tot,smallnum_cooling)
     heat_out=max(heat_tot,smallnum_cooling)
     cool_com=max(cool_com,smallnum_cooling)
     heat_com=max(heat_com,smallnum_cooling)
     cb1s=max(cb1,smallnum_cooling)
     cb2s=max(cb2,smallnum_cooling)
     cb3s=max(cb3,smallnum_cooling)
     ci1s=max(ci1,smallnum_cooling)
     ci2s=max(ci2,smallnum_cooling)
     ci3s=max(ci3,smallnum_cooling)
     cr1s=max(cr1,smallnum_cooling)
     cr2s=max(cr2,smallnum_cooling)
     cr3s=max(cr3,smallnum_cooling)
     cds =max(cd ,smallnum_cooling)
     ce1s=max(ce1,smallnum_cooling)
     ce3s=max(ce3,smallnum_cooling)
     cocs=max(coc,smallnum_cooling)
     cohs=max(coh,smallnum_cooling)
     ch1s=max(ch1,smallnum_cooling)
     ch2s=max(ch2,smallnum_cooling)
     ch3s=max(ch3,smallnum_cooling)
     cohs=max(coh,smallnum_cooling)  
  endif
end subroutine cmp_cooling
!=======================================================================
subroutine cmp_chem_eq(T,n_H,t_rad_spec,n_spec,n_TOT,mu)
!=======================================================================
  implicit none
  real(kind=8)::T,n_H,n_TOT,mu
  real(kind=8),dimension(1:3)::t_rad_spec
  real(kind=8),dimension(1:6)::n_spec
  real(kind=8)::xx,yy
  real(kind=8)::n_HI,n_HII,n_HEI,n_HEII,n_HEIII,n_E
  real(kind=8)::t_rad_HI,t_rad_HEI,t_rad_HEII
  real(kind=8)::t_rec_HI,t_rec_HEI,t_rec_HEII
  real(kind=8)::t_ion_HI,t_ion_HEI,t_ion_HEII
  real(kind=8)::t_ion2_HI,t_ion2_HEI,t_ion2_HEII
  real(kind=8)::x1,err_nE
  
  xx=(1.-Y)
  yy=Y/(1.-Y)/4.
  
  t_rad_HI   = t_rad_spec(HI)
  t_rec_HI   = taux_rec  (HI,T)
  t_ion_HI   = taux_ion  (HI,T)
  
  t_rad_HEI  = t_rad_spec(HEI)
  t_rec_HEI  = taux_rec  (HEI,T)
  t_ion_HEI  = taux_ion  (HEI,T)
  
  t_rad_HEII = t_rad_spec(HEII)
  t_rec_HEII = taux_rec  (HEII,T)
  t_ion_HEII = taux_ion  (HEII,T)
  
  n_E = n_H        
  err_nE = 1.
  
  do while(err_nE > 1.d-8)
     
     t_ion2_HI   = t_ion_HI   + t_rad_HI  /MAX(n_E,1e-15*n_H)
     t_ion2_HEI  = t_ion_HEI  + t_rad_HEI /MAX(n_E,1e-15*n_H)
     t_ion2_HEII = t_ion_HEII + t_rad_HEII/MAX(n_E,1e-15*n_H)
     
     n_HI  = t_rec_HI/(t_ion2_HI+t_rec_HI)*n_H
     n_HII = t_ion2_HI/(t_ion2_HI+t_rec_HI)*n_H
     
     x1 = (t_rec_HEII*t_rec_HEI+t_ion2_HEI*t_rec_HEII+t_ion2_HEII*t_ion2_HEI)
     
     n_HEIII = yy*t_ion2_HEII*t_ion2_HEI/x1*n_H
     n_HEII  = yy*t_ion2_HEI *t_rec_HEII/x1*n_H
     n_HEI   = yy*t_rec_HEII *t_rec_HEI /x1*n_H
     
     err_nE = ABS((n_E - (n_HII + n_HEII + 2.*n_HEIII))/n_H)
     n_E = 0.5*n_E+0.5*(n_HII + n_HEII + 2.*n_HEIII)
     
  end do
    
  n_TOT    =n_E+n_HI+n_HII+n_HEI+n_HEII+n_HEIII
  mu       =n_H/xx/n_TOT
  n_spec(1)=n_E
  n_spec(2)=n_HI
  n_spec(3)=n_HII
  n_spec(4)=n_HEI
  n_spec(5)=n_HEII
  n_spec(6)=n_HEIII
  
end subroutine cmp_chem_eq
!=======================================================================
function cool_bre(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_bre
  if(ispec==HI  )cool_bre = 1.42D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3.D0))
  if(ispec==HEI )cool_bre = 1.42D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3.D0))
  if(ispec==HEII)cool_bre = 5.68D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3.D0))
  return
end function cool_bre
!=======================================================================
function cool_exc(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_exc,T5
  T5=1.d-5*T
  if(ispec==HI  )cool_exc = 7.50D-19/(1.+sqrt(T5))              *exp(-118348.D0/T)
  if(ispec==HEI )cool_exc = 9.10D-27/(1.+sqrt(T5))/(T**0.1687D0)*exp(-13179.D0/T)
  if(ispec==HEII)cool_exc = 5.54D-17/(1.+sqrt(T5))/(T**0.397D0 )*exp(-473638.D0/T)
  return
end function cool_exc
!=======================================================================
function cool_rec(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_rec
  real(kind=8)   ::T3, T6
  T3 = 1.d-03*T
  T6 = 1.d-06*T
  if(ispec==HI  )cool_rec = 8.70D-27*SQRT(T)/T3**(0.2D0)/(1.D0+T6**0.7D0)
  if(ispec==HEI )cool_rec = 1.55D-26*T**0.3647D0
  if(ispec==HEII)cool_rec = 3.48D-26*SQRT(T)/T3**(0.2D0)/(1.D0+T6**0.7D0)
  return
end function cool_rec
!=======================================================================
function cool_die(T)
!=======================================================================
  implicit none
  real(kind=8) :: T,cool_die
  cool_die=1.24D-13*T**(-1.5D0)*exp(-470000.D0/T)*(1.D0+0.3D0*exp(-94000.D0/T))
  return
end function cool_die
!=======================================================================
function taux_rec(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,taux_rec
  real(kind=8)   ::T3, T6
  T3 = 1.d-03*T
  T6 = 1.d-06*T
  if(ispec==HI  )taux_rec = dumfac_rec*8.40e-11/SQRT(T)/T3**(0.2)/(1.+T6**0.7)
  if(ispec==HEI )taux_rec = 1.50e-10/T**0.6353+taux_die(T)
  if(ispec==HEII)taux_rec = 3.36e-10/SQRT(T)/T3**(0.2)/(1.+T6**0.7)
  return
end function taux_rec
!=======================================================================
function taux_die(T)
!=======================================================================
  implicit none
  real(kind=8) :: T,taux_die
  taux_die=1.9D-3*T**(-1.5D0)*exp(-470000.D0/T)*(1.D0+0.3D0*exp(-94000.D0/T))
  return
end function taux_die
!=======================================================================
function cool_ion(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_ion
  real(kind=8)   ::T5
  T5 = 1.d-05*T
  if(ispec==HI  )cool_ion = dumfac_ion*1.27D-21*SQRT(T)/(1.+SQRT(T5))*EXP(-157809.1D0/T)
  if(ispec==HEI )cool_ion = dumfac_ion*9.38D-22*SQRT(T)/(1.+SQRT(T5))*EXP(-285335.4D0/T)
  if(ispec==HEII)cool_ion = dumfac_ion*4.95D-22*SQRT(T)/(1.+SQRT(T5))*EXP(-631515.0D0/T)
  return
end function cool_ion
!=======================================================================
function cool_compton(T,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T,aexp,cool_compton
  cool_compton=5.406D-36*T/aexp**4 
  return
end function cool_compton
!=======================================================================
function heat_compton(T,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T,aexp,heat_compton
  heat_compton=5.406D-36*2.726D0/aexp**5
  return
end function heat_compton
!=======================================================================
function taux_ion(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   :: T,taux_ion
  real(kind=8)   :: T5
  T5 = 1.d-05*T
  if(ispec==HI  )taux_ion = dumfac_ion*5.85D-11*SQRT(T)/(1.+SQRT(T5))*EXP(-157809.1D0/T)
  if(ispec==HEI )taux_ion = dumfac_ion*2.38D-11*SQRT(T)/(1.+SQRT(T5))*EXP(-285335.4D0/T)
  if(ispec==HEII)taux_ion = dumfac_ion*5.68D-12*SQRT(T)/(1.+SQRT(T5))*EXP(-631515.0D0/T)
  return
end function taux_ion
!=======================================================================
function J_nu(e,J0)
!=======================================================================
  implicit none
  real(kind=8) :: e,J_nu,e_L,J0,Jloc
  Jloc = max(J0,J0min) 
  e_L  = 13.598*eV
  J_nu = Jloc*(e_L/e)
  return
end function J_nu
!=======================================================================
function sigma_rad(e,ispec)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::sigma_rad,e,e_i,xxx,alph
  if(ispec==HI  )e_i = 13.598D0*eV
  if(ispec==HEI )e_i = 24.587D0*eV
  if(ispec==HEII)e_i = 54.416D0*eV
  xxx = e/e_i
  alph = sqrt(xxx-1.0d0)
  if(ispec==HI  )sigma_rad = 6.30D-18/xxx**4*exp(4.D0-4.D0*atan(alph)/alph) &
       &                             /(1.D0-exp(-twopi/alph))
  if(ispec==HEI )sigma_rad = 7.42D-18*(1.66D0/xxx**2.05D0-0.66D0/xxx**3.05D0)
  if(ispec==HEII)sigma_rad = 1.58D-18/xxx**4*exp(4.D0-4.D0*atan(alph)/alph) &
       &                             /(1.D0-exp(-twopi/alph))
  return
end function sigma_rad
!=======================================================================
function taux_rad(ispec,J0)
!=======================================================================
  implicit none  
  integer::ispec
  real(kind=8) :: J0,taux_rad,e_i,e,de,error,integ
  if(ispec==HI  )e_i = 13.598D0*eV
  if(ispec==HEI )e_i = 24.587D0*eV
  if(ispec==HEII)e_i = 54.416D0*eV
  integ = 0.0d0
  e = e_i
  de = e/100.D0
  error = 1.D0
  do while(error>1.d-6)
     e = e + de
     de = e/100.D0
     error = 2.0d0*twopi*J_nu(e,J0)*sigma_rad(e,ispec)*de/e
     integ = integ + error
     error = error/abs(integ)
  end do
  taux_rad = integ/hplanck
  return
end function taux_rad
!=======================================================================
function taux_rad_madau(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: z,taux_rad_madau,tt
  if (z < 15.d0) then
     if (ispec==HI  ) taux_rad_madau=normfacJ0*exp(-31.04D0+2.795D0*z-0.5589D0*z**2)
     if (ispec==HEI ) taux_rad_madau=normfacJ0*exp(-31.08D0+2.822D0*z-0.5664D0*z**2)
     if (ispec==HEII) taux_rad_madau=normfacJ0*exp(-34.30D0+1.826D0*z-0.3899D0*z**2)
  else
     taux_rad_madau=0.d0
  endif
  tt=taux_rad_theuns(ispec,J0min)
  if (taux_rad_madau < tt) taux_rad_madau=tt
  return
end function taux_rad_madau
!=======================================================================
function taux_rad_weinbergint(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,taux_rad_weinbergint,hh,tt
  if (z < 8.5d0) then
     if (ispec==HI  ) iweinb=1
     if (ispec==HEI ) iweinb=2
     if (ispec==HEII) iweinb=3
     hh=0.d0
     zz=max(z,1.0d-15)
     do i=1,Norderweinberg
        hh=hh+coefweinberg(i,iweinb)*zz**(i-1)
     enddo
     taux_rad_weinbergint=normfacJ0*exp(hh)
  else
     taux_rad_weinbergint=0.d0
  endif
  tt=taux_rad_theuns(ispec,J0min)
  if (taux_rad_weinbergint < tt) taux_rad_weinbergint=tt
  return
end function taux_rad_weinbergint
!=======================================================================
function taux_rad_theuns(ispec,J0)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: J0,taux_rad_theuns
  if (ispec==HI  ) taux_rad_theuns=1.26D10*J0/(3.D0+alpha)
  if (ispec==HEI ) taux_rad_theuns=1.48D10*J0*0.553D0**alpha &
                     & *(1.66D0/(alpha+2.05D0)-0.66D0/(alpha+3.05D0))
  if (ispec==HEII) taux_rad_theuns=3.34D9*J0*0.249D0**alpha/(3.D0+alpha)
  return
end function taux_rad_theuns
!=======================================================================
function taux_rad_courty(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,taux_rad_courty,hh,tt,hhreion
  if (z < zreioniz) then
     if (ispec==HI  ) iweinb=1
     if (ispec==HEI ) iweinb=2
     if (ispec==HEII) iweinb=3
     hh=0.d0
     zz=max(z,1.0d-15)
     do i=0,Nordercourty
        hh=hh+coefcourty(i,iweinb)*zz**i
     enddo
     hhreion=coef_fit(iweinb)*(zz/zreioniz)**beta_fit(iweinb)
     taux_rad_courty=10.**(hh-hhreion)
  else
     taux_rad_courty=0.d0
  endif
  tt=taux_rad_theuns(ispec,J0min)
  if (taux_rad_courty < tt) taux_rad_courty=tt
  return
end function taux_rad_courty
!=======================================================================
function heat_rad(ispec,J0)
!=======================================================================
  implicit none  
  integer::ispec
  real(kind=8) :: J0,heat_rad,e_i,e,de,error,integ
  if(ispec==HI  )e_i = 13.598D0*eV
  if(ispec==HEI )e_i = 24.587D0*eV
  if(ispec==HEII)e_i = 54.416D0*eV
  integ = 0.0d0
  e = e_i
  de = e/100.D0
  error = 1.D0
  do while(error>1.d-6)
     e = e + de
     de = e/100.D0
     error = 2.0d0*twopi*J_nu(e,J0)*sigma_rad(e,ispec)*(e/e_i-1.D0)*de/e
     integ = integ + error
     error=error/abs(integ)
  end do
  heat_rad = integ/hplanck*e_i
  return
end function heat_rad
!=======================================================================
function heat_rad_madau(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: z,heat_rad_madau,tt
  if (z < 15.d0) then
     if (ispec==HI  ) heat_rad_madau=normfacJ0*exp(-56.62D0+2.788D0*z-0.5594D0*z**2)
     if (ispec==HEI ) heat_rad_madau=normfacJ0*exp(-56.06D0+2.800D0*z-0.5532D0*z**2)
     if (ispec==HEII) heat_rad_madau=normfacJ0*exp(-58.67D0+1.888D0*z-0.3947D0*z**2)
  else
     heat_rad_madau=0.d0
  endif
  tt=heat_rad_theuns(ispec,J0min)
  if (heat_rad_madau < tt) heat_rad_madau=tt
  return
end function heat_rad_madau
!=======================================================================
function heat_rad_weinbergint(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,heat_rad_weinbergint,hh,tt
  if (z < 8.5d0) then
     if (ispec==HI  ) iweinb=4
     if (ispec==HEI ) iweinb=5
     if (ispec==HEII) iweinb=6
     hh=0.d0
     zz=max(z,1.0d-15)
     do i=1,Norderweinberg
        hh=hh+coefweinberg(i,iweinb)*zz**(i-1)
     enddo
     heat_rad_weinbergint=normfacJ0*exp(hh)
  else
     heat_rad_weinbergint=0.d0
  endif
  tt=heat_rad_theuns(ispec,J0min)
  if (heat_rad_weinbergint < tt) heat_rad_weinbergint=tt
  return
end function heat_rad_weinbergint
!=======================================================================
function heat_rad_theuns(ispec,J0)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: J0,heat_rad_theuns
  if (ispec==HI  ) heat_rad_theuns=(2.91D-1*J0/(2.D0+alpha))/(3.D0+alpha)
  if (ispec==HEI ) heat_rad_theuns=5.84D-1*J0*0.553D0**alpha* &
                 & (1.66D0/(alpha+1.05D0)-2.32D0/(alpha+2.05D0)+0.66D0/(alpha+3.05D0))
  if (ispec==HEII) heat_rad_theuns=(2.92D-1*J0*0.249D0**alpha/(2.D0+alpha))/(3.D0+alpha)
  return
end function heat_rad_theuns
!=======================================================================
function heat_rad_courty(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,heat_rad_courty,hh,tt,hhreion
  if (z < zreioniz) then
     if (ispec==HI  ) iweinb=4
     if (ispec==HEI ) iweinb=5
     if (ispec==HEII) iweinb=6
     hh=0.d0
     zz=max(z,1.0d-15)
     do i=0,Nordercourty
        hh=hh+coefcourty(i,iweinb)*zz**i
     enddo
     hhreion=coef_fit(iweinb)*(zz/zreioniz)**beta_fit(iweinb)
     heat_rad_courty=10.**(hh-hhreion)
  else
     heat_rad_courty=0.d0
  endif
  tt=heat_rad_theuns(ispec,J0min)
  if (heat_rad_courty < tt) heat_rad_courty=tt
  return
end function heat_rad_courty
!=======================================================================
function HsurH0(z,omega0,omegaL,OmegaR)
!=======================================================================
  implicit none
  real(kind=8) :: HsurH0,z,omega0,omegaL,omegaR
  HsurH0=sqrt(Omega0*(1.d0+z)**3+OmegaR*(1.d0+z)**2+OmegaL)
end function HsurH0

end module cooling_module

