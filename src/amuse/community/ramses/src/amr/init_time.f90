subroutine init_time
  use amr_commons
  use hydro_commons
  use pm_commons
  use cooling_module
#ifdef RT
  use rt_cooling_module
#endif
  implicit none
  integer::i,Nmodel
  real(kind=8)::T2_sim  

  if(nrestart==0)then
     if(cosmo)then
        ! Get cosmological parameters from input files
        call init_cosmo
     else
        ! Get parameters from input files
        if(initfile(levelmin).ne.' '.and.filetype.eq.'grafic')then
           call init_file
        endif
        t=0.0
        aexp=1.0
     end if
  end if

  if(cosmo)then

     ! Allocate look-up tables
     n_frw=1000
     allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
     allocate(tau_frw(0:n_frw),t_frw(0:n_frw))

     ! Compute Friedman model look up table
     if(myid==1)write(*,*)'Computing Friedman model'
     call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
          & 1.d-6,dble(aexp_ini), &
          & aexp_frw,hexp_frw,tau_frw,t_frw,n_frw)

     ! Compute initial conformal time                                    
     ! Find neighboring expansion factors                                  
     i=1                                                                   
     do while(aexp_frw(i)>aexp.and.i<n_frw)                                
        i=i+1                                                              
     end do                                                                
     ! Interploate time                                                    
     if(nrestart==0)then                                                   
        t=tau_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &   
             & tau_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i)) 
        aexp=aexp_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &     
             & aexp_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))      
        hexp=hexp_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &     
             & hexp_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))      
     end if                                                                
     texp=t_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &           
          & t_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))            
  else                                                                     
     texp=t                                                                
  end if                                                                   

  ! Initialize cooling model
  if(cooling.and..not.(neq_chem.or.rt))then
     if(myid==1)write(*,*)'Computing cooling model'
     Nmodel=-1
     if(.not. haardt_madau)then
        Nmodel=2
     endif
     if(cosmo)then
        ! Reonization redshift has to be later than starting redshift
        z_reion=min(1./(1.1*aexp_ini)-1.,z_reion)
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(1.0),T2_sim)
     endif
  end if

#ifdef RT
  if(neq_chem.or.rt) then
     if(myid==1)write(*,*)'Computing thermochemistry model'
     Nmodel=-1
     if(.not. haardt_madau)then
        Nmodel=2
     endif
     if(cosmo)then
        ! Reonization redshift has to be later than starting redshift
        z_reion=min(1./(1.1*aexp_ini)-1.,z_reion)
        call rt_set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call rt_set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(1.0),T2_sim)
     endif
  end if
#endif

end subroutine init_time

subroutine init_file
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
  !------------------------------------------------------
  ! Read geometrical parameters in the initial condition files.
  ! Initial conditions are supposed to be made by 
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel,nx_loc,ny_loc,nz_loc
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=80)::filename
  logical::ok

  if(verbose)write(*,*)'Entering init_file'

  ! Reading initial conditions parameters only
  nlevelmax_part=levelmin-1
  do ilevel=levelmin,nlevelmax
     if(initfile(ilevel).ne.' ')then
        filename=TRIM(initfile(ilevel))//'/ic_d'
        INQUIRE(file=filename,exist=ok)
        if(.not.ok)then
           if(myid==1)then
              write(*,*)'File '//TRIM(filename)//' does not exist'
           end if
           call clean_stop
        end if
        open(10,file=filename,form='unformatted')
        if(myid==1)write(*,*)'Reading file '//TRIM(filename)
        rewind 10
        read(10)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
             & ,xoff10,xoff20,xoff30 &
             & ,astart0,omega_m0,omega_l0,h00
        close(10)
        dxini(ilevel)=dxini0
        xoff1(ilevel)=xoff10
        xoff2(ilevel)=xoff20
        xoff3(ilevel)=xoff30
        nlevelmax_part=nlevelmax_part+1
     endif
  end do

  ! Check compatibility with run parameters
  nx_loc=icoarse_max-icoarse_min+1
  ny_loc=jcoarse_max-jcoarse_min+1
  nz_loc=kcoarse_max-kcoarse_min+1
  if(         nx_loc.ne.n1(levelmin)/2**levelmin &
       & .or. ny_loc.ne.n2(levelmin)/2**levelmin &
       & .or. nz_loc.ne.n3(levelmin)/2**levelmin) then 
     write(*,*)'coarser grid is not compatible with initial conditions file'
     write(*,*)'Found    n1=',n1(levelmin),&
          &            ' n2=',n2(levelmin),&
          &            ' n3=',n3(levelmin)
     write(*,*)'Expected n1=',nx_loc*2**levelmin &
          &           ,' n2=',ny_loc*2**levelmin &
          &           ,' n3=',nz_loc*2**levelmin
     call clean_stop
  end if

  ! Write initial conditions parameters
  if(myid==1)then
     do ilevel=levelmin,nlevelmax_part
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
             & n1(ilevel),&
             & n2(ilevel),&
             & n3(ilevel)
        write(*,'(" dx=",1pe10.3)')dxini(ilevel)
        write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
             & 1pe10.3)') &
             & xoff1(ilevel),&
             & xoff2(ilevel),&
             & xoff3(ilevel)
     end do
  end if

end subroutine init_file


subroutine init_cosmo
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod

  implicit none
  !------------------------------------------------------
  ! Read cosmological and geometrical parameters
  ! in the initial condition files.
  ! Initial conditions are supposed to be made by 
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=80)::filename
  character(LEN=5)::nchar
  logical::ok
  TYPE(gadgetheadertype) :: gadgetheader 
  integer::i

  if(verbose)write(*,*)'Entering init_cosmo'

  if(initfile(levelmin)==' ')then
     write(*,*)'You need to specifiy at least one level of initial condition'
     call clean_stop
  end if

  SELECT CASE (filetype)
  case ('grafic', 'ascii')
     ! Reading initial conditions parameters only
     aexp=2.0
     nlevelmax_part=levelmin-1
     do ilevel=levelmin,nlevelmax
        if(initfile(ilevel).ne.' ')then
           if(multiple)then
              call title(myid,nchar)
              filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
           else
              filename=TRIM(initfile(ilevel))//'/ic_deltab'
           endif
           INQUIRE(file=filename,exist=ok)
           if(.not.ok)then
              if(myid==1)then
                 write(*,*)'File '//TRIM(filename)//' does not exist'
              end if
              call clean_stop
           end if
           open(10,file=filename,form='unformatted')
           if(myid==1)write(*,*)'Reading file '//TRIM(filename)
           rewind 10
           read(10)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
                & ,xoff10,xoff20,xoff30 &
                & ,astart0,omega_m0,omega_l0,h00
           close(10)
           dxini(ilevel)=dxini0
           xoff1(ilevel)=xoff10
           xoff2(ilevel)=xoff20
           xoff3(ilevel)=xoff30
           astart(ilevel)=astart0
           omega_m=omega_m0
           omega_l=omega_l0
           if(hydro)omega_b=0.045
           !!!if(hydro)omega_b=0.999999*omega_m
           h0=h00
           aexp=MIN(aexp,astart(ilevel))
           nlevelmax_part=nlevelmax_part+1
           ! Compute SPH equivalent mass (initial gas mass resolution)
           mass_sph=omega_b/omega_m*0.5d0**(ndim*ilevel)
           
        endif
     end do

     ! Compute initial expansion factor
     if(aexp_ini.lt.1.0)then
        aexp=aexp_ini
     else
        aexp_ini=aexp
     endif
     
     ! Check compatibility with run parameters
     if(.not. multiple) then
        if(         nx.ne.n1(levelmin)/2**levelmin &
             & .or. ny.ne.n2(levelmin)/2**levelmin &
             & .or. nz.ne.n3(levelmin)/2**levelmin) then 
           write(*,*)'coarser grid is not compatible with initial conditions file'
           write(*,*)'Found    n1=',n1(levelmin),&
                &            ' n2=',n2(levelmin),&
                &            ' n3=',n3(levelmin)
           write(*,*)'Expected n1=',nx*2**levelmin &
                &           ,' n2=',ny*2**levelmin &
                &           ,' n3=',nz*2**levelmin
           call clean_stop
        endif
     end if
     
     ! Compute box length in the initial conditions in units of h-1 Mpc
     boxlen_ini=dble(nx)*2**levelmin*dxini(levelmin)*(h0/100.)
     
  CASE ('gadget')
     if (verbose) write(*,*)'Reading in gadget format from '//TRIM(initfile(levelmin))
     call gadgetreadheader(TRIM(initfile(levelmin)), 0, gadgetheader, ok)
     if(.not.ok) call clean_stop
     do i=1,6
        if (i .ne. 2) then
           if (gadgetheader%nparttotal(i) .ne. 0) then
              write(*,*) 'Non DM particles present in bin ', i
              call clean_stop
           endif
        endif
     enddo
     if (gadgetheader%mass(2) == 0) then
        write(*,*) 'Particles have different masses, not supported'
        call clean_stop
     endif
     omega_m = gadgetheader%omega0
     omega_l = gadgetheader%omegalambda
     h0 = gadgetheader%hubbleparam * 100.d0
     boxlen_ini = gadgetheader%boxsize
     aexp = gadgetheader%time
     aexp_ini = aexp
     ! Compute SPH equivalent mass (initial gas mass resolution)
     mass_sph=omega_b/omega_m*0.5d0**(ndim*levelmin)
     nlevelmax_part = levelmin
     astart(levelmin) = aexp
     xoff1(levelmin)=0
     xoff2(levelmin)=0
     xoff3(levelmin)=0
     dxini(levelmin) = boxlen_ini/(nx*2**levelmin*(h0/100.0))

  CASE DEFAULT
     write(*,*) 'Unsupported input format '//filetype
     call clean_stop
  END SELECT

  ! Write cosmological parameters
  if(myid==1)then
     write(*,'(" Cosmological parameters:")')
     write(*,'(" aexp=",1pe10.3," H0=",1pe10.3," km s-1 Mpc-1")')aexp,h0
     write(*,'(" omega_m=",F7.3," omega_l=",F7.3)')omega_m,omega_l
     write(*,'(" box size=",1pe10.3," h-1 Mpc")')boxlen_ini
  end if
  omega_k=1.d0-omega_l-omega_m
           
  ! Compute linear scaling factor between aexp and astart(ilevel)
  do ilevel=levelmin,nlevelmax_part
     dfact(ilevel)=d1a(aexp)/d1a(astart(ilevel))
     vfact(ilevel)=astart(ilevel)*fpeebl(astart(ilevel)) & ! Same scale factor as in grafic1
          & *sqrt(omega_m/astart(ilevel)+omega_l*astart(ilevel)*astart(ilevel)+omega_k) &
          & /astart(ilevel)*h0
  end do

  ! Write initial conditions parameters
  do ilevel=levelmin,nlevelmax_part
     if(myid==1)then
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" dx=",1pe10.3," h-1 Mpc")')dxini(ilevel)*h0/100.
     endif
     if(.not.multiple)then
        if(myid==1)then
           write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
                & n1(ilevel),&
                & n2(ilevel),&
                & n3(ilevel)
           write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
                & 1pe10.3," h-1 Mpc")') &
                & xoff1(ilevel)*h0/100.,&
                & xoff2(ilevel)*h0/100.,&
                & xoff3(ilevel)*h0/100.
        endif
     else
        write(*,'(" myid=",I4," n1=",I4," n2=",I4," n3=",I4)') &
             & myid,n1(ilevel),n2(ilevel),n3(ilevel)
        write(*,'(" myid=",I4," xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
             & 1pe10.3," h-1 Mpc")') &
             & myid,xoff1(ilevel)*h0/100.,&
             & xoff2(ilevel)*h0/100.,&
             & xoff3(ilevel)*h0/100.
     endif
  end do

  ! Scale displacement in Mpc to code velocity (v=dx/dtau)
  ! in coarse cell units per conformal time
  vfact(1)=aexp*fpeebl(aexp)*sqrt(omega_m/aexp+omega_l*aexp*aexp+omega_k)
  ! This scale factor is different from vfact in grafic by h0/aexp

contains

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fy(a)
    implicit none
    !      Computes the integrand
    real(dp)::fy
    real(dp)::y,a
    
    y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
    fy=1.d0/y**1.5d0
    
    return
  end function fy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function d1a(a)
    implicit none
    real(dp)::d1a
    !     Computes the linear growing mode D1 in a Friedmann-Robertson-Walker
    !     universe. See Peebles LSSU sections 11 and 14.
    real(dp)::a,y12,y,eps
    
    eps=1.0d-6
    if(a .le. 0.0d0)then
       write(*,*)'a=',a
       call clean_stop
    end if
    y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
    if(y .lt. 0.0D0)then
       write(*,*)'y=',y
       call clean_stop
    end if
    y12=y**0.5d0
    d1a=y12/a*rombint(eps,a,eps)
    
    return
  end function d1a
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  function ad1(d1)
!!$    implicit none
!!$    real(dp)::ad1
!!$    real(dp)::a,d1,da
!!$    integer::niter
!!$    ! Inverts the relation d1(a) given by function d1a(a) using 
!!$    ! Newton-Raphson.
!!$    if (d1.eq.0.0) stop 'ad1 undefined for d1=0!'
!!$    ! Initial guess for Newton-Raphson iteration, good for Omega near 1.
!!$    a=1.e-7
!!$    niter=0
!!$10  niter=niter+1
!!$    da=(d1/d1a(a)-1.d0)/fpeebl(a)*a
!!$    a=a+da
!!$    if (abs(da).gt.1.0e-8.and.niter.lt.10) go to 10
!!$    ad1=a
!!$    return
!!$  end function ad1
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fpeebl(a)
    implicit none
    real(dp) :: fpeebl,a
    !     Computes the growth factor f=d\log D1/d\log a.
    real(dp) :: fact,y,eps
    
    eps=1.0d-6
    y=omega_m*(1.d0/a-1.d0) + omega_l*(a*a-1.d0) + 1.d0
    fact=rombint(eps,a,eps)
    fpeebl=(omega_l*a*a-0.5d0*omega_m/a)/y - 1.d0 + a*fy(a)/fact
    return
  end function fpeebl
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rombint(a,b,tol)
    implicit none
    real(dp)::rombint
    !
    !     Rombint returns the integral from a to b of f(x)dx using Romberg 
    !     integration. The method converges provided that f(x) is continuous 
    !     in (a,b). The function f must be double precision and must be 
    !     declared external in the calling routine.  
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fy(a)+fy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !	Calculate next trapezoidal rule approximation to integral.
       
       g0=0.0d0
       do k=1,nint
          g0=g0+fy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0
       
       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    rombint=g0
    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Rombint failed to converge; integral, error=', &
         &    rombint,error
    return
  end function rombint
     
end subroutine init_cosmo

subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable)
  use amr_parameters
  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min
  real(dp),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
     write(*,*)'Error: non-physical cosmological constants'
     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
     write(*,*)'The sum must be equal to 1.0, but '
     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
     call clean_stop
  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

!  write(*,666)-t
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  use amr_parameters
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  use amr_parameters
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt




