MODULE mpiamrvac_interface
    include 'amrvacdef.f'

    CHARACTER(LEN=1024) :: parameters_filename = 'amrvac.par'
    
    INTEGER :: refinement_level = 1
    
    CHARACTER(LEN=1024) :: error_string = ''
    
CONTAINS

    FUNCTION initialize_code()
        use StoppingConditions
        
        IMPLICIT NONE
        INTEGER initialize_code
        logical :: file_exists
        integer :: error
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

        icomm=MPI_COMM_WORLD
        
        print * , "NW IS:", nw
        error = set_support_for_condition(NUMBER_OF_STEPS_DETECTION)
        error = set_support_for_condition(TIMEOUT_DETECTION)
        
        INQUIRE(FILE= TRIM(parameters_filename), EXIST=file_exists)
        
        if (.NOT.file_exists) THEN
            initialize_code = -1
            return
        end if
        inifile = parameters_filename
        
        call readparameters()
        
        time_advance=.true. 
        time_accurate=.true.
        
        eqpar(gamma_) = 5.0d0/3.0d0
        refinement_level = 1
        initialize_code = 0
    END FUNCTION
    
    FUNCTION cleanup_code()
        INTEGER cleanup_code
        cleanup_code = 0
    END FUNCTION
    
    FUNCTION commit_parameters()
        IMPLICIT NONE
        INTEGER commit_parameters
        commit_parameters = 0
        
        commit_parameters = check_method_parameters()
        if(commit_parameters .ne. 0) return
        
        
        call initialize_vars()
        call init_comm_types   
        
        call initglobaldata_usr
        call initglobaldata 
        
        ! form and initialize all grids at level one
        call initlevelone
        
        refinement_level = 1
        
        commit_parameters = 0
    END FUNCTION
    
    FUNCTION recommit_parameters()
        INTEGER recommit_parameters
        recommit_parameters = 0
    END FUNCTION

    FUNCTION initialize_grid()
        IMPLICIT NONE
        INTEGER :: initialize_grid
        
        ! set up and initialize finer level grids, if needed
        !call settree
        if (levmax>levmin) then
            call allocateBflux
        endif
        
        call getbc(t,ixG^LL,pw,pwCoarse,pgeo,&
               pgeoCoarse,.false.)

        initialize_grid = 0
        
    END FUNCTION
    
    function refine_grid(has_advanced)
        IMPLICIT NONE

        ! create and initialize grids on all levels > 1. On entry, all
        ! level=1 grids have been formed and initialized. This subroutine
        ! creates and initializes new level grids

        integer :: igrid, iigrid, levnew, refine_grid
        logical, intent(out) :: has_advanced
        type(walloc) :: pwtmp
        
        has_advanced = .false.
        refine_grid = 0
        
        !----------------------------------------------------------------------------
        ! when only one level allowed, there is nothing to do anymore
        if (mxnest == 1) return
        
        if (refinement_level == mxnest) return
        
        if (refinement_level == 1) then
            call getbc(t,ixG^LL,&
            pw,pwCoarse,pgeo,pgeoCoarse,.false.)
        end if 
        
        levnew = refinement_level + 1
        
        if (errorestimate==1 .or. errorestimate==2) then
            call setdt
            call advance(0)
        end if

        call errest
        
        if (errorestimate==1 .or. errorestimate==2) then
            do iigrid=1,igridstail; igrid=igrids(iigrid);
                pwtmp%w => pwold(igrid)%w
                pwold(igrid)%w => pw(igrid)%w
                pw(igrid)%w => pwtmp%w
            end do
        end if

        call amr_coarsen_refine

        if (levmax/=levnew) return
        
        refinement_level = levmax
        
        has_advanced = .true.
        
    end function refine_grid
    
    
    function get_current_error(outputvalue)
        IMPLICIT NONE
        integer :: get_current_error
        character(len=1024) :: outputvalue
        outputvalue = error_string
        get_current_error = 0   
    end function
    
    function set_gamma(inputvalue)
        IMPLICIT NONE
        integer :: set_gamma
        double precision :: inputvalue
        eqpar(gamma_) = inputvalue
        set_gamma = 0   
    end function

    function get_gamma(outputvalue)
        IMPLICIT NONE
        integer :: get_gamma
        double precision :: outputvalue
        outputvalue = eqpar(gamma_)
        get_gamma = 0   
    end function

    function set_typeentropy(inputvalue)
        IMPLICIT NONE
        integer :: set_typeentropy, iw
        character(len=*),intent(in) :: inputvalue
        
        do iw=1,nw
            typeentropy(iw)=trim(inputvalue)
        end do
        
        typeentropy(1)=trim(inputvalue)
        set_typeentropy = 0   
    end function

    function get_typeentropy(outputvalue)
        IMPLICIT NONE
        integer :: get_typeentropy, iw
        character(len=256), intent(out) :: outputvalue
        character(len=256) :: checkvalue
        
        get_typeentropy = 0 
        
        checkvalue = typeentropy(1)
        
        do iw=1,nw
            if (typeentropy(iw) .ne. checkvalue) then
                get_typeentropy = -1
            end if
        end do
        
        outputvalue = typeentropy(1)
          
    end function

    function set_typefull1(inputvalue)
        IMPLICIT NONE
        integer :: set_typefull1, i
        character(len=*),intent(in) :: inputvalue
        
        do i=1,nlevelshi
            typefull1(i)=trim(inputvalue)
        end do
        
        set_typefull1 = 0   
    end function

    function get_typefull1(outputvalue)
        IMPLICIT NONE
        integer :: get_typefull1, i
        character(len=256), intent(out) :: outputvalue
        character(len=256) :: checkvalue
        
        get_typefull1 = 0 
        
        checkvalue = typefull1(1)
        
        do i=1,nlevelshi
            if (typefull1(i) .ne. checkvalue) then
                get_typefull1 = -1
            end if
        end do
        
        outputvalue = typefull1(1)
          
    end function
    
    
    function set_typepred1(inputvalue)
        IMPLICIT NONE
        integer :: set_typepred1, i
        character(len=*),intent(in) :: inputvalue
        
        do i=1,nlevelshi
            typepred1(i)=trim(inputvalue)
        end do
        
        set_typepred1 = 0   
    end function

    function get_typepred1(outputvalue)
        IMPLICIT NONE
        integer :: get_typepred1, i
        character(len=256), intent(out) :: outputvalue
        character(len=256) :: checkvalue
        
        get_typepred1 = 0 
        
        checkvalue = typepred1(1)
        
        do i=1,nlevelshi
            if (typepred1(i) .ne. checkvalue) then
                get_typepred1 = -1
            end if
        end do
        
        outputvalue = typepred1(1)
          
    end function
    
    function set_dt(inputvalue)
        IMPLICIT NONE
        integer :: set_dt
        double precision :: inputvalue
        dt = inputvalue
        set_dt = 0   
    end function

    function get_dt(outputvalue)
        IMPLICIT NONE
        integer :: get_dt
        double precision :: outputvalue
        outputvalue = dt
        get_dt = 0   
    end function

    function set_nbufferx1(inputvalue)
        IMPLICIT NONE
        integer :: set_nbufferx1
        integer :: inputvalue
        nbufferx1 = inputvalue
        set_nbufferx1 = 0   
    end function

    function get_nbufferx1(outputvalue)
        IMPLICIT NONE
        integer :: get_nbufferx1
        integer :: outputvalue
        outputvalue = nbufferx1
        get_nbufferx1 = 0   
    end function

    function set_nbufferx2(inputvalue)
        IMPLICIT NONE
        integer :: set_nbufferx2
        integer :: inputvalue
{^NOONED
        nbufferx2 = inputvalue
}
        set_nbufferx2 = 0   
    end function

    function get_nbufferx2(outputvalue)
        IMPLICIT NONE
        integer :: get_nbufferx2
        integer :: outputvalue
        outputvalue = 0
{^NOONED
        outputvalue = nbufferx2
}
        get_nbufferx2 = 0   
    end function

    function set_nbufferx3(inputvalue)
        IMPLICIT NONE
        integer :: set_nbufferx3
        integer :: inputvalue
{^IFTHREED
        nbufferx3 = inputvalue
}
        set_nbufferx3 = 0   
    end function

    function get_nbufferx3(outputvalue)
        IMPLICIT NONE
        integer :: get_nbufferx3
        integer :: outputvalue
        outputvalue = 0
{^IFTHREED
        outputvalue = nbufferx3
}
    
        get_nbufferx3 = 0   
    end function

    function set_mxnest(inputvalue)
        IMPLICIT NONE
        integer :: set_mxnest
        integer :: inputvalue
        mxnest = inputvalue
        set_mxnest = 0   
    end function

    function get_mxnest(outputvalue)
        IMPLICIT NONE
        integer :: get_mxnest
        integer :: outputvalue
        outputvalue = mxnest
        get_mxnest = 0   
    end function

    function set_dixb(inputvalue)
        IMPLICIT NONE
        integer :: set_dixb
        integer :: inputvalue
        dixb = inputvalue
        set_dixb = 0   
    end function

    function get_dixb(outputvalue)
        IMPLICIT NONE
        integer :: get_dixb
        integer :: outputvalue
        outputvalue = dixb
        get_dixb = 0   
    end function

    function set_levmin(inputvalue)
        IMPLICIT NONE
        integer :: set_levmin
        integer :: inputvalue
        levmin = inputvalue
        set_levmin = 0   
    end function

    function get_levmin(outputvalue)
        IMPLICIT NONE
        integer :: get_levmin
        integer :: outputvalue
        outputvalue = levmin
        get_levmin = 0   
    end function

    function set_levmax(inputvalue)
        IMPLICIT NONE
        integer :: set_levmax
        integer :: inputvalue
        levmax = inputvalue
        set_levmax = 0   
    end function

    function get_levmax(outputvalue)
        IMPLICIT NONE
        integer :: get_levmax
        integer :: outputvalue
        outputvalue = levmax
        get_levmax = 0   
    end function

    function set_skipfinestep(inputvalue)
        IMPLICIT NONE
        integer :: set_skipfinestep
        logical :: inputvalue
        skipfinestep = inputvalue
        set_skipfinestep = 0   
    end function

    function get_skipfinestep(outputvalue)
        IMPLICIT NONE
        integer :: get_skipfinestep
        logical :: outputvalue
        outputvalue = skipfinestep
        get_skipfinestep = 0   
    end function

    function set_time_advance(inputvalue)
        IMPLICIT NONE
        integer :: set_time_advance
        logical :: inputvalue
        time_advance = inputvalue
        set_time_advance = 0   
    end function

    function get_time_advance(outputvalue)
        IMPLICIT NONE
        integer :: get_time_advance
        logical :: outputvalue
        outputvalue = time_advance
        get_time_advance = 0   
    end function

    function set_courantpar(inputvalue)
        IMPLICIT NONE
        integer :: set_courantpar
        double precision :: inputvalue
        courantpar = inputvalue
        set_courantpar = 0   
    end function

    function get_courantpar(outputvalue)
        IMPLICIT NONE
        integer :: get_courantpar
        double precision :: outputvalue
        outputvalue = courantpar
        get_courantpar = 0   
    end function

    function set_dtpar(inputvalue)
        IMPLICIT NONE
        integer :: set_dtpar
        double precision :: inputvalue
        dtpar = inputvalue
        set_dtpar = 0   
    end function

    function get_dtpar(outputvalue)
        IMPLICIT NONE
        integer :: get_dtpar
        double precision :: outputvalue
        outputvalue = dtpar
        get_dtpar = 0   
    end function

    function set_dtdiffpar(inputvalue)
        IMPLICIT NONE
        integer :: set_dtdiffpar
        double precision :: inputvalue
        dtdiffpar = inputvalue
        set_dtdiffpar = 0   
    end function

    function get_dtdiffpar(outputvalue)
        IMPLICIT NONE
        integer :: get_dtdiffpar
        double precision :: outputvalue
        outputvalue = dtdiffpar
        get_dtdiffpar = 0   
    end function

    function set_t(inputvalue)
        IMPLICIT NONE
        integer :: set_t
        double precision :: inputvalue
        t = inputvalue
        set_t = 0   
    end function

    function get_t(outputvalue)
        IMPLICIT NONE
        integer :: get_t
        double precision :: outputvalue
        outputvalue = t
        get_t = 0   
    end function

    function set_tmax(inputvalue)
        IMPLICIT NONE
        integer :: set_tmax
        double precision :: inputvalue
        tmax = inputvalue
        set_tmax = 0   
    end function

    function get_tmax(outputvalue)
        IMPLICIT NONE
        integer :: get_tmax
        double precision :: outputvalue
        outputvalue = tmax
        get_tmax = 0   
    end function

    function set_dtmin(inputvalue)
        IMPLICIT NONE
        integer :: set_dtmin
        double precision :: inputvalue
        dtmin = inputvalue
        set_dtmin = 0   
    end function

    function get_dtmin(outputvalue)
        IMPLICIT NONE
        integer :: get_dtmin
        double precision :: outputvalue
        outputvalue = dtmin
        get_dtmin = 0   
    end function

    function set_residmin(inputvalue)
        IMPLICIT NONE
        integer :: set_residmin
        double precision :: inputvalue
        residmin = inputvalue
        set_residmin = 0   
    end function

    function get_residmin(outputvalue)
        IMPLICIT NONE
        integer :: get_residmin
        double precision :: outputvalue
        outputvalue = residmin
        get_residmin = 0   
    end function

    function set_residmax(inputvalue)
        IMPLICIT NONE
        integer :: set_residmax
        double precision :: inputvalue
        residmax = inputvalue
        set_residmax = 0   
    end function

    function get_residmax(outputvalue)
        IMPLICIT NONE
        integer :: get_residmax
        double precision :: outputvalue
        outputvalue = residmax
        get_residmax = 0   
    end function

    function set_residual(inputvalue)
        IMPLICIT NONE
        integer :: set_residual
        double precision :: inputvalue
        residual = inputvalue
        set_residual = 0   
    end function

    function get_residual(outputvalue)
        IMPLICIT NONE
        integer :: get_residual
        double precision :: outputvalue
        outputvalue = residual
        get_residual = 0   
    end function

    function set_tfixgrid(inputvalue)
        IMPLICIT NONE
        integer :: set_tfixgrid
        double precision :: inputvalue
        tfixgrid = inputvalue
        set_tfixgrid = 0   
    end function

    function get_tfixgrid(outputvalue)
        IMPLICIT NONE
        integer :: get_tfixgrid
        double precision :: outputvalue
        outputvalue = tfixgrid
        get_tfixgrid = 0   
    end function

    function set_tvdlfeps(inputvalue)
        IMPLICIT NONE
        integer :: set_tvdlfeps
        double precision :: inputvalue
        tvdlfeps = inputvalue
        set_tvdlfeps = 0   
    end function

    function get_tvdlfeps(outputvalue)
        IMPLICIT NONE
        integer :: get_tvdlfeps
        double precision :: outputvalue
        outputvalue = tvdlfeps
        get_tvdlfeps = 0   
    end function

    function set_mcbeta(inputvalue)
        IMPLICIT NONE
        integer :: set_mcbeta
        double precision :: inputvalue
        mcbeta = inputvalue
        set_mcbeta = 0   
    end function

    function get_mcbeta(outputvalue)
        IMPLICIT NONE
        integer :: get_mcbeta
        double precision :: outputvalue
        outputvalue = mcbeta
        get_mcbeta = 0   
    end function

    function set_divbdiff(inputvalue)
        IMPLICIT NONE
        integer :: set_divbdiff
        double precision :: inputvalue
        divbdiff = inputvalue
        set_divbdiff = 0   
    end function

    function get_divbdiff(outputvalue)
        IMPLICIT NONE
        integer :: get_divbdiff
        double precision :: outputvalue
        outputvalue = divbdiff
        get_divbdiff = 0   
    end function

    function set_smallp(inputvalue)
        IMPLICIT NONE
        integer :: set_smallp
        double precision :: inputvalue
        smallp = inputvalue
        set_smallp = 0   
    end function

    function get_smallp(outputvalue)
        IMPLICIT NONE
        integer :: get_smallp
        double precision :: outputvalue
        outputvalue = smallp
        get_smallp = 0   
    end function

    function set_smallrho(inputvalue)
        IMPLICIT NONE
        integer :: set_smallrho
        double precision :: inputvalue
        smallrho = inputvalue
        set_smallrho = 0   
    end function

    function get_smallrho(outputvalue)
        IMPLICIT NONE
        integer :: get_smallrho
        double precision :: outputvalue
        outputvalue = smallrho
        get_smallrho = 0   
    end function

    function set_dmaxvel(inputvalue)
        IMPLICIT NONE
        integer :: set_dmaxvel
        double precision :: inputvalue
        dmaxvel = inputvalue
        set_dmaxvel = 0   
    end function

    function get_dmaxvel(outputvalue)
        IMPLICIT NONE
        integer :: get_dmaxvel
        double precision :: outputvalue
        outputvalue = dmaxvel
        get_dmaxvel = 0   
    end function

    function set_tolernr(inputvalue)
        IMPLICIT NONE
        integer :: set_tolernr
        double precision :: inputvalue
        tolernr = inputvalue
        set_tolernr = 0   
    end function

    function get_tolernr(outputvalue)
        IMPLICIT NONE
        integer :: get_tolernr
        double precision :: outputvalue
        outputvalue = tolernr
        get_tolernr = 0   
    end function

    function set_absaccnr(inputvalue)
        IMPLICIT NONE
        integer :: set_absaccnr
        double precision :: inputvalue
        absaccnr = inputvalue
        set_absaccnr = 0   
    end function

    function get_absaccnr(outputvalue)
        IMPLICIT NONE
        integer :: get_absaccnr
        double precision :: outputvalue
        outputvalue = absaccnr
        get_absaccnr = 0   
    end function

    function set_cfrac(inputvalue)
        IMPLICIT NONE
        integer :: set_cfrac
        double precision :: inputvalue
        cfrac = inputvalue
        set_cfrac = 0   
    end function

    function get_cfrac(outputvalue)
        IMPLICIT NONE
        integer :: get_cfrac
        double precision :: outputvalue
        outputvalue = cfrac
        get_cfrac = 0   
    end function

    function set_x1ptms(inputvalue)
        IMPLICIT NONE
        integer :: set_x1ptms
        double precision :: inputvalue
        x1ptms = inputvalue
        set_x1ptms = 0   
    end function

    function get_x1ptms(outputvalue)
        IMPLICIT NONE
        integer :: get_x1ptms
        double precision :: outputvalue
        outputvalue = x1ptms
        get_x1ptms = 0   
    end function

    function set_x2ptms(inputvalue)
        IMPLICIT NONE
        integer :: set_x2ptms
        double precision :: inputvalue
        x2ptms = inputvalue
        set_x2ptms = 0   
    end function

    function get_x2ptms(outputvalue)
        IMPLICIT NONE
        integer :: get_x2ptms
        double precision :: outputvalue
        outputvalue = x2ptms
        get_x2ptms = 0   
    end function

    function set_x3ptms(inputvalue)
        IMPLICIT NONE
        integer :: set_x3ptms
        double precision :: inputvalue
        x3ptms = inputvalue
        set_x3ptms = 0   
    end function

    function get_x3ptms(outputvalue)
        IMPLICIT NONE
        integer :: get_x3ptms
        double precision :: outputvalue
        outputvalue = x3ptms
        get_x3ptms = 0   
    end function

    function set_ptmass(inputvalue)
        IMPLICIT NONE
        integer :: set_ptmass
        double precision :: inputvalue
        ptmass = inputvalue
        set_ptmass = 0   
    end function

    function get_ptmass(outputvalue)
        IMPLICIT NONE
        integer :: get_ptmass
        double precision :: outputvalue
        outputvalue = ptmass
        get_ptmass = 0   
    end function

    function set_ratebdflux(inputvalue)
        IMPLICIT NONE
        integer :: set_ratebdflux
        double precision :: inputvalue
        ratebdflux = inputvalue
        set_ratebdflux = 0   
    end function

    function get_ratebdflux(outputvalue)
        IMPLICIT NONE
        integer :: get_ratebdflux
        double precision :: outputvalue
        outputvalue = ratebdflux
        get_ratebdflux = 0   
    end function

    function set_normt(inputvalue)
        IMPLICIT NONE
        integer :: set_normt
        double precision :: inputvalue
        normt = inputvalue
        set_normt = 0   
    end function

    function get_normt(outputvalue)
        IMPLICIT NONE
        integer :: get_normt
        double precision :: outputvalue
        outputvalue = normt
        get_normt = 0   
    end function

    function set_time_bc(inputvalue)
        IMPLICIT NONE
        integer :: set_time_bc
        double precision :: inputvalue
        time_bc = inputvalue
        set_time_bc = 0   
    end function

    function get_time_bc(outputvalue)
        IMPLICIT NONE
        integer :: get_time_bc
        double precision :: outputvalue
        outputvalue = time_bc
        get_time_bc = 0   
    end function

    function set_it(inputvalue)
        IMPLICIT NONE
        integer :: set_it
        integer :: inputvalue
        it = inputvalue
        set_it = 0   
    end function

    function get_it(outputvalue)
        IMPLICIT NONE
        integer :: get_it
        integer :: outputvalue
        outputvalue = it
        get_it = 0   
    end function

    function set_itmax(inputvalue)
        IMPLICIT NONE
        integer :: set_itmax
        integer :: inputvalue
        itmax = inputvalue
        set_itmax = 0   
    end function

    function get_itmax(outputvalue)
        IMPLICIT NONE
        integer :: get_itmax
        integer :: outputvalue
        outputvalue = itmax
        get_itmax = 0   
    end function

    function set_itmin(inputvalue)
        IMPLICIT NONE
        integer :: set_itmin
        integer :: inputvalue
        itmin = inputvalue
        set_itmin = 0   
    end function

    function get_itmin(outputvalue)
        IMPLICIT NONE
        integer :: get_itmin
        integer :: outputvalue
        outputvalue = itmin
        get_itmin = 0   
    end function

    function set_slowsteps(inputvalue)
        IMPLICIT NONE
        integer :: set_slowsteps
        integer :: inputvalue
        slowsteps = inputvalue
        set_slowsteps = 0   
    end function

    function get_slowsteps(outputvalue)
        IMPLICIT NONE
        integer :: get_slowsteps
        integer :: outputvalue
        outputvalue = slowsteps
        get_slowsteps = 0   
    end function

    function set_typepario(inputvalue)
        IMPLICIT NONE
        integer :: set_typepario
        integer :: inputvalue
        typepario = inputvalue
        set_typepario = 0   
    end function

    function get_typepario(outputvalue)
        IMPLICIT NONE
        integer :: get_typepario
        integer :: outputvalue
        outputvalue = typepario
        get_typepario = 0   
    end function

    function set_itfixgrid(inputvalue)
        IMPLICIT NONE
        integer :: set_itfixgrid
        integer :: inputvalue
        itfixgrid = inputvalue
        set_itfixgrid = 0   
    end function

    function get_itfixgrid(outputvalue)
        IMPLICIT NONE
        integer :: get_itfixgrid
        integer :: outputvalue
        outputvalue = itfixgrid
        get_itfixgrid = 0   
    end function

    function set_nwauxio(inputvalue)
        IMPLICIT NONE
        integer :: set_nwauxio
        integer :: inputvalue
        nwauxio = inputvalue
        set_nwauxio = 0   
    end function

    function get_nwauxio(outputvalue)
        IMPLICIT NONE
        integer :: get_nwauxio
        integer :: outputvalue
        outputvalue = nwauxio
        get_nwauxio = 0   
    end function

    function set_istep(inputvalue)
        IMPLICIT NONE
        integer :: set_istep
        integer :: inputvalue
        istep = inputvalue
        set_istep = 0   
    end function

    function get_istep(outputvalue)
        IMPLICIT NONE
        integer :: get_istep
        integer :: outputvalue
        outputvalue = istep
        get_istep = 0   
    end function

    function set_nstep(inputvalue)
        IMPLICIT NONE
        integer :: set_nstep
        integer :: inputvalue
        nstep = inputvalue
        set_nstep = 0   
    end function

    function get_nstep(outputvalue)
        IMPLICIT NONE
        integer :: get_nstep
        integer :: outputvalue
        outputvalue = nstep
        get_nstep = 0   
    end function

    function set_errorestimate(inputvalue)
        IMPLICIT NONE
        integer :: set_errorestimate
        integer :: inputvalue
        errorestimate = inputvalue
        set_errorestimate = 0   
    end function

    function get_errorestimate(outputvalue)
        IMPLICIT NONE
        integer :: get_errorestimate
        integer :: outputvalue
        outputvalue = errorestimate
        get_errorestimate = 0   
    end function

    function set_nxdiffusehllc(inputvalue)
        IMPLICIT NONE
        integer :: set_nxdiffusehllc
        integer :: inputvalue
        nxdiffusehllc = inputvalue
        set_nxdiffusehllc = 0   
    end function

    function get_nxdiffusehllc(outputvalue)
        IMPLICIT NONE
        integer :: get_nxdiffusehllc
        integer :: outputvalue
        outputvalue = nxdiffusehllc
        get_nxdiffusehllc = 0   
    end function

    function set_typespherical(inputvalue)
        IMPLICIT NONE
        integer :: set_typespherical
        integer :: inputvalue
        typespherical = inputvalue
        set_typespherical = 0   
    end function

    function get_typespherical(outputvalue)
        IMPLICIT NONE
        integer :: get_typespherical
        integer :: outputvalue
        outputvalue = typespherical
        get_typespherical = 0   
    end function

    function set_maxitnr(inputvalue)
        IMPLICIT NONE
        integer :: set_maxitnr
        integer :: inputvalue
        maxitnr = inputvalue
        set_maxitnr = 0   
    end function

    function get_maxitnr(outputvalue)
        IMPLICIT NONE
        integer :: get_maxitnr
        integer :: outputvalue
        outputvalue = maxitnr
        get_maxitnr = 0   
    end function

    function set_nflatgetaux(inputvalue)
        IMPLICIT NONE
        integer :: set_nflatgetaux
        integer :: inputvalue
        nflatgetaux = inputvalue
        set_nflatgetaux = 0   
    end function

    function get_nflatgetaux(outputvalue)
        IMPLICIT NONE
        integer :: get_nflatgetaux
        integer :: outputvalue
        outputvalue = nflatgetaux
        get_nflatgetaux = 0   
    end function

    function set_level_io(inputvalue)
        IMPLICIT NONE
        integer :: set_level_io
        integer :: inputvalue
        level_io = inputvalue
        set_level_io = 0   
    end function

    function get_level_io(outputvalue)
        IMPLICIT NONE
        integer :: get_level_io
        integer :: outputvalue
        outputvalue = level_io
        get_level_io = 0   
    end function

    function set_ncool(inputvalue)
        IMPLICIT NONE
        integer :: set_ncool
        integer :: inputvalue
        ncool = inputvalue
        set_ncool = 0   
    end function

    function get_ncool(outputvalue)
        IMPLICIT NONE
        integer :: get_ncool
        integer :: outputvalue
        outputvalue = ncool
        get_ncool = 0   
    end function

    function set_cmulti(inputvalue)
        IMPLICIT NONE
        integer :: set_cmulti
        integer :: inputvalue
        cmulti = inputvalue
        set_cmulti = 0   
    end function

    function get_cmulti(outputvalue)
        IMPLICIT NONE
        integer :: get_cmulti
        integer :: outputvalue
        outputvalue = cmulti
        get_cmulti = 0   
    end function

    function set_snapshotini(inputvalue)
        IMPLICIT NONE
        integer :: set_snapshotini
        integer :: inputvalue
        snapshotini = inputvalue
        set_snapshotini = 0   
    end function

    function get_snapshotini(outputvalue)
        IMPLICIT NONE
        integer :: get_snapshotini
        integer :: outputvalue
        outputvalue = snapshotini
        get_snapshotini = 0   
    end function

    function set_ixtest1(inputvalue)
        IMPLICIT NONE
        integer :: set_ixtest1
        integer :: inputvalue
        ixtest1 = inputvalue
        set_ixtest1 = 0   
    end function

    function get_ixtest1(outputvalue)
        IMPLICIT NONE
        integer :: get_ixtest1
        integer :: outputvalue
        outputvalue = ixtest1
        get_ixtest1 = 0   
    end function

    function set_ixtest2(inputvalue)
        IMPLICIT NONE
        integer :: set_ixtest2
        integer :: inputvalue
        ixtest2 = inputvalue
        set_ixtest2 = 0   
    end function

    function get_ixtest2(outputvalue)
        IMPLICIT NONE
        integer :: get_ixtest2
        integer :: outputvalue
        outputvalue = ixtest2
        get_ixtest2 = 0   
    end function

    function set_ixtest3(inputvalue)
        IMPLICIT NONE
        integer :: set_ixtest3
        integer :: inputvalue
        ixtest3 = inputvalue
        set_ixtest3 = 0   
    end function

    function get_ixtest3(outputvalue)
        IMPLICIT NONE
        integer :: get_ixtest3
        integer :: outputvalue
        outputvalue = ixtest3
        get_ixtest3 = 0   
    end function

    function set_iwtest(inputvalue)
        IMPLICIT NONE
        integer :: set_iwtest
        integer :: inputvalue
        iwtest = inputvalue
        set_iwtest = 0   
    end function

    function get_iwtest(outputvalue)
        IMPLICIT NONE
        integer :: get_iwtest
        integer :: outputvalue
        outputvalue = iwtest
        get_iwtest = 0   
    end function

    function set_idimtest(inputvalue)
        IMPLICIT NONE
        integer :: set_idimtest
        integer :: inputvalue
        idimtest = inputvalue
        set_idimtest = 0   
    end function

    function get_idimtest(outputvalue)
        IMPLICIT NONE
        integer :: get_idimtest
        integer :: outputvalue
        outputvalue = idimtest
        get_idimtest = 0   
    end function

    function set_saveigrid(inputvalue)
        IMPLICIT NONE
        integer :: set_saveigrid
        integer :: inputvalue
        saveigrid = inputvalue
        set_saveigrid = 0   
    end function

    function get_saveigrid(outputvalue)
        IMPLICIT NONE
        integer :: get_saveigrid
        integer :: outputvalue
        outputvalue = saveigrid
        get_saveigrid = 0   
    end function

    function set_typecourant(inputvalue)
        IMPLICIT NONE
        integer :: set_typecourant
        character*79 :: inputvalue
        typecourant = inputvalue
        set_typecourant = 0   
    end function

    function get_typecourant(outputvalue)
        IMPLICIT NONE
        integer :: get_typecourant
        character*79 :: outputvalue
        outputvalue = typecourant
        get_typecourant = 0   
    end function

    function set_typeresid(inputvalue)
        IMPLICIT NONE
        integer :: set_typeresid
        character*79 :: inputvalue
        typeresid = inputvalue
        set_typeresid = 0   
    end function

    function get_typeresid(outputvalue)
        IMPLICIT NONE
        integer :: get_typeresid
        character*79 :: outputvalue
        outputvalue = typeresid
        get_typeresid = 0   
    end function

    function set_typeadvance(inputvalue)
        IMPLICIT NONE
        integer :: set_typeadvance
        character*79 :: inputvalue
        typeadvance = inputvalue
        set_typeadvance = 0   
    end function

    function get_typeadvance(outputvalue)
        IMPLICIT NONE
        integer :: get_typeadvance
        character*79 :: outputvalue
        outputvalue = typeadvance
        get_typeadvance = 0   
    end function

    function set_typelimited(inputvalue)
        IMPLICIT NONE
        integer :: set_typelimited
        character*79 :: inputvalue
        typelimited = inputvalue
        set_typelimited = 0   
    end function

    function get_typelimited(outputvalue)
        IMPLICIT NONE
        integer :: get_typelimited
        character*79 :: outputvalue
        outputvalue = typelimited
        get_typelimited = 0   
    end function

    function set_typesourcesplit(inputvalue)
        IMPLICIT NONE
        integer :: set_typesourcesplit
        character*79 :: inputvalue
        typesourcesplit = inputvalue
        set_typesourcesplit = 0   
    end function

    function get_typesourcesplit(outputvalue)
        IMPLICIT NONE
        integer :: get_typesourcesplit
        character*79 :: outputvalue
        outputvalue = typesourcesplit
        get_typesourcesplit = 0   
    end function

    function set_typelimiter(inputvalue)
        IMPLICIT NONE
        integer :: set_typelimiter
        character*79 :: inputvalue
        typelimiter = inputvalue
        set_typelimiter = 0   
    end function

    function get_typelimiter(outputvalue)
        IMPLICIT NONE
        integer :: get_typelimiter
        character*79 :: outputvalue
        outputvalue = typelimiter
        get_typelimiter = 0   
    end function

    function set_typegradlimiter(inputvalue)
        IMPLICIT NONE
        integer :: set_typegradlimiter
        character*79 :: inputvalue
        typegradlimiter = inputvalue
        set_typegradlimiter = 0   
    end function

    function get_typegradlimiter(outputvalue)
        IMPLICIT NONE
        integer :: get_typegradlimiter
        character*79 :: outputvalue
        outputvalue = typegradlimiter
        get_typegradlimiter = 0   
    end function

    function set_typeprolonglimit(inputvalue)
        IMPLICIT NONE
        integer :: set_typeprolonglimit
        character*79 :: inputvalue
        typeprolonglimit = inputvalue
        set_typeprolonglimit = 0   
    end function

    function get_typeprolonglimit(outputvalue)
        IMPLICIT NONE
        integer :: get_typeprolonglimit
        character*79 :: outputvalue
        outputvalue = typeprolonglimit
        get_typeprolonglimit = 0   
    end function

    function set_typetvd(inputvalue)
        IMPLICIT NONE
        integer :: set_typetvd
        character*79 :: inputvalue
        typetvd = inputvalue
        set_typetvd = 0   
    end function

    function get_typetvd(outputvalue)
        IMPLICIT NONE
        integer :: get_typetvd
        character*79 :: outputvalue
        outputvalue = typetvd
        get_typetvd = 0   
    end function

    function set_typetvdlf(inputvalue)
        IMPLICIT NONE
        integer :: set_typetvdlf
        character*79 :: inputvalue
        typetvdlf = inputvalue
        set_typetvdlf = 0   
    end function

    function get_typetvdlf(outputvalue)
        IMPLICIT NONE
        integer :: get_typetvdlf
        character*79 :: outputvalue
        outputvalue = typetvdlf
        get_typetvdlf = 0   
    end function

    function set_typeaverage(inputvalue)
        IMPLICIT NONE
        integer :: set_typeaverage
        character*79 :: inputvalue
        typeaverage = inputvalue
        set_typeaverage = 0   
    end function

    function get_typeaverage(outputvalue)
        IMPLICIT NONE
        integer :: get_typeaverage
        character*79 :: outputvalue
        outputvalue = typeaverage
        get_typeaverage = 0   
    end function

    function set_typedimsplit(inputvalue)
        IMPLICIT NONE
        integer :: set_typedimsplit
        character*79 :: inputvalue
        typedimsplit = inputvalue
        set_typedimsplit = 0   
    end function

    function get_typedimsplit(outputvalue)
        IMPLICIT NONE
        integer :: get_typedimsplit
        character*79 :: outputvalue
        outputvalue = typedimsplit
        get_typedimsplit = 0   
    end function

    function set_typeaxial(inputvalue)
        IMPLICIT NONE
        integer :: set_typeaxial
        character*79 :: inputvalue
        typeaxial = inputvalue
        set_typeaxial = 0   
    end function

    function get_typeaxial(outputvalue)
        IMPLICIT NONE
        integer :: get_typeaxial
        character*79 :: outputvalue
        outputvalue = typeaxial
        get_typeaxial = 0   
    end function

    function set_typepoly(inputvalue)
        IMPLICIT NONE
        integer :: set_typepoly
        character*79 :: inputvalue
        typepoly = inputvalue
        set_typepoly = 0   
    end function

    function get_typepoly(outputvalue)
        IMPLICIT NONE
        integer :: get_typepoly
        character*79 :: outputvalue
        outputvalue = typepoly
        get_typepoly = 0   
    end function

    function set_typedivbdiff(inputvalue)
        IMPLICIT NONE
        integer :: set_typedivbdiff
        character*79 :: inputvalue
        typedivbdiff = inputvalue
        set_typedivbdiff = 0   
    end function

    function get_typedivbdiff(outputvalue)
        IMPLICIT NONE
        integer :: get_typedivbdiff
        character*79 :: outputvalue
        outputvalue = typedivbdiff
        get_typedivbdiff = 0   
    end function

    function set_typedivbfix(inputvalue)
        IMPLICIT NONE
        integer :: set_typedivbfix
        character*79 :: inputvalue
        typedivbfix = inputvalue
        set_typedivbfix = 0   
    end function

    function get_typedivbfix(outputvalue)
        IMPLICIT NONE
        integer :: get_typedivbfix
        character*79 :: outputvalue
        outputvalue = typedivbfix
        get_typedivbfix = 0   
    end function

    function set_typediv(inputvalue)
        IMPLICIT NONE
        integer :: set_typediv
        character*79 :: inputvalue
        typediv = inputvalue
        set_typediv = 0   
    end function

    function get_typediv(outputvalue)
        IMPLICIT NONE
        integer :: get_typediv
        character*79 :: outputvalue
        outputvalue = typediv
        get_typediv = 0   
    end function

    function set_typegrad(inputvalue)
        IMPLICIT NONE
        integer :: set_typegrad
        character*79 :: inputvalue
        typegrad = inputvalue
        set_typegrad = 0   
    end function

    function get_typegrad(outputvalue)
        IMPLICIT NONE
        integer :: get_typegrad
        character*79 :: outputvalue
        outputvalue = typegrad
        get_typegrad = 0   
    end function

    function set_typeglm(inputvalue)
        IMPLICIT NONE
        integer :: set_typeglm
        character*79 :: inputvalue
        typeglm = inputvalue
        set_typeglm = 0   
    end function

    function get_typeglm(outputvalue)
        IMPLICIT NONE
        integer :: get_typeglm
        character*79 :: outputvalue
        outputvalue = typeglm
        get_typeglm = 0   
    end function

    function set_coolcurve(inputvalue)
        IMPLICIT NONE
        integer :: set_coolcurve
        character*79 :: inputvalue
        coolcurve = inputvalue
        set_coolcurve = 0   
    end function

    function get_coolcurve(outputvalue)
        IMPLICIT NONE
        integer :: get_coolcurve
        character*79 :: outputvalue
        outputvalue = coolcurve
        get_coolcurve = 0   
    end function

    function set_coolmethod(inputvalue)
        IMPLICIT NONE
        integer :: set_coolmethod
        character*79 :: inputvalue
        coolmethod = inputvalue
        set_coolmethod = 0   
    end function

    function get_coolmethod(outputvalue)
        IMPLICIT NONE
        integer :: get_coolmethod
        character*79 :: outputvalue
        outputvalue = coolmethod
        get_coolmethod = 0   
    end function

    function set_typeghostfill(inputvalue)
        IMPLICIT NONE
        integer :: set_typeghostfill
        character*79 :: inputvalue
        typeghostfill = inputvalue
        set_typeghostfill = 0   
    end function

    function get_typeghostfill(outputvalue)
        IMPLICIT NONE
        integer :: get_typeghostfill
        character*79 :: outputvalue
        outputvalue = typeghostfill
        get_typeghostfill = 0   
    end function

    function set_typegridfill(inputvalue)
        IMPLICIT NONE
        integer :: set_typegridfill
        character*79 :: inputvalue
        typegridfill = inputvalue
        set_typegridfill = 0   
    end function

    function get_typegridfill(outputvalue)
        IMPLICIT NONE
        integer :: get_typegridfill
        character*79 :: outputvalue
        outputvalue = typegridfill
        get_typegridfill = 0   
    end function

    function set_filenameout(inputvalue)
        IMPLICIT NONE
        integer :: set_filenameout
        character*79 :: inputvalue
        filenameout = inputvalue
        set_filenameout = 0   
    end function

    function get_filenameout(outputvalue)
        IMPLICIT NONE
        integer :: get_filenameout
        character*79 :: outputvalue
        outputvalue = filenameout
        get_filenameout = 0   
    end function

    function set_filenameini(inputvalue)
        IMPLICIT NONE
        integer :: set_filenameini
        character*79 :: inputvalue
        filenameini = inputvalue
        set_filenameini = 0   
    end function

    function get_filenameini(outputvalue)
        IMPLICIT NONE
        integer :: get_filenameini
        character*79 :: outputvalue
        outputvalue = filenameini
        get_filenameini = 0   
    end function

    function set_filenamelog(inputvalue)
        IMPLICIT NONE
        integer :: set_filenamelog
        character*79 :: inputvalue
        filenamelog = inputvalue
        set_filenamelog = 0   
    end function

    function get_filenamelog(outputvalue)
        IMPLICIT NONE
        integer :: get_filenamelog
        character*79 :: outputvalue
        outputvalue = filenamelog
        get_filenamelog = 0   
    end function

    function set_fileheadout(inputvalue)
        IMPLICIT NONE
        integer :: set_fileheadout
        character*79 :: inputvalue
        fileheadout = inputvalue
        set_fileheadout = 0   
    end function

    function get_fileheadout(outputvalue)
        IMPLICIT NONE
        integer :: get_fileheadout
        character*79 :: outputvalue
        outputvalue = fileheadout
        get_fileheadout = 0   
    end function

    function set_wnames(inputvalue)
        IMPLICIT NONE
        integer :: set_wnames
        character*79 :: inputvalue
        wnames = inputvalue
        set_wnames = 0   
    end function

    function get_wnames(outputvalue)
        IMPLICIT NONE
        integer :: get_wnames
        character*79 :: outputvalue
        outputvalue = wnames
        get_wnames = 0   
    end function

    function set_primnames(inputvalue)
        IMPLICIT NONE
        integer :: set_primnames
        character*79 :: inputvalue
        primnames = inputvalue
        set_primnames = 0   
    end function

    function get_primnames(outputvalue)
        IMPLICIT NONE
        integer :: get_primnames
        character*79 :: outputvalue
        outputvalue = primnames
        get_primnames = 0   
    end function

    function set_typefilelog(inputvalue)
        IMPLICIT NONE
        integer :: set_typefilelog
        character*79 :: inputvalue
        typefilelog = inputvalue
        set_typefilelog = 0   
    end function

    function get_typefilelog(outputvalue)
        IMPLICIT NONE
        integer :: get_typefilelog
        character*79 :: outputvalue
        outputvalue = typefilelog
        get_typefilelog = 0   
    end function

    function set_convert_type(inputvalue)
        IMPLICIT NONE
        integer :: set_convert_type
        character*79 :: inputvalue
        convert_type = inputvalue
        set_convert_type = 0   
    end function

    function get_convert_type(outputvalue)
        IMPLICIT NONE
        integer :: get_convert_type
        character*79 :: outputvalue
        outputvalue = convert_type
        get_convert_type = 0   
    end function

    function set_dxfiletype(inputvalue)
        IMPLICIT NONE
        integer :: set_dxfiletype
        character*79 :: inputvalue
        dxfiletype = inputvalue
        set_dxfiletype = 0   
    end function

    function get_dxfiletype(outputvalue)
        IMPLICIT NONE
        integer :: get_dxfiletype
        character*79 :: outputvalue
        outputvalue = dxfiletype
        get_dxfiletype = 0   
    end function

    function set_teststr(inputvalue)
        IMPLICIT NONE
        integer :: set_teststr
        character*79 :: inputvalue
        teststr = inputvalue
        set_teststr = 0   
    end function

    function get_teststr(outputvalue)
        IMPLICIT NONE
        integer :: get_teststr
        character*79 :: outputvalue
        outputvalue = teststr
        get_teststr = 0   
    end function

    function set_time_accurate(inputvalue)
        IMPLICIT NONE
        integer :: set_time_accurate
        logical :: inputvalue
        time_accurate = inputvalue
        set_time_accurate = 0   
    end function

    function get_time_accurate(outputvalue)
        IMPLICIT NONE
        integer :: get_time_accurate
        logical :: outputvalue
        outputvalue = time_accurate
        get_time_accurate = 0   
    end function

    function set_addmpibarrier(inputvalue)
        IMPLICIT NONE
        integer :: set_addmpibarrier
        logical :: inputvalue
        addmpibarrier = inputvalue
        set_addmpibarrier = 0   
    end function

    function get_addmpibarrier(outputvalue)
        IMPLICIT NONE
        integer :: get_addmpibarrier
        logical :: outputvalue
        outputvalue = addmpibarrier
        get_addmpibarrier = 0   
    end function

    function set_tmaxexact(inputvalue)
        IMPLICIT NONE
        integer :: set_tmaxexact
        logical :: inputvalue
        tmaxexact = inputvalue
        set_tmaxexact = 0   
    end function

    function get_tmaxexact(outputvalue)
        IMPLICIT NONE
        integer :: get_tmaxexact
        logical :: outputvalue
        outputvalue = tmaxexact
        get_tmaxexact = 0   
    end function

    function set_treset(inputvalue)
        IMPLICIT NONE
        integer :: set_treset
        logical :: inputvalue
        treset = inputvalue
        set_treset = 0   
    end function

    function get_treset(outputvalue)
        IMPLICIT NONE
        integer :: get_treset
        logical :: outputvalue
        outputvalue = treset
        get_treset = 0   
    end function

    function set_itreset(inputvalue)
        IMPLICIT NONE
        integer :: set_itreset
        logical :: inputvalue
        itreset = inputvalue
        set_itreset = 0   
    end function

    function get_itreset(outputvalue)
        IMPLICIT NONE
        integer :: get_itreset
        logical :: outputvalue
        outputvalue = itreset
        get_itreset = 0   
    end function

    function set_firstprocess(inputvalue)
        IMPLICIT NONE
        integer :: set_firstprocess
        logical :: inputvalue
        firstprocess = inputvalue
        set_firstprocess = 0   
    end function

    function get_firstprocess(outputvalue)
        IMPLICIT NONE
        integer :: get_firstprocess
        logical :: outputvalue
        outputvalue = firstprocess
        get_firstprocess = 0   
    end function

    function set_fixprocess(inputvalue)
        IMPLICIT NONE
        integer :: set_fixprocess
        logical :: inputvalue
        fixprocess = inputvalue
        set_fixprocess = 0   
    end function

    function get_fixprocess(outputvalue)
        IMPLICIT NONE
        integer :: get_fixprocess
        logical :: outputvalue
        outputvalue = fixprocess
        get_fixprocess = 0   
    end function

    function set_flathllc(inputvalue)
        IMPLICIT NONE
        integer :: set_flathllc
        logical :: inputvalue
        flathllc = inputvalue
        set_flathllc = 0   
    end function

    function get_flathllc(outputvalue)
        IMPLICIT NONE
        integer :: get_flathllc
        logical :: outputvalue
        outputvalue = flathllc
        get_flathllc = 0   
    end function

    function set_flatcd(inputvalue)
        IMPLICIT NONE
        integer :: set_flatcd
        logical :: inputvalue
        flatcd = inputvalue
        set_flatcd = 0   
    end function

    function get_flatcd(outputvalue)
        IMPLICIT NONE
        integer :: get_flatcd
        logical :: outputvalue
        outputvalue = flatcd
        get_flatcd = 0   
    end function

    function set_flatsh(inputvalue)
        IMPLICIT NONE
        integer :: set_flatsh
        logical :: inputvalue
        flatsh = inputvalue
        set_flatsh = 0   
    end function

    function get_flatsh(outputvalue)
        IMPLICIT NONE
        integer :: get_flatsh
        logical :: outputvalue
        outputvalue = flatsh
        get_flatsh = 0   
    end function

    function set_flatppm(inputvalue)
        IMPLICIT NONE
        integer :: set_flatppm
        logical :: inputvalue
        flatppm = inputvalue
        set_flatppm = 0   
    end function

    function get_flatppm(outputvalue)
        IMPLICIT NONE
        integer :: get_flatppm
        logical :: outputvalue
        outputvalue = flatppm
        get_flatppm = 0   
    end function

    function set_sourcesplit(inputvalue)
        IMPLICIT NONE
        integer :: set_sourcesplit
        logical :: inputvalue
        sourcesplit = inputvalue
        set_sourcesplit = 0   
    end function

    function get_sourcesplit(outputvalue)
        IMPLICIT NONE
        integer :: get_sourcesplit
        logical :: outputvalue
        outputvalue = sourcesplit
        get_sourcesplit = 0   
    end function

    function set_sourceunsplit(inputvalue)
        IMPLICIT NONE
        integer :: set_sourceunsplit
        logical :: inputvalue
        sourceunsplit = inputvalue
        set_sourceunsplit = 0   
    end function

    function get_sourceunsplit(outputvalue)
        IMPLICIT NONE
        integer :: get_sourceunsplit
        logical :: outputvalue
        outputvalue = sourceunsplit
        get_sourceunsplit = 0   
    end function

    function set_useprimitive(inputvalue)
        IMPLICIT NONE
        integer :: set_useprimitive
        logical :: inputvalue
        useprimitive = inputvalue
        set_useprimitive = 0   
    end function

    function get_useprimitive(outputvalue)
        IMPLICIT NONE
        integer :: get_useprimitive
        logical :: outputvalue
        outputvalue = useprimitive
        get_useprimitive = 0   
    end function

    function set_dimsplit(inputvalue)
        IMPLICIT NONE
        integer :: set_dimsplit
        logical :: inputvalue
        dimsplit = inputvalue
        set_dimsplit = 0   
    end function

    function get_dimsplit(outputvalue)
        IMPLICIT NONE
        integer :: get_dimsplit
        logical :: outputvalue
        outputvalue = dimsplit
        get_dimsplit = 0   
    end function

    function set_restrictprimitive(inputvalue)
        IMPLICIT NONE
        integer :: set_restrictprimitive
        logical :: inputvalue
        restrictprimitive = inputvalue
        set_restrictprimitive = 0   
    end function

    function get_restrictprimitive(outputvalue)
        IMPLICIT NONE
        integer :: get_restrictprimitive
        logical :: outputvalue
        outputvalue = restrictprimitive
        get_restrictprimitive = 0   
    end function

    function set_prolongprimitive(inputvalue)
        IMPLICIT NONE
        integer :: set_prolongprimitive
        logical :: inputvalue
        prolongprimitive = inputvalue
        set_prolongprimitive = 0   
    end function

    function get_prolongprimitive(outputvalue)
        IMPLICIT NONE
        integer :: get_prolongprimitive
        logical :: outputvalue
        outputvalue = prolongprimitive
        get_prolongprimitive = 0   
    end function

    function set_coarsenprimitive(inputvalue)
        IMPLICIT NONE
        integer :: set_coarsenprimitive
        logical :: inputvalue
        coarsenprimitive = inputvalue
        set_coarsenprimitive = 0   
    end function

    function get_coarsenprimitive(outputvalue)
        IMPLICIT NONE
        integer :: get_coarsenprimitive
        logical :: outputvalue
        outputvalue = coarsenprimitive
        get_coarsenprimitive = 0   
    end function

    function set_useprimitiverel(inputvalue)
        IMPLICIT NONE
        integer :: set_useprimitiverel
        logical :: inputvalue
        useprimitiverel = inputvalue
        set_useprimitiverel = 0   
    end function

    function get_useprimitiverel(outputvalue)
        IMPLICIT NONE
        integer :: get_useprimitiverel
        logical :: outputvalue
        outputvalue = useprimitiverel
        get_useprimitiverel = 0   
    end function

    function set_amrentropy(inputvalue)
        IMPLICIT NONE
        integer :: set_amrentropy
        logical :: inputvalue
        amrentropy = inputvalue
        set_amrentropy = 0   
    end function

    function get_amrentropy(outputvalue)
        IMPLICIT NONE
        integer :: get_amrentropy
        logical :: outputvalue
        outputvalue = amrentropy
        get_amrentropy = 0   
    end function

    function set_divbfix(inputvalue)
        IMPLICIT NONE
        integer :: set_divbfix
        logical :: inputvalue
        divbfix = inputvalue
        set_divbfix = 0   
    end function

    function get_divbfix(outputvalue)
        IMPLICIT NONE
        integer :: get_divbfix
        logical :: outputvalue
        outputvalue = divbfix
        get_divbfix = 0   
    end function

    function set_divbwave(inputvalue)
        IMPLICIT NONE
        integer :: set_divbwave
        logical :: inputvalue
        divbwave = inputvalue
        set_divbwave = 0   
    end function

    function get_divbwave(outputvalue)
        IMPLICIT NONE
        integer :: get_divbwave
        logical :: outputvalue
        outputvalue = divbwave
        get_divbwave = 0   
    end function

    function set_compactres(inputvalue)
        IMPLICIT NONE
        integer :: set_compactres
        logical :: inputvalue
        compactres = inputvalue
        set_compactres = 0   
    end function

    function get_compactres(outputvalue)
        IMPLICIT NONE
        integer :: get_compactres
        logical :: outputvalue
        outputvalue = compactres
        get_compactres = 0   
    end function

    function set_bnormlf(inputvalue)
        IMPLICIT NONE
        integer :: set_bnormlf
        logical :: inputvalue
        bnormlf = inputvalue
        set_bnormlf = 0   
    end function

    function get_bnormlf(outputvalue)
        IMPLICIT NONE
        integer :: get_bnormlf
        logical :: outputvalue
        outputvalue = bnormlf
        get_bnormlf = 0   
    end function

    function set_strictnr(inputvalue)
        IMPLICIT NONE
        integer :: set_strictnr
        logical :: inputvalue
        strictnr = inputvalue
        set_strictnr = 0   
    end function

    function get_strictnr(outputvalue)
        IMPLICIT NONE
        integer :: get_strictnr
        logical :: outputvalue
        outputvalue = strictnr
        get_strictnr = 0   
    end function

    function set_strictsmall(inputvalue)
        IMPLICIT NONE
        integer :: set_strictsmall
        logical :: inputvalue
        strictsmall = inputvalue
        set_strictsmall = 0   
    end function

    function get_strictsmall(outputvalue)
        IMPLICIT NONE
        integer :: get_strictsmall
        logical :: outputvalue
        outputvalue = strictsmall
        get_strictsmall = 0   
    end function

    function set_strictzero(inputvalue)
        IMPLICIT NONE
        integer :: set_strictzero
        logical :: inputvalue
        strictzero = inputvalue
        set_strictzero = 0   
    end function

    function get_strictzero(outputvalue)
        IMPLICIT NONE
        integer :: get_strictzero
        logical :: outputvalue
        outputvalue = strictzero
        get_strictzero = 0   
    end function

    function set_strictgetaux(inputvalue)
        IMPLICIT NONE
        integer :: set_strictgetaux
        logical :: inputvalue
        strictgetaux = inputvalue
        set_strictgetaux = 0   
    end function

    function get_strictgetaux(outputvalue)
        IMPLICIT NONE
        integer :: get_strictgetaux
        logical :: outputvalue
        outputvalue = strictgetaux
        get_strictgetaux = 0   
    end function

    function set_usecovariant(inputvalue)
        IMPLICIT NONE
        integer :: set_usecovariant
        logical :: inputvalue
        usecovariant = inputvalue
        set_usecovariant = 0   
    end function

    function get_usecovariant(outputvalue)
        IMPLICIT NONE
        integer :: get_usecovariant
        logical :: outputvalue
        outputvalue = usecovariant
        get_usecovariant = 0   
    end function

    function set_nocartesian(inputvalue)
        IMPLICIT NONE
        integer :: set_nocartesian
        logical :: inputvalue
        nocartesian = inputvalue
        set_nocartesian = 0   
    end function

    function get_nocartesian(outputvalue)
        IMPLICIT NONE
        integer :: get_nocartesian
        logical :: outputvalue
        outputvalue = nocartesian
        get_nocartesian = 0   
    end function

    function set_tfix(inputvalue)
        IMPLICIT NONE
        integer :: set_tfix
        logical :: inputvalue
        tfix = inputvalue
        set_tfix = 0   
    end function

    function get_tfix(outputvalue)
        IMPLICIT NONE
        integer :: get_tfix
        logical :: outputvalue
        outputvalue = tfix
        get_tfix = 0   
    end function

    function set_convert(inputvalue)
        IMPLICIT NONE
        integer :: set_convert
        logical :: inputvalue
        convert = inputvalue
        set_convert = 0   
    end function

    function get_convert(outputvalue)
        IMPLICIT NONE
        integer :: get_convert
        logical :: outputvalue
        outputvalue = convert
        get_convert = 0   
    end function

    function set_saveprim(inputvalue)
        IMPLICIT NONE
        integer :: set_saveprim
        logical :: inputvalue
        saveprim = inputvalue
        set_saveprim = 0   
    end function

    function get_saveprim(outputvalue)
        IMPLICIT NONE
        integer :: get_saveprim
        logical :: outputvalue
        outputvalue = saveprim
        get_saveprim = 0   
    end function

    function set_uselimiter(inputvalue)
        IMPLICIT NONE
        integer :: set_uselimiter
        logical :: inputvalue
        uselimiter = inputvalue
        set_uselimiter = 0   
    end function

    function get_uselimiter(outputvalue)
        IMPLICIT NONE
        integer :: get_uselimiter
        logical :: outputvalue
        outputvalue = uselimiter
        get_uselimiter = 0   
    end function
    
    
    function set_parameters_filename(path)
        IMPLICIT NONE
        integer :: set_parameters_filename
        character(len=1024), intent(in) :: path
        logical :: file_exists
        INQUIRE(FILE= TRIM(path), EXIST=file_exists) 
        if (file_exists) THEN
            parameters_filename = TRIM(path)
            set_parameters_filename = 0 
        else
            set_parameters_filename = -1 
        end if
    end function
    
    function get_parameters_filename(path)
        IMPLICIT NONE
        integer :: get_parameters_filename
        character(len=1024), intent(out) :: path
        path = parameters_filename 
        get_parameters_filename = 0   
    end function
    
    function get_local_index_of_grid(global_index_of_grid)
        use mod_forest
        
        IMPLICIT NONE
        
        integer :: get_local_index_of_grid
        integer :: iigrid, number_of_grids_before, ipe
        integer, intent(in) :: global_index_of_grid
        
        get_local_index_of_grid = 0
        
        if(global_index_of_grid < Morton_start(mype)) then
            get_local_index_of_grid = 0
            return
        end if 
        if(global_index_of_grid > Morton_stop(mype)) then
            get_local_index_of_grid = 0
            return
        end if 
        
        get_local_index_of_grid = sfc_to_igrid(global_index_of_grid)
        
        !number_of_grids_before = 0
        !do ipe = 0, mype-1
        !    do iigrid = 1, ngridshi
        !        if(igrid_inuse(iigrid, ipe)) then
        !           number_of_grids_before = number_of_grids_before + 1
        !         end if
        !    end do
        !end do
        
        !if ( global_index_of_grid .LE. number_of_grids_before ) then
        !    get_local_index_of_grid = 0
        !else if ( global_index_of_grid .GT. number_of_grids_before + igridstail ) then
        !    get_local_index_of_grid = 0
        !else
        !    get_local_index_of_grid = igrids(global_index_of_grid - number_of_grids_before)
        !end if
        
        !print *, "get_local_index_of_grid", global_index_of_grid, get_local_index_of_grid,&
        !&number_of_grids_before,number_of_grids_before+igridstail
                        
    end function
    
    function is_index_valid(i, j, k)
        IMPLICIT NONE
        logical :: is_index_valid
        integer, intent(in) :: i,j,k
        
        is_index_valid = .TRUE.
        
        is_index_valid = is_index_valid .AND. (i .GE. 0)
        is_index_valid = is_index_valid .AND. (j .GE. 0)
        is_index_valid = is_index_valid .AND. (k .GE. 0)
        
        is_index_valid = is_index_valid .AND. (i .LE. (ixMhi1 - ixMlo1))
{^NOONED
        is_index_valid = is_index_valid .AND. (j .LE. (ixMhi2 - ixMlo2))
}
{^IFTHREED
        is_index_valid = is_index_valid .AND. (k .LE. (ixMhi3 - ixMlo3))
}
    end function
    
    
    function get_position_of_index(i, j, k, index_of_grid, x, y, z, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_position_of_index
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        double precision, intent(out), dimension(n) :: x, y, z
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            x(index) = 0.0
            y(index) = 0.0
            z(index) = 0.0
            
            if (local_index_of_grid .EQ. 0) then
                x(index) = 0.0
                y(index) = 0.0
                z(index) = 0.0
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                real1 = ixMlo1 + i(index)
{^NOONED
                real2 = ixMlo2 + j(index)
}
{^IFTHREED
                real3 = ixMlo3 + k(index)
}
                x(index) = px(local_index_of_grid)%x(real^D,1)
{^NOONED
                y(index) = px(local_index_of_grid)%x(real^D,2)
}
                
{^IFTHREED
                z(index) = px(local_index_of_grid)%x(real^D,3)
}
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(x,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(y,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(z,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, x, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, y, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, z, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_position_of_index=0
    end function
    
    function get_index_of_position( x, y, z, index_of_grid, i, j, k, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_index_of_position
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: index_of_grid
        double precision, intent(out), dimension(n) :: i, j, k
        double precision, intent(in), dimension(n) :: x, y, z
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            if (local_index_of_grid .EQ. 0) then
                i(index) = 0
                j(index) = 0
                k(index) = 0
            else
                
                ! TODO TODO TODO
                i(index) = 0
                j(index) = 0
                k(index) = 0
                
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(i,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(j,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(k,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, i, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, j, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, k, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_index_of_position=-2
    end function
    
    function get_mesh_size(nx, ny, nz, index_of_grid)
        IMPLICIT NONE

        integer :: get_mesh_size
        integer, intent(out) :: nx, ny, nz
        integer, intent(in) :: index_of_grid
        
        ny = 1
        nz = 1
        
        nx = ixMhi1 - ixMlo1 + 1
{^NOONED
        ny = ixMhi2 - ixMlo2 + 1
}
{^IFTHREED
        nz = ixMhi3 - ixMlo3 + 1
}
        get_mesh_size = 0
    end function

    function get_number_of_grids(number_of_grids)
        use mod_forest
        
        IMPLICIT NONE
        integer :: ip, igrid
        integer, intent(out) :: number_of_grids
        integer :: get_number_of_grids
        get_number_of_grids = 0
        number_of_grids = 0
        do ip = 0, npe - 1
            do igrid = 1, ngridshi
                if (igrid_inuse(igrid, ip)) then
                    number_of_grids = number_of_grids + 1
                end if
            end do
        end do
    end function
    
    function get_level_of_grid(level, index_of_grid)
        use mod_forest
        
        IMPLICIT NONE
        integer, intent(out) :: level
        integer, intent(in) :: index_of_grid
        integer :: get_level_of_grid
        integer :: local_index_of_grid
        
        if(index_of_grid > ngridshi) then
            level = 0
            get_level_of_grid = -1
            return
        end if
        
        local_index_of_grid = get_local_index_of_grid(index_of_grid)
        if(local_index_of_grid > 0) then
            level = node(plevel_, local_index_of_grid)
        else
            level = 0
        end if
        
        if(mype .GT. 0) then
            call MPI_Reduce(level,  0, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, level, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_level_of_grid = 0
    end function
    
    function get_cell_size_of_grid(dx1, dx2, dx3, index_of_grid)
        use mod_forest
        
        IMPLICIT NONE
        double precision, intent(out) :: dx1, dx2, dx3
        integer, intent(in) :: index_of_grid
        integer :: level
        integer :: get_cell_size_of_grid
        integer :: local_index_of_grid
        
        dx1 = 0
        dx2 = 0
        dx3 = 0
        
        if(index_of_grid > ngridshi) then
            get_cell_size_of_grid = -1
            return
        end if
        
        local_index_of_grid = get_local_index_of_grid(index_of_grid)
        
        if(local_index_of_grid > 0) then
            level = node(plevel_, local_index_of_grid)
        
            if(level > 0) then
                dx1 = 0.0
                dx3 = 0.0
                dx1 = dx(1,level)
{^NOONED
                dx2 = dx(2,level)
}
{^IFTHREED
                dx3 = dx(3,level)
}
            end if
        end if
        
        if(mype .GT. 0) then
            call MPI_Reduce(dx1,  0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(dx2,  0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(dx3,  0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, dx1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, dx2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, dx3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_cell_size_of_grid = 0
        
    end function

    
    function setup_mesh(nmeshx, nmeshy, nmeshz, xlength, ylength, zlength)
        IMPLICIT NONE
        integer, intent(in) :: nmeshx, nmeshy, nmeshz
        double precision, intent(in) :: xlength, ylength, zlength
        integer :: setup_mesh
                
        xprobmin1 = 0.0
        xprobmax1 = xlength
{^NOONED
        xprobmin2 = 0.0
        xprobmax2 = ylength
}
{^IFTHREED
        xprobmin3 = 0.0
        xprobmax3 = zlength
}
        
        if(nmeshx < 1 .OR.  mod(nmeshx,2)/=0) then
            setup_mesh = -1
            return
        end if

{^NOONED
        if(nmeshy < 1 .OR.  mod(nmeshy,2)/=0) then
            setup_mesh = -1
            return
        end if
}
{^IFTHREED
        if(nmeshz < 1 .OR.  mod(nmeshz,2)/=0) then
            setup_mesh = -1
            return
        end if
}
        
        dx(1,1)=(xprobmax1-xprobmin1)/dble(nmeshx)
        
{^NOONED
        dx(2,1)=(xprobmax2-xprobmin2)/dble(nmeshy)
}

{^IFTHREED
        dx(3,1)=(xprobmax3-xprobmin3)/dble(nmeshz)
}
        
        setup_mesh = 0
    end function
    
    function get_grid_state(i, j, k, index_of_grid, rho,  m1, m2, m3, en, n)
    
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_state
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: rho
        double precision, intent(out), dimension(n) :: m1, m2, m3
        double precision, intent(out), dimension(n) :: en
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid .EQ. 0) then
            
                rho(index) = 0.0
                m1(index) = 0.0
                m2(index) = 0.0
                m3(index) = 0.0
                en(index) = 0.0
                
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                real1 = ixMlo1 + i(index)
{^NOONED
                real2 = ixMlo2 + j(index)
}
{^IFTHREED
                real3 = ixMlo3 + k(index)
}
                rho(index) = pw(local_index_of_grid)%w(real^D,rho_)
                m1(index) = pw(local_index_of_grid)%w(real^D,m1_)

{^NOONED
                m2(index) = pw(local_index_of_grid)%w(real^D,m2_)
}
{^IFTHREED
                m3(index) = pw(local_index_of_grid)%w(real^D,m3_)
}
                en(index) = pw(local_index_of_grid)%w(real^D,e_)
            end if
            
        end do
                
        
        if(mype .GT. 0) then
            call MPI_Reduce(rho,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m1,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m2,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m3,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(en,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, rho, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m1, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m2, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m3, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, en, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        
        get_grid_state = 0
    end function
    
    
    function get_grid_density(i, j, k, index_of_grid, rho, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_density
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: rho
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid .EQ. 0) then
            
                rho(index) = 0.0
                
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                                
                real1 = ixMlo1 + i(index)
{^NOONED
                real2 = ixMlo2 + j(index)
}
{^IFTHREED
                real3 = ixMlo3 + k(index)
}              
                rho(index) = pw(local_index_of_grid)%w(real^D,rho_)
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(rho,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, rho, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        get_grid_density = 0
    end function
    
    function get_grid_momentum_density(i, j, k, index_of_grid, m1, m2, m3, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_momentum_density
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: m1, m2, m3
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            m1(index) = 0.0
            m2(index) = 0.0
            m3(index) = 0.0
            if (local_index_of_grid .EQ. 0) then
                m1(index) = 0.0
                m2(index) = 0.0
                m3(index) = 0.0
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                
                real1 = ixMlo1 + i(index)
{^NOONED
                real2 = ixMlo2 + j(index)
}
{^IFTHREED
                real3 = ixMlo3 + k(index)
}              
                
                m1(index) = pw(local_index_of_grid)%w(real^D,m1_)
{^NOONED
                m2(index) = pw(local_index_of_grid)%w(real^D,m2_)
}
{^IFTHREED
                m3(index) = pw(local_index_of_grid)%w(real^D,m3_)
} 
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(m1,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m2,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m3,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, m1, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m2, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m3, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_grid_momentum_density = 0
    end function
    
    function get_grid_energy_density(i, j, k, index_of_grid, en, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_energy_density
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: en
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid .EQ. 0) then
            
                en(index) = 0.0
                
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                real1 = ixMlo1 + i(index)
{^NOONED
                real2 = ixMlo2 + j(index)
}
{^IFTHREED
                real3 = ixMlo3 + k(index)
}              
                
                en(index) = pw(local_index_of_grid)%w(real^D,e_)
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(en,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, en, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_grid_energy_density = 0
    end function
    
    
    function set_grid_state(i, j, k, rho, m1, m2, m3, en, index_of_grid, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_state
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: rho
        double precision, intent(in), dimension(n) :: m1, m2, m3
        double precision, intent(in), dimension(n) :: en
        
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                
                    real1 = ixMlo1 + i(index)
{^NOONED
                    real2 = ixMlo2 + j(index)
}
{^IFTHREED
                    real3 = ixMlo3 + k(index)
}              
                    
                    pw(local_index_of_grid)%w(real^D,rho_) = rho(index)
                    pw(local_index_of_grid)%w(real^D,m1_) = m1(index)
{^NOONED
                    pw(local_index_of_grid)%w(real^D,m2_) = m2(index)
}
{^IFTHREED
                    pw(local_index_of_grid)%w(real^D,m3_) = m3(index)
}              
                    pw(local_index_of_grid)%w(real^D,e_) = en(index)
                
                end if
            end if
            
        end do
        
        set_grid_state = 0
    end function
    
    function set_grid_density(i, j, k, rho, index_of_grid, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_density
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: rho
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                
                    real1 = ixMlo1 + i(index)
{^NOONED
                    real2 = ixMlo2 + j(index)
}
{^IFTHREED
                    real3 = ixMlo3 + k(index)
}              
                    
                    pw(local_index_of_grid)%w(real^D,rho_) = rho(index)
                
                end if
            end if
            
        end do
        
        set_grid_density = 0
    end function

    function set_grid_energy_density(i, j, k, en, index_of_grid, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_energy_density
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: en
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                    real1 = ixMlo1 + i(index)
{^NOONED
                    real2 = ixMlo2 + j(index)
}
{^IFTHREED
                    real3 = ixMlo3 + k(index)
}              
                    
                    pw(local_index_of_grid)%w(real^D,e_) = en(index)
                
                end if
            end if
            
        end do
        
        set_grid_energy_density = 0
    end function
    
    function set_grid_momentum_density(i, j, k, m1, m2, m3, index_of_grid, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_momentum_density
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: m1, m2, m3
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                    real1 = ixMlo1 + i(index)
{^NOONED
                    real2 = ixMlo2 + j(index)
}
{^IFTHREED
                    real3 = ixMlo3 + k(index)
}              
                    
                    pw(local_index_of_grid)%w(real^D,m1_) = m1(index)
{^NOONED
                    pw(local_index_of_grid)%w(real^D,m2_) = m2(index)
}
{^IFTHREED
                    pw(local_index_of_grid)%w(real^D,m3_) = m3(index)
}              
                
                end if
            end if
            
        end do
        
        set_grid_momentum_density = 0
    end function
    
    
    function evolve_model(tend)
        use StoppingConditions
        
        IMPLICIT NONE
        
        integer :: evolve_model
        double precision :: tend

        integer :: is_number_of_steps_detection_enabled, is_timeout_detection_enabled
        integer :: max_number_of_steps
        integer :: level, ifile
        integer :: error
        integer :: localsteps, stopping_index
        double precision :: time_in, timeio0, timeio_tot, timegr0, timegr_tot,&
            timeloop, timeloop0
        double precision :: loop_timeout, time_spent
       
        
        error = reset_stopping_conditions()
        error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, &
                    is_number_of_steps_detection_enabled)
        error = get_stopping_condition_number_of_steps_parameter(max_number_of_steps)
        error = is_stopping_condition_enabled(TIMEOUT_DETECTION, &
                    is_timeout_detection_enabled)
        error = get_stopping_condition_timeout_parameter(loop_timeout)


        time_accurate = .TRUE.
        
        time_in=MPI_WTIME()
        timeio_tot=zero
        timegr_tot=zero
        localsteps=0
        itmin=it
        
        call getbc(t,ixG^LL,pw,pwCoarse,pgeo,&
               pgeoCoarse,.false.)

        !  ------ start of integration loop. ------------------
        
        call MPI_BARRIER(icomm,ierrmpi)
        timeloop0=MPI_WTIME()
        if (mype==0) then
           write(*,'(a,f12.3,a)')'BCs before Advance took    : ',timeloop0&
              -time_in,' sec'
        end if
        
        tmax = tend
        
        time_evol : do
           ! exit time loop criteria
           if (it>=itmax) exit time_evol
           if (time_accurate .and. t>=tmax) exit time_evol
           call setdt
           
           if(fixprocess) call process(it,t)
           if(mype==0.and..false.) print *,'done setdt for it=',it

           timeio0=MPI_WTIME()
           !do ifile=nfile,1,-1
           !   if(timetosave(ifile)) call saveamrfile(ifile)
           !end do
           timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

           if(mype==0.and..false.) print *,'enters advance for it=',it
           call advance(it)
           if(mype==0.and..false.) print *,'done advance for it=',it
            
           if((.not.time_accurate).or.(residmin>smalldouble)) then
              call getresidual(it)
           endif 
           if (time_accurate .and.  (residual<residmin .or. residual&
              >residmax)) exit time_evol
           if (.not.time_accurate .and. (residual<residmin .or. residual&
              >residmax)) exit time_evol

           ! resetting of tree BEFORE IO and setdt
           timegr0=MPI_WTIME()
           !if (mxnest>1 .and. .not.(fixgrid(0))) call resettree
           if (mxnest>1) call resettree
           timegr_tot=timegr_tot+(MPI_WTIME()-timegr0)

           it = it + 1
           if (time_accurate) t = t + dt
           if(addmpibarrier) call MPI_BARRIER(icomm,ierrmpi)
           localsteps = localsteps + 1
           if (is_number_of_steps_detection_enabled.NE.0 .AND. &
              localsteps.GE.max_number_of_steps) then
              stopping_index = next_index_for_stopping_condition()
              error = set_stopping_condition_info(stopping_index,&
                NUMBER_OF_STEPS_DETECTION)
              exit time_evol
           end if    
           if (is_timeout_detection_enabled.NE.0) then
              
              if(.NOT. addmpibarrier) then
                !need to add the barrier
                !if not already done before, to make sure all processes will
                !stop
                call MPI_BARRIER(icomm,ierrmpi)
              end if
              time_spent = MPI_WTIME()-timeloop0
              if (time_spent .GE. loop_timeout) then
                stopping_index = next_index_for_stopping_condition()
                error = set_stopping_condition_info(stopping_index,&
                    TIMEOUT_DETECTION)
                exit time_evol
              end if
           end if                
           if(it>900000) it = slowsteps+10
           
        end do time_evol

        timeloop=MPI_WTIME()-timeloop0

        if (mype==0) then
           write(*,'(a,f12.3,a)')'Total timeloop took        : ',timeloop,' sec'
           write(*,'(a,f12.3,a)')'Time spent on Regrid+Update: ',timegr_tot,' sec'
           write(*,'(a,f9.2)')   '                 Percentage: ',100.0*timegr_tot&
              /timeloop
           write(*,'(a,f12.3,a)')'Time spent on IO in loop   : ',timeio_tot,' sec'
           write(*,'(a,f9.2)')   '                 Percentage: ',100.0*timeio_tot&
              /timeloop
           write(*,'(a,f12.3,a)')'Time spent on BC           : ',time_bc,' sec'
        end if

        timeio0=MPI_WTIME()
        !do ifile=nfile,1,-1
        !   if(itsavelast(ifile)<it)call saveamrfile(ifile)
        !enddo
        if (mype==0) call MPI_FILE_CLOSE(log_fh,ierrmpi)
        timeio_tot=timeio_tot+(MPI_WTIME()-timeio0)

        if (mype==0) then
           write(*,'(a,f12.3,a)')'Total time spent on IO     : ',timeio_tot,' sec'
           write(*,'(a,f12.3,a)')'Total timeintegration took : ',MPI_WTIME()&
              -time_in,' sec'
        end if
    
        evolve_model = 0
    end function

    function get_time(tnow)
        IMPLICIT NONE
        
        integer :: get_time
        double precision :: tnow

        tnow = t
        
        get_time = 0
    end function
    
    function set_boundary(lowx,highx,lowy,highy,lowz,highz)
        IMPLICIT NONE
        
        integer :: set_boundary
        integer :: idim
        character(len=25), intent(in) :: lowx,highx,lowy,highy,lowz,highz
        
        typeB(1:nw,1) = lowx
        typeB(1:nw,2) = highx
{^NOONED
        typeB(1:nw,3) = lowy
        typeB(1:nw,4) = highy
}
{^IFTHREED
        typeB(1:nw,5) = lowz
        typeB(1:nw,6) = highz
}
        do idim=1,ndim
            periodB(idim)=(typeB(1,2*idim-1)=="periodic")
            if ((.not.periodB(idim).and.(any(typeB(:,2*idim-1:2*idim)&
                =="periodic"))) .or.(periodB(idim).and.(any(typeB(:,2*idim&
                -1:2*idim)/="periodic")))) then
                write(unitterm,*)"Each dimension should either have all",&
                " or no variables periodic"
                set_boundary = -2
                return
            end if
        end do
        
        set_boundary = 0
    end function  


    function get_acceleration_grid_size(nx, ny, nz)
        IMPLICIT NONE

        integer :: get_acceleration_grid_size
        integer, intent(out) :: nx, ny, nz
        
        nx = naccel1
        ny = naccel2
        nz = naccel3
    
        get_acceleration_grid_size = 0
    end function
    
    
    function set_acceleration_grid_acceleration(i, j, k, a1, a2, a3, n)
        IMPLICIT NONE
              
        integer :: index
        integer :: set_acceleration_grid_acceleration
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k
        integer :: real1, real2, real3
        
        double precision, intent(in), dimension(n) :: a1, a2, a3
        
        do index = 1,n
            
 !           if (local_index_of_grid /= 0) then
                !if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                !
                !else
                    
                real1 = i(index)
{^NOONED
                real2 = j(index)
}
{^IFTHREED
                real3 = k(index)
}              
                accel(real^D,1) = a1(index)
{^NOONED
                accel(real^D,2) = a2(index)
}
{^IFTHREED
                accel(real^D,3) = a3(index)
}
                    
                
                !end if
!            end if
            
        end do
        
        set_acceleration_grid_acceleration = 0
    end function
    
    function get_acceleration_grid_acceleration(i, j, k, a1, a2, a3, n)
        IMPLICIT NONE
       
        integer :: index
        integer :: get_acceleration_grid_acceleration
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k
        
        double precision, intent(out), dimension(n) :: a1, a2, a3
        integer :: real1, real2, real3
        
        do index = 1,n
            
 !           if (local_index_of_grid /= 0) then
                !if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                !
                !else
                a2(index) = 0.0
                a3(index) = 0.0
                
                real1 = i(index)
{^NOONED
                real2 = j(index)
}
{^IFTHREED
                real3 = k(index)
}              
                a1(index) = accel(real^D,1)
{^NOONED
                a2(index) = accel(real^D,2)
}
{^IFTHREED
                a3(index) = accel(real^D,3)
}
                    
                
                !end if
 !           end if
            
        end do
        
        get_acceleration_grid_acceleration = 0
    end function
    
    function get_acceleration_grid_position_of_index(i, j, k, x, y, z, n)
        IMPLICIT NONE
        
        integer :: index
        integer :: get_acceleration_grid_position_of_index
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k
        
        double precision, intent(out), dimension(n) :: x, y, z
        
        do index = 1,n
            
 !           if (local_index_of_grid /= 0) then
                !if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                !
                !else
                    
                x(index) = xaccel1(i(index))
                y(index) = xaccel2(j(index))
                z(index) = xaccel3(k(index))
                    
                !end if
 !           end if
            
        end do
        
        get_acceleration_grid_position_of_index = 0
    end function
    
    function validation_error(error_description)
        integer :: validation_error
        CHARACTER(LEN=*), INTENT(IN) :: error_description
        PRINT *, error_description
        error_string = error_description
        validation_error = -1
    end function
    
    function check_method_parameters()
        IMPLICIT NONE
        integer :: check_method_parameters
        integer :: level
        
        check_method_parameters = 0
        
        if(compactres)then
         if(typeaxial/='slab') check_method_parameters = validation_error("compactres for MHD only in cartesian case") ; return
        endif

        if(TRIM(wnames)=='default') check_method_parameters = validation_error("Provide wnames and restart code") ; return

        do level=1,nlevelshi
           if(typefull1(level)=='tvdlf1'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" tvdlf1 is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='hll1'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" hll1 is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='hllc1'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" hllc1 is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='hllcd1'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" hllcd1 is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='tvdmu1'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" tvdmu1 is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='tvd'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" tvd is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='tvd1'.and.typeadvance=='twostep') &
              check_method_parameters = validation_error(" tvd1 is onestep method, reset typeadvance=onestep!") ; return
           if(typefull1(level)=='tvd'.or.typefull1(level)=='tvd1')then 
              if(mype==0.and.(.not.dimsplit)) write(unitterm,*) &
                 'Warning: setting dimsplit=T for tvd, as used for level=',level
              dimsplit=.true.
           endif

           if (typepred1(level)=='default') then
              select case (typefull1(level))
              case ('cd')
                 typepred1(level)='cd'
              case ('tvdlf','tvdmu')
                 typepred1(level)='hancock'
              case ('hll')
                 typepred1(level)='hll'
              case ('hllc')
                 typepred1(level)='hllc'
              case ('hllcd')
                 typepred1(level)='hllcd'
              case ('tvdlf1','tvdmu1','tvd1','tvd','hll1','hllc1', &
                    'hllcd1','nul','source')
                 typepred1(level)='nul'
              case default
                 check_method_parameters = validation_error("No default predictor for this full step") ; return
              end select
           end if
        end do

        select case (typeadvance)
            case ("onestep")
               nstep=1
            case ("twostep")
               nstep=2
            case ("fourstep")
               nstep=4
            case default
               check_method_parameters = validation_error("Unknown typeadvance") ; return
        end select


        do level=1,nlevelshi
            if (typelow1(level)=='default') then
                select case (typefull1(level))
                    case ('cd')
                      typelow1(level)='cd'
                    case ('hancock')
                      typelow1(level)='hancock'
                    case ('tvdlf','tvdlf1','tvdmu','tvdmu1','tvd1','tvd')
                      typelow1(level)='tvdlf1'
                    case ('hll','hll1')
                      typelow1(level)='hll1'
                    case ('hllc','hllc1')
                      typelow1(level)='hllc1'
                    case ('hllcd','hllcd1')
                      typelow1(level)='hllcd1'
                    case ('nul')
                      typelow1(level)='nul'
                    case ('source')
                      typelow1(level)='source'
                    case default
                      check_method_parameters = validation_error("No default typelow for this full step") ; return
                end select
            end if
        enddo

        ! Harmonize the parameters for dimensional splitting and source splitting
        if(typedimsplit   =='default'.and.     dimsplit)   typedimsplit='xyyx'
        if(typedimsplit   =='default'.and..not.dimsplit)   typedimsplit='unsplit'
        if(typesourcesplit=='default'.and.     sourcesplit)typesourcesplit='sfs'
        if(typesourcesplit=='default'.and..not.sourcesplit)typesourcesplit='unsplit'
        dimsplit   = typedimsplit   /='unsplit'
        sourcesplit= typesourcesplit/='unsplit'
        if(sourcesplit)sourceunsplit=.false.

        if (typeaxial=="slab") then
           slab=.true.
        else
           slab=.false.
        end if

        if (typeaxial=='spherical') then
           if (dimsplit) then
              if(mype==0)print *,'Warning: spherical symmetry needs dimsplit=F, resetting'
              dimsplit=.false.
           end if
        end if

        if (ndim==1) dimsplit=.false.
        if (.not.dimsplit.and.ndim>1) then
           select case (typeadvance)
           case ("fourstep")
              ! Runge-Kutta needs predictor
              typelimited="predictor"
              if(mype==0)print *,'typelimited to predictor for RK4'
           case ("twostep")
              ! non-conservative Hancock predictor needs the previous typelimited
              !do level=1,nlevelshi
                !!if(typepred1(level)=='hancock')then
                !!!  typelimited="previous"
                !!  if(mype==0)print *,'Warning: typelimited to previous for hancock on level=',level
                !!end if
              !end do
           end select
        end if

        if((divbfix.and.(typephys=='mhd'.or.typephys=='srmhd'.or.typephys=='srmhdeos')).and.&
           (.not.sourcesplit).and.(.not.sourceunsplit))&
           check_method_parameters = validation_error("divbfix=T requires unsplitsource=T or splitsource=T !") ; return
        if((divbdiff>zero.and.(typephys=='mhd'.or.typephys=='srmhd'.or.typephys=='srmhdeos')).and.&
           (.not.sourcesplit).and.(.not.sourceunsplit))&
           check_method_parameters = validation_error("divbdiff>0 requires unsplitsource=T or splitsource=T !") ; return

        if (B0field) then
           if (.not.typephys=='mhd') B0field=.false.
           if(mype==0)print *,'B0+B1 split only for MHD'
        end if

        if (any(typelimiter1(1:nlevelshi)== 'ppm')&
            .and.(flatsh.and.typephys=='rho')) then
            check_method_parameters = validation_error(" PPM with flatsh=.true. can not be used with typephys='rho'!") ; return
        end if
        if (any(typelimiter1(1:nlevelshi)== 'ppm')&
            .and.(flatsh.and.typephys=='hdadiab')) then
             check_method_parameters = validation_error(" PPM with flatsh=.true. can not be used with typephys='hdadiab'!") ; return
        end if
        if (any(typelimiter1(1:nlevelshi)== 'ppm')&
            .and.(flatcd.and.typephys=='hdadiab')) then
             check_method_parameters = validation_error(" PPM with flatcd=.true. can not be used with typephys='hdadiab'!") ; return
        end if
        if (any(typelimiter1(1:nlevelshi)== 'ppm')&
            .and.(flatsh.and..not.useprimitive)) then
             check_method_parameters = validation_error(" PPM with flatsh=.true. needs useprimitive=T!") ; return
        end if
        if (any(typelimiter1(1:nlevelshi)== 'ppm')&
            .and.(flatcd.and..not.useprimitive)) then
             check_method_parameters = validation_error(" PPM with flatcd=.true. needs useprimitive=T!") ; return
        end if

    end function
    
    function check_stop_parameters()
        IMPLICIT NONE
        integer :: check_stop_parameters
        
        if(residmin>=zero) then
            if(residmin<=smalldouble) then
                check_stop_parameters = validation_error("Provide value for residual above smalldouble")
                return
            end if
        end if
        
    end function
    
    
    function get_grid_acceleration(i, j, k, index_of_grid, a1, a2, a3, n)
    
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_acceleration
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid

        double precision, intent(out), dimension(n) :: a1, a2, a3
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            a1(index) = 0.0
            a2(index) = 0.0
            a3(index) = 0.0
            
            if (local_index_of_grid .EQ. 0) then
            
                a1(index) = 0.0
                a2(index) = 0.0
                a3(index) = 0.0
                
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                real1 = ixMlo1 + i(index)
{^NOONED
                real2 = ixMlo2 + j(index)
}
{^IFTHREED
                real3 = ixMlo3 + k(index)
}
                a1(index) = pw(local_index_of_grid)%w(real^D,e_+1)

{^NOONED
                a2(index) = pw(local_index_of_grid)%w(real^D,e_+2)
}
{^IFTHREED
                a3(index) = pw(local_index_of_grid)%w(real^D,e_+3)
}
            end if
            
        end do
                
        
        if(mype .GT. 0) then
            call MPI_Reduce(a1,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(a2,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(a3,  0, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, a1, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, a2, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, a3, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        
        get_grid_acceleration = 0
    end function
    
    
    function set_grid_acceleration(i, j, k, a1, a2, a3, index_of_grid, n)
        IMPLICIT NONE
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_acceleration
        integer :: real1, real2, real3
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: a1, a2, a3
        
        
        local_index_of_grid = 0
        previous_index_of_grid = -1
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
                previous_index_of_grid = index_of_grid(index)
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                
                    real1 = ixMlo1 + i(index)
{^NOONED
                    real2 = ixMlo2 + j(index)
}
{^IFTHREED
                    real3 = ixMlo3 + k(index)
}              
                    
                    pw(local_index_of_grid)%w(real^D,e_+1) = a1(index)
{^NOONED
                    pw(local_index_of_grid)%w(real^D,e_+2) = a2(index)
}
{^IFTHREED
                    pw(local_index_of_grid)%w(real^D,e_+3) = a3(index)
}              
                
                end if
            end if
            
        end do
        
        set_grid_acceleration = 0
    end function
    
    function get_hydro_state_at_point(x1, x2, x3, vx, vy, vz, rho, rhovx, rhovy, rhovz, rhoe,npoints) result(ret)
        IMPLICIT NONE
        integer :: ret, ilist
        integer, intent(in) :: npoints
        integer :: iigrid, igrid, i^D, idims
        double precision :: weighing_factors(2^D&)
        double precision :: dx^D, xleft, xright
        double precision, intent(in) :: x1(npoints), x2(npoints), x3(npoints)
        double precision, intent(in) :: vx(npoints), vy(npoints), vz(npoints)
        double precision, intent(out):: rho(npoints), rhovx(npoints), rhovy(npoints)
        double precision, intent(out):: rhovz(npoints), rhoe(npoints)
        rho=0.0
        rhovx=0.0
        rhovy=0.0
        rhovz=0.0
        rhoe=0.0
        weighing_factors = 1.0
        do ilist = 1,npoints
            weighing_factors = 1.0
            do iigrid=1,igridstail; 
                igrid=igrids(iigrid)
                if({(rnode(rpxmin^D_,igrid).LE.x^D(ilist).and.&
                    rnode(rpxmax^D_,igrid).GT.x^D(ilist))|.and.})then
                    
                    dx^D=rnode(rpdx^D_,igrid);

                    i^D=FLOOR(((x^D(ilist)-rnode(rpxmin^D_,igrid) - (half * dx^D))/dx^D))+dixB+1;
                    do idims=1,ndim
                    select case(idims)
                        {case(^D)
                        xleft = px(igrid)%x(i^DD,^D)
                        xright = px(igrid)%x(i^D+1^D%i^DD,^D)
                        weighing_factors(1^D%1:2^D&) = weighing_factors(1^D%1:2^D&) * ((xright - x^D(ilist)) / dx^D)
                        weighing_factors(2^D%1:2^D&) = weighing_factors(2^D%1:2^D&) * ((x^D(ilist)  - xleft) / dx^D)
                        
                       }
                    end select
                    end do
                
                    rho(ilist) = sum(pw(igrid)%w(i^D:i^D+1, rho_) * weighing_factors)
                    rhoe(ilist) = sum(pw(igrid)%w(i^D:i^D+1, e_) * weighing_factors)
                    rhovx(ilist) = sum(pw(igrid)%w(i^D:i^D+1, m1_) * weighing_factors)
{^NOONED
                    rhovy(ilist) = sum(pw(igrid)%w(i^D:i^D+1, m2_) * weighing_factors)
}
{^IFTHREED
                    rhovz(ilist) = sum(pw(igrid)%w(i^D:i^D+1, m3_) * weighing_factors)
}
                end if
            end do
        end do
        
        
        
        if(mype .GT. 0) then
            call MPI_Reduce(rho,  0, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(rhovx,  0, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(rhovy,  0, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(rhovz,  0, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(rhoe,  0, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, rho, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, rhovx, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, rhovy, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, rhovz, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, rhoe, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        
        
        ret=0  
    end function
END MODULE mpiamrvac_interface

