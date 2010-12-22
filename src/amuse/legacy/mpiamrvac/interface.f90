MODULE mpiamrvac_interface

    CHARACTER(LEN=1024) :: parameters_filename = 'amrvac.par'
    
CONTAINS

    FUNCTION initialize_code()
        include 'amrvacdef.f'
        INTEGER initialize_code
        logical :: file_exists
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierrmpi)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierrmpi)

        icomm=MPI_COMM_WORLD
        
        
        INQUIRE(FILE= TRIM(parameters_filename), EXIST=file_exists)
        
        if (.NOT.file_exists) THEN
            initialize_code = -1
            return
        end if
        
        call readparameters(parameters_filename)
        
        eqpar(gamma_) = 5.0d0/3.0d0
        initialize_code = 0
    END FUNCTION
    
    FUNCTION cleanup_code()
        INTEGER cleanup_code
        cleanup_code = 0
    END FUNCTION
    
    FUNCTION commit_parameters()
        INTEGER commit_parameters
        
        call initialize_vars()
        call init_comm_types   
        
        call initglobaldata_usr
        call initglobaldata 
        
        ! form and initialize all grids at level one
        call initlevelone
        
        commit_parameters = 0
    END FUNCTION
    
    FUNCTION recommit_parameters()
        INTEGER recommit_parameters
        recommit_parameters = 0
    END FUNCTION

    FUNCTION initialize_grid()
        INTEGER :: initialize_grid
        
        ! set up and initialize finer level grids, if needed
        !call settree
   
        if (levmax>levmin) call allocateBflux

        initialize_grid = 0
        
    END FUNCTION
    

    function set_dt(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dt
        double precision :: inputvalue
        dt = inputvalue
        set_dt = 0   
    end function

    function get_dt(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dt
        double precision :: outputvalue
        outputvalue = dt
        get_dt = 0   
    end function

    function set_nbufferx1(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nbufferx1
        integer :: inputvalue
        nbufferx1 = inputvalue
        set_nbufferx1 = 0   
    end function

    function get_nbufferx1(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nbufferx1
        integer :: outputvalue
        outputvalue = nbufferx1
        get_nbufferx1 = 0   
    end function

    function set_nbufferx2(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nbufferx2
        integer :: inputvalue
        nbufferx2 = inputvalue
        set_nbufferx2 = 0   
    end function

    function get_nbufferx2(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nbufferx2
        integer :: outputvalue
        outputvalue = nbufferx2
        get_nbufferx2 = 0   
    end function

    function set_nbufferx3(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nbufferx3
        integer :: inputvalue
        nbufferx3 = inputvalue
        set_nbufferx3 = 0   
    end function

    function get_nbufferx3(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nbufferx3
        integer :: outputvalue
        outputvalue = nbufferx3
        get_nbufferx3 = 0   
    end function

    function set_mxnest(inputvalue)
        include 'amrvacdef.f'
        integer :: set_mxnest
        integer :: inputvalue
        mxnest = inputvalue
        set_mxnest = 0   
    end function

    function get_mxnest(outputvalue)
        include 'amrvacdef.f'
        integer :: get_mxnest
        integer :: outputvalue
        outputvalue = mxnest
        get_mxnest = 0   
    end function

    function set_dixb(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dixb
        integer :: inputvalue
        dixb = inputvalue
        set_dixb = 0   
    end function

    function get_dixb(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dixb
        integer :: outputvalue
        outputvalue = dixb
        get_dixb = 0   
    end function

    function set_levmin(inputvalue)
        include 'amrvacdef.f'
        integer :: set_levmin
        integer :: inputvalue
        levmin = inputvalue
        set_levmin = 0   
    end function

    function get_levmin(outputvalue)
        include 'amrvacdef.f'
        integer :: get_levmin
        integer :: outputvalue
        outputvalue = levmin
        get_levmin = 0   
    end function

    function set_levmax(inputvalue)
        include 'amrvacdef.f'
        integer :: set_levmax
        integer :: inputvalue
        levmax = inputvalue
        set_levmax = 0   
    end function

    function get_levmax(outputvalue)
        include 'amrvacdef.f'
        integer :: get_levmax
        integer :: outputvalue
        outputvalue = levmax
        get_levmax = 0   
    end function

    function set_skipfinestep(inputvalue)
        include 'amrvacdef.f'
        integer :: set_skipfinestep
        logical :: inputvalue
        skipfinestep = inputvalue
        set_skipfinestep = 0   
    end function

    function get_skipfinestep(outputvalue)
        include 'amrvacdef.f'
        integer :: get_skipfinestep
        logical :: outputvalue
        outputvalue = skipfinestep
        get_skipfinestep = 0   
    end function

    function set_time_advance(inputvalue)
        include 'amrvacdef.f'
        integer :: set_time_advance
        logical :: inputvalue
        time_advance = inputvalue
        set_time_advance = 0   
    end function

    function get_time_advance(outputvalue)
        include 'amrvacdef.f'
        integer :: get_time_advance
        logical :: outputvalue
        outputvalue = time_advance
        get_time_advance = 0   
    end function

    function set_courantpar(inputvalue)
        include 'amrvacdef.f'
        integer :: set_courantpar
        double precision :: inputvalue
        courantpar = inputvalue
        set_courantpar = 0   
    end function

    function get_courantpar(outputvalue)
        include 'amrvacdef.f'
        integer :: get_courantpar
        double precision :: outputvalue
        outputvalue = courantpar
        get_courantpar = 0   
    end function

    function set_dtpar(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dtpar
        double precision :: inputvalue
        dtpar = inputvalue
        set_dtpar = 0   
    end function

    function get_dtpar(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dtpar
        double precision :: outputvalue
        outputvalue = dtpar
        get_dtpar = 0   
    end function

    function set_dtdiffpar(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dtdiffpar
        double precision :: inputvalue
        dtdiffpar = inputvalue
        set_dtdiffpar = 0   
    end function

    function get_dtdiffpar(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dtdiffpar
        double precision :: outputvalue
        outputvalue = dtdiffpar
        get_dtdiffpar = 0   
    end function

    function set_t(inputvalue)
        include 'amrvacdef.f'
        integer :: set_t
        double precision :: inputvalue
        t = inputvalue
        set_t = 0   
    end function

    function get_t(outputvalue)
        include 'amrvacdef.f'
        integer :: get_t
        double precision :: outputvalue
        outputvalue = t
        get_t = 0   
    end function

    function set_tmax(inputvalue)
        include 'amrvacdef.f'
        integer :: set_tmax
        double precision :: inputvalue
        tmax = inputvalue
        set_tmax = 0   
    end function

    function get_tmax(outputvalue)
        include 'amrvacdef.f'
        integer :: get_tmax
        double precision :: outputvalue
        outputvalue = tmax
        get_tmax = 0   
    end function

    function set_dtmin(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dtmin
        double precision :: inputvalue
        dtmin = inputvalue
        set_dtmin = 0   
    end function

    function get_dtmin(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dtmin
        double precision :: outputvalue
        outputvalue = dtmin
        get_dtmin = 0   
    end function

    function set_residmin(inputvalue)
        include 'amrvacdef.f'
        integer :: set_residmin
        double precision :: inputvalue
        residmin = inputvalue
        set_residmin = 0   
    end function

    function get_residmin(outputvalue)
        include 'amrvacdef.f'
        integer :: get_residmin
        double precision :: outputvalue
        outputvalue = residmin
        get_residmin = 0   
    end function

    function set_residmax(inputvalue)
        include 'amrvacdef.f'
        integer :: set_residmax
        double precision :: inputvalue
        residmax = inputvalue
        set_residmax = 0   
    end function

    function get_residmax(outputvalue)
        include 'amrvacdef.f'
        integer :: get_residmax
        double precision :: outputvalue
        outputvalue = residmax
        get_residmax = 0   
    end function

    function set_residual(inputvalue)
        include 'amrvacdef.f'
        integer :: set_residual
        double precision :: inputvalue
        residual = inputvalue
        set_residual = 0   
    end function

    function get_residual(outputvalue)
        include 'amrvacdef.f'
        integer :: get_residual
        double precision :: outputvalue
        outputvalue = residual
        get_residual = 0   
    end function

    function set_tfixgrid(inputvalue)
        include 'amrvacdef.f'
        integer :: set_tfixgrid
        double precision :: inputvalue
        tfixgrid = inputvalue
        set_tfixgrid = 0   
    end function

    function get_tfixgrid(outputvalue)
        include 'amrvacdef.f'
        integer :: get_tfixgrid
        double precision :: outputvalue
        outputvalue = tfixgrid
        get_tfixgrid = 0   
    end function

    function set_tvdlfeps(inputvalue)
        include 'amrvacdef.f'
        integer :: set_tvdlfeps
        double precision :: inputvalue
        tvdlfeps = inputvalue
        set_tvdlfeps = 0   
    end function

    function get_tvdlfeps(outputvalue)
        include 'amrvacdef.f'
        integer :: get_tvdlfeps
        double precision :: outputvalue
        outputvalue = tvdlfeps
        get_tvdlfeps = 0   
    end function

    function set_mcbeta(inputvalue)
        include 'amrvacdef.f'
        integer :: set_mcbeta
        double precision :: inputvalue
        mcbeta = inputvalue
        set_mcbeta = 0   
    end function

    function get_mcbeta(outputvalue)
        include 'amrvacdef.f'
        integer :: get_mcbeta
        double precision :: outputvalue
        outputvalue = mcbeta
        get_mcbeta = 0   
    end function

    function set_divbdiff(inputvalue)
        include 'amrvacdef.f'
        integer :: set_divbdiff
        double precision :: inputvalue
        divbdiff = inputvalue
        set_divbdiff = 0   
    end function

    function get_divbdiff(outputvalue)
        include 'amrvacdef.f'
        integer :: get_divbdiff
        double precision :: outputvalue
        outputvalue = divbdiff
        get_divbdiff = 0   
    end function

    function set_smallp(inputvalue)
        include 'amrvacdef.f'
        integer :: set_smallp
        double precision :: inputvalue
        smallp = inputvalue
        set_smallp = 0   
    end function

    function get_smallp(outputvalue)
        include 'amrvacdef.f'
        integer :: get_smallp
        double precision :: outputvalue
        outputvalue = smallp
        get_smallp = 0   
    end function

    function set_smallrho(inputvalue)
        include 'amrvacdef.f'
        integer :: set_smallrho
        double precision :: inputvalue
        smallrho = inputvalue
        set_smallrho = 0   
    end function

    function get_smallrho(outputvalue)
        include 'amrvacdef.f'
        integer :: get_smallrho
        double precision :: outputvalue
        outputvalue = smallrho
        get_smallrho = 0   
    end function

    function set_dmaxvel(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dmaxvel
        double precision :: inputvalue
        dmaxvel = inputvalue
        set_dmaxvel = 0   
    end function

    function get_dmaxvel(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dmaxvel
        double precision :: outputvalue
        outputvalue = dmaxvel
        get_dmaxvel = 0   
    end function

    function set_tolernr(inputvalue)
        include 'amrvacdef.f'
        integer :: set_tolernr
        double precision :: inputvalue
        tolernr = inputvalue
        set_tolernr = 0   
    end function

    function get_tolernr(outputvalue)
        include 'amrvacdef.f'
        integer :: get_tolernr
        double precision :: outputvalue
        outputvalue = tolernr
        get_tolernr = 0   
    end function

    function set_absaccnr(inputvalue)
        include 'amrvacdef.f'
        integer :: set_absaccnr
        double precision :: inputvalue
        absaccnr = inputvalue
        set_absaccnr = 0   
    end function

    function get_absaccnr(outputvalue)
        include 'amrvacdef.f'
        integer :: get_absaccnr
        double precision :: outputvalue
        outputvalue = absaccnr
        get_absaccnr = 0   
    end function

    function set_cfrac(inputvalue)
        include 'amrvacdef.f'
        integer :: set_cfrac
        double precision :: inputvalue
        cfrac = inputvalue
        set_cfrac = 0   
    end function

    function get_cfrac(outputvalue)
        include 'amrvacdef.f'
        integer :: get_cfrac
        double precision :: outputvalue
        outputvalue = cfrac
        get_cfrac = 0   
    end function

    function set_x1ptms(inputvalue)
        include 'amrvacdef.f'
        integer :: set_x1ptms
        double precision :: inputvalue
        x1ptms = inputvalue
        set_x1ptms = 0   
    end function

    function get_x1ptms(outputvalue)
        include 'amrvacdef.f'
        integer :: get_x1ptms
        double precision :: outputvalue
        outputvalue = x1ptms
        get_x1ptms = 0   
    end function

    function set_x2ptms(inputvalue)
        include 'amrvacdef.f'
        integer :: set_x2ptms
        double precision :: inputvalue
        x2ptms = inputvalue
        set_x2ptms = 0   
    end function

    function get_x2ptms(outputvalue)
        include 'amrvacdef.f'
        integer :: get_x2ptms
        double precision :: outputvalue
        outputvalue = x2ptms
        get_x2ptms = 0   
    end function

    function set_x3ptms(inputvalue)
        include 'amrvacdef.f'
        integer :: set_x3ptms
        double precision :: inputvalue
        x3ptms = inputvalue
        set_x3ptms = 0   
    end function

    function get_x3ptms(outputvalue)
        include 'amrvacdef.f'
        integer :: get_x3ptms
        double precision :: outputvalue
        outputvalue = x3ptms
        get_x3ptms = 0   
    end function

    function set_ptmass(inputvalue)
        include 'amrvacdef.f'
        integer :: set_ptmass
        double precision :: inputvalue
        ptmass = inputvalue
        set_ptmass = 0   
    end function

    function get_ptmass(outputvalue)
        include 'amrvacdef.f'
        integer :: get_ptmass
        double precision :: outputvalue
        outputvalue = ptmass
        get_ptmass = 0   
    end function

    function set_ratebdflux(inputvalue)
        include 'amrvacdef.f'
        integer :: set_ratebdflux
        double precision :: inputvalue
        ratebdflux = inputvalue
        set_ratebdflux = 0   
    end function

    function get_ratebdflux(outputvalue)
        include 'amrvacdef.f'
        integer :: get_ratebdflux
        double precision :: outputvalue
        outputvalue = ratebdflux
        get_ratebdflux = 0   
    end function

    function set_normt(inputvalue)
        include 'amrvacdef.f'
        integer :: set_normt
        double precision :: inputvalue
        normt = inputvalue
        set_normt = 0   
    end function

    function get_normt(outputvalue)
        include 'amrvacdef.f'
        integer :: get_normt
        double precision :: outputvalue
        outputvalue = normt
        get_normt = 0   
    end function

    function set_time_bc(inputvalue)
        include 'amrvacdef.f'
        integer :: set_time_bc
        double precision :: inputvalue
        time_bc = inputvalue
        set_time_bc = 0   
    end function

    function get_time_bc(outputvalue)
        include 'amrvacdef.f'
        integer :: get_time_bc
        double precision :: outputvalue
        outputvalue = time_bc
        get_time_bc = 0   
    end function

    function set_it(inputvalue)
        include 'amrvacdef.f'
        integer :: set_it
        integer :: inputvalue
        it = inputvalue
        set_it = 0   
    end function

    function get_it(outputvalue)
        include 'amrvacdef.f'
        integer :: get_it
        integer :: outputvalue
        outputvalue = it
        get_it = 0   
    end function

    function set_itmax(inputvalue)
        include 'amrvacdef.f'
        integer :: set_itmax
        integer :: inputvalue
        itmax = inputvalue
        set_itmax = 0   
    end function

    function get_itmax(outputvalue)
        include 'amrvacdef.f'
        integer :: get_itmax
        integer :: outputvalue
        outputvalue = itmax
        get_itmax = 0   
    end function

    function set_itmin(inputvalue)
        include 'amrvacdef.f'
        integer :: set_itmin
        integer :: inputvalue
        itmin = inputvalue
        set_itmin = 0   
    end function

    function get_itmin(outputvalue)
        include 'amrvacdef.f'
        integer :: get_itmin
        integer :: outputvalue
        outputvalue = itmin
        get_itmin = 0   
    end function

    function set_slowsteps(inputvalue)
        include 'amrvacdef.f'
        integer :: set_slowsteps
        integer :: inputvalue
        slowsteps = inputvalue
        set_slowsteps = 0   
    end function

    function get_slowsteps(outputvalue)
        include 'amrvacdef.f'
        integer :: get_slowsteps
        integer :: outputvalue
        outputvalue = slowsteps
        get_slowsteps = 0   
    end function

    function set_typepario(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typepario
        integer :: inputvalue
        typepario = inputvalue
        set_typepario = 0   
    end function

    function get_typepario(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typepario
        integer :: outputvalue
        outputvalue = typepario
        get_typepario = 0   
    end function

    function set_itfixgrid(inputvalue)
        include 'amrvacdef.f'
        integer :: set_itfixgrid
        integer :: inputvalue
        itfixgrid = inputvalue
        set_itfixgrid = 0   
    end function

    function get_itfixgrid(outputvalue)
        include 'amrvacdef.f'
        integer :: get_itfixgrid
        integer :: outputvalue
        outputvalue = itfixgrid
        get_itfixgrid = 0   
    end function

    function set_nwauxio(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nwauxio
        integer :: inputvalue
        nwauxio = inputvalue
        set_nwauxio = 0   
    end function

    function get_nwauxio(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nwauxio
        integer :: outputvalue
        outputvalue = nwauxio
        get_nwauxio = 0   
    end function

    function set_istep(inputvalue)
        include 'amrvacdef.f'
        integer :: set_istep
        integer :: inputvalue
        istep = inputvalue
        set_istep = 0   
    end function

    function get_istep(outputvalue)
        include 'amrvacdef.f'
        integer :: get_istep
        integer :: outputvalue
        outputvalue = istep
        get_istep = 0   
    end function

    function set_nstep(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nstep
        integer :: inputvalue
        nstep = inputvalue
        set_nstep = 0   
    end function

    function get_nstep(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nstep
        integer :: outputvalue
        outputvalue = nstep
        get_nstep = 0   
    end function

    function set_errorestimate(inputvalue)
        include 'amrvacdef.f'
        integer :: set_errorestimate
        integer :: inputvalue
        errorestimate = inputvalue
        set_errorestimate = 0   
    end function

    function get_errorestimate(outputvalue)
        include 'amrvacdef.f'
        integer :: get_errorestimate
        integer :: outputvalue
        outputvalue = errorestimate
        get_errorestimate = 0   
    end function

    function set_nxdiffusehllc(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nxdiffusehllc
        integer :: inputvalue
        nxdiffusehllc = inputvalue
        set_nxdiffusehllc = 0   
    end function

    function get_nxdiffusehllc(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nxdiffusehllc
        integer :: outputvalue
        outputvalue = nxdiffusehllc
        get_nxdiffusehllc = 0   
    end function

    function set_typespherical(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typespherical
        integer :: inputvalue
        typespherical = inputvalue
        set_typespherical = 0   
    end function

    function get_typespherical(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typespherical
        integer :: outputvalue
        outputvalue = typespherical
        get_typespherical = 0   
    end function

    function set_maxitnr(inputvalue)
        include 'amrvacdef.f'
        integer :: set_maxitnr
        integer :: inputvalue
        maxitnr = inputvalue
        set_maxitnr = 0   
    end function

    function get_maxitnr(outputvalue)
        include 'amrvacdef.f'
        integer :: get_maxitnr
        integer :: outputvalue
        outputvalue = maxitnr
        get_maxitnr = 0   
    end function

    function set_nflatgetaux(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nflatgetaux
        integer :: inputvalue
        nflatgetaux = inputvalue
        set_nflatgetaux = 0   
    end function

    function get_nflatgetaux(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nflatgetaux
        integer :: outputvalue
        outputvalue = nflatgetaux
        get_nflatgetaux = 0   
    end function

    function set_level_io(inputvalue)
        include 'amrvacdef.f'
        integer :: set_level_io
        integer :: inputvalue
        level_io = inputvalue
        set_level_io = 0   
    end function

    function get_level_io(outputvalue)
        include 'amrvacdef.f'
        integer :: get_level_io
        integer :: outputvalue
        outputvalue = level_io
        get_level_io = 0   
    end function

    function set_ncool(inputvalue)
        include 'amrvacdef.f'
        integer :: set_ncool
        integer :: inputvalue
        ncool = inputvalue
        set_ncool = 0   
    end function

    function get_ncool(outputvalue)
        include 'amrvacdef.f'
        integer :: get_ncool
        integer :: outputvalue
        outputvalue = ncool
        get_ncool = 0   
    end function

    function set_cmulti(inputvalue)
        include 'amrvacdef.f'
        integer :: set_cmulti
        integer :: inputvalue
        cmulti = inputvalue
        set_cmulti = 0   
    end function

    function get_cmulti(outputvalue)
        include 'amrvacdef.f'
        integer :: get_cmulti
        integer :: outputvalue
        outputvalue = cmulti
        get_cmulti = 0   
    end function

    function set_snapshotini(inputvalue)
        include 'amrvacdef.f'
        integer :: set_snapshotini
        integer :: inputvalue
        snapshotini = inputvalue
        set_snapshotini = 0   
    end function

    function get_snapshotini(outputvalue)
        include 'amrvacdef.f'
        integer :: get_snapshotini
        integer :: outputvalue
        outputvalue = snapshotini
        get_snapshotini = 0   
    end function

    function set_ixtest1(inputvalue)
        include 'amrvacdef.f'
        integer :: set_ixtest1
        integer :: inputvalue
        ixtest1 = inputvalue
        set_ixtest1 = 0   
    end function

    function get_ixtest1(outputvalue)
        include 'amrvacdef.f'
        integer :: get_ixtest1
        integer :: outputvalue
        outputvalue = ixtest1
        get_ixtest1 = 0   
    end function

    function set_ixtest2(inputvalue)
        include 'amrvacdef.f'
        integer :: set_ixtest2
        integer :: inputvalue
        ixtest2 = inputvalue
        set_ixtest2 = 0   
    end function

    function get_ixtest2(outputvalue)
        include 'amrvacdef.f'
        integer :: get_ixtest2
        integer :: outputvalue
        outputvalue = ixtest2
        get_ixtest2 = 0   
    end function

    function set_ixtest3(inputvalue)
        include 'amrvacdef.f'
        integer :: set_ixtest3
        integer :: inputvalue
        ixtest3 = inputvalue
        set_ixtest3 = 0   
    end function

    function get_ixtest3(outputvalue)
        include 'amrvacdef.f'
        integer :: get_ixtest3
        integer :: outputvalue
        outputvalue = ixtest3
        get_ixtest3 = 0   
    end function

    function set_iwtest(inputvalue)
        include 'amrvacdef.f'
        integer :: set_iwtest
        integer :: inputvalue
        iwtest = inputvalue
        set_iwtest = 0   
    end function

    function get_iwtest(outputvalue)
        include 'amrvacdef.f'
        integer :: get_iwtest
        integer :: outputvalue
        outputvalue = iwtest
        get_iwtest = 0   
    end function

    function set_idimtest(inputvalue)
        include 'amrvacdef.f'
        integer :: set_idimtest
        integer :: inputvalue
        idimtest = inputvalue
        set_idimtest = 0   
    end function

    function get_idimtest(outputvalue)
        include 'amrvacdef.f'
        integer :: get_idimtest
        integer :: outputvalue
        outputvalue = idimtest
        get_idimtest = 0   
    end function

    function set_saveigrid(inputvalue)
        include 'amrvacdef.f'
        integer :: set_saveigrid
        integer :: inputvalue
        saveigrid = inputvalue
        set_saveigrid = 0   
    end function

    function get_saveigrid(outputvalue)
        include 'amrvacdef.f'
        integer :: get_saveigrid
        integer :: outputvalue
        outputvalue = saveigrid
        get_saveigrid = 0   
    end function

    function set_typecourant(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typecourant
        character*79 :: inputvalue
        typecourant = inputvalue
        set_typecourant = 0   
    end function

    function get_typecourant(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typecourant
        character*79 :: outputvalue
        outputvalue = typecourant
        get_typecourant = 0   
    end function

    function set_typeresid(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeresid
        character*79 :: inputvalue
        typeresid = inputvalue
        set_typeresid = 0   
    end function

    function get_typeresid(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeresid
        character*79 :: outputvalue
        outputvalue = typeresid
        get_typeresid = 0   
    end function

    function set_typeadvance(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeadvance
        character*79 :: inputvalue
        typeadvance = inputvalue
        set_typeadvance = 0   
    end function

    function get_typeadvance(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeadvance
        character*79 :: outputvalue
        outputvalue = typeadvance
        get_typeadvance = 0   
    end function

    function set_typelimited(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typelimited
        character*79 :: inputvalue
        typelimited = inputvalue
        set_typelimited = 0   
    end function

    function get_typelimited(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typelimited
        character*79 :: outputvalue
        outputvalue = typelimited
        get_typelimited = 0   
    end function

    function set_typesourcesplit(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typesourcesplit
        character*79 :: inputvalue
        typesourcesplit = inputvalue
        set_typesourcesplit = 0   
    end function

    function get_typesourcesplit(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typesourcesplit
        character*79 :: outputvalue
        outputvalue = typesourcesplit
        get_typesourcesplit = 0   
    end function

    function set_typelimiter(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typelimiter
        character*79 :: inputvalue
        typelimiter = inputvalue
        set_typelimiter = 0   
    end function

    function get_typelimiter(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typelimiter
        character*79 :: outputvalue
        outputvalue = typelimiter
        get_typelimiter = 0   
    end function

    function set_typegradlimiter(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typegradlimiter
        character*79 :: inputvalue
        typegradlimiter = inputvalue
        set_typegradlimiter = 0   
    end function

    function get_typegradlimiter(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typegradlimiter
        character*79 :: outputvalue
        outputvalue = typegradlimiter
        get_typegradlimiter = 0   
    end function

    function set_typeprolonglimit(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeprolonglimit
        character*79 :: inputvalue
        typeprolonglimit = inputvalue
        set_typeprolonglimit = 0   
    end function

    function get_typeprolonglimit(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeprolonglimit
        character*79 :: outputvalue
        outputvalue = typeprolonglimit
        get_typeprolonglimit = 0   
    end function

    function set_typetvd(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typetvd
        character*79 :: inputvalue
        typetvd = inputvalue
        set_typetvd = 0   
    end function

    function get_typetvd(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typetvd
        character*79 :: outputvalue
        outputvalue = typetvd
        get_typetvd = 0   
    end function

    function set_typetvdlf(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typetvdlf
        character*79 :: inputvalue
        typetvdlf = inputvalue
        set_typetvdlf = 0   
    end function

    function get_typetvdlf(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typetvdlf
        character*79 :: outputvalue
        outputvalue = typetvdlf
        get_typetvdlf = 0   
    end function

    function set_typeaverage(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeaverage
        character*79 :: inputvalue
        typeaverage = inputvalue
        set_typeaverage = 0   
    end function

    function get_typeaverage(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeaverage
        character*79 :: outputvalue
        outputvalue = typeaverage
        get_typeaverage = 0   
    end function

    function set_typedimsplit(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typedimsplit
        character*79 :: inputvalue
        typedimsplit = inputvalue
        set_typedimsplit = 0   
    end function

    function get_typedimsplit(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typedimsplit
        character*79 :: outputvalue
        outputvalue = typedimsplit
        get_typedimsplit = 0   
    end function

    function set_typeaxial(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeaxial
        character*79 :: inputvalue
        typeaxial = inputvalue
        set_typeaxial = 0   
    end function

    function get_typeaxial(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeaxial
        character*79 :: outputvalue
        outputvalue = typeaxial
        get_typeaxial = 0   
    end function

    function set_typepoly(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typepoly
        character*79 :: inputvalue
        typepoly = inputvalue
        set_typepoly = 0   
    end function

    function get_typepoly(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typepoly
        character*79 :: outputvalue
        outputvalue = typepoly
        get_typepoly = 0   
    end function

    function set_typedivbdiff(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typedivbdiff
        character*79 :: inputvalue
        typedivbdiff = inputvalue
        set_typedivbdiff = 0   
    end function

    function get_typedivbdiff(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typedivbdiff
        character*79 :: outputvalue
        outputvalue = typedivbdiff
        get_typedivbdiff = 0   
    end function

    function set_typedivbfix(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typedivbfix
        character*79 :: inputvalue
        typedivbfix = inputvalue
        set_typedivbfix = 0   
    end function

    function get_typedivbfix(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typedivbfix
        character*79 :: outputvalue
        outputvalue = typedivbfix
        get_typedivbfix = 0   
    end function

    function set_typediv(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typediv
        character*79 :: inputvalue
        typediv = inputvalue
        set_typediv = 0   
    end function

    function get_typediv(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typediv
        character*79 :: outputvalue
        outputvalue = typediv
        get_typediv = 0   
    end function

    function set_typegrad(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typegrad
        character*79 :: inputvalue
        typegrad = inputvalue
        set_typegrad = 0   
    end function

    function get_typegrad(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typegrad
        character*79 :: outputvalue
        outputvalue = typegrad
        get_typegrad = 0   
    end function

    function set_typeglm(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeglm
        character*79 :: inputvalue
        typeglm = inputvalue
        set_typeglm = 0   
    end function

    function get_typeglm(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeglm
        character*79 :: outputvalue
        outputvalue = typeglm
        get_typeglm = 0   
    end function

    function set_coolcurve(inputvalue)
        include 'amrvacdef.f'
        integer :: set_coolcurve
        character*79 :: inputvalue
        coolcurve = inputvalue
        set_coolcurve = 0   
    end function

    function get_coolcurve(outputvalue)
        include 'amrvacdef.f'
        integer :: get_coolcurve
        character*79 :: outputvalue
        outputvalue = coolcurve
        get_coolcurve = 0   
    end function

    function set_coolmethod(inputvalue)
        include 'amrvacdef.f'
        integer :: set_coolmethod
        character*79 :: inputvalue
        coolmethod = inputvalue
        set_coolmethod = 0   
    end function

    function get_coolmethod(outputvalue)
        include 'amrvacdef.f'
        integer :: get_coolmethod
        character*79 :: outputvalue
        outputvalue = coolmethod
        get_coolmethod = 0   
    end function

    function set_typeghostfill(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typeghostfill
        character*79 :: inputvalue
        typeghostfill = inputvalue
        set_typeghostfill = 0   
    end function

    function get_typeghostfill(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typeghostfill
        character*79 :: outputvalue
        outputvalue = typeghostfill
        get_typeghostfill = 0   
    end function

    function set_typegridfill(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typegridfill
        character*79 :: inputvalue
        typegridfill = inputvalue
        set_typegridfill = 0   
    end function

    function get_typegridfill(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typegridfill
        character*79 :: outputvalue
        outputvalue = typegridfill
        get_typegridfill = 0   
    end function

    function set_filenameout(inputvalue)
        include 'amrvacdef.f'
        integer :: set_filenameout
        character*79 :: inputvalue
        filenameout = inputvalue
        set_filenameout = 0   
    end function

    function get_filenameout(outputvalue)
        include 'amrvacdef.f'
        integer :: get_filenameout
        character*79 :: outputvalue
        outputvalue = filenameout
        get_filenameout = 0   
    end function

    function set_filenameini(inputvalue)
        include 'amrvacdef.f'
        integer :: set_filenameini
        character*79 :: inputvalue
        filenameini = inputvalue
        set_filenameini = 0   
    end function

    function get_filenameini(outputvalue)
        include 'amrvacdef.f'
        integer :: get_filenameini
        character*79 :: outputvalue
        outputvalue = filenameini
        get_filenameini = 0   
    end function

    function set_filenamelog(inputvalue)
        include 'amrvacdef.f'
        integer :: set_filenamelog
        character*79 :: inputvalue
        filenamelog = inputvalue
        set_filenamelog = 0   
    end function

    function get_filenamelog(outputvalue)
        include 'amrvacdef.f'
        integer :: get_filenamelog
        character*79 :: outputvalue
        outputvalue = filenamelog
        get_filenamelog = 0   
    end function

    function set_fileheadout(inputvalue)
        include 'amrvacdef.f'
        integer :: set_fileheadout
        character*79 :: inputvalue
        fileheadout = inputvalue
        set_fileheadout = 0   
    end function

    function get_fileheadout(outputvalue)
        include 'amrvacdef.f'
        integer :: get_fileheadout
        character*79 :: outputvalue
        outputvalue = fileheadout
        get_fileheadout = 0   
    end function

    function set_wnames(inputvalue)
        include 'amrvacdef.f'
        integer :: set_wnames
        character*79 :: inputvalue
        wnames = inputvalue
        set_wnames = 0   
    end function

    function get_wnames(outputvalue)
        include 'amrvacdef.f'
        integer :: get_wnames
        character*79 :: outputvalue
        outputvalue = wnames
        get_wnames = 0   
    end function

    function set_primnames(inputvalue)
        include 'amrvacdef.f'
        integer :: set_primnames
        character*79 :: inputvalue
        primnames = inputvalue
        set_primnames = 0   
    end function

    function get_primnames(outputvalue)
        include 'amrvacdef.f'
        integer :: get_primnames
        character*79 :: outputvalue
        outputvalue = primnames
        get_primnames = 0   
    end function

    function set_typefilelog(inputvalue)
        include 'amrvacdef.f'
        integer :: set_typefilelog
        character*79 :: inputvalue
        typefilelog = inputvalue
        set_typefilelog = 0   
    end function

    function get_typefilelog(outputvalue)
        include 'amrvacdef.f'
        integer :: get_typefilelog
        character*79 :: outputvalue
        outputvalue = typefilelog
        get_typefilelog = 0   
    end function

    function set_convert_type(inputvalue)
        include 'amrvacdef.f'
        integer :: set_convert_type
        character*79 :: inputvalue
        convert_type = inputvalue
        set_convert_type = 0   
    end function

    function get_convert_type(outputvalue)
        include 'amrvacdef.f'
        integer :: get_convert_type
        character*79 :: outputvalue
        outputvalue = convert_type
        get_convert_type = 0   
    end function

    function set_dxfiletype(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dxfiletype
        character*79 :: inputvalue
        dxfiletype = inputvalue
        set_dxfiletype = 0   
    end function

    function get_dxfiletype(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dxfiletype
        character*79 :: outputvalue
        outputvalue = dxfiletype
        get_dxfiletype = 0   
    end function

    function set_teststr(inputvalue)
        include 'amrvacdef.f'
        integer :: set_teststr
        character*79 :: inputvalue
        teststr = inputvalue
        set_teststr = 0   
    end function

    function get_teststr(outputvalue)
        include 'amrvacdef.f'
        integer :: get_teststr
        character*79 :: outputvalue
        outputvalue = teststr
        get_teststr = 0   
    end function

    function set_time_accurate(inputvalue)
        include 'amrvacdef.f'
        integer :: set_time_accurate
        logical :: inputvalue
        time_accurate = inputvalue
        set_time_accurate = 0   
    end function

    function get_time_accurate(outputvalue)
        include 'amrvacdef.f'
        integer :: get_time_accurate
        logical :: outputvalue
        outputvalue = time_accurate
        get_time_accurate = 0   
    end function

    function set_addmpibarrier(inputvalue)
        include 'amrvacdef.f'
        integer :: set_addmpibarrier
        logical :: inputvalue
        addmpibarrier = inputvalue
        set_addmpibarrier = 0   
    end function

    function get_addmpibarrier(outputvalue)
        include 'amrvacdef.f'
        integer :: get_addmpibarrier
        logical :: outputvalue
        outputvalue = addmpibarrier
        get_addmpibarrier = 0   
    end function

    function set_tmaxexact(inputvalue)
        include 'amrvacdef.f'
        integer :: set_tmaxexact
        logical :: inputvalue
        tmaxexact = inputvalue
        set_tmaxexact = 0   
    end function

    function get_tmaxexact(outputvalue)
        include 'amrvacdef.f'
        integer :: get_tmaxexact
        logical :: outputvalue
        outputvalue = tmaxexact
        get_tmaxexact = 0   
    end function

    function set_treset(inputvalue)
        include 'amrvacdef.f'
        integer :: set_treset
        logical :: inputvalue
        treset = inputvalue
        set_treset = 0   
    end function

    function get_treset(outputvalue)
        include 'amrvacdef.f'
        integer :: get_treset
        logical :: outputvalue
        outputvalue = treset
        get_treset = 0   
    end function

    function set_itreset(inputvalue)
        include 'amrvacdef.f'
        integer :: set_itreset
        logical :: inputvalue
        itreset = inputvalue
        set_itreset = 0   
    end function

    function get_itreset(outputvalue)
        include 'amrvacdef.f'
        integer :: get_itreset
        logical :: outputvalue
        outputvalue = itreset
        get_itreset = 0   
    end function

    function set_firstprocess(inputvalue)
        include 'amrvacdef.f'
        integer :: set_firstprocess
        logical :: inputvalue
        firstprocess = inputvalue
        set_firstprocess = 0   
    end function

    function get_firstprocess(outputvalue)
        include 'amrvacdef.f'
        integer :: get_firstprocess
        logical :: outputvalue
        outputvalue = firstprocess
        get_firstprocess = 0   
    end function

    function set_fixprocess(inputvalue)
        include 'amrvacdef.f'
        integer :: set_fixprocess
        logical :: inputvalue
        fixprocess = inputvalue
        set_fixprocess = 0   
    end function

    function get_fixprocess(outputvalue)
        include 'amrvacdef.f'
        integer :: get_fixprocess
        logical :: outputvalue
        outputvalue = fixprocess
        get_fixprocess = 0   
    end function

    function set_flathllc(inputvalue)
        include 'amrvacdef.f'
        integer :: set_flathllc
        logical :: inputvalue
        flathllc = inputvalue
        set_flathllc = 0   
    end function

    function get_flathllc(outputvalue)
        include 'amrvacdef.f'
        integer :: get_flathllc
        logical :: outputvalue
        outputvalue = flathllc
        get_flathllc = 0   
    end function

    function set_flatcd(inputvalue)
        include 'amrvacdef.f'
        integer :: set_flatcd
        logical :: inputvalue
        flatcd = inputvalue
        set_flatcd = 0   
    end function

    function get_flatcd(outputvalue)
        include 'amrvacdef.f'
        integer :: get_flatcd
        logical :: outputvalue
        outputvalue = flatcd
        get_flatcd = 0   
    end function

    function set_flatsh(inputvalue)
        include 'amrvacdef.f'
        integer :: set_flatsh
        logical :: inputvalue
        flatsh = inputvalue
        set_flatsh = 0   
    end function

    function get_flatsh(outputvalue)
        include 'amrvacdef.f'
        integer :: get_flatsh
        logical :: outputvalue
        outputvalue = flatsh
        get_flatsh = 0   
    end function

    function set_flatppm(inputvalue)
        include 'amrvacdef.f'
        integer :: set_flatppm
        logical :: inputvalue
        flatppm = inputvalue
        set_flatppm = 0   
    end function

    function get_flatppm(outputvalue)
        include 'amrvacdef.f'
        integer :: get_flatppm
        logical :: outputvalue
        outputvalue = flatppm
        get_flatppm = 0   
    end function

    function set_sourcesplit(inputvalue)
        include 'amrvacdef.f'
        integer :: set_sourcesplit
        logical :: inputvalue
        sourcesplit = inputvalue
        set_sourcesplit = 0   
    end function

    function get_sourcesplit(outputvalue)
        include 'amrvacdef.f'
        integer :: get_sourcesplit
        logical :: outputvalue
        outputvalue = sourcesplit
        get_sourcesplit = 0   
    end function

    function set_sourceunsplit(inputvalue)
        include 'amrvacdef.f'
        integer :: set_sourceunsplit
        logical :: inputvalue
        sourceunsplit = inputvalue
        set_sourceunsplit = 0   
    end function

    function get_sourceunsplit(outputvalue)
        include 'amrvacdef.f'
        integer :: get_sourceunsplit
        logical :: outputvalue
        outputvalue = sourceunsplit
        get_sourceunsplit = 0   
    end function

    function set_useprimitive(inputvalue)
        include 'amrvacdef.f'
        integer :: set_useprimitive
        logical :: inputvalue
        useprimitive = inputvalue
        set_useprimitive = 0   
    end function

    function get_useprimitive(outputvalue)
        include 'amrvacdef.f'
        integer :: get_useprimitive
        logical :: outputvalue
        outputvalue = useprimitive
        get_useprimitive = 0   
    end function

    function set_dimsplit(inputvalue)
        include 'amrvacdef.f'
        integer :: set_dimsplit
        logical :: inputvalue
        dimsplit = inputvalue
        set_dimsplit = 0   
    end function

    function get_dimsplit(outputvalue)
        include 'amrvacdef.f'
        integer :: get_dimsplit
        logical :: outputvalue
        outputvalue = dimsplit
        get_dimsplit = 0   
    end function

    function set_restrictprimitive(inputvalue)
        include 'amrvacdef.f'
        integer :: set_restrictprimitive
        logical :: inputvalue
        restrictprimitive = inputvalue
        set_restrictprimitive = 0   
    end function

    function get_restrictprimitive(outputvalue)
        include 'amrvacdef.f'
        integer :: get_restrictprimitive
        logical :: outputvalue
        outputvalue = restrictprimitive
        get_restrictprimitive = 0   
    end function

    function set_prolongprimitive(inputvalue)
        include 'amrvacdef.f'
        integer :: set_prolongprimitive
        logical :: inputvalue
        prolongprimitive = inputvalue
        set_prolongprimitive = 0   
    end function

    function get_prolongprimitive(outputvalue)
        include 'amrvacdef.f'
        integer :: get_prolongprimitive
        logical :: outputvalue
        outputvalue = prolongprimitive
        get_prolongprimitive = 0   
    end function

    function set_coarsenprimitive(inputvalue)
        include 'amrvacdef.f'
        integer :: set_coarsenprimitive
        logical :: inputvalue
        coarsenprimitive = inputvalue
        set_coarsenprimitive = 0   
    end function

    function get_coarsenprimitive(outputvalue)
        include 'amrvacdef.f'
        integer :: get_coarsenprimitive
        logical :: outputvalue
        outputvalue = coarsenprimitive
        get_coarsenprimitive = 0   
    end function

    function set_useprimitiverel(inputvalue)
        include 'amrvacdef.f'
        integer :: set_useprimitiverel
        logical :: inputvalue
        useprimitiverel = inputvalue
        set_useprimitiverel = 0   
    end function

    function get_useprimitiverel(outputvalue)
        include 'amrvacdef.f'
        integer :: get_useprimitiverel
        logical :: outputvalue
        outputvalue = useprimitiverel
        get_useprimitiverel = 0   
    end function

    function set_amrentropy(inputvalue)
        include 'amrvacdef.f'
        integer :: set_amrentropy
        logical :: inputvalue
        amrentropy = inputvalue
        set_amrentropy = 0   
    end function

    function get_amrentropy(outputvalue)
        include 'amrvacdef.f'
        integer :: get_amrentropy
        logical :: outputvalue
        outputvalue = amrentropy
        get_amrentropy = 0   
    end function

    function set_divbfix(inputvalue)
        include 'amrvacdef.f'
        integer :: set_divbfix
        logical :: inputvalue
        divbfix = inputvalue
        set_divbfix = 0   
    end function

    function get_divbfix(outputvalue)
        include 'amrvacdef.f'
        integer :: get_divbfix
        logical :: outputvalue
        outputvalue = divbfix
        get_divbfix = 0   
    end function

    function set_divbwave(inputvalue)
        include 'amrvacdef.f'
        integer :: set_divbwave
        logical :: inputvalue
        divbwave = inputvalue
        set_divbwave = 0   
    end function

    function get_divbwave(outputvalue)
        include 'amrvacdef.f'
        integer :: get_divbwave
        logical :: outputvalue
        outputvalue = divbwave
        get_divbwave = 0   
    end function

    function set_compactres(inputvalue)
        include 'amrvacdef.f'
        integer :: set_compactres
        logical :: inputvalue
        compactres = inputvalue
        set_compactres = 0   
    end function

    function get_compactres(outputvalue)
        include 'amrvacdef.f'
        integer :: get_compactres
        logical :: outputvalue
        outputvalue = compactres
        get_compactres = 0   
    end function

    function set_bnormlf(inputvalue)
        include 'amrvacdef.f'
        integer :: set_bnormlf
        logical :: inputvalue
        bnormlf = inputvalue
        set_bnormlf = 0   
    end function

    function get_bnormlf(outputvalue)
        include 'amrvacdef.f'
        integer :: get_bnormlf
        logical :: outputvalue
        outputvalue = bnormlf
        get_bnormlf = 0   
    end function

    function set_strictnr(inputvalue)
        include 'amrvacdef.f'
        integer :: set_strictnr
        logical :: inputvalue
        strictnr = inputvalue
        set_strictnr = 0   
    end function

    function get_strictnr(outputvalue)
        include 'amrvacdef.f'
        integer :: get_strictnr
        logical :: outputvalue
        outputvalue = strictnr
        get_strictnr = 0   
    end function

    function set_strictsmall(inputvalue)
        include 'amrvacdef.f'
        integer :: set_strictsmall
        logical :: inputvalue
        strictsmall = inputvalue
        set_strictsmall = 0   
    end function

    function get_strictsmall(outputvalue)
        include 'amrvacdef.f'
        integer :: get_strictsmall
        logical :: outputvalue
        outputvalue = strictsmall
        get_strictsmall = 0   
    end function

    function set_strictzero(inputvalue)
        include 'amrvacdef.f'
        integer :: set_strictzero
        logical :: inputvalue
        strictzero = inputvalue
        set_strictzero = 0   
    end function

    function get_strictzero(outputvalue)
        include 'amrvacdef.f'
        integer :: get_strictzero
        logical :: outputvalue
        outputvalue = strictzero
        get_strictzero = 0   
    end function

    function set_strictgetaux(inputvalue)
        include 'amrvacdef.f'
        integer :: set_strictgetaux
        logical :: inputvalue
        strictgetaux = inputvalue
        set_strictgetaux = 0   
    end function

    function get_strictgetaux(outputvalue)
        include 'amrvacdef.f'
        integer :: get_strictgetaux
        logical :: outputvalue
        outputvalue = strictgetaux
        get_strictgetaux = 0   
    end function

    function set_usecovariant(inputvalue)
        include 'amrvacdef.f'
        integer :: set_usecovariant
        logical :: inputvalue
        usecovariant = inputvalue
        set_usecovariant = 0   
    end function

    function get_usecovariant(outputvalue)
        include 'amrvacdef.f'
        integer :: get_usecovariant
        logical :: outputvalue
        outputvalue = usecovariant
        get_usecovariant = 0   
    end function

    function set_nocartesian(inputvalue)
        include 'amrvacdef.f'
        integer :: set_nocartesian
        logical :: inputvalue
        nocartesian = inputvalue
        set_nocartesian = 0   
    end function

    function get_nocartesian(outputvalue)
        include 'amrvacdef.f'
        integer :: get_nocartesian
        logical :: outputvalue
        outputvalue = nocartesian
        get_nocartesian = 0   
    end function

    function set_tfix(inputvalue)
        include 'amrvacdef.f'
        integer :: set_tfix
        logical :: inputvalue
        tfix = inputvalue
        set_tfix = 0   
    end function

    function get_tfix(outputvalue)
        include 'amrvacdef.f'
        integer :: get_tfix
        logical :: outputvalue
        outputvalue = tfix
        get_tfix = 0   
    end function

    function set_convert(inputvalue)
        include 'amrvacdef.f'
        integer :: set_convert
        logical :: inputvalue
        convert = inputvalue
        set_convert = 0   
    end function

    function get_convert(outputvalue)
        include 'amrvacdef.f'
        integer :: get_convert
        logical :: outputvalue
        outputvalue = convert
        get_convert = 0   
    end function

    function set_saveprim(inputvalue)
        include 'amrvacdef.f'
        integer :: set_saveprim
        logical :: inputvalue
        saveprim = inputvalue
        set_saveprim = 0   
    end function

    function get_saveprim(outputvalue)
        include 'amrvacdef.f'
        integer :: get_saveprim
        logical :: outputvalue
        outputvalue = saveprim
        get_saveprim = 0   
    end function

    function set_uselimiter(inputvalue)
        include 'amrvacdef.f'
        integer :: set_uselimiter
        logical :: inputvalue
        uselimiter = inputvalue
        set_uselimiter = 0   
    end function

    function get_uselimiter(outputvalue)
        include 'amrvacdef.f'
        integer :: get_uselimiter
        logical :: outputvalue
        outputvalue = uselimiter
        get_uselimiter = 0   
    end function
    
    
    function set_parameters_filename(path)
        include 'amrvacdef.f'
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
        include 'amrvacdef.f'
        integer :: get_parameters_filename
        character(len=1024), intent(out) :: path
        path = parameters_filename 
        get_parameters_filename = 0   
    end function
    
    function get_local_index_of_grid(global_index_of_grid)
        use mod_forest
        
        include 'amrvacdef.f'
        
        integer :: get_local_index_of_grid
        integer :: iigrid
        integer, intent(in) :: global_index_of_grid
        
        if(igrid_inuse(global_index_of_grid, mype)) then
            do iigrid=1,igridstail
               if (igrids(iigrid) .EQ. global_index_of_grid) then
                  get_local_index_of_grid = iigrid
                  return 
               end if
            end do
        else
            get_local_index_of_grid = 0
            return
        end if
                
    end function
    
    function is_index_valid(i, j, k)
        include 'amrvacdef.f'
        logical :: is_index_valid
        integer, intent(in) :: i,j,k
        
        is_index_valid = .TRUE.
        
        is_index_valid = is_index_valid .AND. (i .GE. 0)
        is_index_valid = is_index_valid .AND. (i .LE. (ixMhi1 - ixMlo1))
        
        is_index_valid = is_index_valid .AND. (j .GE. 0)
        is_index_valid = is_index_valid .AND. (j .LE. (ixMhi2 - ixMlo2))
        
        is_index_valid = is_index_valid .AND. (k .GE. 0)
        is_index_valid = is_index_valid .AND. (k .LE. (ixMhi3 - ixMlo3))
    end function
    
    
    function get_position_of_index(i, j, k, index_of_grid, x, y, z, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_position_of_index
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        double precision, intent(out), dimension(n) :: x, y, z
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid .EQ. 0) then
                x(index) = 0.0
                y(index) = 0.0
                z(index) = 0.0
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                reali = ixMlo1 + i(index)
                realj = ixMlo2 + j(index)
                realk = ixMlo3 + k(index)
                
                x(index) = px(local_index_of_grid)%x(reali,realj,realk,1)
                y(index) = px(local_index_of_grid)%x(reali,realj,realk,2)
                z(index) = px(local_index_of_grid)%x(reali,realj,realk,3)
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(x,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(y,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(z,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, x, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, y, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, z, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
    end function
    
    function get_mesh_size(nx, ny, nz, index_of_grid)
        include 'amrvacdef.f'

        integer :: get_mesh_size
        integer, intent(out) :: nx, ny, nz
        integer, intent(in) :: index_of_grid
        
        nx = ixMhi1 - ixMlo1 + 1
        ny = ixMhi2 - ixMlo2 + 1
        nz = ixMhi3 - ixMlo3 + 1
    
        get_mesh_size = 0
    end function

    function get_number_of_grids(number_of_grids)
        use mod_forest
        
        include 'amrvacdef.f'
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
        
        include 'amrvacdef.f'
        integer, intent(out) :: level
        integer, intent(in) :: index_of_grid
        integer :: get_level_of_grid
        
        if(index_of_grid > ngridshi) then
            level = 0
            get_level_of_grid = -1
            return
        end if
        
        level = node(plevel_, index_of_grid)
        get_level_of_grid = 0
    end function
    
    function get_cell_size_of_grid(dx1, dx2, dx3, index_of_grid)
        use mod_forest
        
        include 'amrvacdef.f'
        double precision, intent(out) :: dx1, dx2, dx3
        integer, intent(in) :: index_of_grid
        integer :: level
        integer :: get_cell_size_of_grid
        
        dx1 = 0
        dx2 = 0
        dx3 = 0
        if(index_of_grid > ngridshi) then
            get_cell_size_of_grid = -1
            return
        end if
        
        level = node(plevel_, index_of_grid)
        
        if(level > 0) then
            dx1 = dx(1,level)
            dx2 = dx(2,level)
            dx3 = dx(3,level)
            get_cell_size_of_grid = 0
        else
            get_cell_size_of_grid = -2
        end if
        
    end function

    
    function setup_mesh(nmeshx, nmeshy, nmeshz, xlength, ylength, zlength)
        include 'amrvacdef.f'
        integer, intent(in) :: nmeshx, nmeshy, nmeshz
        double precision, intent(in) :: xlength, ylength, zlength
        integer :: setup_mesh
                
        xprobmin1 = 0.0
        xprobmin2 = 0.0
        xprobmin3 = 0.0
        xprobmax1 = xlength
        xprobmax2 = ylength
        xprobmax3 = zlength
        
        if(nmeshx < 1 .OR.  mod(nmeshx,2)/=0) then
            setup_mesh = -1
            return
        end if
        if(nmeshy < 1 .OR.  mod(nmeshy,2)/=0) then
            setup_mesh = -1
            return
        end if
        if(nmeshz < 1 .OR.  mod(nmeshz,2)/=0) then
            setup_mesh = -1
            return
        end if
        
        
        dx(1,1)=(xprobmax1-xprobmin1)/dble(nmeshx)
        dx(2,1)=(xprobmax2-xprobmin2)/dble(nmeshy)
        dx(3,1)=(xprobmax3-xprobmin3)/dble(nmeshz)
        
        setup_mesh = 0
    end function
    
    
    function get_grid_density(i, j, k, index_of_grid, rho, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_density
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: rho
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid .EQ. 0) then
            
                rho(index) = 0.0
                
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                reali = ixMlo1 + i(index)
                realj = ixMlo2 + j(index)
                realk = ixMlo3 + k(index)
                
                rho(index) = pw(local_index_of_grid)%w(reali,realj,realk,1)
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(rho,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, rho, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        get_grid_density = 0
    end function
    
    function get_grid_momentum_density(i, j, k, index_of_grid, m1, m2, m3, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_momentum_density
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: m1, m2, m3
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid .EQ. 0) then
                m1(index) = 0.0
                m2(index) = 0.0
                m3(index) = 0.0
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                reali = ixMlo1 + i(index)
                realj = ixMlo2 + j(index)
                realk = ixMlo3 + k(index)
                
                m1(index) = pw(local_index_of_grid)%w(reali,realj,realk,2)
                m2(index) = pw(local_index_of_grid)%w(reali,realj,realk,3)
                m3(index) = pw(local_index_of_grid)%w(reali,realj,realk,4)
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(m1,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m2,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(m3,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, m1, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m2, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
            call MPI_Reduce(MPI_IN_PLACE, m3, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_grid_momentum_density = 0
    end function
    
    function get_grid_energy_density(i, j, k, index_of_grid, en, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: get_grid_energy_density
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(out), dimension(n) :: en
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid .EQ. 0) then
            
                en(index) = 0.0
                
            else
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                end if
                reali = ixMlo1 + i(index)
                realj = ixMlo2 + j(index)
                realk = ixMlo3 + k(index)
                
                en(index) = pw(local_index_of_grid)%w(reali,realj,realk,5)
            end if
            
        end do
                
        if(mype .GT. 0) then
            call MPI_Reduce(en,  0, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        else
            call MPI_Reduce(MPI_IN_PLACE, en, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierrmpi)
        end if
        
        get_grid_energy_density = 0
    end function
    
    
    
    function set_grid_density(i, j, k, rho, index_of_grid, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_density
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: rho
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                    reali = ixMlo1 + i(index)
                    realj = ixMlo2 + j(index)
                    realk = ixMlo3 + k(index)
                    
                    pw(local_index_of_grid)%w(reali,realj,realk,1) = rho(index)
                
                end if
            end if
            
        end do
        
        set_grid_density = 0
    end function

    function set_grid_energy_density(i, j, k, en, index_of_grid, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_energy_density
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: en
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                    reali = ixMlo1 + i(index)
                    realj = ixMlo2 + j(index)
                    realk = ixMlo3 + k(index)
                    
                    pw(local_index_of_grid)%w(reali,realj,realk,5) = en(index)
                
                end if
            end if
            
        end do
        
        set_grid_energy_density = 0
    end function
    
    function set_grid_momentum_density(i, j, k, m1, m2, m3, index_of_grid, n)
        include 'amrvacdef.f'
        
        integer :: local_index_of_grid, index, previous_index_of_grid = -1
        integer :: set_grid_momentum_density
        integer :: reali, realj, realk
        
        integer, intent(in) :: n
        integer, intent(in), dimension(n) :: i, j, k, index_of_grid
        
        double precision, intent(in), dimension(n) :: m1, m2, m3
        
        local_index_of_grid = 0
        
        do index = 1,n
            if ( previous_index_of_grid .NE. index_of_grid(index)) then
                local_index_of_grid = get_local_index_of_grid(index_of_grid(index))
            end if 
            
            if (local_index_of_grid /= 0) then
                if(.NOT. is_index_valid(i(index), j(index), k(index))) then
                
                else
                    
                    reali = ixMlo1 + i(index)
                    realj = ixMlo2 + j(index)
                    realk = ixMlo3 + k(index)
                    
                    pw(local_index_of_grid)%w(reali,realj,realk,2) = m1(index)
                    pw(local_index_of_grid)%w(reali,realj,realk,3) = m2(index)
                    pw(local_index_of_grid)%w(reali,realj,realk,4) = m3(index)
                
                end if
            end if
            
        end do
        
        set_grid_momentum_density = 0
    end function
    
    
    function evolve(tend) 
        include 'amrvacdef.f'
        
        integer :: evolve
        double precision :: tend

        integer :: level, ifile
        double precision :: time_in, timeio0, timeio_tot, timegr0, timegr_tot,&
            timeloop, timeloop0
       
        
        time_advance=.true. 
        time_accurate=.true.
        
        
        time_in=MPI_WTIME()
        timeio_tot=zero
        timegr_tot=zero

        

        itmin=it

        call getbc(t,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,pw,pwCoarse,pgeo,&
           pgeoCoarse,.false.)

        !  ------ start of integration loop. ------------------
        timeloop0=MPI_WTIME()
        if (mype==0) then
           write(*,'(a,f12.3,a)')'BCs before Advance took    : ',timeloop0&
              -time_in,' sec'
        end if

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
    
        evolve = 0
    end function

    function get_time(tnow)
        include 'amrvacdef.f'
        
        integer :: get_time
        double precision :: tnow

        tnow = t
        
        get_time = 0
    end function
    
    function set_boundary(lowx,highx,lowy,highy,lowz,highz)
        include 'amrvacdef.f'
        
        integer :: set_boundary
        integer :: idim
        character(len=25), intent(in) :: lowx,highx,lowy,highy,lowz,highz
        
        typeB(1:nw,1) = lowx
        typeB(1:nw,2) = highx
        typeB(1:nw,3) = lowy
        typeB(1:nw,4) = highy
        typeB(1:nw,5) = lowz
        typeB(1:nw,6) = highz
        
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


    
END MODULE mpiamrvac_interface

