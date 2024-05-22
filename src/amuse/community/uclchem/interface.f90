  function initialize_code() result(ret)
    use uclchemhelper
    integer :: ret
    ret=chem_initialize()
  end function
  
  function cleanup_code() result(ret)
    use uclchemhelper
    integer :: ret
    ret=chem_end()
  end function
  
  ! function commit_particles() result(ret)
  !   use uclchemhelper
  !   integer :: ret
  !   ret=chem_commit_particles()
  ! end function
  
  ! function recommit_particles() result(ret)
  !   use uclchemhelper
  !   integer :: ret
  !   ret=chem_commit_particles()
  ! end function
  
  ! function commit_parameters() result(ret)
  !   use uclchemhelper
  !   integer :: ret
  !   ret=chem_commit_parameters()
  ! end function
  
  ! function recommit_parameters() result(ret)
  !   use uclchemhelper
  !   integer :: ret
  !   ret=chem_commit_parameters()
  ! end function
  
  ! function get_number_of_particles(n) result(ret)
  !   use uclchemhelper
  !   integer n,ret
  !   n=nparticle
  !   ret=0
  ! end function
  
  function new_particle(id,dens,temperature,ionrate) result(ret)
    use uclchemhelper
    integer :: ret,id
    double precision :: dens,temperature,ionrate 
    print *, 'fortran'
    ret=add_particle(id,dens,temperature,ionrate)
  end function
  
  function set_state(id,dens,temperature,ionrate) result(ret)   
    use uclchemhelper
    integer :: ret,id
    double precision :: dens,temperature,ionrate
    ret=set_particle_state(id,dens,temperature,ionrate)
  end function
  
  function get_state(id,dens,temperature,ionrate) result(ret)   
    use uclchemhelper
    integer :: ret,id
    double precision :: dens,temperature,ionrate
    ret=get_particle_state(id,dens,temperature,ionrate)
  end function
  
  function get_abundance(id,aid,x) result(ret)
    use uclchemhelper
    integer :: ret,id,aid
    double precision :: x
    ret=get_particle_abundance(id,aid,x)
  end function
  
  function set_abundance(id,aid,x) result(ret)
    use uclchemhelper
    integer ret,id,aid
    double precision x
    ret=set_particle_abundance(id,aid,x)
  end function
  
  function get_firstlast_abundance(first,last) result(ret)
   use network
   integer :: ret,first,last
   first=1
   last=nSpec 
   ret=0
  end function
  
  ! !function get_name_of_species(index,s) result(ret)
  ! !  use uclchemhelper
  ! !  integer ret,index
  ! !  character*16 :: ss(krome_nmols),s
  ! !  if(index.LT.1.OR.index.GT.krome_nmols) then
  ! !    ret=-1
  ! !    return
  ! !  endif
  ! !  ss=krome_get_names()
  ! !  s=ss(index)
  ! !  ret=0
  ! !end function
  
  ! !function get_index_of_species(s,index) result(ret)
  ! !  use uclchemhelper
  ! !  integer ret,index
  ! !  character*16 s
  ! !  i=krome_get_index(s)
  ! !  ret=0
  ! !end function
  
  ! function evolve_model(tend) result(ret)
  !   use uclchemhelper
  !   integer :: ret
  !   double precision :: tend  
  !   ret=evolve_chem_model(tend)
  ! end function
  
  ! function get_time(outtime) result(ret)
  !   use uclchemhelper
  !   integer :: ret
  !   double precision :: outtime  
  !   ret=chem_model_time(outtime)
  ! end function
  
  function delete_particle(id) result(ret)
    use uclchemhelper
    integer :: ret,id
    ret=remove_particle(id)
  end function


  function sim_cloud(outSpeciesIn, dictionary, abundance_out) result(successFlag)
    use uclchemhelper
    DOUBLE PRECISION :: abundance_out
    character(len=20) :: outSpeciesIn
    character(len=20) :: dictionary
    INTEGER :: successFlag
    !abundance(:) = 0.00
    !outSpeciesIn = 'H H2'
    !dictionary = "{'outspecies': 2}"
    !abundance_out(:) = 0.0_dp
    print *, 'check1'
    !do i = 1, 500
    !  abundance_out(i) = 0.0_dp
    !end do
    !if (allocated(abundance_out)) then
    !  deallocate(abundance_out)
    !end if

    print *,'check2'
    !allocate(abundance_out(500),source=0.0_dp)
    abundance_out=14
    print *, 'check3'
    !call test_cloud(dictionary, outSpeciesIn,abundance_out,successFlag)
    write(*,*) abundance_out
    successFlag=0
  end function
  