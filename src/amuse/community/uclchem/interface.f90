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
  
  function commit_particles() result(ret)
    use uclchemhelper
    integer :: ret
    ret=chem_commit_particles()
  end function
  
  function recommit_particles() result(ret)
    use uclchemhelper
    integer :: ret
    ret=chem_commit_particles()
  end function
  
  function commit_parameters() result(ret)
    use uclchemhelper
    integer :: ret
    ret=chem_commit_parameters()
  end function
  
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
  
  function get_name_of_species(index,s) result(ret)
   use network
   integer ret,index
   character*16 :: ss(nSpec),s
   index = index + 1
   if(index.LT.1.OR.index.GT.nSpec) then
     ret=-1
     return
   endif
   ss=specname
   s=ss(index)
   ret=0
  end function
  
  function get_index_of_species(s,index) result(ret)
   use network
   integer ret,index
   character*16 s
   if (any(specname == s)) then
    index = FINDLOC(specname, s, dim=1) - 1
    ret = 0
   else 
    ret = 1
   end if 
  end function
  
  function run_model(dictionary, out_species) result(ret)
    use uclchemhelper
    integer :: ret 
    character(len=*) :: out_species
    character(len=*) :: dictionary(nparticle)
    ret=simple_evolution(dictionary, out_species)
  end function
  
  function get_time(time) result(ret)
    use uclchemhelper
    integer :: ret
    double precision :: time  
    ret=get_current_time(time)
  end function
  
  function delete_particle(id) result(ret)
    use uclchemhelper
    integer :: ret,id
    ret=remove_particle(id)
  end function
