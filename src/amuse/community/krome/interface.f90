function initialize_code() result(ret)
  use chem_mod
  integer :: ret
  ret=chem_initialize()
end function

function cleanup_code() result(ret)
  use chem_mod
  integer :: ret
  ret=chem_end()
end function

function commit_particles() result(ret)
  use chem_mod
  integer :: ret
  ret=chem_commit_particles()
end function

function recommit_particles() result(ret)
  use chem_mod
  integer :: ret
  ret=chem_commit_particles()
end function

function commit_parameters() result(ret)
  use chem_mod
  integer :: ret
  ret=chem_commit_parameters()
end function

function recommit_parameters() result(ret)
  use chem_mod
  integer :: ret
  ret=chem_commit_parameters()
end function

function get_number_of_particles(n) result(ret)
  use chem_mod
  integer n,ret
  n=nparticle
  ret=0
end function

function new_particle(id,density,temperature,ionrate) result(ret)
  use chem_mod
  integer :: ret,id
  double precision :: density,temperature,ionrate 
  ret=add_particle(id,density,temperature,ionrate)
end function

function set_state(id,density,temperature,ionrate) result(ret)   
  use chem_mod
  integer :: ret,id
  double precision :: density,temperature,ionrate
  ret=set_particle_state(id,density,temperature,ionrate)
end function

function get_state(id,density,temperature,ionrate) result(ret)   
  use chem_mod
  integer :: ret,id
  double precision :: density,temperature,ionrate
  ret=get_particle_state(id,density,temperature,ionrate)
end function

function get_abundance(id,aid,x) result(ret)
  use chem_mod
  integer ret,id,aid
  double precision x
  ret=get_particle_abundance(id,aid,x)
end function

function set_abundance(id,aid,x) result(ret)
  use chem_mod
  integer ret,id,aid
  double precision x
  ret=set_particle_abundance(id,aid,x)
end function

function get_firstlast_abundance(first,last) result(ret)
  use chem_mod
  integer ret,first,last
  first=1
  last=krome_nmols ! this is the last species defined in amuse_helpers.f90
  ret=0
end function

function get_name_of_species(i,s) result(ret)
  use chem_mod
  integer ret,i
  character*16 :: ss(krome_nmols),s
  if(i.LT.1.OR.i.GT.krome_nmols) then
    ret=-1
    return
  endif
  ss=krome_get_names()
  s=ss(i)
  ret=0
end function

function get_index_of_species(s,i) result(ret)
  use chem_mod
  integer ret,i
  character*16 s
  i=krome_get_index(s)
  ret=0
end function

function evolve_model(tend) result(ret)
  use chem_mod
  integer :: ret
  double precision :: tend  
  ret=evolve_chem_model(tend)
end function

function get_time(outtime) result(ret)
  use chem_mod
  integer :: ret
  double precision :: outtime  
  ret=chem_model_time(outtime)
end function

function delete_particle(id) result(ret)
  use chem_mod
  integer :: ret,id
  ret=remove_particle(id)
end function
