module amuseinterface
  use map_helpers
  
contains

function initialize_code() result(ret)
  integer :: ret
  ret=0
end function

function generate_projection() result(ret)
  integer :: ret
  call map_generate_projection()
  ret=0
end function

function init_map() result(ret)
  integer :: ret
  ret=map_init()
end function

function reset_map() result(ret)
  integer :: ret
  call map_reset()
  ret=0
end function

function erase_map() result(ret)
  integer :: ret
  call map_erase()
  ret=0
end function    

function new_particle(id,m,x,y,z,r,o,i) result(ret)
  integer :: ret,i,id
  real :: m,x,y,z,r,o
  ret=add_particle(id,m,x,y,z,r,o,i)
end function

function set_state(id,m,x,y,z,r,o) result(ret)
  integer :: ret,id
  real :: m,x,y,z,r,o
  ret=set_particle_state(id,m,x,y,z,r,o)
end function

function get_state(id,m,x,y,z,r,o) result(ret)
  integer :: ret,id
  real :: m,x,y,z,r,o
  ret=get_particle_state(id,m,x,y,z,r,o)
end function


function delete_particle(id) result(ret)
  integer ret,id
  ret=map_remove_particle(id)
end function  

function set_random_seed(x) result(ret)
  integer :: x
  integer :: ret
  ret=map_set_random_seed(x)
end function
function get_random_seed(x) result(ret)
  integer :: x
  integer :: ret
  ret=map_get_random_seed(x)
end function

function set_minimum_distance(x) result(ret)
  real :: x
  integer :: ret
  call map_set_zm(x)
  ret=0
end function
function get_minimum_distance(x) result(ret)
  real :: x
  integer :: ret
  call map_get_zm(x)
  ret=0
end function

function set_extinction_flag(x) result(ret)
  integer :: x
  integer :: ret
  call map_set_ext(x)
  ret=0
end function
function get_extinction_flag(x) result(ret)
  integer :: x
  integer :: ret
  call map_get_ext(x)
  ret=0
end function

function set_image_angle(x) result(ret)
  real :: x
  integer :: ret
  call map_set_angle(x)
  ret=0
end function
function get_image_angle(x) result(ret)
  real :: x
  integer :: ret
  call map_get_angle(x)
  ret=0
end function

function set_image_width(x) result(ret)
  real :: x
  integer :: ret
  call map_set_width(x)
  ret=0
end function
function get_image_width(x) result(ret)
  real :: x
  integer :: ret
  call map_get_width(x)
  ret=0
end function

function set_image_pixel_size(nx,ny) result(ret)
  integer :: nx,ny
  integer :: ret
  call map_set_imsize(nx,ny)
  ret=0
end function
function get_image_pixel_size(nx,ny) result(ret)
  integer :: nx,ny
  integer :: ret
  call map_get_imsize(nx,ny)
  ret=0
end function

function set_image_target(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_set_focus(x,y,z)
  ret=0
end function
function get_image_target(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_get_focus(x,y,z)
  ret=0
end function

function set_viewpoint(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_set_viewpoint(x,y,z)
  ret=0
end function
function get_viewpoint(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_get_viewpoint(x,y,z)
  ret=0
end function

function set_upvector(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_set_upvector(x,y,z)
  ret=0
end function
function get_upvector(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_get_upvector(x,y,z)
  ret=0
end function

function set_projection_direction(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_set_direction(x,y,z)
  ret=0
end function
function get_projection_direction(x,y,z) result(ret)
  real :: x,y,z
  integer :: ret
  call map_get_direction(x,y,z)
  ret=0
end function

function set_projection_mode(x) result(ret)
  integer :: ret
  character(len=15), intent(in) :: x
  if(x.NE."parallel".AND.x.NE."perspective") then
    ret=-1
    return
  endif
  call map_set_projection_mode(x)
  ret=0
end function
function get_projection_mode(x) result(ret)
  integer :: ret
  character(len=15), intent(out) :: x
  call map_get_projection_mode(x)
  ret=0
end function



function get_image(i,j,pic,n) result(ret)
  integer n,ret,i(n),j(n)
  real :: pic(n)
  ret=map_get_pic(i,j,pic,n)
end function  

function get_opdepth_map(i,j,pic,n) result(ret)
  integer n,ret,i(n),j(n)
  real :: pic(n)
  ret=map_get_opdepth(i,j,pic,n)
end function

end module

