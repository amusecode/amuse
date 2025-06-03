! This file contains routines to make a visualization of a nbody simulation
! it is self contained and could be called from a running simulation.
! the data is passed through common blocks
 
!


module view_modifier

! This module provides facilities to modify the view in an OpenGL window.
! The mouse buttons and keyboard arrow keys can be used to zoom, pan,
! rotate and change the scale.  A menu or submenu can be used to select which
! buttons perform which function and to reset the view to the initial settings.
! This is limited to one window.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! April, 1998

! To use this module:
!
! 1) put a USE view_modifier statement in any program unit that calls a
!    procedure in this module
!
! 2) set the initial operation assignments, view and scale below the
!    "Initial configuration" comment below
!
! 3) call view_modifier_init after glutCreateWindow
!    This is a function that returns integer(kind=glcint) menuid.  The menuid
!    is the ID returned by glutCreateMenu.  You can either use the view_modifier
!    menu as your menu by calling glutAttachMenu immediately after
!    view_modifier_init, as in
!       menuid = view_modifier_init()
!       call glutAttachMenu(GLUT_RIGHT_BUTTON)
!    or by using the menuid to attach a submenu to your own menu, as in
!       call glutAddSubMenu("View Modifier",menuid)
!
! 4) in any callback functions that update the display, put
!       call reset_view
!    as the first executable statement
!
! Note that view_modifier_init sets the callback functions for glutMouseFunc,
! glutMotionFunc and glutSpecialFunc, so don't call these yourself
!
! The menu allows you to select what operation is attached to the left and
! middle mouse buttons and arrow keys, reset to the initial view, and quit.
! The right mouse button should be used for the menu.

use opengl_gl
use opengl_glu
use opengl_glut
implicit none
private
public :: view_modifier_init, reset_view
private :: ZOOM, PAN, ROTATE, SCALEX, SCALEY, SCALEZ, RESET, QUIT, &
           PI, &
           left_button_func, middle_button_func, arrow_key_func, &
           init_lookat, init_lookfrom, &
           init_xscale_factor, init_yscale_factor, init_zscale_factor, &
           angle, shift, xscale_factor, yscale_factor, zscale_factor, &
           moving_left, moving_middle, begin_left, begin_middle, &
           cart2sphere, sphere2cart, cart3D_plus_cart3D, cart3D_minus_cart3D, &
           reset_to_init, mouse, motion, arrows, &
           menu_handler, set_left_button, set_middle_button, set_arrow_keys

integer(kind=glcint), parameter :: ZOOM = 1, PAN = 2, ROTATE = 3, SCALEX = 4, &
                      SCALEY = 5, SCALEZ = 6
integer(kind=glcint), parameter :: RESET = 10, QUIT = 11
real(kind=gldouble), parameter :: PI = 3.141592653589793_gldouble

type, private :: cart2D ! 2D cartesian coordinates
   real(kind=gldouble) :: x, y
end type cart2D

type, private :: cart3D ! 3D cartesian coordinates
   real(kind=gldouble) :: x, y, z
end type cart3D

type, private :: sphere3D ! 3D spherical coordinates
   real(kind=gldouble) :: theta, phi, rho
end type sphere3D

type(cart2D), save :: angle
type(cart3D), save :: shift
real(kind=gldouble), save :: xscale_factor, yscale_factor, zscale_factor
logical, save :: moving_left, moving_middle
type(cart2D), save :: begin_left, begin_middle

interface operator(+)
   module procedure cart3D_plus_cart3D
end interface
interface operator(-)
   module procedure cart3D_minus_cart3D
end interface

! ------- Initial configuration -------

! Set the initial operation performed by each button and the arrow keys.
! The operations are ZOOM, PAN, ROTATE, SCALEX, SCALEY, and SCALEZ

integer, save ::   left_button_func = ROTATE, &
                 middle_button_func = ZOOM, &
                     arrow_key_func = PAN

! Set the initial view as the point you are looking at, the point you are
! looking from, and the scale factors

type(cart3D), parameter :: &
     init_lookat = cart3D(0.0_gldouble, 0.0_gldouble, 0.0_gldouble), &
   init_lookfrom = cart3D(6.0_gldouble, -12.0_gldouble, 3.0_gldouble)

real(kind=gldouble), save:: &
   init_xscale_factor = 1.0_gldouble, &
   init_yscale_factor = 1.0_gldouble, &
   init_zscale_factor = 1.0_gldouble

! -------- end of Initial configuration ------

contains

!          -------------
subroutine reset_to_init
!          -------------

! This resets the view to the initial configuration

type(sphere3D) :: slookfrom

slookfrom = cart2sphere(init_lookfrom-init_lookat)
angle%x = -180.0_gldouble*slookfrom%theta/PI - 90.0_gldouble
angle%y = -180.0_gldouble*slookfrom%phi/PI
shift%x = 0.0_gldouble
shift%y = 0.0_gldouble
shift%z = -slookfrom%rho
xscale_factor = init_xscale_factor
yscale_factor = init_yscale_factor
zscale_factor = init_zscale_factor

 call glutPostRedisplay

return
end subroutine reset_to_init

!          ----------
subroutine reset_view
!          ----------

! This routine resets the view to the current orientation and scale

 call glMatrixMode(GL_MODELVIEW)
 call glPopMatrix
 call glPushMatrix
 call glTranslated(shift%x, shift%y, shift%z)
 call glRotated(angle%x, 0.0_gldouble, 0.0_gldouble, 1.0_gldouble)
 call glRotated(angle%y, cos(PI*angle%x/180.0_gldouble), &
               -sin(PI*angle%x/180.0_gldouble), 0.0_gldouble)
 call glTranslated(-init_lookat%x, -init_lookat%y, -init_lookat%z)
 call glScaled(xscale_factor,yscale_factor,zscale_factor)

return
end subroutine reset_view

!          -----
subroutine mouse(button, state, x, y)
!          -----
integer(kind=glcint), intent(in out) :: button, state, x, y

! This gets called when a mouse button changes
 
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_DOWN) then
    moving_left = .true.
    begin_left = cart2D(x,y)
  endif
  if (button == GLUT_LEFT_BUTTON .and. state == GLUT_UP) then
    moving_left = .false.
  endif
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_DOWN) then
    moving_middle = .true.
    begin_middle = cart2D(x,y)
  endif
  if (button == GLUT_MIDDLE_BUTTON .and. state == GLUT_UP) then
    moving_middle = .false.
  endif
end subroutine mouse

!          ------
subroutine motion(x, y)
!          ------
integer(kind=glcint), intent(in out) :: x, y

! This gets called when the mouse moves

integer :: button_function
type(cart2D) :: begin
real(kind=gldouble) :: factor

! Determine and apply the button function

if (moving_left) then
   button_function = left_button_func
   begin = begin_left
else if(moving_middle) then
   button_function = middle_button_func
   begin = begin_middle
end if

select case(button_function)
case (ZOOM)
   if (y < begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(begin%y-y))
   else if (y > begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(y-begin%y)
   else
      factor = 1.0_gldouble
   end if
   shift%z = factor*shift%z
case (PAN)
   shift%x = shift%x + .01*(x - begin%x)
   shift%y = shift%y - .01*(y - begin%y)
case (ROTATE)
   angle%x = angle%x + (x - begin%x)
   angle%y = angle%y + (y - begin%y)
case (SCALEX)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   xscale_factor = xscale_factor * factor
case (SCALEY)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   yscale_factor = yscale_factor * factor
case (SCALEZ)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   zscale_factor = zscale_factor * factor
end select

! update private variables and redisplay

if (moving_left) then
   begin_left = cart2D(x,y)
else if(moving_middle) then
   begin_middle = cart2D(x,y)
endif

if (moving_left .or. moving_middle) then
   call glutPostRedisplay
endif

return
end subroutine motion

!          ------
subroutine arrows(key, x, y)
!          ------
integer(glcint), intent(in out) :: key, x, y

! This routine handles the arrow key operations

real(kind=gldouble) :: factor

select case(arrow_key_func)
case(ZOOM)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble + .02_gldouble
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case default
      factor = 1.0_gldouble
   end select
   shift%z = factor*shift%z
case(PAN)
   select case(key)
   case(GLUT_KEY_LEFT)
      shift%x = shift%x - .02
   case(GLUT_KEY_RIGHT)
      shift%x = shift%x + .02
   case(GLUT_KEY_DOWN)
      shift%y = shift%y - .02
   case(GLUT_KEY_UP)
      shift%y = shift%y + .02
   end select
case(ROTATE)
   select case(key)
   case(GLUT_KEY_LEFT)
      angle%x = angle%x - 1.0_gldouble
   case(GLUT_KEY_RIGHT)
      angle%x = angle%x + 1.0_gldouble
   case(GLUT_KEY_DOWN)
      angle%y = angle%y + 1.0_gldouble
   case(GLUT_KEY_UP)
      angle%y = angle%y - 1.0_gldouble
   end select
case(SCALEX)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   xscale_factor = xscale_factor * factor
case(SCALEY)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   yscale_factor = yscale_factor * factor
case(SCALEZ)
   select case(key)
   case(GLUT_KEY_DOWN)
      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
   case(GLUT_KEY_UP)
      factor = 1.0_gldouble + .02_gldouble
   case default
      factor = 1.0_gldouble
   end select
   zscale_factor = zscale_factor * factor

end select
   
call glutPostRedisplay

return
end subroutine arrows

!          ------------
subroutine menu_handler(value)
!          ------------
integer(kind=glcint), intent(in out) :: value

! This routine handles the first level entries in the menu

select case(value)

case(RESET)
   call reset_to_init
case(QUIT)
   call glutLeaveMainLoop()

end select

return
end subroutine menu_handler

!          ---------------
subroutine set_left_button(value)
!          ---------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the left button as given by menu selection

left_button_func = value

return
end subroutine set_left_button

!          -----------------
subroutine set_middle_button(value)
!          -----------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the middle button as given by menu selection

middle_button_func = value

return
end subroutine set_middle_button

!          --------------
subroutine set_arrow_keys(value)
!          --------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the arrow keys as given by menu selection

arrow_key_func = value

return
end subroutine set_arrow_keys

!        ------------------
function view_modifier_init(scale) result(menuid)
!        ------------------

real(kind=gldouble), optional,intent(in) :: scale

integer(kind=glcint) :: menuid

! This initializes the view modifier variables and sets initial view.
! It should be called immediately after glutCreateWindow

integer(kind=glcint) :: button_left, button_middle, arrow_keys

if(present(scale)) then
init_xscale_factor=scale
init_yscale_factor=scale
init_zscale_factor=scale
endif

! set the callback functions

call glutMouseFunc(mouse)
call glutMotionFunc(motion)
call glutSpecialFunc(arrows)

! create the menu

button_left = glutCreateMenu(set_left_button)
call glutAddMenuEntry("rotate",ROTATE)
call glutAddMenuEntry("zoom",ZOOM)
call glutAddMenuEntry("pan",PAN)
call glutAddMenuEntry("scale x",SCALEX)
call glutAddMenuEntry("scale y",SCALEY)
call glutAddMenuEntry("scale z", SCALEZ)
button_middle = glutCreateMenu(set_middle_button)
call glutAddMenuEntry("rotate",ROTATE)
call glutAddMenuEntry("zoom",ZOOM)
call glutAddMenuEntry("pan",PAN)
call glutAddMenuEntry("scale x",SCALEX)
call glutAddMenuEntry("scale y",SCALEY)
call glutAddMenuEntry("scale z", SCALEZ)
arrow_keys = glutCreateMenu(set_arrow_keys)
 call glutAddMenuEntry("rotate",ROTATE)
 call glutAddMenuEntry("zoom",ZOOM)
 call glutAddMenuEntry("pan",PAN)
 call glutAddMenuEntry("scale x",SCALEX)
 call glutAddMenuEntry("scale y",SCALEY)
 call glutAddMenuEntry("scale z", SCALEZ)
 menuid = glutCreateMenu(menu_handler)
 call glutAddSubMenu("left mouse button",button_left)
 call glutAddSubMenu("middle mouse button",button_middle)
 call glutAddSubMenu("arrow keys",arrow_keys)
 call glutAddMenuEntry("reset to initial view",RESET)
 call glutAddMenuEntry("quit",QUIT)

! set the perspective

 call glMatrixMode(GL_PROJECTION)
 call gluPerspective(30.0_gldouble, 1.0_gldouble, 0.01_gldouble, 200.0_gldouble)

! set the initial view

 call glPushMatrix
 call reset_to_init

return
end function view_modifier_init

!        -----------
function sphere2cart(spoint) result(cpoint)
!        -----------
type(sphere3D), intent(in) :: spoint
type(cart3D) :: cpoint

! This converts a 3D point from spherical to cartesean coordinates

real(kind=gldouble) :: t,p,r

t=spoint%theta
p=spoint%phi
r=spoint%rho

cpoint%x = r*cos(t)*sin(p)
cpoint%y = r*sin(t)*sin(p)
cpoint%z = r*cos(p)

return
end function sphere2cart

!        -----------
function cart2sphere(cpoint) result(spoint)
!        -----------
type(cart3D), intent(in) :: cpoint
type(sphere3D) :: spoint

! This converts a 3D point from cartesean to spherical coordinates

real(kind=gldouble) :: x,y,z

x=cpoint%x
y=cpoint%y
z=cpoint%z

spoint%rho = sqrt(x*x+y*y+z*z)
if (x==0.0_gldouble .and. y==0.0_gldouble) then
   spoint%theta = 0.0_gldouble
else
   spoint%theta = atan2(y,x)
end if
if (spoint%rho == 0.0_gldouble) then
   spoint%phi = 0.0_gldouble
else
   spoint%phi = acos(z/spoint%rho)
endif

return
end function cart2sphere

!        ------------------
function cart3D_plus_cart3D(cart1,cart2) result(cart3)
!        ------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the sum of two 3D cartesean points

cart3%x = cart1%x + cart2%x
cart3%y = cart1%y + cart2%y
cart3%z = cart1%z + cart2%z

return
end function cart3D_plus_cart3D

!        -------------------
function cart3D_minus_cart3D(cart1,cart2) result(cart3)
!        -------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the difference of two 3D cartesean points

 cart3%x = cart1%x - cart2%x
 cart3%y = cart1%y - cart2%y
 cart3%z = cart1%z - cart2%z

return
end function cart3D_minus_cart3D

end module view_modifier

!********************************************************************

module nbody_viewer
use opengl_gl
use opengl_glu
use opengl_glut
use view_modifier

implicit none

integer, parameter :: NMAX=2000000 

real(kind=gldouble), parameter :: M_PI = 3.14159265358979323846_gldouble

type star
	real(kind=glfloat) x(0:1),y(0:1),z(0:1)
!	real(glfloat),dimension(3) :: &
!	color=(/1.0_glfloat,1.0_glfloat,1.0_glfloat,1.0_glfloat/)
end type star

integer(GLenum) doubleBuffer

logical,save :: show_gas=.TRUE., show_stars=.TRUE., show_halo=.TRUE., &
	 update_flag=.TRUE.,show_age=.FALSE.,show_temp=.FALSE.,     &
	 show_fuvheat=.FALSE.,show_div=.FALSE.,show_hsm=.FALSE.
	 
integer, save :: gas_color,star_color,halo_color	  

integer(glint), save ::framerate=100;

integer(glint), parameter:: FASTUPDATE=100,SLOWUPDATE=5000

real(kind=glfloat),save :: pointsize_gas=3,pointsize_stars=1.5,pointsize_halo=10.

real(kind=glfloat), dimension(3) :: pointpar=(/0.15,0.075,0.026/)

real(kind=gldouble), save :: initialbox
real(kind=gldouble) :: wwidth, wheight

 contains

recursive subroutine starreset(dummy)
	include 'globals.h'
	integer(glint),intent(inout) :: dummy
        real(kind=glfloat) :: xt,yt,zt,temp,red,green,blue,alpha_c
	integer :: i,step
	real(kind=glfloat) :: starscale=10.   ! scaling of brightness of stars with age
	real(kind=glfloat) :: lhmin,lhmax

 step=max(nbodies/NMAX,1)
 call gldeletelists(1,3)
	
 lhmin=log(minval(hsmooth(1:nsph))) 
 lhmax=log(maxval(hsmooth(1:nsph))) 
  
if(show_gas) then

	call glNewList(1,GL_COMPILE)	
	call glPointSize( pointsize_gas)
	if(show_fuvheat.OR.show_temp) then
	call glColor4f(.9_glfloat, .6_glfloat, .6_glfloat,.6_glfloat)
	else
	call glColor4f(1._glfloat, .8_glfloat, .8_glfloat,.6_glfloat)
	endif
	call glBegin(GL_POINTS)
	do i=1,nsph,step
	
	xt=pos(i,1)
	yt=pos(i,2)
	zt=pos(i,3)

	if(show_temp) then
	
	  temp=(alog10(csound(i)/.015))/2.
	  red=.4+.8*temp
	  blue=1.-2*temp*(1-temp)
	  green=.2+.6*temp**2
          alpha_c=.6
!	  alpha_c=.2+.5*temp
	  call glColor4f(red, green, blue,alpha_c)
	
	end if
	
	if(show_fuvheat) then
	
	  if(fuvheat(i).ge.0) then  
	  temp=log(fuvheat(i))-log(7.)
	  red=.2+temp/3
	  blue=1.+temp/3
	  green=.2+temp/3
	  alpha_c=.6	  
	  call glColor4f(red, green, blue,alpha_c)
	  endif
	end if

	if(show_div) then
	
	  red=0.
	  blue=0.
	  green=0.
	  if(hsmdivv(i).GT.0) then
          red=1.
          else
          blue=1.
          endif
          alpha_c=.6
	  call glColor4f(red, green, blue,alpha_c)
	
	end if

	if(show_hsm) then
	
	  green=(log(hsmooth(i))-lhmin)/(lhmax-lhmin)
	  blue=0.
	  red=1-green
    alpha_c=.7-.7*green
	  call glColor4f(red, green, blue,alpha_c)
	
	end if

	
	call glVertex3f(xt,yt,zt)
	
	
	end do
        call glEnd()
	call glEndList()

end if

if(show_stars) then

	call glNewList(2,GL_COMPILE)
	call glPointSize( pointsize_stars)
	call glColor4f(.8_glfloat, .8_glfloat, .6_glfloat,.6_glfloat)
	call glBegin(GL_POINTS)
	do i=nbodies-nstar+1,nbodies,step
	
	xt=pos(i,1)
	yt=pos(i,2)
	zt=pos(i,3)
	
	if(show_age) then
	
	 temp=exp((tform(i)-tnow)/starscale)
	
	red=.7+.3*temp;green=.6+.4*temp;blue=.3+.7*temp
		 
	call glColor4f(red, green, blue,.6_glfloat)
!	call glPointSize( (1._glfloat+3._glfloat*temp))
	end if
		
	call glVertex3f(xt,yt,zt)	
	end do
        call glEnd()
	call glEndList()

end if

! halo
if(show_halo) then

	call glNewList(3,GL_COMPILE)
	call glPointSize( pointsize_halo)
	call glColor4f(.25_glfloat, .25_glfloat, 1._glfloat,.9_glfloat)
	call glBegin(GL_POINTS)
	do i=nsph+1,nbodies-nstar,step
	
	xt=pos(i,1)
	yt=pos(i,2)
	zt=pos(i,3)
	
!	red=.75;blue=.75;green=.75
!	call glColor4f(red, green, blue,.8_glfloat)
		
	
	call glVertex3f(xt,yt,zt)
	
	end do
        call glEnd()
	call glEndList()

end if

 call glutPostRedisplay

 if(update_flag) call glutTimerFunc( framerate, starreset, 0)

end subroutine starreset

function starscale() result(scale_factor)

	include 'globals.h'
	real(kind=gldouble) :: maxcor,maxcorg,maxcors,scale_factor
! schaal:
	maxcor=maxval(abs(pos(1:nbodies,1:3)))	
	if(maxcor.le.0) maxcor=1.
        if(nsph.gt.0) then
	maxcorg=maxval(abs(pos(1:nsph,1:3)))	
	endif
	if(maxcorg.le.0) maxcorg=maxcor	
        if(nstar.gt.0) then
	maxcors=maxval(abs(pos(nbodies-nstar+1:nbodies,1:3)))	
	endif
	if(maxcors.le.0) maxcors=maxcor	
	maxcor=(maxcor*maxcorg*maxcors)**.33333
	scale_factor=2._gldouble/maxcor
	initialbox=2*maxcor
	return
end function starscale

subroutine outputtext(x,y,s)
  real(kind=gldouble) :: x,y
  character :: s*(*)
  character :: c
  integer :: i,lenc
  
  call glrasterpos2d(x,y)
  lenc = len(s)
  do i=1,lenc
    c = s(i:i)
    call glutbitmapcharacter(GLUT_BITMAP_HELVETICA_18, &
          ichar(c))
  end do
end subroutine outputtext

subroutine text
 include 'globals.h'
  character*8 rs
  call glmatrixmode(GL_PROJECTION)
  call glPushMatrix
  call glLoadIdentity
  call glortho(0._gldouble, wwidth, 0._gldouble, wheight,-1._gldouble,1._gldouble)
  call glmatrixmode(GL_MODELVIEW)
  call glPushMatrix
  call glLoadIdentity
  call glColor4f(0.0_glfloat,.5_glfloat,1.0_glfloat,0.5_glfloat)
  if(wheight.gt.120) then
  if(wwidth.gt.310) then
  write(rs,'(f8.2)') tnow*timescale/year/1.e6
  call outputtext(3._gldouble,3._gldouble,'time (Myr):'//rs)
  write(rs,'(f8.2)') initialbox*unitl_in_kpc
  call outputtext(wwidth-156._gldouble,3._gldouble,'box (kpc):'//rs)
  endif

  if(wwidth.gt.200) then
!  write(rs,'(I10)') nbodies-nstar-nsph
  call itos(nbodies-nsph-nstar,8,rs)
  call outputtext(3._gldouble,wheight-20._gldouble,'halo:')
  call outputtext(70._gldouble,wheight-20._gldouble,trim(rs))
  call itos(nsph,8,rs)
  call outputtext(3._gldouble,wheight-40._gldouble,'gas:')
  call outputtext(70._gldouble,wheight-40._gldouble,rs)
  call itos(nstar,8,rs)
  call outputtext(3._gldouble,wheight-60._gldouble,'stars:')
  call outputtext(70._gldouble,wheight-60._gldouble,rs)
  endif
  endif
  call glPopMatrix
  call glmatrixmode(GL_PROJECTION)
  call glPopMatrix
  call glmatrixmode(GL_MODELVIEW)
end subroutine text

subroutine reshape( w,h)
 integer(glint), intent(inout) :: w,h
 real(gldouble) ratio
 wwidth=w
 if(h.GT.1) ratio=w/(1.*h)
 wheight=MAX(h,1)
 call glViewPort(0,0,w,h)
 call glMatrixMode(GL_PROJECTION)
 call glLoadIdentity
 call gluPerspective(30.0_gldouble, ratio, 0.01_gldouble, 200.0_gldouble)

end subroutine reshape

subroutine display

 call reset_view
 call glClearColor(0.02_glfloat,0.02_glfloat,0.02_glfloat, 0.0_glfloat)
 call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
 
 call glBlendFunc(GL_SRC_ALPHA, GL_ONE)
 call glDepthFunc(GL_ALWAYS)
 call glCallList(2)
 call glCallList(1)
 call glCallList(3)
 call glBlendFunc(GL_SRC_ALPHA, GL_ZERO)
 call glDepthFunc(GL_LESS)
 call glColor4f(.3_glfloat, .3_glfloat, .9_glfloat,.9_glfloat)
 call glutWireCube(initialbox)
 call glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA)
! call glDisable(GL_BLEND)
 call glDisable(GL_DEPTH_TEST)
 call text
! call glEnable(GL_BLEND)
 call glEnable(GL_DEPTH_TEST)
 call glutSwapBuffers

end subroutine display

subroutine axes

 call glNewList(3,GL_COMPILE)

! Draw axes so we know the orientation

 call glBegin(GL_LINES)
 call glVertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
 call glVertex3f(2.0_glfloat,0.0_glfloat,0.0_glfloat)
 call glVertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
 call glVertex3f(0.0_glfloat,2.0_glfloat,0.0_glfloat)
 call glVertex3f(0.0_glfloat,0.0_glfloat,0.0_glfloat)
 call glVertex3f(0.0_glfloat,0.0_glfloat,2.0_glfloat)

! Draw crude x, y and z to label the axes

 call glVertex3f(2.1_glfloat,-0.1_glfloat,0.1_glfloat) ! X
 call glVertex3f(2.1_glfloat,0.1_glfloat,-0.1_glfloat)
 call glVertex3f(2.1_glfloat,-0.1_glfloat,-0.1_glfloat)
 call glVertex3f(2.1_glfloat,0.1_glfloat,0.1_glfloat)
 call glVertex3f(0.1_glfloat,2.1_glfloat,0.1_glfloat) ! Y
 call glVertex3f(0.0_glfloat,2.1_glfloat,0.0_glfloat)
 call glVertex3f(-0.1_glfloat,2.1_glfloat,0.1_glfloat)
 call glVertex3f(0.1_glfloat,2.1_glfloat,-0.1_glfloat)
 call glVertex3f(-0.1_glfloat,0.1_glfloat,2.1_glfloat) ! Z
 call glVertex3f(0.1_glfloat,0.1_glfloat,2.1_glfloat)
 call glVertex3f(0.1_glfloat,0.1_glfloat,2.1_glfloat)
 call glVertex3f(-0.1_glfloat,-0.1_glfloat,2.1_glfloat)
 call glVertex3f(-0.1_glfloat,-0.1_glfloat,2.1_glfloat)
 call glVertex3f(0.1_glfloat,-0.1_glfloat,2.1_glfloat)
 call glEnd

 call glEndList

end subroutine axes

subroutine visible(state)
integer(glcint),intent(inout) :: state
  if (state == GLUT_VISIBLE) then
!   if(.not.update_flag) call glutTimerFunc( framerate, starreset,0)
   else
   update_flag=.FALSE. 
  endif
end subroutine visible


subroutine update_menu(optie)
integer(glcint),intent(inout) ::optie
logical ::flag
integer(glint) :: dummy

select case(optie)

case(1)
flag=update_flag
update_flag=.FALSE.
if(.NOT.update_flag) call starreset(dummy)
update_flag=flag
case(2)
framerate=SLOWUPDATE
if(.not.update_flag) then
  update_flag=.TRUE.
  call starreset(dummy)
end if
case(3)
framerate=FASTUPDATE
if(.not.update_flag) then
  update_flag=.TRUE.
  call starreset(dummy)
end if
case(4)
update_flag=.FALSE.
case default 

end select

end subroutine update_menu



subroutine show_menu(optie)
integer(glcint),intent(inout) :: optie
logical :: flag
integer(glint) :: dummy

select case (optie)
case(1)
show_stars=.NOT.show_stars
case(2)
show_gas=.NOT.show_gas
case(3)
show_halo=.NOT.show_halo
case(4)
 if(pointsize_gas.LT.10) then
 pointsize_gas=pointsize_gas*1.5
 pointsize_stars=pointsize_stars*1.5
 pointsize_halo=pointsize_halo*1.5
 endif
case(5)
 if(pointsize_gas.GT.1) then
 pointsize_gas=pointsize_gas/1.5
 pointsize_stars=pointsize_stars/1.5
 pointsize_halo=pointsize_halo/1.5
 endif
case default

end select

flag=update_flag
update_flag=.FALSE.
call starreset(dummy)
update_flag=flag

end subroutine show_menu

subroutine prop_menu(optie)
integer(glcint), intent(inout)::optie
logical :: flag
integer(glint) :: dummy

select case (optie)
case(1)
show_temp=.FALSE.
show_age=.FALSE.
show_fuvheat=.FALSE.
show_div=.FALSE.
show_hsm=.FALSE.
case(2)
show_temp=.NOT.show_temp
case(3)
show_age=.NOT.show_age
case(4)
show_fuvheat=.NOT.show_fuvheat
case(5)
show_div=.NOT.show_div
case(6)
show_hsm=.NOT.show_hsm
case default

end select

flag=update_flag
update_flag=.FALSE.
call starreset(dummy)
update_flag=flag

end subroutine prop_menu

end module nbody_viewer

subroutine viewbodies
 use nbody_viewer
 use view_modifier
 use opengl_gl
 use opengl_glut
 use opengl_glu

implicit none

integer(glcint) :: i,j,k,l
integer(glcint) :: win,menu
 
 if(glutGet(GLUT_ELAPSED_TIME).NE.0) then
!  print*,'already viewer present?'
!  return
 endif
 call glutInit
 call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGBA,GLUT_DEPTH)))
 call glutInitWindowPosition(100,100)
 call glutInitWindowSize(600,350)


win=glutCreateWindow("Stars3D Viewer")

i=glutCreateMenu(update_menu)
 call glutaddmenuentry('refresh view', 1)
 call glutaddmenuentry('slow update', 2)
 call glutaddmenuentry('fast update', 3)
 call glutaddmenuentry('no update', 4)

j=glutCreateMenu(show_menu)
 call glutaddmenuentry('toggle stars', 1)
 call glutaddmenuentry('toggle gas', 2)
 call glutaddmenuentry('toggle halo', 3)
 call glutaddmenuentry('increase pointsize', 4)
 call glutaddmenuentry('decrease pointsize', 5)

k=glutCreateMenu(prop_menu)
 call glutaddmenuentry('plain',1)
 call glutaddmenuentry('toggle temperature', 2)
 call glutaddmenuentry('toggle age', 3)
 call glutaddmenuentry('toggle fuvheat', 4)
 call glutaddmenuentry('toggle divv', 5)
 call glutaddmenuentry('toggle hsm', 6)

menu= view_modifier_init(starscale())

 call glutAttachMenu(GLUT_RIGHT_BUTTON)

 call glutaddsubmenu('update options', i)
 call glutaddsubmenu('view options', j)
 call glutaddsubmenu('show properties',k)

 call glutDisplayFunc(display)
 call glutReshapeFunc(reshape)

 call glClearColor(0.0_glclampf, 0.0_glclampf, 0.0_glclampf, 1.0_glclampf)

 call glEnable(GL_DEPTH_TEST)
 call glEnable(GL_BLEND)
 call glBlendFunc(GL_SRC_ALPHA, GL_ONE)
 call glEnable(GL_LINE_SMOOTH)
 call glEnable(GL_POINT_SMOOTH)
 call glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
 call glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST)
! CALL glPointParameterfvEXT(GL_DISTANCE_ATTENUATION_EXT,pointpar)
! call glPointParameterfEXT(GL_POINT_SIZE_MIN_EXT, 1._glfloat) 
! call glPointParameterfEXT(GL_POINT_FADE_THRESHOLD_SIZE_EX, 0._glfloat) 

 call glutTimerFunc(100, starreset,0)
 call glutVisibilityFunc(visible)

 call glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS)
 call glutMainLoop()

end subroutine viewbodies


subroutine pt_setup
	
end subroutine pt_setup
