!> \file viewer.F90

!> \brief SPHRAY OpenGL viewer routines. 
!! 
!<

module viewermod
use particle_system_mod
use ray_mod
use raylist_mod
use opengl_gl
use opengl_glu
use opengl_glut
implicit none

private
save
public :: initViewer
public :: reset_to_init      ! initialize viewer
public :: NormalKey
public :: Idlefunc
public :: MouseClick         ! public to allow application
public :: MouseClickMotion
public :: NormalKeyUp        ! to 'take over'
public :: freeflightQ
public :: leftbuttonQ
public :: rightbuttonQ
public :: slewQ
public :: outputtext   
                                                     
integer, parameter :: MouseMovementScale=200
integer, parameter :: MouseMirror=1
real, parameter :: minSpotDistance=.001

real(kind=gldouble), parameter :: curveScale=0.02
real(kind=gldouble), parameter :: freeFlightSpeed=0.005
real(kind=gldouble), parameter :: M_PI=3.14159265358979323846_gldouble
real(kind=gldouble), dimension(3), parameter :: zero=(/0.,0.,0./)

integer(kind=glcint), parameter :: RESET = 10
integer(kind=glcint), parameter :: QUIT = 11

real, dimension(2) :: mousexyzero
integer :: npart, leftbutton=1,rightbutton=1,free_flight=0
logical :: slew=.TRUE.

real(kind=gldouble) :: scaleFactor=1.
real(kind=gldouble), dimension(3) :: cam,spot,upvector=(/0.,0.,1./),up,ri,dummy

real(kind=gldouble) :: viewangle=45.,ratio,angle
real(kind=gldouble):: nearplane=0.01,farplane=10000.
real(kind=gldouble) :: mx,my
real(kind=gldouble) :: fi=0.,theta=0.,r=1.
real(kind=gldouble),dimension(3) :: da=(/0.,0.03, 0.1/) 



real, public, parameter :: kpcunit=1.0

 contains

! state query functions
integer function freeflightQ()
  freeflightQ=free_flight
end function freeflightQ

integer function leftbuttonQ()
  leftbuttonQ=leftbutton
end function leftbuttonQ

integer function rightbuttonQ()
  rightbuttonQ=rightbutton
end function rightbuttonQ

logical function slewQ()
  slewQ=slew
end function slewQ




subroutine rotate(moving,fixed,r,fi,theta)
  real(kind=gldouble), dimension(3), intent(inout) :: moving
  real(kind=gldouble), dimension(3), intent(in) :: fixed
  real(kind=gldouble), intent(in) :: r, fi, theta
 
  moving(1)=fixed(1)+r*sin(fi)*cos(theta)
  moving(2)=fixed(2)+r*cos(fi)*cos(theta)
  moving(3)=fixed(3)+r*sin(theta)
end subroutine rotate


subroutine translate(one,two,trans)
  real(kind=gldouble), dimension(3), intent(inout) :: one,two
  real(kind=gldouble), dimension(3), intent(in) :: trans
  
  one(1)=one(1)+trans(1);two(1)=two(1)+trans(1)
  one(2)=one(2)+trans(2);two(2)=two(2)+trans(2)
  one(3)=one(3)+trans(3);two(3)=two(3)+trans(3)
end subroutine translate


function distance2(o, t)
  real(kind=gldouble), dimension(3), intent(in) :: o,t
  real(kind=gldouble) :: distance2
  distance2=(o(1)-t(1))**2+(o(2)-t(2))**2+(o(3)-t(3))**2
  return 
end function distance2


subroutine pointCamera()
  call glMatrixMode(GL_MODELVIEW)
  call glLoadIdentity()
  call gluLookAt( cam(1),      cam(2),      cam(3),     &
                  spot(1),     spot(2),     spot(3),    &
                  upvector(1), upvector(2), upvector(3) )
end subroutine pointCamera


subroutine ReshapeWindow(width, height) bind(c,name="ReshapeWindow")
  integer(kind=glint), value :: width, height
 
  if(height.LT.1) height=1_glint
  ratio=width/real(height)
  call glViewPort(0,0,width,height)
  call glMatrixMode(GL_PROJECTION)
  call glLoadIdentity()
  call gluPerspective(viewangle, ratio, nearplane, farplane)
  call pointCamera()
 
end subroutine ReshapeWindow


subroutine MouseClick(button, state, x, y) bind(c)
  integer(kind=glcint), value :: button,state,x,y 
   
  if(button.EQ.GLUT_LEFT_BUTTON) then
     if(state.EQ.GLUT_DOWN) leftbutton=0
     if(state.EQ.GLUT_UP) leftbutton=1
     mousexyzero(1)=x;mousexyzero(2)=y
     mx=0.
     my=0.
  endif
  if(button.EQ.GLUT_RIGHT_BUTTON) then
     if(state.EQ.GLUT_DOWN) rightbutton=0
     if(state.EQ.GLUT_UP) rightbutton=1
     mousexyzero(1)=x;mousexyzero(2)=y
  endif
  
end subroutine MouseClick


subroutine MouseClickMotion(x,y) bind(c)
  integer(kind=glcint), value :: x,y 
  real(kind=gldouble) :: rscl
  real(kind=gldouble),dimension(3) :: trans=(/0.,0.,0./)

  if((.NOT.(leftbutton.EQ.0)).AND.(.NOT.(rightbutton.EQ.0))) return

  if((.NOT.(leftbutton.EQ.0)).AND.(rightbutton.EQ.0)) then
     if(.NOT.slew) then
        r=r*2**((y-mousexyzero(2))/MouseMovementScale)
        if(r.LT.minSpotDistance) r=minSpotDistance
        if(r.gt.1000*scaleFactor) r=1000*scaleFactor
        call rotate(cam, spot,r,fi,theta)
     endif
     if(slew) then
        call rotate(trans,zero, &
             scaleFactor*((y-mousexyzero(2))/MouseMovementScale/2),fi,theta)
        call translate(cam,spot,trans)
     endif
  endif

  if((leftbutton.EQ.0).AND.(.NOT.(rightbutton.EQ.0))) then
     if(free_flight.EQ.1) then
        mx=mx+(x-mousexyzero(1))/MouseMovementScale
        my=my+(y-mousexyzero(2))/MouseMovementScale
        if(abs(mx).LT.0.02) mx=0.
        if(abs(my).LT.0.02) my=0.  
     endif
     if(free_flight.EQ.0) then
        fi=fi+(x-mousexyzero(1))/MouseMovementScale
        theta=theta+(y-mousexyzero(2))/MouseMovementScale
        if(theta.GT.M_PI/2.) theta=M_PI/2.-0.001
        if(theta.LT.-M_PI/2.) theta=-M_PI/2.+0.001
        if(slew) then
           call rotate(spot, cam,r,fi+M_PI,-theta)
        else
           call rotate(cam, spot,r,fi,theta)
        endif
     endif
  endif
  
  if((leftbutton.EQ.0).AND.(rightbutton.EQ.0)) then
     viewangle=viewangle*(1+(y-mousexyzero(2))/MouseMovementScale)
     if(viewangle.LT.1) viewangle=1
     if(viewangle.GT.120) viewangle=120
     call glMatrixMode(GL_PROJECTION)
     call glLoadIdentity()
     call gluPerspective(viewangle, ratio, nearplane, farplane)
  endif
  mousexyzero(1)=x;mousexyzero(2)=y;
  call pointCamera
  call glutPostRedisplay
  
end subroutine MouseClickMotion


subroutine NormalKey( key, x, y) bind(c)
  integer(kind=glcint), value :: key
  integer(kind=glcint), value :: x,y
  mx=0.;my=0.
  if(key.EQ.32) then 
     free_flight=1 
     slew=.TRUE.
  endif
  if(key.eq.113) stop
  if(key.eq.99) then
     print*,'----- view parameters -----'
     write(*,'("spot: ",3f8.2)') spot  
     write(*,'("cam: ",3f8.2)') cam
     write(*,'("dir: ",3f8.2)') spot-cam  
  endif
end subroutine NormalKey


subroutine NormalKeyUp( key, x, y) bind(c)
  integer(kind=glint), value :: key
  integer(kind=glcint), value :: x,y
  mx=0.;my=0.
  if(key.EQ.32) free_flight=0 
end subroutine NormalKeyUp


subroutine Idlefunc() bind(c)
  real(kind=gldouble),dimension(3) :: trans=(/0.,0.,0./)
end subroutine Idlefunc


recursive subroutine frame(dummy) bind(c)
  integer(kind=glint), value :: dummy
  real(kind=gldouble),dimension(3) :: trans=(/0.,0.,0./)

  if(free_flight.EQ.1) then
     call rotate(trans,zero,-freeFlightSpeed*scaleFactor,fi,theta)
     call translate(cam,spot,trans)
     if(leftbutton.EQ.0) then
        fi=fi+mx*curveScale
        theta=theta+my*curveScale
        if(theta.GT.M_PI/2.) theta=M_PI/2.-0.001
        if(theta.LT.-M_PI/2.) theta=-M_PI/2.+0.001
        call rotate(spot,cam,r,fi+M_PI,-theta)
     endif
     call pointCamera
     call glutPostRedisplay
  endif
  call glutTimerFunc(50, frame, 0)
end subroutine frame


subroutine MouseMoves( x, y)
  integer(glint) :: x,y
  
  if(free_flight.EQ.1) then
     mx=mx+(x-mousexyzero(1))/MouseMovementScale
     my=my+(y-mousexyzero(2))/MouseMovementScale
     if(abs(mx).LT.0.02) mx=0.
     if(abs(my).LT.0.02) my=0.
  endif
end subroutine MouseMoves


subroutine menu_handler(value) bind(c)
  integer(kind=glcint), value :: value
  select case(value)
  case(RESET)
     call reset_to_init
  case(QUIT)
     stop
  end select
  return
end subroutine menu_handler


subroutine reset_to_init(scale,centre)
  real,optional :: scale,centre(3)

  if(present(scale)) scaleFactor=scale
  spot(1)=0.;spot(2)=0. ;spot(3)=0.
  if(present(centre)) spot=centre

  slew=.FALSE.
  r=scaleFactor;fi=0.57;theta=0.3
  call rotate( cam,spot,r,fi,theta);
  call pointCamera
  call glutPostRedisplay 
end subroutine reset_to_init


function initViewer(scale,centre) result(menuid)
  real,optional :: scale,centre(3)
  integer(kind=glcint) :: menuid
 
  menuid = glutCreateMenu(menu_handler)
  call glutAddMenuEntry("reset to initial view",RESET)
  call glutAddMenuEntry("quit",QUIT)  
  call glutReshapeFunc( ReshapeWindow)
  call glutMousefunc( MouseClick)
  call glutMotionFunc( MouseClickMotion)
  call glutIgnoreKeyRepeat(1)
  call glutKeyboardFunc( NormalKey )
  call glutKeyboardUpFunc( NormalKeyUp)
  call glutIdleFunc( Idlefunc )
  call glutTimerFunc(50, frame, 0)
  call reset_to_init(scale,centre)
end function initViewer


subroutine outputtext(x,y,s)
  real(kind=gldouble) :: x,y
  character :: s*(*)
  character :: c
  integer :: i,lenc
  
  call glrasterpos2d(x,y)
  lenc = len(s)
  do i=1,lenc
     c = s(i:i)
     call glutbitmapcharacter(GLUT_BITMAP_HELVETICA_18,ichar(c))
  end do
end subroutine outputtext


end module viewermod








module snap_viewer
use particle_system_mod
use ray_mod
use raylist_mod
use opengl_gl
use opengl_glut
use opengl_glu

save
private
public :: initnbody
public :: display
public :: attachmenu
public :: createsubmenus
public :: IdleFunc2
public :: NormalKey2
public :: starscale
public :: starcentre
public :: starreset
public :: visible
public :: da
 
integer :: current_frame
integer, parameter :: NMAX=3000000
integer, parameter :: GAS=1,STAR=2,RAY=3,AXES=4,BOX=5
integer :: ngas,nstar,nray
real(kind=gldouble) :: boxsize,kpcboxsize,snaptime
real(kind=glfloat),dimension(3) :: da=(/0.,0., 1./) 
real(kind=glfloat),dimension(3) :: no=(/1.,0., 0./) 
!  real(kind=glfloat),dimension(3) :: da=(/.15,.075, .026/) 
!  real(kind=glfloat),dimension(3) :: no=(/1.,0., .0/) 
  
logical :: show_gas=.TRUE.,show_stars=.TRUE.,SHOW_ray=.FALSE.
logical :: show_axes=.FALSE.,show_box=.TRUE.,show_data=.TRUE.
logical :: show_brightness=.FALSE.
integer(glint) :: update_flag=0 ! 1+2+4 for gas+src+ray
integer :: show_temp=1
    
integer(glint) :: framerate=5000;
  
integer(glint), parameter:: FASTUPDATE=50,SLOWUPDATE=5000
  
real(kind=gldouble) :: top(3),bot(3)
  
character(len=200) :: string
    
contains

  
subroutine visible(state) bind(c)
  integer(glcint), value :: state
  if (state == GLUT_VISIBLE) then
     !    if(update_flag.EQ.0) call starreset(0_glint)
  else
     update_flag=0 
  endif
end subroutine visible
  
  
subroutine gaslist
use global_mod, only: psys,GV

    
real(kind=glfloat) :: xt,yt,zt,temp,dT,Tmin,Tmax,Tmid
real(kind=glfloat) :: red,green,blue,alpha_c
integer :: i,step
    
ngas=0
if (.NOT. allocated(psys%par)) return
ngas=size(psys%par)
if(.not.show_gas) return

call glNewList(1,GL_COMPILE)	
call glPointSize( 0.5_glfloat)
call glBegin(GL_POINTS)
call glColor4f(5._glfloat, .3_glfloat, .3_glfloat, 0.1_glfloat)


Tmax = 6.0
Tmin = 4.0
Tmid = 4.0
dT = Tmax-Tmin

do i=1,size(psys%par)
   
   xt=psys%par(i)%pos(1)
   yt=psys%par(i)%pos(2)
   zt=psys%par(i)%pos(3)
   
   ! temperature
   !-------------------
   if(show_temp.EQ.2) then

      temp = log10(psys%par(i)%T)
      if (temp <= Tmin) then
         temp = 0.0
      else if (temp >= Tmax) then
         temp = 1.0
      else
         temp = (temp - Tmin) / dT
      end if

      red=temp
      blue=(1.0_glfloat - temp)
      green=0.00
      alpha_c=0.1

      call glColor4f(red, green, blue, alpha_c)
   end if

   ! ionization
   !---------------
   if(show_temp.EQ.3) then
      temp = psys%par(i)%xHI 
      green = temp
      red = 1.0-temp
      blue=0.05
      if (temp < 0.5) then
         alpha_c=0.8
      else
         alpha_c=0.1
      endif
      call glColor4f(red, green, blue, alpha_c)
   endif


   ! only neutral
   !---------------
   if(show_temp.EQ.4) then
      temp = psys%par(i)%xHII 
      if (temp <= 0.5) then
         red = 0.5
         green = 0.2
         blue=0.5
         alpha_c=0.1
      else
         cycle
      endif
      call glColor4f(red, green, blue, alpha_c)
   endif


   ! only ionized
   !---------------
   if(show_temp.EQ.5) then
      temp = psys%par(i)%xHI 
      if (temp <= 0.5) then
         red = 5.0
         green = 0.3
         blue=0.3
         alpha_c=0.1
      else
         cycle
      endif
      call glColor4f(red, green, blue, alpha_c)
   endif


   call glVertex3f(xt,yt,zt)

end do

call glEnd()
call glEndList()

end subroutine gaslist
  

subroutine srclist
use global_mod, only: psys
real(kind=glfloat) :: xt,yt,zt,temp,red,green,blue,alpha_c
integer :: i,step
    
nstar=0
if(.NOT.allocated(psys%src)) return
nstar=size(psys%src)
if(.not.show_stars) return
    
call glNewList(2,GL_COMPILE)
call glPointSize( 4.5_glfloat)
call glBegin(GL_POINTS)
call glColor4f(.6_glfloat, 1._glfloat, .6_glfloat,1._glfloat)
nstar=size(psys%src)
do i=1,size(psys%src)
   xt=psys%src(i)%pos(1)
   yt=psys%src(i)%pos(2)
   zt=psys%src(i)%pos(3)
   if(show_brightness) then
      !	call glColor4f(red, green, blue,.6_glfloat)
      !	call glPointSize( (1._glfloat+3._glfloat*temp))
   end if
   call glVertex3f(xt,yt,zt)	
end do
call glEnd()
call glEndList()
    
end subroutine srclist

  
subroutine rlist
use global_mod, only: psys,globalraylist
integer :: i,pi
real(kind=glfloat) :: xt,yt,zt
    
if(.not.show_ray) return
if(.NOT.allocated(psys%par)) return
if(.NOT.allocated(globalraylist%intersection)) return
if(globalraylist%searchcell.NE.0) return

nray=globalraylist%nnb
    
call glNewList(3,GL_COMPILE)
call glPointSize( 3._glfloat)
call glBegin(GL_POINTS)
call glColor4f(1._glfloat, 1._glfloat, .5_glfloat,.8_glfloat)
do i=1,globalraylist%nnb
   pi=globalraylist%intersection(i)%pindx
   if(pi.LE.0.OR.pi.GT.ngas) exit
   xt=psys%par(pi)%pos(1)
   yt=psys%par(pi)%pos(2)
   zt=psys%par(pi)%pos(3)
   call glVertex3f(xt,yt,zt)
end do
call glEnd()
call glEndList()

end subroutine rlist

  
subroutine boxlist
    
  call glNewList(BOX,GL_COMPILE)
  call glBegin(GL_LINES)
  call glColor4f(.4_glfloat, .5_glfloat, .9_glfloat,1._glfloat)
    
  call glvertex3d(top(1),top(2),top(3))
  call glvertex3d(top(1),top(2),bot(3))
    
  call glvertex3d(top(1),top(2),top(3))
  call glvertex3d(top(1),bot(2),top(3))
    
  call glvertex3d(top(1),top(2),top(3))
  call glvertex3d(bot(1),top(2),top(3))
    
  call glvertex3d(bot(1),top(2),top(3))
  call glvertex3d(bot(1),top(2),bot(3))
    
  call glvertex3d(bot(1),top(2),top(3))
  call glvertex3d(bot(1),bot(2),top(3))
    
  call glvertex3d(top(1),bot(2),top(3))
  call glvertex3d(top(1),bot(2),bot(3))
    
  call glvertex3d(top(1),bot(2),top(3))
  call glvertex3d(bot(1),bot(2),top(3))
    
  call glvertex3d(top(1),top(2),bot(3))
  call glvertex3d(top(1),bot(2),bot(3))
  
  call glvertex3d(top(1),top(2),bot(3))
  call glvertex3d(bot(1),top(2),bot(3))
  
  call glvertex3d(bot(1),bot(2),bot(3))
  call glvertex3d(top(1),bot(2),bot(3))
  
  call glvertex3d(bot(1),bot(2),bot(3))
  call glvertex3d(bot(1),top(2),bot(3))
  
  call glvertex3d(bot(1),bot(2),bot(3))
  call glvertex3d(bot(1),bot(2),top(3))
  
  call glEnd()
  call glEndList()
  
end subroutine boxlist
  
  
function starscale() result(scale_factor)
  use global_mod, only: psys
  real :: maxcor,maxcor1,maxcor2,scale_factor
  !  real :: top(3), bot(3)
  ! schaal:
  top=psys%box%tops
  bot=psys%box%bots
  maxcor1=maxval(top)	
  maxcor2=maxval(bot)	
  if(maxcor1.le.0) maxcor1=1.
  if(maxcor2.le.0) maxcor2=1.
  maxcor=max(maxcor1,maxcor2)
  scale_factor=2*maxcor
  return
end function starscale
  

function starcentre() result(centre)
  use global_mod, only: psys
  real :: centre(3)
  top=psys%box%tops
  bot=psys%box%bots
  centre=top/2+bot/2
  return
end function starcentre


recursive subroutine starreset(dummy) bind(c)
  use global_mod
  integer(glint), value :: dummy
  call prerender()
  call glutPostRedisplay
  if(update_flag.NE.0) call glutTimerFunc( framerate, starreset, 1_glint)
end subroutine starreset
  
  
subroutine show_menu(optie) bind(c)
  integer(kind=glint), value :: optie
    
  select case (optie)
  case(1)
     show_stars=.NOT.show_stars
  case(2)
     show_gas=.NOT.show_gas
  case(3)
     show_ray=.NOT.show_ray
  case(4)
     !       call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,no)
  case(5) 
     !       call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,da)
  case(6)
     show_data=.NOT.show_data
  case(7)
     show_axes=.NOT.show_axes 
  case default
  end select
    
  call glutPostRedisplay
    
end subroutine show_menu

  
subroutine prop_menu(optie) bind(c)
  integer(glcint), value :: optie
  integer(glint) :: dummy,flag
    
  select case (optie)
  case(1)
     show_temp=1
  case(2)
     show_temp=2
  case(3)
     show_temp=3
  case(4)
     show_temp=4
  case(5)
     show_temp=5
  case(6)
     show_brightness=.NOT.show_brightness
  case default
     
  end select
    
  flag=update_flag
  update_flag=0
  call starreset(flag)
  update_flag=flag
  
end subroutine prop_menu

  
subroutine update_menu(optie) bind(c)
  use viewermod
  integer(glcint), value :: optie
  integer(glint) :: dummy,flag
  real :: rsize
  
  select case(optie)
     
  case(0)
     update_flag=0
  case(1)
     flag=update_flag
     update_flag=0
     if(flag.EQ.0) call starreset(0_glint)
     update_flag=flag
  case(2)
     framerate=SLOWUPDATE
     if(update_flag.EQ.0) then 
        update_flag=7  
        call starreset(0_glint)
     endif
  case(3)
     framerate=FASTUPDATE
     if(update_flag.EQ.0) then 
        update_flag=7  
        call starreset(0_glint)
     endif
  case(4)
  case(5)
  case(6)
     rsize=starscale()
     call reset_to_init(rsize,starcentre())
     da(3)=1./rsize**2*16.
     da(2)=da(3)*rsize/2.
     
     ! i added this in the hopes that it would update the new box for a new psys
     boxsize=starscale()
     kpcboxsize=boxsize*kpcunit
     da(3)=1./boxsize**2*16.
     da(2)=da(3)*boxsize/2.
     
     flag=update_flag
     update_flag=0
     call starreset(0_glint)
     update_flag=flag
  case default 
     
  end select
  
end subroutine update_menu

  
subroutine createsubmenus(i,j,k)
  integer(glcint) :: i,j,k
  
  i=glutCreateMenu(update_menu)
  call glutaddmenuentry('update scale', 6)
  call glutaddmenuentry('no updates', 0)
  call glutaddmenuentry('refresh view', 1)
  call glutaddmenuentry('slow update', 2)
  call glutaddmenuentry('fast update', 3)
  call glutaddmenuentry('update rays', 4)
  call glutaddmenuentry('update all', 5)
  
  j=glutCreateMenu(show_menu)
  call glutaddmenuentry('toggle sources', 1)
  call glutaddmenuentry('toggle gas', 2)
  call glutaddmenuentry('toggle rays', 3)
  call glutaddmenuentry('normal points', 4)
  call glutaddmenuentry('EXT points', 5)
  call glutaddmenuentry('toggle data', 6)
  call glutaddmenuentry('toggle axes', 7)
  
  k=glutCreateMenu(prop_menu)
  call glutaddmenuentry('plain', 1)
  call glutaddmenuentry('show temperature', 2)
  call glutaddmenuentry('show ionization', 3)
  call glutaddmenuentry('show only neutral', 4)
  call glutaddmenuentry('show only ionized', 5)
  call glutaddmenuentry('show luminosity', 6)
  
end subroutine createsubmenus


subroutine attachmenu(i,j,k)
  integer(glcint), intent(in) :: i,j,k
  call glutaddsubmenu('update options', i)
  call glutaddsubmenu('view options', j)
  call glutaddsubmenu('show options',k)
end subroutine attachmenu


subroutine prerender()
  integer(glint) :: flag
  call gaslist
  call srclist
  call rlist
  call boxlist
  call axeslist
end subroutine prerender


subroutine NormalKey2( key, x, y) bind(c)
  use viewermod
  integer(kind=glcint), value :: key
  integer(kind=glcint), value :: x,y
  if(key.EQ.110) then 
     ! 'N'
  endif
  if(key.EQ.112) then
     ! 'P'
  endif
  call Normalkey( key,x,y)
end subroutine NormalKey2

subroutine display bind(c)  
  call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
  call glDepthFunc(GL_ALWAYS)
  if( show_stars) call renderstars
  if( show_ray) call renderray
  if( show_gas) call rendergas
  
  
  call glDepthFunc(GL_LESS)
  call glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
  call glCallList(BOX)
  
  call glDepthFunc(GL_GEQUAL)
  call glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_ONE)
  call glCallList(BOX)
  
  call glDepthFunc(GL_LESS)
  call glBlendFunc(GL_ONE,GL_ZERO)
  call drawtext
  
  
  
  !  call glDepthFunc(GL_LESS)
  !  call glDepthMask(.TRUE.)
  !  call glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA)
  !  if( show_axes) call glCallList(AXES)
  !  call glCallList(BOX)
  !  if( show_axes) call glCallList(AXES)
  
  !  call glDepthMask(.FALSE.)
  !  call glDepthFunc(GL_ALWAYS)
  !  call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
  !  if( show_stars) call renderstars
  !  if( show_ray) call renderray
  !  if( show_gas) call rendergas
  
  !  call glDepthMask(.TRUE.)
  !  call glDepthFunc(GL_LESS)
  !  call glBlendFunc(GL_ONE,GL_ZERO)
  !  call drawtext
  call glutSwapBuffers();
end subroutine display

  
subroutine IdleFunc2  
end subroutine Idlefunc2
  

subroutine rendergas
  call glCallList(GAS)
end subroutine rendergas

  
subroutine renderstars
  call glCallList(STAR)
end subroutine renderstars

  
subroutine renderray
  call glCallList(RAY)
end subroutine renderray


subroutine axeslist
  
  call glNewList(AXES,GL_COMPILE)
  call glPushMatrix
  call glScaled(boxsize/3,boxsize/3,boxsize/3)
  call glColor4f(.2_glfloat, 1._glfloat, .2_glfloat,.8_glfloat)
  call glBegin(GL_LINES)
  call glVertex3f(-2.0_glfloat,0.0_glfloat,0.0_glfloat)
  call glVertex3f(2.0_glfloat,0.0_glfloat,0.0_glfloat)
  call glVertex3f(0.0_glfloat,-2.0_glfloat,0.0_glfloat)
  call glVertex3f(0.0_glfloat,2.0_glfloat,0.0_glfloat)
  call glVertex3f(0.0_glfloat,0.0_glfloat,-2.0_glfloat)
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
  
  call glPopMatrix
  
  call glEndList
  
end subroutine axeslist

  
subroutine itos(i,n,s)
  integer,intent(in) :: i,n
  character(len=n), intent(inout) :: s
  character(len=11) :: nstring
  data nstring/'0123456789X'/
  integer :: j,k,l
  
  if(i.LT.0.OR.i.GE.10**n) then
     do k=1,n
        s(k:k)=nstring(11:11)
     enddo
     return
  endif
  j=1
  do k=1,n
     l=1+mod(i,10*j)/j
     s(n-k+1:n-k+1)=nstring(l:l)
     j=j*10
  enddo
  
end subroutine itos


subroutine drawtext
  use global_mod, only: GV, PLAN
  use viewermod
  character*3 rs3
  character*8 rs
  character*10 rs10
  integer(kind=glint) :: viewport(4),rayn
  real(kind=gldouble) :: xp1,xp2,yp1,yp2
  
  if(.not.show_data) return
  call glGetIntegerv(GL_VIEWPORT,viewport)
  xp1=viewport(1)
  yp1=viewport(2)
  xp2=viewport(3)
  yp2=viewport(4) 
  call glmatrixmode(GL_PROJECTION)  
  call glPushMatrix
  call glLoadIdentity
  call glortho(xp1,xp2,yp1,yp2,-1._gldouble,1._gldouble) 
  call glmatrixmode(GL_MODELVIEW)
  call glPushMatrix
  call glLoadIdentity
  call glColor4f(0.0_glfloat,.5_glfloat,1.0_glfloat,0.5_glfloat)
  if(yp2-yp1.gt.120) then
     if(xp2-xp1.gt.320) then
        write(rs,'(f8.3)') GV%time_elapsed_myr
        call outputtext(xp1+3._gldouble,yp1+3._gldouble,'time elapsed (Myr):'//rs)
        
        write(rs,'(f8.2)') top(1)-bot(1)
        call outputtext(xp2-156._gldouble,yp1+3._gldouble,'box (Mpc):'//rs)
     endif
     
     if(xp2-xp1.gt.360) then
        !  write(rs,'(I10)') nbodies-nstar-nsph
        call itos(nray,8,rs)
        call outputtext(xp1+3._gldouble,yp2-20._gldouble,'in ray:')
        if(.not.show_ray) call glColor4f(0.4_glfloat,0.4_glfloat,0.4_glfloat,0.5_glfloat)
        call outputtext(xp1+70._gldouble,yp2-20._gldouble,rs)
        call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)
        
        call itos(ngas,8,rs)
        call outputtext(xp1+3._gldouble,yp2-40._gldouble,'gas:')
        if(.not.show_gas) call glColor4f(0.4_glfloat,0.4_glfloat,0.4_glfloat,0.5_glfloat)
        call outputtext(xp1+70._gldouble,yp2-40._gldouble,rs)
        call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)
        
        call itos(nstar,8,rs)
        call outputtext(xp1+3._gldouble,yp2-60._gldouble,'srcs:')
        if(.not.show_stars) call glColor4f(0.4_glfloat,.4_glfloat,.4_glfloat,0.5_glfloat)
        call outputtext(xp1+70._gldouble,yp2-60._gldouble,rs)
        
        rayn=GV%rayn
        call itos(rayn,10,rs10)
        !    write(rs10,'(i10)') GV%rayn
        call outputtext(xp2-204._gldouble,yp2-20._gldouble,rs10//' rays traced')
        !    call outputtext(xp1+108._gldouble,yp1-20._gldouble,rs10)
        write(rs3,'(i3)') int(100.*GV%rayn/real(PLAN%snap(GV%CurSnapNum)%SrcRays))
        call outputtext(xp2-54._gldouble,yp2-40._gldouble,rs3//'%')
        
        
     endif
  endif
  call glPopMatrix
  call glmatrixmode(GL_PROJECTION)
  call glPopMatrix
  call glmatrixmode(GL_MODELVIEW)
end subroutine drawtext


function initnbody()
  use global_mod, only: psys
  real :: initnbody
  real, parameter :: kpcunit=1.0
  
  ! while wait
  do while(.not.allocated(psys%par))
  enddo
  
  
  boxsize=starscale()
  kpcboxsize=boxsize*kpcunit
  initnbody=boxsize
  da(3)=1./initnbody**2*16.
  da(2)=da(3)*initnbody/2.
  
  call prerender()
  
end function initnbody

end module snap_viewer











subroutine pt_setup
end subroutine pt_setup



subroutine viewbodies
use ray_mod
use raylist_mod
  use viewermod
  use opengl_gl
  use opengl_glut
  use opengl_glu
  use snap_viewer
  implicit none
  
  real :: rsize
  integer(glcint) :: i,j,k,flag
  integer(glcint) :: win,menu
  character(len=200) :: string
  
  call glutInit()
  call glutInitDisplayMode(ior(GLUT_DEPTH,ior(GLUT_RGBA,GLUT_DOUBLE)))
  !  call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))
  
  call glutInitWindowPosition(100,100)
  call glutInitWindowSize(600,600)
  
  string="please, don't close me.."
  win=glutCreateWindow(string)
  
  
  call createsubmenus(i,j,k)
  
  rsize=initnbody()
  
  menu=initViewer(rsize,starcentre()) 
  ! NB
  call glutKeyboardFunc( NormalKey2 )
  ! NB (in initnbody?)
  call glutAttachMenu(GLUT_MIDDLE_BUTTON)
  call attachmenu(i,j,k)
  call glutDisplayFunc(display)
  
  call glEnable(GL_DEPTH_TEST)
  call glEnable(GL_BLEND)
  call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
  call glEnable(GL_POINT_SMOOTH)
  call glEnable(GL_LINE_SMOOTH)
  call glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
  call glHint(GL_POINT_SMOOTH_HINT, GL_NICEST)
  call glLineWidth(1.55_glfloat)
  !  call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,da)
  !  call glPointParameterfARB(GL_POINT_SIZE_MIN_ARB, 1._glfloat) 
  !  call glPointParameterfARB(GL_POINT_FADE_THRESHOLD_SIZE_ARB, 0._glfloat) 
  call glutVisibilityFunc(visible)
  
  call glutTimerFunc(100, starreset,0_glint)
  
  call glutMainLoop()
  
end subroutine viewbodies

