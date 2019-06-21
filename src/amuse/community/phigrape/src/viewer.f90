module viewermod

use opengl_gl
use opengl_glu
use opengl_glut

implicit none
 private
 save
 public :: initViewer,reset_to_init            ! initialize viewer
 public :: NormalKey,MouseClick, &             ! public to allow application
           MouseClickMotion,NormalKeyUp        ! to 'take over'
 public :: freeflightQ, leftbuttonQ, rightbuttonQ,slewQ,outputtext   
                                                     
integer, parameter :: MouseMovementScale=100
real(kind=gldouble), parameter :: curveScale=0.02
real(kind=gldouble), parameter :: freeFlightSpeed=0.0005
integer, parameter :: MouseMirror=1
real(kind=gldouble), parameter :: M_PI=3.14159265358979323846_gldouble
real(kind=gldouble),dimension(3), parameter :: zero=(/0.,0.,0./)
real, parameter :: minSpotDistance=.001
integer(kind=glcint), parameter :: RESET = 10, QUIT = 11
real(kind=gldouble) :: scaleFactor=1.
real(kind=gldouble), dimension(3) :: cam,spot,upvector=(/0.,0.,1./),up,ri,dummy
real, dimension(2) :: mousexyzero
integer :: npart, leftbutton=1,rightbutton=1,free_flight=0
real(kind=gldouble) :: viewangle=45.,ratio,angle
real(kind=gldouble):: nearplane=0.01,farplane=10000.
real(kind=gldouble) :: mx,my
real(kind=gldouble) :: fi=0.,theta=0.,r=1.
real(kind=gldouble),dimension(3) :: da=(/0.,0.03, 0.1/) 
logical :: slew=.TRUE.


 contains

! state query functions
integer function freeflightQ();freeflightQ=free_flight;end function
integer function leftbuttonQ();leftbuttonQ=leftbutton;end function leftbuttonQ
integer function rightbuttonQ();rightbuttonQ=rightbutton;end function rightbuttonQ
logical function slewQ();slewQ=slew;end function slewQ

subroutine rotate(moving,fixed,r,fi,theta)
  real(kind=gldouble), dimension(3), intent(inout) :: moving
  real(kind=gldouble), dimension(3), intent(in) :: fixed
  real(kind=gldouble),intent(in) :: r, fi, theta
 
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

function distance2( o, t)
  real(kind=gldouble), dimension(3), intent(in) :: o,t
  real(kind=gldouble) :: distance2
 distance2=(o(1)-t(1))**2+(o(2)-t(2))**2+(o(3)-t(3))**2
 return 

end function distance2

subroutine pointCamera()

 call glMatrixMode(GL_MODELVIEW)
 call glLoadIdentity()
 call gluLookAt( cam(1), cam(2), cam(3), &
                 spot(1), spot(2), spot(3), &
		 upvector(1), upvector(2), upvector(3))

end subroutine pointCamera


subroutine ReshapeWindow( width, height)
  integer(kind=glint), intent(inout) :: width,height
 
 if(height.LT.1) height=1
 ratio=width/real(height)
 call glViewPort(0,0,width,height)
 call glMatrixMode(GL_PROJECTION)
 call glLoadIdentity()
 call gluPerspective(viewangle, ratio, nearplane, farplane);
 call pointCamera()
 
end subroutine ReshapeWindow


subroutine MouseClick(button, state, x, y)
  integer(kind=glcint),intent(in out) :: button,state,x,y 
   
 if(button.EQ.GLUT_LEFT_BUTTON) then
  if(state.EQ.GLUT_DOWN) leftbutton=0
  if(state.EQ.GLUT_UP) leftbutton=1
 mousexyzero(1)=x;mousexyzero(2)=y
 mx=0.;my=0.
 endif
 if(button.EQ.GLUT_MIDDLE_BUTTON) then
  if(state.EQ.GLUT_DOWN) rightbutton=0
  if(state.EQ.GLUT_UP) rightbutton=1
 mousexyzero(1)=x;mousexyzero(2)=y
 endif
  
end subroutine MouseClick

subroutine MouseClickMotion(x,y)
  integer(kind=glcint),intent(in out) :: x,y 
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

subroutine NormalKey( key, x, y)
  integer(kind=glcint),intent(in out) :: key
  integer(kind=glcint),intent(in out) :: x,y
 mx=0.;my=0.
 if(key.EQ.32) then 
  free_flight=1 
  slew=.TRUE.
 endif
 if(key.eq.113) call glutLeaveMainLoop()
 if(key.eq.99) then
  print*,'----- view parameters -----'
  write(*,'("spot: ",3f8.2)') spot  
  write(*,'("cam: ",3f8.2)') cam
  write(*,'("dir: ",3f8.2)') spot-cam  
 endif
end subroutine NormalKey

subroutine NormalKeyUp( key, x, y)
  integer(kind=glint),intent(in out) :: key
  integer(kind=glcint),intent(in out) :: x,y
 mx=0.;my=0.
 if(key.EQ.32) free_flight=0 
end subroutine NormalKeyUp

recursive subroutine frame(dummy)
  real(kind=gldouble),dimension(3) :: trans=(/0.,0.,0./)
integer(kind=glint), intent(in out) :: dummy
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

subroutine menu_handler(value)
  integer(kind=glcint), intent(in out) :: value
 select case(value)
 case(RESET)
   call reset_to_init
 case(QUIT)
   call glutLeaveMainLoop()
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
    call glutbitmapcharacter(GLUT_BITMAP_HELVETICA_18, &
          ichar(c))
  end do
end subroutine outputtext


end module viewermod

module snap_viewer
 
  use opengl_gl
  use opengl_glut
  use opengl_glu
 save
 private
 public :: initnbody,display,attachmenu,createsubmenus, &
  NormalKey2,starscale,starcentre,starreset,visible
 public :: da
 
  integer, parameter :: NMAX=5000000
  integer, parameter :: GAS=1,STAR=2,HALO=3,AXES=4,BOX=5
  real(kind=gldouble) :: boxsize
  real(kind=glfloat),dimension(3) :: da=(/0.,0., 1./) 
  real(kind=glfloat),dimension(3) :: no=(/1.,0., 0./) 
!  real(kind=glfloat),dimension(3) :: da=(/.15,.075, .026/) 
!  real(kind=glfloat),dimension(3) :: no=(/1.,0., .0/) 
 
  logical :: show_gas=.FALSE.,show_stars=.TRUE.,show_halo=.FALSE.
  logical :: show_axes=.FALSE.,show_box=.TRUE.,show_data=.TRUE.
  integer(glint) :: update_flag=0 
    
  integer(glint) ::framerate=5000;
  integer(glint), parameter:: FASTUPDATE=50,SLOWUPDATE=5000
  
  real(kind=glfloat),save :: pointsize_gas=3.,pointsize_stars=1.5,pointsize_halo=1.
  
  real(kind=gldouble) :: top(3),bot(3)
    
 contains

subroutine visible(state)
integer(glcint),intent(inout) :: state
  if (state == GLUT_VISIBLE) then
!    if(update_flag.EQ.0) call starreset(0_glint)
   else
   update_flag=0 
  endif
end subroutine visible


subroutine starslist
#include "general.inc"
  real(kind=glfloat) :: xt,yt,zt,temp,red,green,blue,alpha_c
  integer :: i,step

  step=max(N/NMAX,1)
  if(.not.show_stars) return

  call glNewList(2,GL_COMPILE)
  call glPointSize( pointsize_stars)
  call glColor4f(.8_glfloat, .8_glfloat, .6_glfloat,.6_glfloat)
  call glBegin(GL_POINTS)
  do i=1,N,step

    xt=x(1,i)
    yt=x(2,i)
    zt=x(3,i)

    call glVertex3f(xt,yt,zt)	
  end do
  call glEnd()
  call glEndList()
end subroutine starslist

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

end subroutine


function starscale() result(scale_factor)
#include "general.inc"
	real(kind=gldouble) :: maxcor,maxcorg,maxcors,scale_factor
! schaal:
	maxcor=maxval(abs(x(1:3,1:N)))	
	if(maxcor.le.0) maxcor=1.
	scale_factor=2._gldouble*maxcor
	initialbox=2*maxcor
        top=initialbox/2
        bot=-initialbox/2
	return
end function starscale

function starcentre() result(centre)
  real :: centre(3)
        centre=top/2+bot/2
end function starcentre

recursive subroutine starreset(dummy)
 integer(glint) :: dummy
 call prerender()
 call glutPostRedisplay
 if(update_flag.NE.0) call glutTimerFunc( framerate, starreset, 1_glint)
end subroutine starreset

subroutine show_menu(optie)
  integer(kind=glint), intent(inout) :: optie
  
 select case (optie)
 case(1)
  show_stars=.NOT.show_stars
 case(2)
  show_gas=.NOT.show_gas
 case(3)
  show_halo=.NOT.show_halo
 case(4)
  call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,no)
 case(5) 
  call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,da)
 case(6)
  show_data=.NOT.show_data
 case(7)
  show_axes=.NOT.show_axes 
 case(8)
 if(pointsize_gas.LT.10) then
  pointsize_gas=pointsize_gas*1.5
  pointsize_stars=pointsize_stars*1.5
  pointsize_halo=pointsize_halo*1.5
 endif
 case(9)
 if(pointsize_gas.GT.1) then
  pointsize_gas=pointsize_gas/1.5
  pointsize_stars=pointsize_stars/1.5
  pointsize_halo=pointsize_halo/1.5
 endif
 case default
 end select

 call glutPostRedisplay
 
end subroutine show_menu

subroutine prop_menu(optie)
integer(glcint), intent(inout)::optie
integer(glint) :: dummy,flag

select case (optie)
 case(1)
 show_temp=0
 case(2)
 show_temp=1
 case(3)
 show_temp=2
 case(4)
 case default

end select

 flag=update_flag
 update_flag=0
 call starreset(flag)
 update_flag=flag

end subroutine prop_menu

subroutine update_menu(optie)
 use viewermod
integer(glcint),intent(inout) ::optie
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
!  boxsize=starscale()
!  da(3)=1./boxsize**2*16.
!  da(2)=da(3)*boxsize/2.

  flag=update_flag
  update_flag=0
  call starreset(0_glint)
  update_flag=flag
 case default 

end select

end subroutine update_menu

subroutine createsubmenus(i,j,k)
  integer(glcint),intent(inout) :: i,j,k
 
 i=glutCreateMenu(update_menu)
 call glutaddmenuentry('update scale', 6)
 call glutaddmenuentry('no updates', 0)
 call glutaddmenuentry('refresh view', 1)
 call glutaddmenuentry('slow update', 2)
 call glutaddmenuentry('fast update', 3)

 j=glutCreateMenu(show_menu)
 call glutaddmenuentry('toggle stars', 1)
 call glutaddmenuentry('toggle gas', 2)
 call glutaddmenuentry('toggle halo', 3)
 call glutaddmenuentry('normal points', 4)
 call glutaddmenuentry('EXT points', 5)
 call glutaddmenuentry('toggle data', 6)
 call glutaddmenuentry('toggle axes', 7)
 call glutaddmenuentry('increase pointsize', 8)
 call glutaddmenuentry('decrease pointsize', 9)

 k=glutCreateMenu(prop_menu)
 call glutaddmenuentry('plain', 1)

end subroutine createsubmenus

subroutine attachmenu(i,j,k)
 integer(glcint),intent(in) :: i,j,k
 call glutaddsubmenu('update options', i)
 call glutaddsubmenu('view options', j)
 call glutaddsubmenu('show options',k)

end subroutine attachmenu

subroutine prerender()
 integer(glint) :: flag
! call gaslist
 call starslist
! call rlist
 call boxlist
 call axeslist
end subroutine prerender

subroutine NormalKey2( key, x, y)
 use viewermod
  integer(kind=glcint),intent(in out) :: key
  integer(kind=glcint),intent(in out) :: x,y
 if(key.EQ.110) then 
! 'N'
 endif
  if(key.EQ.112) then
! 'P'
 endif
 call Normalkey( key,x,y)
end subroutine NormalKey2

subroutine display  
 call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
  call glDepthFunc(GL_ALWAYS)
  if( show_stars) call renderstars
  if( show_halo) call renderhalo
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
!  if( show_halo) call renderhalo
!  if( show_gas) call rendergas

!  call glDepthMask(.TRUE.)
!  call glDepthFunc(GL_LESS)
!  call glBlendFunc(GL_ONE,GL_ZERO)
!  call drawtext
 call glutSwapBuffers();
end subroutine display

subroutine rendergas
 call glCallList(GAS)
end subroutine rendergas

subroutine renderstars
 call glCallList(STAR)
end subroutine renderstars

subroutine renderhalo
  call glCallList(HALO)
end subroutine renderhalo

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

subroutine drawtext
 use viewermod
#include "general.inc"
  character*3 rs3
  character*8 rs
  character*10 rs10
  integer(kind=glint) :: viewport(4)
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
    write(rs,'(f8.2)') time_cur
    call outputtext(xp1+3._gldouble,yp1+3._gldouble,'time:'//rs)
    write(rs,'(f8.2)') top(1)-bot(1)
    call outputtext(xp2-120._gldouble,yp1+3._gldouble,'box:'//rs)
   endif

   if(xp2-xp1.gt.200) then
    write(rs,'(i8.8)') 0
    call outputtext(xp1+3._gldouble,yp2-20._gldouble,'halo:')
    if(.not.show_halo) call glColor4f(0.4_glfloat,0.4_glfloat,0.4_glfloat,0.5_glfloat)
    call outputtext(xp1+70._gldouble,yp2-20._gldouble,rs)
    call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)

    write(rs,'(i8.8)') 0
    call outputtext(xp1+3._gldouble,yp2-40._gldouble,'gas:')
    if(.not.show_gas) call glColor4f(0.4_glfloat,0.4_glfloat,0.4_glfloat,0.5_glfloat)
    call outputtext(xp1+70._gldouble,yp2-40._gldouble,rs)
    call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)

    write(rs,'(i8.8)') N
    call outputtext(xp1+3._gldouble,yp2-60._gldouble,'stars:')
    if(.not.show_stars) call glColor4f(0.4_glfloat,.4_glfloat,.4_glfloat,0.5_glfloat)
    call outputtext(xp1+70._gldouble,yp2-60._gldouble,rs)
  endif
  endif
  call glPopMatrix
  call glmatrixmode(GL_PROJECTION)
  call glPopMatrix
  call glmatrixmode(GL_MODELVIEW)
end subroutine drawtext

function initnbody()
#include "general.inc"
  real :: initnbody
  
 boxsize=starscale()
 initnbody=boxsize
 da(3)=1./initnbody**2*16.
 da(2)=da(3)*initnbody/2.

 call prerender()
 
end function initnbody

end module snap_viewer


subroutine pt_setup
	
end subroutine pt_setup

subroutine viewbodies
 use viewermod
 use snap_viewer
 use opengl_gl
 use opengl_glut
 use opengl_glu
 implicit none
  real :: rsize
  integer(glcint) :: i,j,k,flag
  integer(glcint) :: win,menu
 
 if(glutGet(GLUT_ELAPSED_TIME).NE.0) then
  print*,'already viewer present'
!  return
 endif
  
 call glutInit
 call glutInitWindowPosition(100,100)
 call glutInitWindowSize(700,400)
 call glutInitDisplayMode(ior(GLUT_DEPTH,ior(GLUT_RGBA,GLUT_DOUBLE)))
 
 win=glutCreateWindow("3D Viewer")
 call createsubmenus(i,j,k)
 
 rsize=initnbody()
 
 menu=initViewer(rsize,starcentre()) 
! NB
 call glutKeyboardFunc( NormalKey2 )
! NB (in initnbody?)
 call glutAttachMenu(GLUT_RIGHT_BUTTON)
 call attachmenu(i,j,k)
 call glutDisplayFunc(display)

 call glClearColor(0.02_glclampf, 0.02_glclampf, 0.02_glclampf, 1.0_glclampf)

 call glEnable(GL_DEPTH_TEST)
 call glEnable(GL_BLEND)
 call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
 call glEnable(GL_POINT_SMOOTH)
 call glEnable(GL_LINE_SMOOTH)
 call glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
 call glHint(GL_POINT_SMOOTH_HINT, GL_NICEST)
 call glLineWidth(1.55_glfloat)
 call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,da)
 call glPointParameterfARB(GL_POINT_SIZE_MIN_ARB, 1._glfloat) 
 call glPointParameterfARB(GL_POINT_FADE_THRESHOLD_SIZE_ARB, 0._glfloat) 


 call glutTimerFunc(100, starreset,0_glint)
 call glutVisibilityFunc(visible)

 call glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS)
 call glutMainLoop()

end subroutine viewbodies

