module viewer

use opengl_gl
use opengl_glu
use opengl_glut

implicit none
 private
 save
 public :: initViewer                               ! initialize viewer
 public :: NormalKey, Idlefunc,MouseClick, &        ! public to allow application
           MouseClickMotion,NormalKeyUp             !  to  'take over'
 public :: freeflightQ, leftbuttonQ, rightbuttonQ,slewQ,outputtext   
                                                     
integer, parameter :: MouseMovementScale=200
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
 if(button.EQ.GLUT_RIGHT_BUTTON) then
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

subroutine Idlefunc()
  real(kind=gldouble),dimension(3) :: trans=(/0.,0.,0./)

end subroutine Idlefunc

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
 call glutTimerFunc(20, frame, 0)
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

subroutine reset_to_init
 slew=.FALSE.
 r=scaleFactor;fi=0.;theta=0.
 spot(1)=0.;spot(2)=0. ;spot(3)=0.
 call rotate( cam,spot,r,fi,theta);
 call pointCamera
 call glutPostRedisplay 
end subroutine reset_to_init

function initViewer(starScale) result(menuid)
  real,intent(in) :: starScale
  integer(kind=glcint) :: menuid
 
 
 if(starScale.GT.0) scaleFactor=starScale

 
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
 call glutTimerFunc(40, frame, 0)

 call reset_to_init
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


end module viewer

module snap_viewer
 
  use opengl_gl
  use opengl_glut
  use opengl_glu
 save
 private
 public :: initnbody,display,attachmenu,createsubmenus,IdleFunc2,NormalKey2,da
 
  integer :: current_frame
  integer, parameter :: NMAX=3000000
  integer, parameter :: GAS=1,STAR=2,HALO=3,AXES=4,BOX=5
  integer :: ngasp,nstarp,nhalop,ntot,ngas,nstars  
  real(kind=gldouble) :: boxsize,kpcboxsize,snaptime
  real(kind=glfloat),dimension(3) :: da=(/0.,0., 1./) 
  real(kind=glfloat),dimension(3) :: no=(/1.,0., 0./) 
!  real(kind=glfloat),dimension(3) :: da=(/.15,.075, .026/) 
!  real(kind=glfloat),dimension(3) :: no=(/1.,0., .0/) 
 
  logical :: show_gas=.TRUE.,show_stars=.TRUE.,SHOW_halo=.FALSE.
  logical :: show_axes=.TRUE.,show_box=.TRUE.,show_data=.TRUE.
  integer :: show_temp=0,show_age=.FALSE.

  integer :: ok=0
  integer :: step
  real :: tfac
    
 contains

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
 case default
 end select

 call glutPostRedisplay
 
end subroutine show_menu

subroutine prop_menu(optie)
  integer(kind=glcint), intent(inout)::optie

 select case (optie)
 case(1)
  show_temp=0
  show_age=.FALSE.
 case(2)
  if(show_temp.eq.1) then 
   show_temp=0
  else
   show_temp=1
  endif
 case(3)
  show_age=.NOT.show_age
 case(4)
  if(show_temp.eq.2) then
   show_temp=0
  else 
   show_temp=2
   endif
 case(5)
  if(show_temp.eq.3) then 
   show_temp=0
  else
   show_temp=3
  endif
 case(6)
  if(show_temp.eq.4) then 
   show_temp=0
  else
   show_temp=4
  endif

 case default
end select
 call prerender
 call glutPostRedisplay
end subroutine prop_menu

subroutine createsubmenus(i,j,k)
  integer(glcint),intent(inout) :: i,j,k
 
 i=glutCreateMenu(file_menu)
 call glutaddmenuentry('initial ', -1000)
 call glutaddmenuentry('rewind 5 ', -5)
 call glutaddmenuentry('rewind', -1)
 call glutaddmenuentry('forward', 1)
 call glutaddmenuentry('forward 5', 5)
 call glutaddmenuentry('forward 10', 10)
 
 j=glutCreateMenu(show_menu)
 call glutaddmenuentry('toggle stars', 1)
 call glutaddmenuentry('toggle gas', 2)
 call glutaddmenuentry('toggle halo', 3)
 call glutaddmenuentry('normal points', 4)
 call glutaddmenuentry('EXT points', 5)
 call glutaddmenuentry('toggle data', 6)
 call glutaddmenuentry('toggle axes', 7)
 

 k=glutCreateMenu(prop_menu)
 call glutaddmenuentry('plain',1)
 call glutaddmenuentry('toggle temperature', 2)
 call glutaddmenuentry('toggle fuvheat', 4)
 call glutaddmenuentry('toggle luminosity', 3)
 call glutaddmenuentry('toggle unstable', 5)
 call glutaddmenuentry('toggle ionization', 6)

end subroutine createsubmenus

subroutine attachmenu(i,j,k)
 integer(glcint),intent(in) :: i,j,k
 call glutaddsubmenu('load file', i)
 call glutaddsubmenu('view options', j)
 call glutaddsubmenu('show properties',k)

end subroutine attachmenu

function starscale() result(box)

	include 'globals.h'
	real(kind=gldouble) :: maxcor,maxcorg,maxcors,box
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
	box=2*maxcor
	return
end function starscale


function initnbody()
  use ElementsMod
  use starsMod
  include 'globals.h'
  real :: initnbody
  call initmem(nbodsmax,nsphmax,ncells)
 
 current_frame=0
 
 CALL set_parameters(0)
 CALL readinput
 CALL initpars
  call InitElements(metallicity)
  call InitStars(datadir,zQ())
 CALL setbox('sph ')
 
 boxsize=starscale()
 kpcboxsize=boxsize*unitl_in_kpc
 initnbody=1.5*boxsize

 call selectparticles
 call prerender
 
end function initnbody

subroutine readinput
 include 'globals.h'
 integer,parameter :: ndigits=6
 character(len=ndigits) :: nstring
 character(len=24) :: filenaam

 if(current_frame.lt.0) current_frame=0
 call itos(current_frame,ndigits,nstring)
 filenaam=TRIM(outputfile)//'.'//nstring
 call readbods(filenaam)

 ntot=nbodies
 ngas=nsph
 nstars=nstar
 snaptime=tnow*timescale/year/1.e6

end subroutine readinput

subroutine selectparticles
  include 'globals.h'
  integer :: p,i
  real :: temp
 step=nbodies/NMAX+1
 tfac=meanmwt*mhboltz*(gamma-1.)
 ok=1
end subroutine

subroutine copyfromnbody
 
 call readinput
 call selectparticles
 call prerender
 call glutPostRedisplay

end subroutine copyfromnbody

subroutine prerender

 if(ok.NE.1) return
 call gaslist
 call starlist
 call darklist
 call drawaxes
end subroutine prerender

subroutine file_menu(i)
 integer(kind=glcint),intent(inout) :: i

 current_frame=current_frame+i
 call copyfromnbody
  
end subroutine file_menu

subroutine NormalKey2( key, x, y)
 use viewer
  integer(kind=glcint),intent(in out) :: key
  integer(kind=glcint),intent(in out) :: x,y
 if(key.EQ.110) then 
  current_frame=current_frame+1
  call copyfromnbody
 endif
  if(key.EQ.112) then
  current_frame=current_frame-1
  call copyfromnbody
 endif
 call Normalkey( key,x,y)
end subroutine NormalKey2

subroutine display
 call glClearColor(0.02_glfloat,0.02_glfloat,0.02_glfloat, 0.0_glfloat)
 call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
  call glDepthFunc(GL_ALWAYS)
  call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
  if( show_stars) call renderstars
  if( show_halo) call renderhalo
  if( show_gas) call rendergas
  call glDepthFunc(GL_LESS)
  if( show_axes) call glCallList(AXES)
  call glCallList(BOX)
  call glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA)
  call drawtext
 call glutSwapBuffers();
end subroutine display

subroutine IdleFunc2

end subroutine Idlefunc2


subroutine gaslist
 include 'globals.h'
 integer i
 real(kind=glfloat) :: xt,yt,zt,red,green,blue,alpha_c,temp

 call glNewList(GAS,GL_COMPILE)
 call glPointSize( 3._glfloat)
 call glBegin(GL_POINTS)
 call glColor4f(1._glfloat, .8_glfloat, .8_glfloat,.6_glfloat)

 do i=1,nsph,step
        xt=pos(i,1)
	yt=pos(i,2)
	zt=pos(i,3)
	
	if(show_temp.EQ.1) then
	
	  temp=(log10(temperat(i))-1.)/8.
	  red=-0.125+2*temp
	  blue=MAX(1-8/3.*temp,-3./8.+temp)+0.2
	  green=-1.+3*temp+0.2
          alpha_c=.6
!	  alpha_c=.2+.5*temp
	  call glColor4f(red, green, blue,alpha_c)
	
	end if
	
	if(show_temp.eq.2) then
	
	  temp=log(fuvheat(i))-log(7.)
	  red=.5+temp/4
	  blue=Max(0.6,1.+temp/4)
	  green=.5+temp/4
	  alpha_c=.8
	
	  call glColor4f(red, green, blue,alpha_c)
	
	end if

	if(show_temp.eq.3) then
          red=0.4;green=0.4;blue=0.4;alpha_c=0.6
	 if(rho(i)*pi/6.*((gamma*gamma1*temperat(i)/tfac)*pi/rho(i))**1.5.LE.masscrit)then
	   red=0.;green=1.;blue=0.
	 endif 
	 call glColor4f(red, green, blue,alpha_c)
	endif

	if(show_temp.eq.4) then
	
	  temp=elecfrac(i)
	  if(temp.GT.1) temp=1.
	  temp=3+log10(temp)
	  if(temp.LE.0) temp=0.
	  temp=temp/3.
	  red=temp
	  green=1.-temp
	  blue=0.
	  alpha_c=.6
	  
	  call glColor4f(red, green, blue,alpha_c)
	
	end if


        call glVertex3f(xt,yt,zt)



 enddo

 call glEnd()
 call glEndList()


end subroutine gaslist

subroutine starlist
 use starsMod
 include 'globals.h'
  integer i,nbands
  real(kind=glfloat) :: xt,yt,zt,red,green,blue,alpha_c,temp
  real :: age,band(1:20)
 nbands=nbandsQ()

 call glNewList(STAR,GL_COMPILE)
 call glPointSize( 1.5_glfloat)
 if(show_age) call glPointSize( 1.5_glfloat)
 call glColor4f(1._glfloat, 1._glfloat, .8_glfloat,.75_glfloat)
 call glBegin(GL_POINTS)

 do i=nbodies-nstar+1,nbodies,step
 
  xt=pos(i,1)
  yt=pos(i,2)
  zt=pos(i,3)
  if(show_age) then
   age=(tnow-tform(i))*timescale/year
   band=mbands(age)/mbands(0.)
   red=0.9+alog10(band(5))/3
   green=0.9+alog10(band(2))/3
   blue=0.9+alog10(band(1))/3
   alpha_c=0.5+green/2
   call glColor4f(red, green, blue,alpha_c)
 !  call glPointSize( (1._glfloat+3._glfloat*temp))
  end if
		
  call glVertex3f(xt,yt,zt)
 
 enddo
 call glEnd()
 call glEndList()


end subroutine starlist

subroutine darklist
 include 'globals.h'
  integer i
  real(kind=glfloat) :: xt,yt,zt,red,green,blue,alpha_c,temp

 call glNewList(HALO,GL_COMPILE)
 call glPointSize( 2._glfloat)

 call glColor4f(.75_glfloat, .75_glfloat, .75_glfloat,.7_glfloat)
 call glBegin(GL_POINTS)

 do i=nsph+1,nbodies-nstar
  xt=pos(i,1)
  yt=pos(i,2)
  zt=pos(i,3)
  call glVertex3f(xt,yt,zt)
 enddo
 call glEnd()
 call glEndList()
end subroutine darklist


subroutine rendergas
 call glCallList(GAS)
end subroutine rendergas

subroutine renderstars
 call glCallList(STAR)
end subroutine renderstars

subroutine renderhalo
  call glCallList(HALO)
end subroutine renderhalo

subroutine drawaxes

 call glNewList(BOX,GL_COMPILE)
! Draw axes so we know the orientation
 call glPushMatrix
 call glScaled(boxsize/3,boxsize/3,boxsize/3)
 call glColor4f(.3_glfloat,.3_glfloat,1._glfloat,1._glfloat);
 call glutWireCube(3._gldouble);
 call glPopMatrix
 call glEndList

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

end subroutine drawaxes


subroutine drawtext
 use viewer
  character*8 rs
  integer(kind=glint) :: viewport(4)
  real(kind=gldouble) :: xp1,xp2,yp1,yp2
  
  if(.not.show_data) return
  call glGetIntegerv(GL_VIEWPORT,viewport) 
  call glmatrixmode(GL_PROJECTION)
  call glPushMatrix
  call glLoadIdentity
  xp1=viewport(1)
  yp1=viewport(2)
  xp2=viewport(3)
  yp2=viewport(4)
  call glortho(xp1,xp2,yp1,yp2,-1._gldouble,1._gldouble)
  call glmatrixmode(GL_MODELVIEW)
  call glPushMatrix
  call glLoadIdentity
  call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)
  if(yp2-yp1.gt.120) then
   if(xp2-xp1.gt.320) then
    write(rs,'(f8.2)') snaptime
    call outputtext(xp1+3._gldouble,yp1+3._gldouble,'time (Myr):'//rs)
    write(rs,'(f8.2)') kpcboxsize
    call outputtext(xp2-156._gldouble,yp1+3._gldouble,'box (kpc):'//rs)
   endif

   if(xp2-xp1.gt.200) then
    call itos(ntot-ngas-nstars,8,rs)
    call outputtext(xp1+3._gldouble,yp2-20._gldouble,'halo:')
    if(.not.show_halo) call glColor4f(0.4_glfloat,0.4_glfloat,0.4_glfloat,0.5_glfloat)
    call outputtext(xp1+70._gldouble,yp2-20._gldouble,rs)
    call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)

    call itos(ngas,8,rs)
    call outputtext(xp1+3._gldouble,yp2-40._gldouble,'gas:')
    if(.not.show_gas) call glColor4f(0.4_glfloat,0.4_glfloat,0.4_glfloat,0.5_glfloat)
    call outputtext(xp1+70._gldouble,yp2-40._gldouble,rs)
    call glColor4f(0._glfloat,.6_glfloat,1.0_glfloat,0.8_glfloat)

    call itos(nstars,8,rs)
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


end module snap_viewer


program test
 use viewer
 use opengl_gl
 use opengl_glut
 use opengl_glu
 use snap_viewer
 implicit none
  real :: rsize
  integer(glcint) :: i,j,k
  integer(glcint) :: win,menu
 
 call glutInit
 call glutInitWindowPosition(100,100)
 call glutInitWindowSize(700,400)
 call glutInitDisplayMode(ior(GLUT_DEPTH,ior(GLUT_RGBA,GLUT_DOUBLE)))
 
 win=glutCreateWindow("SNAP viewer")
 call createsubmenus(i,j,k)
 
 rsize=initnbody()
 da(3)=1./rsize**2*16.
 da(2)=da(3)*rsize/2.

 menu=initViewer(rsize) 
! NB
 call glutKeyboardFunc( NormalKey2 )
! NB (in initnbody?)
 call glutAttachMenu(GLUT_MIDDLE_BUTTON)
 call attachmenu(i,j,k)
 call glutDisplayFunc( display)

 call glEnable(GL_DEPTH_TEST)
 call glEnable(GL_BLEND)
 call glBlendFunc(GL_SRC_ALPHA,GL_ONE)
 call glEnable(GL_POINT_SMOOTH)
 call glEnable(GL_LINE_SMOOTH)
 call glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
 call glHint(GL_POINT_SMOOTH_HINT, GL_NICEST)
 call glLineWidth(1.25_glfloat)
 call glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB,da)
 call glPointParameterfARB(GL_POINT_SIZE_MIN_ARB, 2._glfloat) 
 call glPointParameterfARB(GL_POINT_FADE_THRESHOLD_SIZE_ARB, 0._glfloat) 

 call glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS)
 call glutMainLoop()

end program test

