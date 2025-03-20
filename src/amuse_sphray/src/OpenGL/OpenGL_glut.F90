MODULE OpenGL_glut
USE OpenGL_kinds
IMPLICIT NONE
PRIVATE

!  Copyright (c) Mark J. Kilgard, 1994, 1995, 1996, 1998. 

!     provided without guarantee or warrantee expressed or  implied. This

!  Cut to provide just basic OpenGL interface by A. Donev 


!  Define GLUTAPIENTRY and GLUTCALLBACK to nothing if we aren't on Win32. 

!   GLUT API revision history:
!  
!   GLUT_API_VERSION is updated to reflect incompatible GLUT
!   API changes (interface changes, semantic changes, deletions,
!   or additions).
!  
!   GLUT_API_VERSION=1  First public release of GLUT.  11/29/94
!  
!   GLUT_API_VERSION=2  Added support for OpenGL/GLX multisampling,
!   extension.  Supports new input devices like tablet, dial and button
!   box, and Spaceball.  Easy to query OpenGL extensions.
!  
!   GLUT_API_VERSION=3  glutMenuStatus added.
!  
!   GLUT_API_VERSION=4  glutInitDisplayString, glutWarpPointer,
!   glutBitmapLength, glutStrokeLength, glutWindowStatusFunc, dynamic
!   video resize subAPI, glutPostWindowRedisplay, glutKeyboardUpFunc,
!   glutSpecialUpFunc, glutIgnoreKeyRepeat, glutSetKeyRepeat,
!   glutJoystickFunc, glutForceJoystickFunc (NOT FINALIZED!).
!  
!   GLUT_API_VERSION=5  glutGetProcAddress (added by BrianP)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_API_VERSION = 5

!   GLUT implementation revision history:
!  
!   GLUT_XLIB_IMPLEMENTATION is updated to reflect both GLUT
!   API revisions and implementation revisions (ie, bug fixes).
!  
!   GLUT_XLIB_IMPLEMENTATION=1  mjk's first public release of
!   GLUT Xlib-based implementation.  11/29/94
!  
!   GLUT_XLIB_IMPLEMENTATION=2  mjk's second public release of
!   GLUT Xlib-based implementation providing GLUT version 2
!   interfaces.
!  
!   GLUT_XLIB_IMPLEMENTATION=3  mjk's GLUT 2.2 images. 4/17/95
!  
!   GLUT_XLIB_IMPLEMENTATION=4  mjk's GLUT 2.3 images. 6/?/95
!  
!   GLUT_XLIB_IMPLEMENTATION=5  mjk's GLUT 3.0 images. 10/?/95
!  
!   GLUT_XLIB_IMPLEMENTATION=7  mjk's GLUT 3.1+ with glutWarpPoitner.  7/24/96
!  
!   GLUT_XLIB_IMPLEMENTATION=8  mjk's GLUT 3.1+ with glutWarpPoitner
!   and video resize.  1/3/97
!  
!   GLUT_XLIB_IMPLEMENTATION=9 mjk's GLUT 3.4 release with early GLUT 4 routines.
!  
!   GLUT_XLIB_IMPLEMENTATION=11 Mesa 2.5's GLUT 3.6 release.
!  
!   GLUT_XLIB_IMPLEMENTATION=12 mjk's GLUT 3.6 release with early GLUT 4 routines + signal handling.
!  
!   GLUT_XLIB_IMPLEMENTATION=13 mjk's GLUT 3.7 beta with GameGLUT support.
!  
!   GLUT_XLIB_IMPLEMENTATION=14 mjk's GLUT 3.7 beta with f90gl friend interface.
!  
!   GLUT_XLIB_IMPLEMENTATION=15 mjk's GLUT 3.7 beta sync'ed with Mesa <GL/glut.h>
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_XLIB_IMPLEMENTATION = 15

!  Display mode bit masks. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RGB = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RGBA = GLUT_RGB
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INDEX = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SINGLE = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DOUBLE = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACCUM = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ALPHA = 8
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEPTH = 16
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_STENCIL = 32
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MULTISAMPLE = 128
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_STEREO = 256
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LUMINANCE = 512

!  Mouse buttons. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LEFT_BUTTON = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MIDDLE_BUTTON = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RIGHT_BUTTON = 2

!  Mouse button  state. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DOWN = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_UP = 1

!  function keys 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F1 = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F2 = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F3 = 3
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F4 = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F5 = 5
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F6 = 6
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F7 = 7
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F8 = 8
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F9 = 9
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F10 = 10
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F11 = 11
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F12 = 12
!  directional keys 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_LEFT = 100
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_UP = 101
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_RIGHT = 102
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_DOWN = 103
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_PAGE_UP = 104
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_PAGE_DOWN = 105
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_HOME = 106
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_END = 107
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_INSERT = 108

!  Entry/exit  state. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LEFT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ENTERED = 1

!  Menu usage  state. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_NOT_IN_USE = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_IN_USE = 1

!  Visibility  state. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NOT_VISIBLE = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VISIBLE = 1

!  Window status  state. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HIDDEN = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULLY_RETAINED = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_PARTIALLY_RETAINED = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULLY_COVERED = 3

!  Color index component selection values. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RED = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GREEN = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_BLUE = 2

!  Layers for use. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NORMAL = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY = 1

!  Stroke font opaque addresses (use constants instead in source code). 

!  Stroke font constants (use these in GLUT program). 
! ???, , PUBLIC :: GLUT_STROKE_ROMAN=(&glutStrokeRoman)
! ???, , PUBLIC :: GLUT_STROKE_MONO_ROMAN=(&glutStrokeMonoRoman)

!  Bitmap font opaque addresses (use constants instead in source code). 

!  Bitmap font constants (use these in GLUT program). 
! ???, , PUBLIC :: GLUT_BITMAP_9_BY_15=(&glutBitmap9By15)
! ???, , PUBLIC :: GLUT_BITMAP_8_BY_13=(&glutBitmap8By13)
! ???, , PUBLIC :: GLUT_BITMAP_TIMES_ROMAN_10=(&glutBitmapTimesRoman10)
! ???, , PUBLIC :: GLUT_BITMAP_TIMES_ROMAN_24=(&glutBitmapTimesRoman24)
! ???, , PUBLIC :: GLUT_BITMAP_HELVETICA_10=(&glutBitmapHelvetica10)
! ???, , PUBLIC :: GLUT_BITMAP_HELVETICA_12=(&glutBitmapHelvetica12)
! ???, , PUBLIC :: GLUT_BITMAP_HELVETICA_18=(&glutBitmapHelvetica18)

!  glutGet parameters. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_X = 100
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_Y = 101
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_WIDTH = 102
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_HEIGHT = 103
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BUFFER_SIZE = 104
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_STENCIL_SIZE = 105
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_DEPTH_SIZE = 106
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_RED_SIZE = 107
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_GREEN_SIZE = 108
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BLUE_SIZE = 109
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ALPHA_SIZE = 110
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_RED_SIZE = 111
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_GREEN_SIZE = 112
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_BLUE_SIZE = 113
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_ALPHA_SIZE = 114
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_DOUBLEBUFFER = 115
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_RGBA = 116
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_PARENT = 117
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_NUM_CHILDREN = 118
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_COLORMAP_SIZE = 119
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_NUM_SAMPLES = 120
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_STEREO = 121
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_CURSOR = 122
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_WIDTH = 200
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_HEIGHT = 201
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_WIDTH_MM = 202
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_HEIGHT_MM = 203
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_NUM_ITEMS = 300
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DISPLAY_MODE_POSSIBLE = 400
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_X = 500
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_Y = 501
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_WIDTH = 502
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_HEIGHT = 503
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_DISPLAY_MODE = 504
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ELAPSED_TIME = 700
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_FORMAT_ID = 123

!  glutDeviceGet parameters. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_KEYBOARD = 600
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_MOUSE = 601
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_SPACEBALL = 602
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_DIAL_AND_BUTTON_BOX = 603
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_TABLET = 604
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_MOUSE_BUTTONS = 605
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_SPACEBALL_BUTTONS = 606
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_BUTTON_BOX_BUTTONS = 607
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_DIALS = 608
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_TABLET_BUTTONS = 609
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEVICE_IGNORE_KEY_REPEAT = 610
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEVICE_KEY_REPEAT = 611
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_JOYSTICK = 612
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OWNS_JOYSTICK = 613
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTONS = 614
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_AXES = 615
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_POLL_RATE = 616

!  glutLayerGet parameters. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY_POSSIBLE = 800
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LAYER_IN_USE = 801
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_OVERLAY = 802
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_TRANSPARENT_INDEX = 803
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NORMAL_DAMAGED = 804
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY_DAMAGED = 805

!  glutVideoResizeGet parameters. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_POSSIBLE = 900
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_IN_USE = 901
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_X_DELTA = 902
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_Y_DELTA = 903
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_WIDTH_DELTA = 904
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_HEIGHT_DELTA = 905
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_X = 906
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_Y = 907
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_WIDTH = 908
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_HEIGHT = 909

!  glutGetModifiers return mask. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_SHIFT = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_CTRL = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_ALT = 4

!  glutSetCursor parameters. 
!  Basic arrows. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_RIGHT_ARROW = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_ARROW = 1
!  Symbolic cursor shapes. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_INFO = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_DESTROY = 3
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_HELP = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_CYCLE = 5
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_SPRAY = 6
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_WAIT = 7
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TEXT = 8
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_CROSSHAIR = 9
!  Directional cursors. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_UP_DOWN = 10
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_RIGHT = 11
!  Sizing cursors. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_SIDE = 12
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_SIDE = 13
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_SIDE = 14
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_RIGHT_SIDE = 15
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_LEFT_CORNER = 16
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_RIGHT_CORNER = 17
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_RIGHT_CORNER = 18
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_LEFT_CORNER = 19
!  Inherit from parent window. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_INHERIT = 100
!  Blank cursor. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_NONE = 101
!  Fullscreen crosshair (if available). 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_FULL_CROSSHAIR = 102

!  GLUT initialization sub-API. 
!  void glutInit(int *argcp, char **argv)
PUBLIC glutInit
INTERFACE glutInit
MODULE PROCEDURE glutInit_f03
END INTERFACE glutInit
INTERFACE
SUBROUTINE glutInit_gl(argcp, argv) BIND(C,NAME="glutInit")
IMPORT
! INTEGER(GLint), DIMENSION(*) :: argcp
INTEGER(GLint) :: argcp
TYPE(C_PTR), INTENT(IN) :: argv
END SUBROUTINE glutInit_gl
END INTERFACE

!  void glutInitDisplayMode(unsigned int mode)
PUBLIC glutInitDisplayMode
INTERFACE
SUBROUTINE glutInitDisplayMode(mode) BIND(C,NAME="glutInitDisplayMode")
IMPORT
INTEGER(GLuint), VALUE :: mode
END SUBROUTINE glutInitDisplayMode
END INTERFACE

!  void glutInitDisplayString(const char *string)
PUBLIC glutInitDisplayString
INTERFACE
SUBROUTINE glutInitDisplayString(string) BIND(C,NAME="glutInitDisplayString")
IMPORT
! CHARACTER, INTENT(IN) :: string
CHARACTER, DIMENSION(*), INTENT(IN) :: string
END SUBROUTINE glutInitDisplayString
END INTERFACE

!  void glutInitWindowPosition(int x, int y)
PUBLIC glutInitWindowPosition
INTERFACE
SUBROUTINE glutInitWindowPosition(x, y) BIND(C,NAME="glutInitWindowPosition")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutInitWindowPosition
END INTERFACE

!  void glutInitWindowSize(int width, int height)
PUBLIC glutInitWindowSize
INTERFACE
SUBROUTINE glutInitWindowSize(width, height) BIND(C,NAME="glutInitWindowSize")
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE glutInitWindowSize
END INTERFACE

!  void glutMainLoop(void)
PUBLIC glutMainLoop
INTERFACE
SUBROUTINE glutMainLoop() BIND(C,NAME="glutMainLoop")
IMPORT
END SUBROUTINE glutMainLoop
END INTERFACE


!  GLUT window sub-API. 
!  int glutCreateWindow(const char *title)
PUBLIC glutCreateWindow
INTERFACE
FUNCTION glutCreateWindow(title) BIND(C,NAME="glutCreateWindow")
IMPORT
INTEGER(GLint) :: glutCreateWindow
! CHARACTER, INTENT(IN) :: title
CHARACTER, DIMENSION(*), INTENT(IN) :: title
END FUNCTION glutCreateWindow
END INTERFACE

!  int glutCreateSubWindow(int win, int x, int y, int width, int height)
PUBLIC glutCreateSubWindow
INTERFACE
FUNCTION glutCreateSubWindow(win, x, y, width, height) BIND(C,NAME="glutCreateSubWindow")
IMPORT
INTEGER(GLint) :: glutCreateSubWindow
INTEGER(GLint), VALUE :: win, x, y, width, height
END FUNCTION glutCreateSubWindow
END INTERFACE

!  void glutDestroyWindow(int win)
PUBLIC glutDestroyWindow
INTERFACE
SUBROUTINE glutDestroyWindow(win) BIND(C,NAME="glutDestroyWindow")
IMPORT
INTEGER(GLint), VALUE :: win
END SUBROUTINE glutDestroyWindow
END INTERFACE

!  void glutPostRedisplay(void)
PUBLIC glutPostRedisplay
INTERFACE
SUBROUTINE glutPostRedisplay() BIND(C,NAME="glutPostRedisplay")
IMPORT
END SUBROUTINE glutPostRedisplay
END INTERFACE

!  void glutPostWindowRedisplay(int win)
PUBLIC glutPostWindowRedisplay
INTERFACE
SUBROUTINE glutPostWindowRedisplay(win) BIND(C,NAME="glutPostWindowRedisplay")
IMPORT
INTEGER(GLint), VALUE :: win
END SUBROUTINE glutPostWindowRedisplay
END INTERFACE

!  void glutSwapBuffers(void)
PUBLIC glutSwapBuffers
INTERFACE
SUBROUTINE glutSwapBuffers() BIND(C,NAME="glutSwapBuffers")
IMPORT
END SUBROUTINE glutSwapBuffers
END INTERFACE

!  int glutGetWindow(void)
PUBLIC glutGetWindow
INTERFACE
FUNCTION glutGetWindow() BIND(C,NAME="glutGetWindow")
IMPORT
INTEGER(GLint) :: glutGetWindow
END FUNCTION glutGetWindow
END INTERFACE

!  void glutSetWindow(int win)
PUBLIC glutSetWindow
INTERFACE
SUBROUTINE glutSetWindow(win) BIND(C,NAME="glutSetWindow")
IMPORT
INTEGER(GLint), VALUE :: win
END SUBROUTINE glutSetWindow
END INTERFACE

!  void glutSetWindowTitle(const char *title)
PUBLIC glutSetWindowTitle
INTERFACE
SUBROUTINE glutSetWindowTitle(title) BIND(C,NAME="glutSetWindowTitle")
IMPORT
! CHARACTER, INTENT(IN) :: title
CHARACTER, DIMENSION(*), INTENT(IN) :: title
END SUBROUTINE glutSetWindowTitle
END INTERFACE

!  void glutSetIconTitle(const char *title)
PUBLIC glutSetIconTitle
INTERFACE
SUBROUTINE glutSetIconTitle(title) BIND(C,NAME="glutSetIconTitle")
IMPORT
! CHARACTER, INTENT(IN) :: title
CHARACTER, DIMENSION(*), INTENT(IN) :: title
END SUBROUTINE glutSetIconTitle
END INTERFACE

!  void glutPositionWindow(int x, int y)
PUBLIC glutPositionWindow
INTERFACE
SUBROUTINE glutPositionWindow(x, y) BIND(C,NAME="glutPositionWindow")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutPositionWindow
END INTERFACE

!  void glutReshapeWindow(int width, int height)
PUBLIC glutReshapeWindow
INTERFACE
SUBROUTINE glutReshapeWindow(width, height) BIND(C,NAME="glutReshapeWindow")
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE glutReshapeWindow
END INTERFACE

!  void glutPopWindow(void)
PUBLIC glutPopWindow
INTERFACE
SUBROUTINE glutPopWindow() BIND(C,NAME="glutPopWindow")
IMPORT
END SUBROUTINE glutPopWindow
END INTERFACE

!  void glutPushWindow(void)
PUBLIC glutPushWindow
INTERFACE
SUBROUTINE glutPushWindow() BIND(C,NAME="glutPushWindow")
IMPORT
END SUBROUTINE glutPushWindow
END INTERFACE

!  void glutIconifyWindow(void)
PUBLIC glutIconifyWindow
INTERFACE
SUBROUTINE glutIconifyWindow() BIND(C,NAME="glutIconifyWindow")
IMPORT
END SUBROUTINE glutIconifyWindow
END INTERFACE

!  void glutShowWindow(void)
PUBLIC glutShowWindow
INTERFACE
SUBROUTINE glutShowWindow() BIND(C,NAME="glutShowWindow")
IMPORT
END SUBROUTINE glutShowWindow
END INTERFACE

!  void glutHideWindow(void)
PUBLIC glutHideWindow
INTERFACE
SUBROUTINE glutHideWindow() BIND(C,NAME="glutHideWindow")
IMPORT
END SUBROUTINE glutHideWindow
END INTERFACE

!  void glutFullScreen(void)
PUBLIC glutFullScreen
INTERFACE
SUBROUTINE glutFullScreen() BIND(C,NAME="glutFullScreen")
IMPORT
END SUBROUTINE glutFullScreen
END INTERFACE

!  void glutSetCursor(int cursor)
PUBLIC glutSetCursor
INTERFACE
SUBROUTINE glutSetCursor(cursor) BIND(C,NAME="glutSetCursor")
IMPORT
INTEGER(GLint), VALUE :: cursor
END SUBROUTINE glutSetCursor
END INTERFACE

!  void glutWarpPointer(int x, int y)
PUBLIC glutWarpPointer
INTERFACE
SUBROUTINE glutWarpPointer(x, y) BIND(C,NAME="glutWarpPointer")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutWarpPointer
END INTERFACE


!  GLUT overlay sub-API. 
!  void glutEstablishOverlay(void)
PUBLIC glutEstablishOverlay
INTERFACE
SUBROUTINE glutEstablishOverlay() BIND(C,NAME="glutEstablishOverlay")
IMPORT
END SUBROUTINE glutEstablishOverlay
END INTERFACE

!  void glutRemoveOverlay(void)
PUBLIC glutRemoveOverlay
INTERFACE
SUBROUTINE glutRemoveOverlay() BIND(C,NAME="glutRemoveOverlay")
IMPORT
END SUBROUTINE glutRemoveOverlay
END INTERFACE

!  void glutUseLayer(GLenum layer)
PUBLIC glutUseLayer
INTERFACE
SUBROUTINE glutUseLayer(layer) BIND(C,NAME="glutUseLayer")
IMPORT
INTEGER(GLenum), VALUE :: layer
END SUBROUTINE glutUseLayer
END INTERFACE

!  void glutPostOverlayRedisplay(void)
PUBLIC glutPostOverlayRedisplay
INTERFACE
SUBROUTINE glutPostOverlayRedisplay() BIND(C,NAME="glutPostOverlayRedisplay")
IMPORT
END SUBROUTINE glutPostOverlayRedisplay
END INTERFACE

!  void glutPostWindowOverlayRedisplay(int win)
PUBLIC glutPostWindowOverlayRedisplay
INTERFACE
SUBROUTINE glutPostWindowOverlayRedisplay(win) BIND(C,NAME="glutPostWindowOverlayRedisplay")
IMPORT
INTEGER(GLint), VALUE :: win
END SUBROUTINE glutPostWindowOverlayRedisplay
END INTERFACE

!  void glutShowOverlay(void)
PUBLIC glutShowOverlay
INTERFACE
SUBROUTINE glutShowOverlay() BIND(C,NAME="glutShowOverlay")
IMPORT
END SUBROUTINE glutShowOverlay
END INTERFACE

!  void glutHideOverlay(void)
PUBLIC glutHideOverlay
INTERFACE
SUBROUTINE glutHideOverlay() BIND(C,NAME="glutHideOverlay")
IMPORT
END SUBROUTINE glutHideOverlay
END INTERFACE


!  GLUT menu sub-API. 
!  int glutCreateMenu(void (GLUTCALLBACK *func)(int item))
PUBLIC glutCreateMenu
INTERFACE
FUNCTION glutCreateMenu(func) BIND(C,NAME="glutCreateMenu")
IMPORT
INTEGER(GLint) :: glutCreateMenu
!  void func(int item)
INTERFACE
SUBROUTINE func(item) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: item
END SUBROUTINE func
END INTERFACE
END FUNCTION glutCreateMenu
END INTERFACE

!  void glutDestroyMenu(int menu)
PUBLIC glutDestroyMenu
INTERFACE
SUBROUTINE glutDestroyMenu(menu) BIND(C,NAME="glutDestroyMenu")
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE glutDestroyMenu
END INTERFACE

!  int glutGetMenu(void)
PUBLIC glutGetMenu
INTERFACE
FUNCTION glutGetMenu() BIND(C,NAME="glutGetMenu")
IMPORT
INTEGER(GLint) :: glutGetMenu
END FUNCTION glutGetMenu
END INTERFACE

!  void glutSetMenu(int menu)
PUBLIC glutSetMenu
INTERFACE
SUBROUTINE glutSetMenu(menu) BIND(C,NAME="glutSetMenu")
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE glutSetMenu
END INTERFACE

!  void glutAddMenuEntry(const char *label, int value)
PUBLIC glutAddMenuEntry
INTERFACE
SUBROUTINE glutAddMenuEntry(label, value) BIND(C,NAME="glutAddMenuEntry")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER, DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: value
END SUBROUTINE glutAddMenuEntry
END INTERFACE

!  void glutAddSubMenu(const char *label, int submenu)
PUBLIC glutAddSubMenu
INTERFACE
SUBROUTINE glutAddSubMenu(label, submenu) BIND(C,NAME="glutAddSubMenu")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER, DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: submenu
END SUBROUTINE glutAddSubMenu
END INTERFACE

!  void glutChangeToMenuEntry(int item, const char *label, int value)
PUBLIC glutChangeToMenuEntry
INTERFACE
SUBROUTINE glutChangeToMenuEntry(item, label, value) BIND(C,NAME="glutChangeToMenuEntry")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER, DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: item, value
END SUBROUTINE glutChangeToMenuEntry
END INTERFACE

!  void glutChangeToSubMenu(int item, const char *label, int submenu)
PUBLIC glutChangeToSubMenu
INTERFACE
SUBROUTINE glutChangeToSubMenu(item, label, submenu) BIND(C,NAME="glutChangeToSubMenu")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER, DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: item, submenu
END SUBROUTINE glutChangeToSubMenu
END INTERFACE

!  void glutRemoveMenuItem(int item)
PUBLIC glutRemoveMenuItem
INTERFACE
SUBROUTINE glutRemoveMenuItem(item) BIND(C,NAME="glutRemoveMenuItem")
IMPORT
INTEGER(GLint), VALUE :: item
END SUBROUTINE glutRemoveMenuItem
END INTERFACE

!  void glutAttachMenu(int button)
PUBLIC glutAttachMenu
INTERFACE
SUBROUTINE glutAttachMenu(button) BIND(C,NAME="glutAttachMenu")
IMPORT
INTEGER(GLint), VALUE :: button
END SUBROUTINE glutAttachMenu
END INTERFACE

!  void glutDetachMenu(int button)
PUBLIC glutDetachMenu
INTERFACE
SUBROUTINE glutDetachMenu(button) BIND(C,NAME="glutDetachMenu")
IMPORT
INTEGER(GLint), VALUE :: button
END SUBROUTINE glutDetachMenu
END INTERFACE


!  GLUT window callback sub-API. 
!  void glutDisplayFunc(void (GLUTCALLBACK *func)(void))
PUBLIC glutDisplayFunc
INTERFACE glutDisplayFunc
MODULE PROCEDURE glutDisplayFunc_f03
END INTERFACE glutDisplayFunc
INTERFACE
SUBROUTINE glutDisplayFunc_gl(func) BIND(C,NAME="glutDisplayFunc")
IMPORT
!  void func(void)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutDisplayFunc_gl
END INTERFACE

!  void glutReshapeFunc(void (GLUTCALLBACK *func)(int width, int height))
PUBLIC glutReshapeFunc
INTERFACE glutReshapeFunc
MODULE PROCEDURE glutReshapeFunc_f03
END INTERFACE glutReshapeFunc
INTERFACE
SUBROUTINE glutReshapeFunc_gl(func) BIND(C,NAME="glutReshapeFunc")
IMPORT
!  void func(int width, int height)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutReshapeFunc_gl
END INTERFACE

!  void glutKeyboardFunc(void (GLUTCALLBACK *func)(unsigned char key, int x, int y))
PUBLIC glutKeyboardFunc
INTERFACE glutKeyboardFunc
MODULE PROCEDURE glutKeyboardFunc_f03
END INTERFACE glutKeyboardFunc
INTERFACE
SUBROUTINE glutKeyboardFunc_gl(func) BIND(C,NAME="glutKeyboardFunc")
IMPORT
!  void func(unsigned char key, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutKeyboardFunc_gl
END INTERFACE

!  void glutMouseFunc(void (GLUTCALLBACK *func)(int button, int state, int x, int y))
PUBLIC glutMouseFunc
INTERFACE glutMouseFunc
MODULE PROCEDURE glutMouseFunc_f03
END INTERFACE glutMouseFunc
INTERFACE
SUBROUTINE glutMouseFunc_gl(func) BIND(C,NAME="glutMouseFunc")
IMPORT
!  void func(int button, int state, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutMouseFunc_gl
END INTERFACE

!  void glutMotionFunc(void (GLUTCALLBACK *func)(int x, int y))
PUBLIC glutMotionFunc
INTERFACE glutMotionFunc
MODULE PROCEDURE glutMotionFunc_f03
END INTERFACE glutMotionFunc
INTERFACE
SUBROUTINE glutMotionFunc_gl(func) BIND(C,NAME="glutMotionFunc")
IMPORT
!  void func(int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutMotionFunc_gl
END INTERFACE

!  void glutPassiveMotionFunc(void (GLUTCALLBACK *func)(int x, int y))
PUBLIC glutPassiveMotionFunc
INTERFACE glutPassiveMotionFunc
MODULE PROCEDURE glutPassiveMotionFunc_f03
END INTERFACE glutPassiveMotionFunc
INTERFACE
SUBROUTINE glutPassiveMotionFunc_gl(func) BIND(C,NAME="glutPassiveMotionFunc")
IMPORT
!  void func(int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutPassiveMotionFunc_gl
END INTERFACE

!  void glutEntryFunc(void (GLUTCALLBACK *func)(int state))
PUBLIC glutEntryFunc
INTERFACE glutEntryFunc
MODULE PROCEDURE glutEntryFunc_f03
END INTERFACE glutEntryFunc
INTERFACE
SUBROUTINE glutEntryFunc_gl(func) BIND(C,NAME="glutEntryFunc")
IMPORT
!  void func(int state)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutEntryFunc_gl
END INTERFACE

!  void glutVisibilityFunc(void (GLUTCALLBACK *func)(int state))
PUBLIC glutVisibilityFunc
INTERFACE glutVisibilityFunc
MODULE PROCEDURE glutVisibilityFunc_f03
END INTERFACE glutVisibilityFunc
INTERFACE
SUBROUTINE glutVisibilityFunc_gl(func) BIND(C,NAME="glutVisibilityFunc")
IMPORT
!  void func(int state)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutVisibilityFunc_gl
END INTERFACE

!  void glutIdleFunc(void (GLUTCALLBACK *func)(void))
PUBLIC glutIdleFunc
INTERFACE glutIdleFunc
MODULE PROCEDURE glutIdleFunc_f03
END INTERFACE glutIdleFunc
INTERFACE
SUBROUTINE glutIdleFunc_gl(func) BIND(C,NAME="glutIdleFunc")
IMPORT
!  void func(void)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutIdleFunc_gl
END INTERFACE

!  void glutTimerFunc(unsigned int millis, void (GLUTCALLBACK *func)(int value), int value)
PUBLIC glutTimerFunc
INTERFACE glutTimerFunc
MODULE PROCEDURE glutTimerFunc_f03
END INTERFACE glutTimerFunc
INTERFACE
SUBROUTINE glutTimerFunc_gl(millis, func, value) BIND(C,NAME="glutTimerFunc")
IMPORT
INTEGER(GLint), VALUE :: value
INTEGER(GLuint), VALUE :: millis
!  void func(int value)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutTimerFunc_gl
END INTERFACE

!  void glutMenuStateFunc(void (GLUTCALLBACK *func)(int state))
PUBLIC glutMenuStateFunc
INTERFACE glutMenuStateFunc
MODULE PROCEDURE glutMenuStateFunc_f03
END INTERFACE glutMenuStateFunc
INTERFACE
SUBROUTINE glutMenuStateFunc_gl(func) BIND(C,NAME="glutMenuStateFunc")
IMPORT
!  void func(int state)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutMenuStateFunc_gl
END INTERFACE

!  void glutSpecialFunc(void (GLUTCALLBACK *func)(int key, int x, int y))
PUBLIC glutSpecialFunc
INTERFACE glutSpecialFunc
MODULE PROCEDURE glutSpecialFunc_f03
END INTERFACE glutSpecialFunc
INTERFACE
SUBROUTINE glutSpecialFunc_gl(func) BIND(C,NAME="glutSpecialFunc")
IMPORT
!  void func(int key, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutSpecialFunc_gl
END INTERFACE

!  void glutSpaceballMotionFunc(void (GLUTCALLBACK *func)(int x, int y, int z))
PUBLIC glutSpaceballMotionFunc
INTERFACE glutSpaceballMotionFunc
MODULE PROCEDURE glutSpaceballMotionFunc_f03
END INTERFACE glutSpaceballMotionFunc
INTERFACE
SUBROUTINE glutSpaceballMotionFunc_gl(func) BIND(C,NAME="glutSpaceballMotionFunc")
IMPORT
!  void func(int x, int y, int z)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutSpaceballMotionFunc_gl
END INTERFACE

!  void glutSpaceballRotateFunc(void (GLUTCALLBACK *func)(int x, int y, int z))
PUBLIC glutSpaceballRotateFunc
INTERFACE glutSpaceballRotateFunc
MODULE PROCEDURE glutSpaceballRotateFunc_f03
END INTERFACE glutSpaceballRotateFunc
INTERFACE
SUBROUTINE glutSpaceballRotateFunc_gl(func) BIND(C,NAME="glutSpaceballRotateFunc")
IMPORT
!  void func(int x, int y, int z)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutSpaceballRotateFunc_gl
END INTERFACE

!  void glutSpaceballButtonFunc(void (GLUTCALLBACK *func)(int button, int state))
PUBLIC glutSpaceballButtonFunc
INTERFACE glutSpaceballButtonFunc
MODULE PROCEDURE glutSpaceballButtonFunc_f03
END INTERFACE glutSpaceballButtonFunc
INTERFACE
SUBROUTINE glutSpaceballButtonFunc_gl(func) BIND(C,NAME="glutSpaceballButtonFunc")
IMPORT
!  void func(int button, int state)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutSpaceballButtonFunc_gl
END INTERFACE

!  void glutButtonBoxFunc(void (GLUTCALLBACK *func)(int button, int state))
PUBLIC glutButtonBoxFunc
INTERFACE glutButtonBoxFunc
MODULE PROCEDURE glutButtonBoxFunc_f03
END INTERFACE glutButtonBoxFunc
INTERFACE
SUBROUTINE glutButtonBoxFunc_gl(func) BIND(C,NAME="glutButtonBoxFunc")
IMPORT
!  void func(int button, int state)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutButtonBoxFunc_gl
END INTERFACE

!  void glutDialsFunc(void (GLUTCALLBACK *func)(int dial, int value))
PUBLIC glutDialsFunc
INTERFACE glutDialsFunc
MODULE PROCEDURE glutDialsFunc_f03
END INTERFACE glutDialsFunc
INTERFACE
SUBROUTINE glutDialsFunc_gl(func) BIND(C,NAME="glutDialsFunc")
IMPORT
!  void func(int dial, int value)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutDialsFunc_gl
END INTERFACE

!  void glutTabletMotionFunc(void (GLUTCALLBACK *func)(int x, int y))
PUBLIC glutTabletMotionFunc
INTERFACE glutTabletMotionFunc
MODULE PROCEDURE glutTabletMotionFunc_f03
END INTERFACE glutTabletMotionFunc
INTERFACE
SUBROUTINE glutTabletMotionFunc_gl(func) BIND(C,NAME="glutTabletMotionFunc")
IMPORT
!  void func(int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutTabletMotionFunc_gl
END INTERFACE

!  void glutTabletButtonFunc(void (GLUTCALLBACK *func)(int button, int state, int x, int y))
PUBLIC glutTabletButtonFunc
INTERFACE glutTabletButtonFunc
MODULE PROCEDURE glutTabletButtonFunc_f03
END INTERFACE glutTabletButtonFunc
INTERFACE
SUBROUTINE glutTabletButtonFunc_gl(func) BIND(C,NAME="glutTabletButtonFunc")
IMPORT
!  void func(int button, int state, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutTabletButtonFunc_gl
END INTERFACE

!  void glutMenuStatusFunc(void (GLUTCALLBACK *func)(int status, int x, int y))
PUBLIC glutMenuStatusFunc
INTERFACE glutMenuStatusFunc
MODULE PROCEDURE glutMenuStatusFunc_f03
END INTERFACE glutMenuStatusFunc
INTERFACE
SUBROUTINE glutMenuStatusFunc_gl(func) BIND(C,NAME="glutMenuStatusFunc")
IMPORT
!  void func(int status, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutMenuStatusFunc_gl
END INTERFACE

!  void glutOverlayDisplayFunc(void (GLUTCALLBACK *func)(void))
PUBLIC glutOverlayDisplayFunc
INTERFACE glutOverlayDisplayFunc
MODULE PROCEDURE glutOverlayDisplayFunc_f03
END INTERFACE glutOverlayDisplayFunc
INTERFACE
SUBROUTINE glutOverlayDisplayFunc_gl(func) BIND(C,NAME="glutOverlayDisplayFunc")
IMPORT
!  void func(void)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutOverlayDisplayFunc_gl
END INTERFACE

!  void glutWindowStatusFunc(void (GLUTCALLBACK *func)(int state))
PUBLIC glutWindowStatusFunc
INTERFACE glutWindowStatusFunc
MODULE PROCEDURE glutWindowStatusFunc_f03
END INTERFACE glutWindowStatusFunc
INTERFACE
SUBROUTINE glutWindowStatusFunc_gl(func) BIND(C,NAME="glutWindowStatusFunc")
IMPORT
!  void func(int state)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutWindowStatusFunc_gl
END INTERFACE

!  void glutKeyboardUpFunc(void (GLUTCALLBACK *func)(unsigned char key, int x, int y))
PUBLIC glutKeyboardUpFunc
INTERFACE glutKeyboardUpFunc
MODULE PROCEDURE glutKeyboardUpFunc_f03
END INTERFACE glutKeyboardUpFunc
INTERFACE
SUBROUTINE glutKeyboardUpFunc_gl(func) BIND(C,NAME="glutKeyboardUpFunc")
IMPORT
!  void func(unsigned char key, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutKeyboardUpFunc_gl
END INTERFACE

!  void glutSpecialUpFunc(void (GLUTCALLBACK *func)(int key, int x, int y))
PUBLIC glutSpecialUpFunc
INTERFACE glutSpecialUpFunc
MODULE PROCEDURE glutSpecialUpFunc_f03
END INTERFACE glutSpecialUpFunc
INTERFACE
SUBROUTINE glutSpecialUpFunc_gl(func) BIND(C,NAME="glutSpecialUpFunc")
IMPORT
!  void func(int key, int x, int y)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutSpecialUpFunc_gl
END INTERFACE

!  void glutJoystickFunc(void (GLUTCALLBACK *func)(unsigned int buttonMask, int x, int y, int z), int pollInterval)
PUBLIC glutJoystickFunc
INTERFACE glutJoystickFunc
MODULE PROCEDURE glutJoystickFunc_f03
END INTERFACE glutJoystickFunc
INTERFACE
SUBROUTINE glutJoystickFunc_gl(func, pollInterval) BIND(C,NAME="glutJoystickFunc")
IMPORT
INTEGER(GLint), VALUE :: pollInterval
!  void func(unsigned int buttonMask, int x, int y, int z)
TYPE(C_FUNPTR), VALUE :: func
END SUBROUTINE glutJoystickFunc_gl
END INTERFACE


!  GLUT color index sub-API. 
!  void glutSetColor(int ndx, GLfloat red, GLfloat green, GLfloat blue)
PUBLIC glutSetColor
INTERFACE
SUBROUTINE glutSetColor(ndx, red, green, blue) BIND(C,NAME="glutSetColor")
IMPORT
INTEGER(GLint), VALUE :: ndx
REAL(GLfloat), VALUE :: red, green, blue
END SUBROUTINE glutSetColor
END INTERFACE

!  GLfloat glutGetColor(int ndx, int component)
PUBLIC glutGetColor
INTERFACE
FUNCTION glutGetColor(ndx, component) BIND(C,NAME="glutGetColor")
IMPORT
REAL(GLfloat) :: glutGetColor
INTEGER(GLint), VALUE :: ndx, component
END FUNCTION glutGetColor
END INTERFACE

!  void glutCopyColormap(int win)
PUBLIC glutCopyColormap
INTERFACE
SUBROUTINE glutCopyColormap(win) BIND(C,NAME="glutCopyColormap")
IMPORT
INTEGER(GLint), VALUE :: win
END SUBROUTINE glutCopyColormap
END INTERFACE


!  GLUT state retrieval sub-API. 
!  int glutGet(GLenum type)
PUBLIC glutGet
INTERFACE
FUNCTION glutGet(type) BIND(C,NAME="glutGet")
IMPORT
INTEGER(GLint) :: glutGet
INTEGER(GLenum), VALUE :: type
END FUNCTION glutGet
END INTERFACE

!  int glutDeviceGet(GLenum type)
PUBLIC glutDeviceGet
INTERFACE
FUNCTION glutDeviceGet(type) BIND(C,NAME="glutDeviceGet")
IMPORT
INTEGER(GLint) :: glutDeviceGet
INTEGER(GLenum), VALUE :: type
END FUNCTION glutDeviceGet
END INTERFACE

!  GLUT extension support sub-API 
!  int glutExtensionSupported(const char *name)
PUBLIC glutExtensionSupported
INTERFACE
FUNCTION glutExtensionSupported(name) BIND(C,NAME="glutExtensionSupported")
IMPORT
INTEGER(GLint) :: glutExtensionSupported
! CHARACTER, INTENT(IN) :: name
CHARACTER, DIMENSION(*), INTENT(IN) :: name
END FUNCTION glutExtensionSupported
END INTERFACE

!  int glutGetModifiers(void)
PUBLIC glutGetModifiers
INTERFACE
FUNCTION glutGetModifiers() BIND(C,NAME="glutGetModifiers")
IMPORT
INTEGER(GLint) :: glutGetModifiers
END FUNCTION glutGetModifiers
END INTERFACE

!  int glutLayerGet(GLenum type)
PUBLIC glutLayerGet
INTERFACE
FUNCTION glutLayerGet(type) BIND(C,NAME="glutLayerGet")
IMPORT
INTEGER(GLint) :: glutLayerGet
INTEGER(GLenum), VALUE :: type
END FUNCTION glutLayerGet
END INTERFACE

!  void * glutGetProcAddress(const char *procName)
PUBLIC glutGetProcAddress
INTERFACE
FUNCTION glutGetProcAddress(procName) BIND(C,NAME="glutGetProcAddress")
IMPORT
TYPE(C_PTR) :: glutGetProcAddress
! CHARACTER, INTENT(IN) :: procName
CHARACTER, DIMENSION(*), INTENT(IN) :: procName
END FUNCTION glutGetProcAddress
END INTERFACE


!  GLUT font sub-API 
!  void glutBitmapCharacter(void *font, int character)
PUBLIC glutBitmapCharacter
INTERFACE
SUBROUTINE glutBitmapCharacter(font, character) BIND(C,NAME="glutBitmapCharacter")
IMPORT
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutBitmapCharacter
END INTERFACE

!  int glutBitmapWidth(void *font, int character)
PUBLIC glutBitmapWidth
INTERFACE
FUNCTION glutBitmapWidth(font, character) BIND(C,NAME="glutBitmapWidth")
IMPORT
INTEGER(GLint) :: glutBitmapWidth
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapWidth
END INTERFACE

!  void glutStrokeCharacter(void *font, int character)
PUBLIC glutStrokeCharacter
INTERFACE
SUBROUTINE glutStrokeCharacter(font, character) BIND(C,NAME="glutStrokeCharacter")
IMPORT
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutStrokeCharacter
END INTERFACE

!  int glutStrokeWidth(void *font, int character)
PUBLIC glutStrokeWidth
INTERFACE
FUNCTION glutStrokeWidth(font, character) BIND(C,NAME="glutStrokeWidth")
IMPORT
INTEGER(GLint) :: glutStrokeWidth
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeWidth
END INTERFACE

!  int glutBitmapLength(void *font, const unsigned char *string)
PUBLIC glutBitmapLength
INTERFACE
FUNCTION glutBitmapLength(font, string) BIND(C,NAME="glutBitmapLength")
IMPORT
INTEGER(GLint) :: glutBitmapLength
! CHARACTER, INTENT(IN) :: string
CHARACTER, DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapLength
END INTERFACE

!  int glutStrokeLength(void *font, const unsigned char *string)
PUBLIC glutStrokeLength
INTERFACE
FUNCTION glutStrokeLength(font, string) BIND(C,NAME="glutStrokeLength")
IMPORT
INTEGER(GLint) :: glutStrokeLength
! CHARACTER, INTENT(IN) :: string
CHARACTER, DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeLength
END INTERFACE


!  GLUT pre-built models sub-API 
!  void glutWireSphere(GLdouble radius, GLint slices, GLint stacks)
PUBLIC glutWireSphere
INTERFACE
SUBROUTINE glutWireSphere(radius, slices, stacks) BIND(C,NAME="glutWireSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
END SUBROUTINE glutWireSphere
END INTERFACE

!  void glutSolidSphere(GLdouble radius, GLint slices, GLint stacks)
PUBLIC glutSolidSphere
INTERFACE
SUBROUTINE glutSolidSphere(radius, slices, stacks) BIND(C,NAME="glutSolidSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
END SUBROUTINE glutSolidSphere
END INTERFACE

!  void glutWireCone(GLdouble base, GLdouble height, GLint slices, GLint stacks)
PUBLIC glutWireCone
INTERFACE
SUBROUTINE glutWireCone(base, height, slices, stacks) BIND(C,NAME="glutWireCone")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, height
END SUBROUTINE glutWireCone
END INTERFACE

!  void glutSolidCone(GLdouble base, GLdouble height, GLint slices, GLint stacks)
PUBLIC glutSolidCone
INTERFACE
SUBROUTINE glutSolidCone(base, height, slices, stacks) BIND(C,NAME="glutSolidCone")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, height
END SUBROUTINE glutSolidCone
END INTERFACE

!  void glutWireCube(GLdouble size)
PUBLIC glutWireCube
INTERFACE
SUBROUTINE glutWireCube(size) BIND(C,NAME="glutWireCube")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutWireCube
END INTERFACE

!  void glutSolidCube(GLdouble size)
PUBLIC glutSolidCube
INTERFACE
SUBROUTINE glutSolidCube(size) BIND(C,NAME="glutSolidCube")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutSolidCube
END INTERFACE

!  void glutWireTorus(GLdouble innerRadius, GLdouble outerRadius, GLint sides, GLint rings)
PUBLIC glutWireTorus
INTERFACE
SUBROUTINE glutWireTorus(innerRadius, outerRadius, sides, rings) BIND(C,NAME="glutWireTorus")
IMPORT
INTEGER(GLint), VALUE :: sides, rings
REAL(GLdouble), VALUE :: innerRadius, outerRadius
END SUBROUTINE glutWireTorus
END INTERFACE

!  void glutSolidTorus(GLdouble innerRadius, GLdouble outerRadius, GLint sides, GLint rings)
PUBLIC glutSolidTorus
INTERFACE
SUBROUTINE glutSolidTorus(innerRadius, outerRadius, sides, rings) BIND(C,NAME="glutSolidTorus")
IMPORT
INTEGER(GLint), VALUE :: sides, rings
REAL(GLdouble), VALUE :: innerRadius, outerRadius
END SUBROUTINE glutSolidTorus
END INTERFACE

!  void glutWireDodecahedron(void)
PUBLIC glutWireDodecahedron
INTERFACE
SUBROUTINE glutWireDodecahedron() BIND(C,NAME="glutWireDodecahedron")
IMPORT
END SUBROUTINE glutWireDodecahedron
END INTERFACE

!  void glutSolidDodecahedron(void)
PUBLIC glutSolidDodecahedron
INTERFACE
SUBROUTINE glutSolidDodecahedron() BIND(C,NAME="glutSolidDodecahedron")
IMPORT
END SUBROUTINE glutSolidDodecahedron
END INTERFACE

!  void glutWireTeapot(GLdouble size)
PUBLIC glutWireTeapot
INTERFACE
SUBROUTINE glutWireTeapot(size) BIND(C,NAME="glutWireTeapot")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutWireTeapot
END INTERFACE

!  void glutSolidTeapot(GLdouble size)
PUBLIC glutSolidTeapot
INTERFACE
SUBROUTINE glutSolidTeapot(size) BIND(C,NAME="glutSolidTeapot")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutSolidTeapot
END INTERFACE

!  void glutWireOctahedron(void)
PUBLIC glutWireOctahedron
INTERFACE
SUBROUTINE glutWireOctahedron() BIND(C,NAME="glutWireOctahedron")
IMPORT
END SUBROUTINE glutWireOctahedron
END INTERFACE

!  void glutSolidOctahedron(void)
PUBLIC glutSolidOctahedron
INTERFACE
SUBROUTINE glutSolidOctahedron() BIND(C,NAME="glutSolidOctahedron")
IMPORT
END SUBROUTINE glutSolidOctahedron
END INTERFACE

!  void glutWireTetrahedron(void)
PUBLIC glutWireTetrahedron
INTERFACE
SUBROUTINE glutWireTetrahedron() BIND(C,NAME="glutWireTetrahedron")
IMPORT
END SUBROUTINE glutWireTetrahedron
END INTERFACE

!  void glutSolidTetrahedron(void)
PUBLIC glutSolidTetrahedron
INTERFACE
SUBROUTINE glutSolidTetrahedron() BIND(C,NAME="glutSolidTetrahedron")
IMPORT
END SUBROUTINE glutSolidTetrahedron
END INTERFACE

!  void glutWireIcosahedron(void)
PUBLIC glutWireIcosahedron
INTERFACE
SUBROUTINE glutWireIcosahedron() BIND(C,NAME="glutWireIcosahedron")
IMPORT
END SUBROUTINE glutWireIcosahedron
END INTERFACE

!  void glutSolidIcosahedron(void)
PUBLIC glutSolidIcosahedron
INTERFACE
SUBROUTINE glutSolidIcosahedron() BIND(C,NAME="glutSolidIcosahedron")
IMPORT
END SUBROUTINE glutSolidIcosahedron
END INTERFACE


!  GLUT video resize sub-API. 
!  int glutVideoResizeGet(GLenum param)
PUBLIC glutVideoResizeGet
INTERFACE
FUNCTION glutVideoResizeGet(param) BIND(C,NAME="glutVideoResizeGet")
IMPORT
INTEGER(GLint) :: glutVideoResizeGet
INTEGER(GLenum), VALUE :: param
END FUNCTION glutVideoResizeGet
END INTERFACE

!  void glutSetupVideoResizing(void)
PUBLIC glutSetupVideoResizing
INTERFACE
SUBROUTINE glutSetupVideoResizing() BIND(C,NAME="glutSetupVideoResizing")
IMPORT
END SUBROUTINE glutSetupVideoResizing
END INTERFACE

!  void glutStopVideoResizing(void)
PUBLIC glutStopVideoResizing
INTERFACE
SUBROUTINE glutStopVideoResizing() BIND(C,NAME="glutStopVideoResizing")
IMPORT
END SUBROUTINE glutStopVideoResizing
END INTERFACE

!  void glutVideoResize(int x, int y, int width, int height)
PUBLIC glutVideoResize
INTERFACE
SUBROUTINE glutVideoResize(x, y, width, height) BIND(C,NAME="glutVideoResize")
IMPORT
INTEGER(GLint), VALUE :: x, y, width, height
END SUBROUTINE glutVideoResize
END INTERFACE

!  void glutVideoPan(int x, int y, int width, int height)
PUBLIC glutVideoPan
INTERFACE
SUBROUTINE glutVideoPan(x, y, width, height) BIND(C,NAME="glutVideoPan")
IMPORT
INTEGER(GLint), VALUE :: x, y, width, height
END SUBROUTINE glutVideoPan
END INTERFACE


!  GLUT debugging sub-API. 
!  void glutReportErrors(void)
PUBLIC glutReportErrors
INTERFACE
SUBROUTINE glutReportErrors() BIND(C,NAME="glutReportErrors")
IMPORT
END SUBROUTINE glutReportErrors
END INTERFACE


!  GLUT device control sub-API. 
!  glutSetKeyRepeat modes. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_OFF = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_ON = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_DEFAULT = 2

!  Joystick button masks. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_A = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_B = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_C = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_D = 8

!  void glutIgnoreKeyRepeat(int ignore)
PUBLIC glutIgnoreKeyRepeat
INTERFACE
SUBROUTINE glutIgnoreKeyRepeat(ignore) BIND(C,NAME="glutIgnoreKeyRepeat")
IMPORT
INTEGER(GLint), VALUE :: ignore
END SUBROUTINE glutIgnoreKeyRepeat
END INTERFACE

!  void glutSetKeyRepeat(int repeatMode)
PUBLIC glutSetKeyRepeat
INTERFACE
SUBROUTINE glutSetKeyRepeat(repeatMode) BIND(C,NAME="glutSetKeyRepeat")
IMPORT
INTEGER(GLint), VALUE :: repeatMode
END SUBROUTINE glutSetKeyRepeat
END INTERFACE

!  void glutForceJoystickFunc(void)
PUBLIC glutForceJoystickFunc
INTERFACE
SUBROUTINE glutForceJoystickFunc() BIND(C,NAME="glutForceJoystickFunc")
IMPORT
END SUBROUTINE glutForceJoystickFunc
END INTERFACE


!  GLUT game mode sub-API. 
!  glutGameModeGet. 
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_ACTIVE = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_POSSIBLE = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_WIDTH = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_HEIGHT = 3
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_PIXEL_DEPTH = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_REFRESH_RATE = 5
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_DISPLAY_CHANGED = 6

!  void glutGameModeString(const char *string)
PUBLIC glutGameModeString
INTERFACE
SUBROUTINE glutGameModeString(string) BIND(C,NAME="glutGameModeString")
IMPORT
! CHARACTER, INTENT(IN) :: string
CHARACTER, DIMENSION(*), INTENT(IN) :: string
END SUBROUTINE glutGameModeString
END INTERFACE

!  int glutEnterGameMode(void)
PUBLIC glutEnterGameMode
INTERFACE
FUNCTION glutEnterGameMode() BIND(C,NAME="glutEnterGameMode")
IMPORT
INTEGER(GLint) :: glutEnterGameMode
END FUNCTION glutEnterGameMode
END INTERFACE

!  void glutLeaveGameMode(void)
PUBLIC glutLeaveGameMode
INTERFACE
SUBROUTINE glutLeaveGameMode() BIND(C,NAME="glutLeaveGameMode")
IMPORT
END SUBROUTINE glutLeaveGameMode
END INTERFACE

!  int glutGameModeGet(GLenum mode)
PUBLIC glutGameModeGet
INTERFACE
FUNCTION glutGameModeGet(mode) BIND(C,NAME="glutGameModeGet")
IMPORT
INTEGER(GLint) :: glutGameModeGet
INTEGER(GLenum), VALUE :: mode
END FUNCTION glutGameModeGet
END INTERFACE



! Font variables in GLUT_fonts.c
TYPE(C_PTR), BIND(C), PUBLIC, PROTECTED :: GLUT_STROKE_ROMAN,         &
    GLUT_STROKE_MONO_ROMAN, GLUT_BITMAP_9_BY_15, GLUT_BITMAP_8_BY_13, &
    GLUT_BITMAP_TIMES_ROMAN_10, GLUT_BITMAP_TIMES_ROMAN_24,           &
    GLUT_BITMAP_HELVETICA_10, GLUT_BITMAP_HELVETICA_12,               &
    GLUT_BITMAP_HELVETICA_18

! A special callback function for compatibility with f90gl
TYPE(C_FUNPTR), PUBLIC, SAVE, BIND(C) :: GLUT_NULL_FUNC=C_NULL_FUNPTR

CONTAINS

SUBROUTINE glutInit_f03()
  INTEGER(C_INT) :: argcp=1
  TYPE(C_PTR), DIMENSION(1), TARGET :: argv=C_NULL_PTR
  CHARACTER(C_CHAR), DIMENSION(1), TARGET :: empty_string=C_NULL_CHAR

  ! A hack
  INTERFACE
   SUBROUTINE SetNullFunc() BIND(C,NAME='SetNullFunc')
   END SUBROUTINE
  END INTERFACE  

  argv(1)=C_LOC(empty_string)
  CALL glutInit_gl(argcp, C_LOC(argv))
  
END SUBROUTINE

SUBROUTINE glutDisplayFunc_f03(func)
INTERFACE
SUBROUTINE func() BIND(C)
IMPORT
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutDisplayFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutDisplayFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutDisplayFunc_f03
SUBROUTINE glutReshapeFunc_f03(func)
INTERFACE
SUBROUTINE func(width, height) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutReshapeFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutReshapeFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutReshapeFunc_f03
SUBROUTINE glutKeyboardFunc_f03(func)
INTERFACE
SUBROUTINE func(key, x, y) BIND(C)
IMPORT
INTEGER(GLbyte), VALUE :: key
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutKeyboardFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutKeyboardFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutKeyboardFunc_f03
SUBROUTINE glutMouseFunc_f03(func)
INTERFACE
SUBROUTINE func(button, state, x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: button, state, x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutMouseFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutMouseFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMouseFunc_f03
SUBROUTINE glutMotionFunc_f03(func)
INTERFACE
SUBROUTINE func(x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutMotionFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMotionFunc_f03
SUBROUTINE glutPassiveMotionFunc_f03(func)
INTERFACE
SUBROUTINE func(x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutPassiveMotionFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutPassiveMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutPassiveMotionFunc_f03
SUBROUTINE glutEntryFunc_f03(func)
INTERFACE
SUBROUTINE func(state) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: state
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutEntryFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutEntryFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutEntryFunc_f03
SUBROUTINE glutVisibilityFunc_f03(func)
INTERFACE
SUBROUTINE func(state) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: state
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutVisibilityFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutVisibilityFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutVisibilityFunc_f03
SUBROUTINE glutIdleFunc_f03(func)
INTERFACE
SUBROUTINE func() BIND(C)
IMPORT
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutIdleFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutIdleFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutIdleFunc_f03
SUBROUTINE glutTimerFunc_f03(millis, func, value)
INTEGER(GLint), VALUE :: value
INTEGER(GLuint), VALUE :: millis
INTERFACE
SUBROUTINE func(value) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: value
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutTimerFunc_gl(millis, C_FUNLOC(func), value)
ELSE
   CALL glutTimerFunc_gl(millis, C_NULL_FUNPTR, value)
END IF
END SUBROUTINE glutTimerFunc_f03
SUBROUTINE glutMenuStateFunc_f03(func)
INTERFACE
SUBROUTINE func(state) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: state
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutMenuStateFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutMenuStateFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMenuStateFunc_f03
SUBROUTINE glutSpecialFunc_f03(func)
INTERFACE
SUBROUTINE func(key, x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: key, x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutSpecialFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutSpecialFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpecialFunc_f03
SUBROUTINE glutSpaceballMotionFunc_f03(func)
INTERFACE
SUBROUTINE func(x, y, z) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: x, y, z
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutSpaceballMotionFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutSpaceballMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpaceballMotionFunc_f03
SUBROUTINE glutSpaceballRotateFunc_f03(func)
INTERFACE
SUBROUTINE func(x, y, z) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: x, y, z
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutSpaceballRotateFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutSpaceballRotateFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpaceballRotateFunc_f03
SUBROUTINE glutSpaceballButtonFunc_f03(func)
INTERFACE
SUBROUTINE func(button, state) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: button, state
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutSpaceballButtonFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutSpaceballButtonFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpaceballButtonFunc_f03
SUBROUTINE glutButtonBoxFunc_f03(func)
INTERFACE
SUBROUTINE func(button, state) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: button, state
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutButtonBoxFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutButtonBoxFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutButtonBoxFunc_f03
SUBROUTINE glutDialsFunc_f03(func)
INTERFACE
SUBROUTINE func(dial, value) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: dial, value
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutDialsFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutDialsFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutDialsFunc_f03
SUBROUTINE glutTabletMotionFunc_f03(func)
INTERFACE
SUBROUTINE func(x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutTabletMotionFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutTabletMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutTabletMotionFunc_f03
SUBROUTINE glutTabletButtonFunc_f03(func)
INTERFACE
SUBROUTINE func(button, state, x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: button, state, x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutTabletButtonFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutTabletButtonFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutTabletButtonFunc_f03
SUBROUTINE glutMenuStatusFunc_f03(func)
INTERFACE
SUBROUTINE func(status, x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: status, x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutMenuStatusFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutMenuStatusFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMenuStatusFunc_f03
SUBROUTINE glutOverlayDisplayFunc_f03(func)
INTERFACE
SUBROUTINE func() BIND(C)
IMPORT
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutOverlayDisplayFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutOverlayDisplayFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutOverlayDisplayFunc_f03
SUBROUTINE glutWindowStatusFunc_f03(func)
INTERFACE
SUBROUTINE func(state) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: state
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutWindowStatusFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutWindowStatusFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutWindowStatusFunc_f03
SUBROUTINE glutKeyboardUpFunc_f03(func)
INTERFACE
SUBROUTINE func(key, x, y) BIND(C)
IMPORT
INTEGER(GLbyte), VALUE :: key
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutKeyboardUpFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutKeyboardUpFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutKeyboardUpFunc_f03
SUBROUTINE glutSpecialUpFunc_f03(func)
INTERFACE
SUBROUTINE func(key, x, y) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: key, x, y
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutSpecialUpFunc_gl(C_FUNLOC(func))
ELSE
   CALL glutSpecialUpFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpecialUpFunc_f03
SUBROUTINE glutJoystickFunc_f03(func, pollInterval)
INTEGER(GLint), VALUE :: pollInterval
INTERFACE
SUBROUTINE func(buttonMask, x, y, z) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: x, y, z
INTEGER(GLuint), VALUE :: buttonMask
END SUBROUTINE func
END INTERFACE
OPTIONAL :: func
IF(PRESENT(func)) THEN
   CALL glutJoystickFunc_gl(C_FUNLOC(func), pollInterval)
ELSE
   CALL glutJoystickFunc_gl(C_NULL_FUNPTR, pollInterval)
END IF
END SUBROUTINE glutJoystickFunc_f03


END MODULE OpenGL_glut

