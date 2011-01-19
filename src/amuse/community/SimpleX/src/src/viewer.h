/* gcc -lm -I/software/local/include -L/software/local/lib -lMesaGL 
-lMesaGLU -lglut -L/usr/X11R6/lib -lXaw -lXt -lXmu -lXi -lXext -lX11 
Untstarsed2.cpp
*/

#ifndef VIEWER_H
#define VIEWER_H

#include <stdio.h>           /* Standard C/C++ Input-Output */
#include <stdlib.h>          /* Standard C/C++ Library */
#include "GL/freeglut.h"         /* The GL Utility Toolkit (GLUT) Header */
//#include <GL/gl.h>           /* open gl lib */
//#include <GL/glu.h>          /* GL utilities */
//#include <GL/glext.h>
#include <math.h>            /* mathematical functions */

//for tiff output
#include <tiffio.h>

#include "SimpleX.h"

using namespace std;

#define MouseMirror +1     /* +1 or -1 */
#define MouseMovementScale 100
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ABS(x)   (x>0?x:-x)


void *viewer(void *simplex);


#endif
