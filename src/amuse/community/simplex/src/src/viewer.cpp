#include "viewer.h"

//to write a png: look at libpng.org
//and also at http://www.opengl.org/discussion_boards/ubbthreads.php?ubb=showflat&Number=267760
//for tiff's, see:
// http://www.opengl.org/resources/code/samples/mjktips/libtiff/index.html
using namespace std;

//c function
extern "C" int writetiff(char *filename, char *description, int x, int y, int width, int height, int compression);

#ifndef __VIEWER
#define SIMPLEX SimpleX
SIMPLEX *SimpleXGrid;
#else
#define SIMPLEX AMUSE_SimpleX
extern SIMPLEX *SimpleXGrid;
#endif

float fi=.9*M_PI/4.,theta=M_PI/7,r=2.25;
float cam[3],dummy[3];
float spot[3];
float upvector[3]={0.,0.,1.};
float mousexyzero[2];
int npart,leftbutton=1,rightbutton=1,free_flight;
float freeFlightSpeed=0.05;
bool slew=false;
float scaleFactor=1.;
float curveScale=0.02;
float minSpotDistance=0.001;
float viewangle=45.;
float ratio,angle;
float mx,my;
float up[3],ri[3];
float da[]={ 0.0,0.03,0.1};
GLuint edges;
GLuint simplices;
int framerate=1000;
int updateflag=0;

void set_lists(int dummy){

  int plotValue = 0;

  glNewList(edges,GL_COMPILE);
  //plot Delaunay edges
  if(plotValue == 0){
    glBegin(GL_LINES);
    if ( (*SimpleXGrid).sites.end()-(*SimpleXGrid).sites.begin() > 0){
      for(SITE_ITERATOR it=(*SimpleXGrid).sites.begin();it<(*SimpleXGrid).sites.end();it++){
	float x1,x2,y1,y2,z1,z2,d1,d2;
	for( unsigned int j=0; j<it->get_numNeigh(); j++ ){
	  if( !it->get_border() && !(*SimpleXGrid).sites[ it->get_neighId(j) ].get_border() ){
	    if( it->get_neighId(j) < it->get_site_id()){
	      x1 = it->get_x();
	      x2 = (*SimpleXGrid).sites[ it->get_neighId(j) ].get_x();
	      y1 = it->get_y();
	      y2 = (*SimpleXGrid).sites[ it->get_neighId(j) ].get_y();
	      z1 = it->get_z();
	      z2 = (*SimpleXGrid).sites[ it->get_neighId(j) ].get_z();
	      d1 = it->get_ionised_fraction();
	      d2 = (*SimpleXGrid).sites[ it->get_neighId(j) ].get_ionised_fraction();
	      //        d1 = it->get_number_density();
	      //        d2 = (*SimpleXGrid).sites[ it->get_neighId(j) ].get_number_density();
	      //        printf("%g %g %g %g %g %g\n",x1,y1,z1,x2,y2,z2);
	      //        d1=(log10(d1)+8.)/5.;
	      //        d2=(log10(d2)+8.)/5.;
	      if(d1>0.05 || d2>0.05){
		glColor4f(d1,0.2,.5-d1/2,d1);
		glVertex3d(x1,y1,z1);
		glColor4f(d2,0.2,.5-d2/2,d2);
		glVertex3d(x2,y2,z2);
	      }
	    }
	  }
	}
      }
    }
  } else if( plotValue == 1 ){
    // plot simplices

    glBegin(GL_QUADS);
    if ( (*SimpleXGrid).simplices.size() > 0){
      for(SIMPL_ITERATOR it=(*SimpleXGrid).simplices.begin();it<(*SimpleXGrid).simplices.end();it++){
	if( !(*SimpleXGrid).sites[it->get_id1()].get_border() || !(*SimpleXGrid).sites[it->get_id2()].get_border() || 
	    !(*SimpleXGrid).sites[it->get_id3()].get_border() || !(*SimpleXGrid).sites[it->get_id4()].get_border() ){

	  float x1 = (*SimpleXGrid).sites[it->get_id1()].get_x();
	  float y1 = (*SimpleXGrid).sites[it->get_id1()].get_y();
	  float z1 = (*SimpleXGrid).sites[it->get_id1()].get_z();
	  float d1 = (*SimpleXGrid).sites[it->get_id1()].get_ionised_fraction();

	  float x2 = (*SimpleXGrid).sites[it->get_id2()].get_x();
	  float y2 = (*SimpleXGrid).sites[it->get_id2()].get_y();
	  float z2 = (*SimpleXGrid).sites[it->get_id2()].get_z();
	  float d2 = (*SimpleXGrid).sites[it->get_id2()].get_ionised_fraction();

	  float x3 = (*SimpleXGrid).sites[it->get_id3()].get_x();
	  float y3 = (*SimpleXGrid).sites[it->get_id3()].get_y();
	  float z3 = (*SimpleXGrid).sites[it->get_id3()].get_z();
	  float d3 = (*SimpleXGrid).sites[it->get_id3()].get_ionised_fraction();

	  float x4 = (*SimpleXGrid).sites[it->get_id4()].get_x();
	  float y4 = (*SimpleXGrid).sites[it->get_id4()].get_y();
	  float z4 = (*SimpleXGrid).sites[it->get_id4()].get_z();
	  float d4 = (*SimpleXGrid).sites[it->get_id4()].get_ionised_fraction();

	  //        d1 = it->get_number_density();
	  //        d2 = (*SimpleXGrid).sites[ it->get_neighId(j) ].get_number_density();
	  //        printf("%g %g %g %g %g %g\n",x1,y1,z1,x2,y2,z2);
	  //        d1=(log10(d1)+8.)/5.;
	  //        d2=(log10(d2)+8.)/5.;
	  if(d1>0.05 || d2>0.05 || d3>0.05 || d4>0.05){

	    float factor = 0.1;

	    //glColor4f(d1,0.2,.5-d1/2,d1);
	    glColor4f(d1,0.2,.5-d1/2,factor*d1);
	    glVertex3d(x1,y1,z1);

	    //glColor4f(d2,0.2,.5-d2/2,d2);
	    glColor4f(d2,0.2,.5-d2/2,factor*d2);
	    glVertex3d(x2,y2,z2);

	    //glColor4f(d3,0.2,.5-d3/2,d3);
	    glColor4f(d3,0.2,.5-d3/2,factor*d3);
	    glVertex3d(x3,y3,z3);

	    //glColor4f(d4,0.2,.5-d4/2,d4);
	    glColor4f(d4,0.2,.5-d4/2,factor*d4);
	    glVertex3d(x4,y4,z4);
	  }
	}
      }
    }

  }else{
    cerr << " Error in viewer, don't know what to plot!" << endl;
  }


  glEnd();
  glEndList();

  glutPostRedisplay();

  if(updateflag!=0) glutTimerFunc(framerate, set_lists,0);

}



void rotate(float *moving,float fixed[3],float r,float fi, float theta)
{
 *moving= fixed[0]+r*sin(fi)*cos(theta);
 *(moving+1)= fixed[1]+r*cos(fi)*cos(theta);
 *(moving+2)= fixed[2]+r*sin(theta);
}

void translate(float *one,float *two, float trans[3])
{
 *one+=trans[0];      *two+=trans[0];
 *(one+1)+=trans[1];  *(two+1)+=trans[1];
 *(one+2)+=trans[2];  *(two+2)+=trans[2];
}

float distance2(float o[3],float t[3])
{
 float d2;
 d2=((o[1]-t[1])*(o[1]-t[1])+(o[0]-t[0])*(o[0]-t[0])+(o[2]-t[2])*(o[2]-t[2]));
 return d2;
}

void display(){
  glClearColor(0.02,0.02,0.02, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();

  glTranslatef(-.5,-.5,-.5);
//  glColor4f(0.,0.,1.,1.);
  glLineWidth(1.25);
  glDepthFunc(GL_ALWAYS);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE);
  glCallList(edges);
  glPopMatrix();
  //almost white, thin box
  //glColor4f(238./255.,233./255.,233./255.,1.);
  glColor4f(220./255.,220./255.,220./255.,1.);
  glLineWidth(1.0);
  //red, thick box
  //glColor4f(1.,0.,0.,1.);
  //glLineWidth(2.25);
  glDepthFunc(GL_LESS);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  //comment out if no bounding box is desired
  glutWireCube(1.);


  glutSwapBuffers();
}

void pointCamera(){

 glMatrixMode(GL_MODELVIEW);
 glLoadIdentity();

 gluLookAt( cam[0], cam[1], cam[2],
            spot[0], spot[1], spot[2],
            upvector[0], upvector[1], upvector[2]);
 }


 void ReshapeWindow(int width,int height)
{

 if(height<1) height=1;
 ratio=width/((float) height);

 glViewport(0, 0, width, height);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 gluPerspective(viewangle,ratio,.001,100.); /* nearplane, farplane */
 pointCamera();

}

void reset_to_init()
{
 slew=false;
 //fi=.9*M_PI/4.,theta=M_PI/7,r=2.25*scaleFactor;
 fi=.9*M_PI,theta=0.1*M_PI,r=2.25*scaleFactor;
 spot[0]=0.;spot[1]=0.,spot[2]=0.;
 rotate(cam,spot,r,fi,theta);
 pointCamera();
 glutPostRedisplay();
}

bool myinit(){

  glEnable(GL_DEPTH_TEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE);
  glLineWidth(1.25);

  scaleFactor=1;

  reset_to_init();
  edges=glGenLists(1);
  //simplices=glGenLists(1);
  set_lists(0);

  return true;

}



void MouseClick(int button,int state, int x, int y)
{
  if(button==GLUT_LEFT_BUTTON)
    {
      if(state==GLUT_DOWN) {
	// printf("left mouse button down at %d %d\n",x,y);
	leftbutton=0;
      }
      if(state==GLUT_UP) {
	// printf("left mouse button up at %d %d\n",x,y);
	leftbutton=1;
      }
      mousexyzero[0]=x; mousexyzero[1]=y;
      mx=0.;my=0.;
    }

  if(button==GLUT_RIGHT_BUTTON)
    {
      if(state==GLUT_DOWN) {
	rightbutton=0;
      }
      if(state==GLUT_UP) {
	rightbutton=1;
      }
      mousexyzero[0]=x; mousexyzero[1]=y;
    }


}

void MouseClickMotion(int x, int y)
{

float rot,kantel,abs1,abs2;
// printf("mouse position: %d %d \n",x,y);

 if(!(leftbutton==0) && !(rightbutton==0)) return;

 if(!(leftbutton==0) && rightbutton==0){
   if(!slew) {
     r*=pow(2.,(y-mousexyzero[1])/MouseMovementScale);
     if(r<minSpotDistance) r=minSpotDistance;
     if(r>1000*scaleFactor) r=1000*scaleFactor;
     rotate(cam,spot,r,fi,theta);
   }
   if(slew){
     float zero[3]={0.,0.,0.},trans[3]={0.,0.,0.};
     rotate(trans,zero,scaleFactor*(y-mousexyzero[1])/MouseMovementScale,fi,theta);
     translate(cam,spot,trans);}
 }
 
 if(leftbutton==0 && !(rightbutton==0)) {
  if(free_flight==1){
  mx+=(x-mousexyzero[0])/MouseMovementScale;
  my+=(y-mousexyzero[1])/MouseMovementScale;
  if(ABS(mx)<0.02) mx=0.;
  if(ABS(my)<0.02) my=0.;
  }
  if(free_flight==0){
   fi+=(x-mousexyzero[0])/MouseMovementScale;
   theta+=(y-mousexyzero[1])/MouseMovementScale;
   if(theta>M_PI/2.) theta=M_PI/2-0.001;
   if(theta<-M_PI/2) theta=-M_PI/2+0.001;
   if(slew){
      rotate(spot,cam,r,fi+M_PI,-theta);
   } else {
      rotate(cam, spot,r,fi,theta);
   }
  }
  }

 if(rightbutton==0 && leftbutton==0){
  viewangle*=1+(y-mousexyzero[1])/MouseMovementScale;
  if(viewangle<1.) viewangle=1;
  if(viewangle>120.) viewangle=120;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //glViewport(0, 0, width, height);
  gluPerspective(viewangle,ratio,.001,100.);
 }
 mousexyzero[0]=x; mousexyzero[1]=y;
 pointCamera();
 glutPostRedisplay();
}

void NormalKey( unsigned char key,int x , int y)
{  mx=0.;my=0;
if(key==' ') {free_flight=1;slew=true;} 
if(key=='q') glutLeaveMainLoop();
if(key=='i') {
 printf("---- view parameters ----\n");
 printf("  spot: %g %g %g \n",spot[0],spot[1],spot[2]);
 printf("  cam: %g %g %g \n",cam[0],cam[1],cam[2]);
 printf("  dir: %g %g %g \n",spot[0]-cam[0],spot[1]-cam[0],spot[2]-cam[0]);
}
}

void NormalKeyUp(unsigned char key,int x , int y)
{  mx=0.;my=0;
if(key='f')free_flight=0; }

void MouseMoves(int x, int y)
{
 if(free_flight==1){
 mx=(x-mousexyzero[0])/MouseMovementScale/1.;
if(fabs(mx)<0.02) mx=0.;
 my=(y-mousexyzero[1])/MouseMovementScale/1.;}
if(fabs(my)<0.02) my=0.;
}


void frame(int i)
{
  float zero[3]={0.,0.,0.};
  float trans[3]={0.,0.,0.};
  if(free_flight==1){
    rotate(trans,zero,-freeFlightSpeed*scaleFactor,fi,theta);
    translate(cam,spot,trans);
    if(leftbutton==0) {
      fi+=mx*curveScale;
      theta+=my*curveScale;
      if(theta>M_PI/2) theta=M_PI/2-0.001;
      if(theta<-M_PI/2) theta=-M_PI/2+0.001;
      rotate(spot,cam,r,fi+M_PI,-theta);
    }
    pointCamera();
    glutPostRedisplay();
  }
  glutTimerFunc(50, frame,0);
}

void show_menu(int optie)
{
 switch(optie)
 {
 case 0:
  if(updateflag==0) set_lists(0);
  break;
 case 1:
  updateflag=0;
  break;
 case 2:
  framerate=2000;
  if(updateflag==0){
    updateflag=1;
    set_lists(0);
    }
  break;
 case 3:
  framerate=200;
  if(updateflag==0){
    updateflag=1;
    set_lists(0);
    }
  break;
 }
}

void prop_menu(int optie)
{
 switch(optie)
 {
 case 1:
  printf("prop optie 1\n");
  break;
 }
}

void write_img_menu(int optie){

  int win_width;
  int win_height;

  switch(optie){
  case 0:
    cerr << "Writing SimpleX_output.tif to disk" << endl;
    win_width  = glutGet(GLUT_WINDOW_WIDTH);
    win_height = glutGet(GLUT_WINDOW_HEIGHT);
    writetiff("SimpleX_output.tif", "OpenGL-rendered SimpleX output", 0, 0, win_width, win_height, COMPRESSION_NONE );
    break;
  case 1:
    cerr << "No implementation to write png's yet. Why don't YOU write it?" << endl;
    break;
  }

}

int create_show_menu(){
  int m;
  m=glutCreateMenu(show_menu);
  glutAddMenuEntry("update now", 0) ;
  glutAddMenuEntry("no update", 1) ;
  glutAddMenuEntry("slow update", 2) ;
  glutAddMenuEntry("fast update", 3) ;
  return m;
}

int create_prop_menu(){
  int m;
  m=glutCreateMenu(prop_menu);
  glutAddMenuEntry("color dummy", 1) ;
  return m;
}

int create_write_img_menu(){
  int m;
  m=glutCreateMenu(write_img_menu);
  glutAddMenuEntry("write tiff", 0) ;
  glutAddMenuEntry("write png", 1) ;
  return m;
}

void menu_handler(int value)
{
  switch(value)
    {
    case 1:
      reset_to_init();
      break;
    case 2:
      glutLeaveMainLoop();  
      break;
    }
}

//function written by Mark J. Kilgard, taken from
// http://www.opengl.org/resources/code/samples/mjktips/libtiff/writetiff.c
int writetiff(char *filename, char *description,
	      int x, int y, int width, int height, int compression){

  TIFF *file;
  GLubyte *image, *p;
  int i;

  file = TIFFOpen(filename, "w");
  if (file == NULL) {
    return 1;
  }
  image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

  /* OpenGL's default 4 byte pack alignment would leave extra bytes at the
     end of each image row so that each full row contained a number of bytes
     divisible by 4.  Ie, an RGB row with 3 pixels and 8-bit componets would
     be laid out like "RGBRGBRGBxxx" where the last three "xxx" bytes exist
     just to pad the row out to 12 bytes (12 is divisible by 4). To make sure
     the rows are packed as tight as possible (no row padding), set the pack
     alignment to 1. */
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
  p = image;
  for (i = height - 1; i >= 0; i--) {
    if (TIFFWriteScanline(file, p, i, 0) < 0) {
      free(image);
      TIFFClose(file);
      return 1;
    }
    p += width * sizeof(GLubyte) * 3;
  }
  TIFFClose(file);
  return 0;

}


void *viewer(void* simplex)
{
    int __argc;
    char ** __argv;
    int sm1,sm2,sm3,mm;

    SimpleXGrid=(SIMPLEX*) simplex;

      __argc = 1;
      __argv = (char **) malloc((__argc+1)*sizeof(char));
      __argv[0] = "noname";
      __argv[1] = NULL;


    glutInit( &__argc, __argv);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(700,400);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
    glutCreateWindow(" SimpleX2 ");
    

    sm1=create_show_menu();
    sm2=create_prop_menu();
    sm3=create_write_img_menu();
    mm=glutCreateMenu(menu_handler);

    glutAddMenuEntry("Reset View",1);
    glutAddMenuEntry("Quit",2);
    glutAttachMenu(GLUT_MIDDLE_BUTTON);
    

    glutAddSubMenu("Show",sm1);
    glutAddSubMenu("Properties",sm2);
    glutAddSubMenu("Write Image",sm3);


    glutReshapeFunc( ReshapeWindow );
    glutDisplayFunc( display );
    glutMouseFunc( MouseClick);
    glutMotionFunc( MouseClickMotion );

    glutIgnoreKeyRepeat(1);
    glutKeyboardFunc( NormalKey );
    glutKeyboardUpFunc( NormalKeyUp);

    myinit();

    glutTimerFunc(50, frame, 0);
   
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);

    glutMainLoop();

    return 0;

}
