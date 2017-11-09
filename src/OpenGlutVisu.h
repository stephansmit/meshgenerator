#ifndef OPENGLUTVISU_H
#define OPENGLUTVISU_H

//#define WITH_OPENGL

#ifdef WITH_OPENGL
#include <GL/glu.h>
#include <GL/glut.h>

#include "point.h"

#include "ArcBall.h"

bool flagMakeMesh;
bool bTransform;
bool bCurrTransform;
POINT cent,eyeVec,lookPos, trans;
double eyeDistance;
int mouseDownX, mouseDownY, modifierKey;
bool rightButton, leftButton, middleButton, ctrKey;
int width, height;



// User Defined Variables
//GLUquadricObj *quadratic;                     // Used For Our Quadric

const float PI2 = 2.0*3.1415926535f;                // PI Squared

Matrix4fT   Transform   = {  1.0f,  0.0f,  0.0f,  0.0f,       // NEW: Final Transform
                             0.0f,  1.0f,  0.0f,  0.0f,
                             0.0f,  0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  0.0f,  1.0f };

Matrix3fT   LastRot     = {  1.0f,  0.0f,  0.0f,          // NEW: Last Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

Matrix3fT   ThisRot     = {  1.0f,  0.0f,  0.0f,          // NEW: This Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

ArcBallT    ArcBall(800.0f, 500.0f);                        // NEW: ArcBall Instance
Point2fT    MousePt;                        // NEW: Current Mouse Point
bool isDragging;





enum key{REDRAW,TRANSFORM_TO_XYZ,TRANSFORM_TO_MPS};

#define clamp(x) x = x > 360.0 ? x-360.0 : x < -360.0 ? x+=360.0 : x


void reshape (int w, int h)
{
  width = w;
  height = h;
  glViewport (0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(eyeDistance, (GLfloat)w / (GLfloat)h, 0.1, 1000.0);
//  glOrtho(-width/2., width/2., -height/2., height/2., 0.1, 100.0);
  glMatrixMode (GL_MODELVIEW);

  ArcBall.setBounds((float)width, (float)height);
}

void keyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
  case 'r':
    flagMakeMesh = true;
    break;

  case 't':
    bTransform = bTransform ? false : true;
    break;

  case '1':
    eyeDistance = 50.0;
    trans.x = trans.y = trans.z = 0.0;
    eyeVec  = POINT(0.0,0.0,-eyeDistance);
    lookPos = POINT(0.0,0.0,0.0);

    LastRot.M[0] = 1.0;  LastRot.M[1] = 0.0;  LastRot.M[2] = 0.0;
    LastRot.M[3] = 0.0;  LastRot.M[4] = 1.0;  LastRot.M[5] = 0.0;
    LastRot.M[6] = 0.0;  LastRot.M[7] = 0.0;  LastRot.M[8] = 1.0;

    for (int i=0; i<9; i++) ThisRot.M[i] = LastRot.M[i];

    Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);

    glutPostRedisplay();
    break;

  case 'z':
    LastRot.M[0] = 1.0;  LastRot.M[1] = 0.0;  LastRot.M[2] = 0.0;
    LastRot.M[3] = 0.0;  LastRot.M[4] = 1.0;  LastRot.M[5] = 0.0;
    LastRot.M[6] = 0.0;  LastRot.M[7] = 0.0;  LastRot.M[8] = 1.0;

    for (int i=0; i<9; i++) ThisRot.M[i] = LastRot.M[i];

    Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);

    glutPostRedisplay();
    break;

  case 'y':
    LastRot.M[0] = 0.0;  LastRot.M[1] = 0.0;  LastRot.M[2] = -1.0;
    LastRot.M[3] = 0.0;  LastRot.M[4] = 1.0;  LastRot.M[5] =  0.0;
    LastRot.M[6] = 1.0;  LastRot.M[7] = 0.0;  LastRot.M[8] =  0.0;

    for (int i=0; i<9; i++) ThisRot.M[i] = LastRot.M[i];

    Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);

    glutPostRedisplay();
    break;
  }
}


void mouseButtonEvent(int button, int state, int x, int y)
{
  mouseDownX = x; mouseDownY = y;

  leftButton   = ((button == GLUT_LEFT_BUTTON)   && (state == GLUT_DOWN));
  rightButton  = ((button == GLUT_RIGHT_BUTTON)  && (state == GLUT_DOWN));
  middleButton = ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN));

  isDragging = true;                      // Prepare For Dragging
  LastRot = ThisRot;                      // Set Last Static Rotation To Last Dynamic One
  MousePt.s.X = x;   MousePt.s.Y = y;
  ArcBall.click(&MousePt);                // Update Start Vector And Prepare For Dragging

  modifierKey = glutGetModifiers();
}

void mouseButtonEventMove(int x, int y)
{
  if (rightButton && (modifierKey != GLUT_ACTIVE_CTRL))
  {
    trans.x += (double)(x-mouseDownX)*eyeDistance/50000.0;
    trans.y += (double)(mouseDownY-y)*eyeDistance/50000.0;
  } // translate

  if (leftButton && (modifierKey != GLUT_ACTIVE_CTRL))
  {
    Quat4fT     ThisQuat;
    MousePt.s.X = x;   MousePt.s.Y = y;
    ArcBall.drag(&MousePt, &ThisQuat);                        // Update End Vector And Get Rotation As Quaternion
    Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);       // Convert Quaternion Into Matrix3fT
    Matrix3fMulMatrix3f(&ThisRot, &LastRot);                  // Accumulate Last Rotation Into This One
    Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);    // Set Our Final Transform's Rotation From This One
  } // rotate

  if (middleButton && (modifierKey != GLUT_ACTIVE_CTRL))
  {
    eyeDistance += (double)(-(mouseDownX-x)+(mouseDownY-y))/1000.0*eyeDistance;
    eyeDistance = max(eyeDistance, 0.001);
  } // scale

  mouseDownX = x;   mouseDownY = y;

  glutPostRedisplay();
}

void initOpenGl(int argc, char *argv[])
{
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(800, 500);
  glutInitWindowPosition (100, 100);
  glutCreateWindow("mesh visu");


  glEnable( GL_POINT_SMOOTH );
//  glEnable (GL_LINE_SMOOTH);
//  glEnable (GL_BLEND);
//  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//  glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glPointSize(6.0);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glClearDepth(1.0f);


  rightButton = false;
  leftButton = false;
  middleButton = false;
  ctrKey = false;
  modifierKey = 0;

  eyeDistance = 50.0;

  trans.x = trans.y = trans.z = 0.0;
  eyeVec  = POINT(0.0,0.0,-eyeDistance);
  lookPos = POINT(0.0,0.0,0.0);

  flagMakeMesh = true;

  glutReshapeFunc (reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouseButtonEvent);
  glutMotionFunc(mouseButtonEventMove);
}




class VISUAL
{
public:
  // these guys are copies, pointer didn't work
  // because object loose their scope after make is finished
  // also, this is the reason why I need to have LINE, BEZIER and SPLINE, urggghg
  deque<POINT> pts;
  deque<LINE> lines;
  deque<SPLINE> splines;
  deque<BEZIER> bezierlines;
  deque<UNSTRUCTMESH> umeshes;
  deque<STRUCTMESH> smeshes;

  POINT urc, llc; // UpperRightcCorner and LowerLeftCorner

public:
  VISUAL() : urc(-1.0e20, -1.0e20, -1.0e20), llc(1.0e20, 1.0e20, 1.0e20) {}
  ~VISUAL() {}

  void printHelp()
  {
    printf("\n\nKEY SHORTCUTS:\n");
    printf("  'r' ... re-read input file and re-mesh\n");
    printf("  't' ... toggle transformations\n");
    printf("  '1' ... reset view\n");
    printf("  'x' ... view along x-axis\n");
    printf("  'y' ... view along y-axis\n");
    printf("  'z' ... view along z-axis\n");
  }

  void resetMinMax()
  {
    urc = POINT(-1.0e20, -1.0e20, -1.0e20);
    llc = POINT(1.0e20, 1.0e20, 1.0e20);
  }

  void clearVisu()
  {
    pts.clear();
    lines.clear();
    splines.clear();
    bezierlines.clear();
    umeshes.clear();
    smeshes.clear();
    resetMinMax();
  }

  void setMinMax(const POINT &pt)
  {
    urc.x = max(urc.x, pt.x);
    urc.y = max(urc.y, pt.y);
    urc.z = max(urc.z, pt.z);
    llc.x = min(llc.x, pt.x);
    llc.y = min(llc.y, pt.y);
    llc.z = min(llc.z, pt.z);
  }

  void recomputeMinMax()
  {
    resetMinMax();

    for (int i=0; i<pts.size(); i++)              setMinMax(pts[i]);

    for (int i=0; i<lines.size(); i++)
      for (int j=0; j<lines[i].size(); j++)       setMinMax(lines[i][j]);

    for (int i=0; i<splines.size(); i++)
      for (int j=0; j<splines[i].size(); j++)     setMinMax(splines[i][j]);
    for (int i=0; i<bezierlines.size(); i++)
      for (int j=0; j<bezierlines[i].size(); j++)     setMinMax(bezierlines[i][j]);


    for (int u=0; u<umeshes.size(); u++)
      for (int i=0; i<umeshes[u].nodes.size(); i++)  setMinMax(umeshes[u].nodes[i].pt);

    for (int s=0; s<smeshes.size(); s++)
      for (int i=0; i<smeshes[s].imax; i++)
        for (int j=0; j<smeshes[s].jmax; j++)        setMinMax(smeshes[s].mesh[i][j]);
  }

  void addToDisplay(POINT pt)
  {
    pts.push_back(pt);
    setMinMax(pt);
  }

  void addToDisplay(deque<POINT> &pts)
  {
    for (int i=0; i<pts.size(); i++)
      addToDisplay(pts[i]);
  }

  void addToDisplay(LINE &line)
  {
    lines.push_back(line);
    for (int i=0; i<line.size(); i++)
      setMinMax(line[i]);
  }

  void addToDisplay(SPLINE &line)
  {
    splines.push_back(line);
    for (int i=0; i<line.size(); i++)
      setMinMax(line[i]);
  }
  void addToDisplay(BEZIER &line)
  {
    bezierlines.push_back(line);
    for (int i=0; i<line.size(); i++)
      setMinMax(line[i]);
  }
  void addToDisplay(UNSTRUCTMESH &mesh)
  {
    umeshes.push_back(mesh);
    for (int i=0; i<mesh.nodes.size(); i++)
      setMinMax(mesh.nodes[i].pt);
  }

  void addToDisplay(STRUCTMESH &mesh)
  {
    smeshes.push_back(mesh);
    for (int i=0; i<mesh.imax; i++)
    for (int j=0; j<mesh.jmax; j++)
      setMinMax(mesh.mesh[i][j]);
  }

  void draw()
  {
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(eyeDistance, (GLfloat)width / (GLfloat)height, 0.1, 100.0);
//    glOrtho(-width/2., width/2., -height/2., height/2., 0.1, 100.0);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();


    cent = 0.5*(urc+llc);
    double maxLen = mag(urc-llc);
    cent = cent/maxLen;

    gluLookAt (0.0, 0.0, 0.0+1.0,   //eyeDistance,  // eye
               0.0, 0.0, 0.0,       // lookat
               0.0, 1.0, 0.0);      // up vector
    glTranslated(trans.x, trans.y, trans.z);
    glMultMatrixf(Transform.M);
    glTranslated(-cent.x, -cent.y, -cent.z);

    glScaled(1.0/maxLen, 1.0/maxLen, 1.0/maxLen);


    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for (int i=0; i<umeshes.size(); i++)  umeshes[i].drawSolid();
    for (int i=0; i<smeshes.size(); i++)  smeshes[i].drawSolid();
    glDisable(GL_POLYGON_OFFSET_FILL);

    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(-1.0f, -1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (int i=0; i<umeshes.size(); i++)  umeshes[i].drawMesh();
    for (int i=0; i<smeshes.size(); i++)  smeshes[i].drawMesh();
    for (int i=0; i<splines.size(); i++)  splines[i].drawLine();
    for (int i=0; i<bezierlines.size(); i++)  bezierlines[i].drawLine();
    for (int i=0; i<lines.size(); i++)    lines[i].drawLine();
    glDisable(GL_POLYGON_OFFSET_LINE);

    glEnable(GL_POLYGON_OFFSET_POINT);
    glPolygonOffset(-2.0f, -1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    for (int i=0; i<pts.size(); i++)      pts[i].drawPoint();
    glDisable(GL_POLYGON_OFFSET_LINE);

    glutSwapBuffers();
    glFlush ();
  }
};


#else

bool flagMakeMesh;


class VISUAL
{
public:
  void clearVisu(){};
  void draw(){};
  void printHelp(){};
};

#endif



#endif





