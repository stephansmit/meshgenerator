#ifndef GLDUMMY_H
#define GLDUMMY_H

#define initOpenGl(argc, argv)  {flagMakeMesh = true;}
#define glutDisplayFunc(A)(A())
#define glutPostRedisplay()
#define glutIdleFunc(A)(A())
#define glutMainLoop()
#define glBegin(a)
#define glColor3d(a, b, c)
#define glVertex3d(x, y, z)
#define glEnd()
#define glLineWidth(A)

#endif


