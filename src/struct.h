#ifndef STRUCTMESH_H_
#define STRUCTMESH_H_



#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "point.h"
#include "lines.h"
#include "functions.h"

#include <iostream>
#include <deque>
#include <algorithm>
using namespace std;

#ifdef WITH_OPENGL
  #include <GL/glu.h>
  #include <GL/glut.h>
#else
  #include "gldummy.h"
#endif


class STRUCTMESH
{
public:
  deque< deque<POINT> > mesh;
  int imax, jmax;

public:
  STRUCTMESH() : imax(0), jmax(0), mesh() {}//mesh(NULL) {}

  void init(int imax_, int jmax_)
  {
    imax = imax_;
    jmax = jmax_;

    mesh.resize(imax);
    for (int i=0; i<imax; i++)
      mesh[i].resize(jmax);
  }

  STRUCTMESH(int imax_, int jmax_)   { init(imax_, jmax_); }

  STRUCTMESH(const char *name)
  {
    FILE *fp;
    if ((fp = fopen(name, "rt")) == NULL)
    {
      printf("couldn't open %s for writing\n", name);
      return;
    }

    int nBl, imax, jmax, kmax;
    fscanf(fp, "%i\n%i%i%i\n",&nBl, &imax, &jmax, &kmax);

    printf("%d\t%d\t%d\t%d\n", nBl, imax, jmax, kmax);

    init(imax+1, jmax+1);

    double dummy;

    for (int i=0; i<=imax; i++)
    for (int j=0; j<=jmax; j++)
    for (int k=0; k<=kmax; k++)
    {
      fscanf(fp, "%lf\t%lf\t%lf\n", &mesh[i][j].x, &mesh[i][j].y, &dummy);
//      printf("%d\t%d\t%d\t%lf\t%lf\t%lf\n", i,j,k,mesh[i][j].x, mesh[i][j].y, dummy);
    }
    fclose(fp);

//    for (int i=0; i<imax; i++)
//    for (int j=0; j<jmax; j++)
//      printf("%lf\t%lf\n", mesh[i][j].x, mesh[i][j].y);


  //
  //  fprintf(fp, "ZONE I = %d, J = %d, DATAPACKING=POINT\n", imax, jmax);
  //  for (int j=0; j<jmax; j++)
  //  for (int i=0; i<imax; i++)
  //    fprintf(fp, "%.8lf\t%.8lf\t%.8lf\n", mesh[i][j].x, mesh[i][j].y, mesh[i][j].z);
  //  fclose(fp);
  }

  void writeMeshTec2D(const char *name);

  STRUCTMESH structBlockMesh2D_splines(SPLINE &bottom, SPLINE &top, const int imax, const int jmax, const double firstLength,
                                    const int interp=0, const double par1=1.0, const double par2=1.0);

  void drawMesh();


  void drawSolid();



};

#endif
