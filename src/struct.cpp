/*
 * struct.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: enrico
 */

#include "point.h"
#include "lines.h"
#include "functions.h"
#include "struct.h"

void STRUCTMESH::writeMeshTec2D(const char *name)
{
	FILE *fp;
	if ((fp = fopen(name, "wt")) == NULL)
	{
		printf("couldn't open %s for writing\n", name);
		return;
	}
	fprintf(fp, "VARIABLES = \"x\" \"y\" \"z\" \n");
	fprintf(fp, "ZONE I = %d, J = %d, DATAPACKING=POINT\n", imax, jmax);
	for (int j=0; j<jmax; j++)
	for (int i=0; i<imax; i++)
		fprintf(fp, "%.8lf\t%.8lf\t%.8lf\n", mesh[i][j].x, mesh[i][j].y, mesh[i][j].z);
	fclose(fp);
}

STRUCTMESH STRUCTMESH::structBlockMesh2D_splines(SPLINE &bottom, SPLINE &top, const int imax, const int jmax, const double firstLength,
                                  const int interp, const double par1, const double par2)
{
  /*
   * This function takes two spline and generates a mesh. Currently use to generate the boundary layer mesh
   * Has a distribution along the splines (imax) and along the normal direction between the splines
   * The connecting line between the splines can be done into way:
   *              1. Using same length distance in both splines
   *              2. Using the normal
   * Currently it is assume the option 1, because the normal for very differnt spline will give not so good meshes!
   * To use option 2 just change true to false inside the for loop from 0 to imax.
   */

  //Output mesh
  STRUCTMESH block2D(imax,jmax);

  double thickness=1;
  double t_i[imax];       // Distribution along the spline
  double t_j[jmax];       // Distribution between the splines

  //Generating the point distribution along the spline
  for (int i=0; i<imax; i++)
    t_i[i] = (double)i/(double)(imax-1);

  if (interp==1) // tanh
    for (int i=0; i<imax; i++)
      t_i[i] = stretchingTanh(t_i[i], par1, par2);
  else if (interp==2) // atanh
    for (int i=0; i<imax; i++)
      t_i[i] = stretchingAtanh(t_i[i], par2, par1);
  else if (interp==3) // first-last length
    firstLastLengthDistr(par1/bottom.length, par2/bottom.length, imax, t_i);


  for (int i=0; i<imax; i++)
  {
    POINT p1, p2;
    //Points at the same lenght distance in the spline
    p1 = bottom.calcPoint(t_i[i]);
    p2 = top.controlPt[i];

    //Connecting the splines
    if(true) //Using same length distance in both splines
    {

      //Generating the point distribution between the spline
      firstLengthDistr(  firstLength/mag(p2-p1), jmax, t_j);
      for (int j=0; j<jmax; j++)
      {
        POINT pt = p1 + t_j[j]*(p2 - p1);
        block2D.mesh[i][j] = pt;
      }
    }
    else //Using the normal
    {
      //Finding the vector normal to the bottom spline
      POINT norm = bottom.calcNorm2D(t_i[i]);
      LINE line1(p1, p1 + thickness*(10)*norm);

      //Finding the intersection between normal and the top spline
      POINT p3 = top.intersectXY(line1, false);

      //Generating the point distribution between the spline
      firstLengthDistr(  firstLength/mag(p3-p1), jmax, t_j);

      for (int j=0; j<jmax; j++)
      {
        POINT pt = p1 + t_j[j]*(p3 - p1);
        block2D.mesh[i][j] = pt;
      }
    }

  }

  return (block2D);
}


void STRUCTMESH::drawMesh()
{
	for (int i=0; i<imax; i++)
	{
		glBegin(GL_LINE_STRIP);
		glColor3d(0.0, 0.0, 1.0);
		for (int j=0; j<jmax; j++)
			glVertex3d(mesh[i][j].x, mesh[i][j].y, mesh[i][j].z);
		glEnd();
	}

	for (int j=0; j<jmax; j++)
	{
		glBegin(GL_LINE_STRIP);
		glColor3d(0.0, 0.0, 1.0);
		for (int i=0; i<imax; i++)
			glVertex3d(mesh[i][j].x, mesh[i][j].y, mesh[i][j].z);
		glEnd();
	}
}

void STRUCTMESH::drawSolid()
{
	glBegin(GL_QUADS);
	glColor3d(0.4, 0.4, 0.4);
	for (int i=0; i<imax-1; i++)
	for (int j=0; j<jmax-1; j++)
	{
		glVertex3d(mesh[i][j].x, mesh[i][j].y, mesh[i][j].z);
		glVertex3d(mesh[i+1][j].x, mesh[i+1][j].y, mesh[i+1][j].z);
		glVertex3d(mesh[i+1][j+1].x, mesh[i+1][j+1].y, mesh[i+1][j+1].z);
		glVertex3d(mesh[i][j+1].x, mesh[i][j+1].y, mesh[i][j+1].z);
	}
	glEnd();
}
