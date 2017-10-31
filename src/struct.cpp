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
