/*
 * transform.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: enrico
 */
#include "lines.h"
#include "struct.h"
#include "unstruct.h"
#include "transform.h"

void TRANSFORM::initTransform(const LINE &hub, const LINE &shr)
{
	hubXYZ = hub;
	transformHubShr(hubXYZ, hubMPS, 0.0);

	shrXYZ = shr;
	transformHubShr(shrXYZ, shrMPS, 1.0);
}

void TRANSFORM::transformRZPHI_to_XYZ(deque<POINT> &pts)
{
	for (int i=0; i<pts.size(); i++)
	{
		double rad = pts[i].x;
		double z   = pts[i].y;
		double phi = pts[i].z;
		pts[i].x = rad*sin(phi);
		pts[i].y = rad*cos(phi);
		pts[i].z = z;
	}
}

void TRANSFORM::transformHubShr(const LINE &lineXYZ, LINE &lineMPS, double span)
{
	int NP = lineXYZ.controlPt.size();
	lineMPS.controlPt.resize(NP);

	// calculate new x coordinate for hub and shroud
	lineMPS[0].x = 0.0;
	for (int i=1; i<NP; i++)
	{
		double dr    = lineXYZ[i].radZ() - lineXYZ[i-1].radZ();
		double dz    = lineXYZ[i].z - lineXYZ[i-1].z;
		double drphi = (lineXYZ[i].phiRadZ() - lineXYZ[i-1].phiRadZ())/300.0;
		lineMPS[i].x = lineMPS[i-1].x + sqrt(dr*dr + dz*dz + drphi*drphi);
	}

	double offset = lineMPS[1].x;
	double bladeLength = lineMPS[NP-2].x - lineMPS[1].x;
	for (int i=0; i<NP; i++)
	{
		lineMPS[i].x -= offset;      // set offset to leading edge
		lineMPS[i].x /= bladeLength; // normalize
	}

	// calculate new y coordinate for hub and shroud: angle divided by 2
	for (int i=0; i<NP; i++)
		lineMPS[i].y = lineXYZ[i].phiRadZ()/2.0;

	// calculate new z coordinate for hub and shroud: span-wise location
	for (int i=0; i<NP; i++)
		lineMPS[i].z = span;
}

double TRANSFORM::calcSpan(POINT &pts)
{
  // calculate span in RZ plane
  deque<POINT> hubRZ, shrRZ;

  for (int i=0; i<hubXYZ.size(); i++)
  {
    hubRZ.push_back(POINT(hubXYZ[i].radZ(), hubXYZ[i].z, 0.0));
    shrRZ.push_back(POINT(shrXYZ[i].radZ(), shrXYZ[i].z, 0.0));
  }

  POINT ptsRZ = POINT(pts.radZ(), pts.z, 0.0);

  int j=1;
  while ((ptsRZ.y>hubRZ[j].y+(shrRZ[j].y-hubRZ[j].y)/(shrRZ[j].x-hubRZ[j].x)*(ptsRZ.x-hubRZ[j].x)) && (j<hubRZ.size())) j++;

  SPLINE hubRZspl(hubRZ, j-1, j), shrRZspl(shrRZ, j-1, j);
  double dist = 1e20;
  int np = 200;
  int ind = -1;
  deque<POINT> proPts, hubPts, shrPts;
  for (int i=0; i<np; i++)
  {
    double t = (double)i/(double)(np-1);
    hubPts.push_back(hubRZspl.calcPoint(t));
    shrPts.push_back(shrRZspl.calcPoint(t));

    deque<POINT> tmpPts;
    tmpPts.push_back(hubPts[i]);
    tmpPts.push_back(shrPts[i]);
    LINE tmpLine(tmpPts);
    proPts.push_back(tmpLine.project2D(ptsRZ));

    if (mag(proPts[i]-ptsRZ)<dist)
      {dist = mag(proPts[i]-ptsRZ); ind = i;}
  }
  return (mag(proPts[ind]-hubPts[ind])/mag(shrPts[ind]-hubPts[ind]));
}

void TRANSFORM::transform(POINT &pts, double span)
{

  POINT p(pts.z, pts.radZ(), 0.0);

  pts.y = pts.phiRadZ()/2.0; // ATTENTION  ATTENTION  ATTENTION  ATTENTION
  pts.z = span;

  int j = 0;
  double dist = 10e10;
  do
  {
    POINT phub1(hubXYZ[j].z,   hubXYZ[j].radZ(),   0.0);
    POINT phub2(hubXYZ[j+1].z, hubXYZ[j+1].radZ(), 0.0);

    POINT pshr1(shrXYZ[j].z,   shrXYZ[j].radZ(),   0.0);
    POINT pshr2(shrXYZ[j+1].z, shrXYZ[j+1].radZ(), 0.0);

    POINT pSpan1 = phub1 + span*(pshr1-phub1);
    POINT pSpan2 = phub2 + span*(pshr2-phub2);

    POINT svec(p     -pSpan1);
    POINT tang(pSpan2-pSpan1);

    dist = svec.dot(tang)/mag(tang)/mag(tang);

    double m1 = hubMPS[j  ].x + span*(shrMPS[j  ].x - hubMPS[j  ].x);
    double m2 = hubMPS[j+1].x + span*(shrMPS[j+1].x - hubMPS[j+1].x);
    pts.x = m1 + dist*(m2-m1);

    j++;

  }while((dist > 1.0) && (j < hubXYZ.controlPt.size()-1));
}

void TRANSFORM::transform(deque<POINT> &pts, double span)
{
  for (int i=0; i<pts.size(); i++)
    transform(pts[i], span);
}

void TRANSFORM::transform(STRUCTMESH &mesh, double span)
{
  for (int i=0; i<mesh.imax; i++)
    for (int j=0; j<mesh.jmax; j++)
      transform(mesh.mesh[i][j], span);
}

void TRANSFORM::transform(UNSTRUCTMESH &unstr, double span)
{
  for (int i=0; i<unstr.nodes.size(); i++)
    transform(unstr.nodes[i].pt, span);
}

void TRANSFORM::transform(POINT &pts)
{
	double span = calcSpan(pts);
	transform(pts, span);
}

void TRANSFORM::transform(deque<POINT> &pts)
{
	for (int i=0; i<pts.size(); i++)
		transform(pts[i]);
}

void TRANSFORM::transform(STRUCTMESH &mesh)
{
  for (int i=0; i<mesh.imax; i++)
    for (int j=0; j<mesh.jmax; j++)
      transform(mesh.mesh[i][j]);
}

void TRANSFORM::transform(UNSTRUCTMESH &unstr)
{
  for (int i=0; i<unstr.nodes.size(); i++)
    transform(unstr.nodes[i].pt);
}



void TRANSFORM::transformBack(POINT &pts)
{
	double m   = pts.x;
	double phi = pts.y*2.0; // ATTENTION  ATTENTION  ATTENTION  ATTENTION
	double z   = pts.z;

	int j=-1;
	double mspan = 0.0;
	do
	{
		j++;
		double mh = hubMPS[j].x;
		double ms = shrMPS[j].x;
		mspan = mh + z*(ms-mh);

	}while((mspan<=m) && (j<hubMPS.controlPt.size()-1));

	if (j==0) j++;

	double mspan1 = hubMPS[j-1].x + z*(shrMPS[j-1].x-hubMPS[j-1].x);
	double mspan2 = hubMPS[j  ].x + z*(shrMPS[j  ].x-hubMPS[j  ].x);

	double fact = (m-mspan1)/(mspan2-mspan1);

	double z1 = hubXYZ[j-1].z + z*(shrXYZ[j-1].z-hubXYZ[j-1].z);
	double z2 = hubXYZ[j  ].z + z*(shrXYZ[j  ].z-hubXYZ[j  ].z);

	double r1 = hubXYZ[j-1].radZ() + z*(shrXYZ[j-1].radZ()-hubXYZ[j-1].radZ());
	double r2 = hubXYZ[j  ].radZ() + z*(shrXYZ[j  ].radZ()-hubXYZ[j  ].radZ());

	double znew = z1 + fact*(z2-z1);
	double rnew = r1 + fact*(r2-r1);

	double x = rnew*sin(phi);
	double y = rnew*cos(phi);
//	z = znew;

	pts.x = x;
	pts.y = y;
	pts.z = znew;
	pts.dummy = z;
}

void TRANSFORM::transformBack(deque<POINT> &pts)
{
	for (int i=0; i<pts.size(); i++)
		transformBack(pts[i]);
}

void TRANSFORM::transformBack(deque<NODE> &nodes)
{
	for (int i=0; i<nodes.size(); i++)
		transformBack(nodes[i].pt);
}

void TRANSFORM::transformBack(STRUCTMESH &mesh)
{
	for (int i=0; i<mesh.imax; i++)
		for (int j=0; j<mesh.jmax; j++)
			transformBack(mesh.mesh[i][j]);
}

void TRANSFORM::transformBack(UNSTRUCTMESH &unstr)
{
	for (int i=0; i<unstr.nodes.size(); i++)
		transformBack(unstr.nodes[i].pt);
}
