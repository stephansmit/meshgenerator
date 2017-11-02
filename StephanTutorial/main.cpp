
/*
 * main.cpp
 *
 *  Created on: Jan 27, 2012
 *      Author: renep
 */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include <iostream>
#include <deque>
#include <algorithm>
using namespace std;

#include "Param.h"

#include "point.h"
#include "lines.h"
#include "transform.h"
#include "struct.h"
#include "unstruct.h"

#include "OpenGlutVisu.h"

#include "meshTools.h"


#define ROTOR_MESH 1;
#define ROTOR_BL   0;
#define ROTOR_UNST 0;
#define ROTOR_3D   1;


class RADIALCOMP : public MESHTOOLS, public VISUAL
{
public:
  RADIALCOMP(const char *name) : MESHTOOLS(name) {}

public:

  void coordTransformForVisualization()
  {
    if (bCurrTransform == bTransform) return;

    bCurrTransform = bTransform;

    if (bTransform == true)
    {
      printf("perform transform to MPS for VISU \n");

      for (int i=0; i<pts.size(); i++)       transform(pts[i], pts[i].dummy);
      for (int i=0; i<lines.size(); i++)     transform(lines[i].controlPt, lines[i][0].dummy);
      for (int i=0; i<splines.size(); i++)   transform(splines[i].controlPt, splines[i][0].dummy);
      for (int u=0; u<umeshes.size(); u++)
        for (int n=0; n<umeshes[u].nodes.size(); n++)  transform(umeshes[u].nodes[n].pt, umeshes[u].nodes[n].pt.dummy);

      for (int s=0; s<smeshes.size(); s++)
        for (int i=0; i<smeshes[s].imax; i++)
          for (int j=0; j<smeshes[s].jmax; j++)        transform(smeshes[s].mesh[i][j], smeshes[s].mesh[i][j].dummy);

      recomputeMinMax();
    }
    else
    {
      printf("perform transform to XYZ for VISU \n");

      transformBack(pts);
      for (int i=0; i<lines.size(); i++)     transformBack(lines[i].controlPt);
      for (int i=0; i<splines.size(); i++)   transformBack(splines[i].controlPt);
      for (int u=0; u<umeshes.size(); u++)   transformBack(umeshes[u]);
      for (int s=0; s<smeshes.size(); s++)   transformBack(smeshes[s]);

      recomputeMinMax();
    }
    printHelp();
  }


  SPLINE makeRadius(POINT startpoint, double radius, double angle)
  {
	  int n = 1000;
	  deque<POINT> points;
	  for (int i=0; i<=n; i++)
	  {
		POINT newpoint(rotateDegZ(angle-angle/(double)n*i, POINT(0.0, 0.0), startpoint));
		points.push_back(newpoint);
	  }

	  SPLINE transline(points);
	  return(transline);
  }
  deque<POINT> addDeques(deque< deque<POINT> > deques)
  {
	  deque<POINT> points;
	  for (int i=0; i<deques.size(); i++)
	  {
		  for (int j=0; j<deques[i].size(); j++)
		  	  {
		  		  points.push_front(deques[i][j]);
		  	  }
	  }
	  return points;
  }


  deque<POINT> discretizeLine(SPLINE line, int npoints, bool startend, bool reverse)
  {
	  deque<POINT> points;
	  int k, l;
	  double s;
	  if (startend) {
		  k=0;
		  l=1;
	  } else {
		  k=1;
		  l=0;
	  }
	  for(int i=k; i< (npoints+l); i++)
	  {
		  s = (double)i/(double)npoints;
		  if(reverse){
			  points.push_back(line.calcPoint(s));
		  }
		  else{
			  points.push_front(line.calcPoint(s));
		  }
	  }
	  return(points);
  }
  deque<POINT> discretizeLineS(SPLINE line, int npoints, bool startend, bool reverse, double sstart, double send)
    {
  	  deque<POINT> points;
  	  int k, l;
  	  double s;
  	  if (startend) {
  		  k=0;
  		  l=1;
  	  } else {
  		  k=1;
  		  l=0;
  	  }
  	  for(int i=k; i< (npoints+l); i++)
  	  {
  		  s = (send-sstart)*(double)i/(double)npoints+sstart;
  		  if(reverse){
  			  points.push_back(line.calcPoint(s));
  		  }
  		  else{
  			  points.push_front(line.calcPoint(s));
  		  }
  	  }
  	  return(points);
    }
  deque<POINT> discretizeLine2(BEZIER line, int npoints, bool startend, bool reverse)
  {
	  deque<POINT> points;
	  int k, l;
	  double s;
	  if (startend) {
		  k=0;
		  l=1;
	  } else {
		  k=1;
		  l=0;
	  }
	  for(int i=k; i< (npoints+l); i++)
	  {
		  s = (double)i/(double)npoints;
		  if(reverse){
			  points.push_back(line.calcPoint(s));
		  }
		  else{
			  points.push_front(line.calcPoint(s));
		  }
	  }
	  return(points);
  }
  deque<POINT> translateLine(SPLINE line, POINT point)
  {
	  deque<POINT> points;
	  POINT newpoint;
	  for (int i=line.iStart;i<=line.iEnd; i++)
	  {
		  points.push_back(line.controlPt[i]+point);
	  }
	  return points;
  }


  void makeMesh()
  {
    // make sure to re-read the file and to delete copies in visu
    clearParamMap();
    addParamsFromFile("param.dat");
    clearVisu();


//--------------------------------------------------------

    double Rin  = getDoubleParam("Rin");
    double Rout  = getDoubleParam("Rout");
    double angle_deg  = getDoubleParam("angle_deg");
    double stagger_deg  = getDoubleParam("stagger_deg");
    double Rtrailing = getDoubleParam("Rtrailing");
    double Rleading = getDoubleParam("Rleading");
    double Radtrailing = getDoubleParam("Radtrailing");
    double Radleading = getDoubleParam("Radleading");
    double alphalead = getDoubleParam("alphalead");
    double alphatrail = getDoubleParam("alphatrail");
    double thickness = getDoubleParam("thickness");
    double thicknesp = getDoubleParam("thicknesp");
    double posthickness = getDoubleParam("posthickness");
    int ntrail = getIntParam("ntrail");
    int nlead = getIntParam("nlead");
    int nsuction = getIntParam("nsuction");
    int npressure = getIntParam("npressure");
    int nboundtop=getIntParam("nboundtop");
    int nboundbot=getIntParam("nboundbot");
    int nboundinner=getIntParam("nboundinner");
    int nboundouter=getIntParam("nboundouter");
    int nbladeBL = getIntParam("nbladeBL");
    double thickBL =getDoubleParam("thickBL");//thickness of the  boundary layer
    double first_l =getDoubleParam("first_l");//2.5e-4//thickness of the  boundary layer


    double angle_rad = (angle_deg/180.)*M_PI;
    double alphalead_rad=(alphalead/180.)*M_PI;
    double alphatrail_rad=(alphatrail/180.)*M_PI;

    // ----------make the leading and trailing edge
    //leading
    POINT rlead(Rleading,0.0,0.0);
    POINT p(0.0,-Radleading,0.0);
    SPLINE leadingedge = makeRadius(p, Radleading,180.0);
    leadingedge.rotateDegZ(-alphalead);
    deque<POINT> leadingedge2 = translateLine(leadingedge,rlead);
    SPLINE leadingedge3(leadingedge2);
    deque<POINT> leadingedge_dis = discretizeLine(leadingedge3, nlead, false, true);

    //trailing
    POINT rtrail(Rtrailing,0.0,0.0);
    rtrail.rotateDegZ(stagger_deg);
    POINT p1(0.0,-Radtrailing,0.0);
    SPLINE trailingedge = makeRadius(p1, Radtrailing,180.0);
    trailingedge.rotateDegZ(180+alphatrail);
    deque<POINT> trailingedge2 = translateLine(trailingedge,rtrail);
    SPLINE trailingedge3(trailingedge2);
    deque<POINT> trailingedge_dis = discretizeLine(trailingedge3, ntrail, false, true);


    //--------------------make the Bezierline
    //calculate the intersection
    double A = .2;
    POINT bezierleaddirectionpoint(-A*cos(alphalead_rad),A*sin(alphalead_rad),0.0);
    POINT bezierleadpoint= rlead + bezierleaddirectionpoint;
    deque<POINT> bezierleadpoints;
    bezierleadpoints.push_back(rlead); bezierleadpoints.push_back(bezierleadpoint);
    LINE bezierleadline(bezierleadpoints);
    POINT beziertraildirectionpoint(A*cos(alphatrail_rad),A*sin(alphatrail_rad),0.0);
    POINT beziertrailpoint= rtrail + beziertraildirectionpoint;
    deque<POINT> beziertrailpoints;
    beziertrailpoints.push_back(rtrail); beziertrailpoints.push_back(beziertrailpoint);
    LINE beziertrailline(beziertrailpoints);
    POINT intersect = beziertrailline.intersectXY(bezierleadline, true);
    //construct the bezierline
    deque<POINT> bezierpoints;
    bezierpoints.push_back(rlead);bezierpoints.push_back(intersect);bezierpoints.push_back(rtrail);
    BEZIER beziercurve(bezierpoints);
    deque<POINT> bezierpointdis = discretizeLine2(beziercurve,1000, true, false);
    SPLINE bezierline(bezierpointdis);
    POINT rmiddle = bezierline.calcPoint(posthickness);
    //bezierline.rotateDegZ(-stagger_deg/2);
    SPLINE test = bezierline;
    test.rotateDegZ(-stagger_deg/2);

    //------------------------calculate the thickness
    POINT s = bezierline.calcNorm2D(posthickness);;
    POINT middlesuction = rmiddle + thickness* (s-rmiddle);
    POINT middlepressure = rmiddle - thicknesp* (s-rmiddle);


    //-----------------------construct the blade
    //suctionside
    deque<POINT> pointsuction;
    pointsuction.push_back(trailingedge3.calcPoint(.99));
    pointsuction.push_back(trailingedge3.calcPoint(1.0));
    pointsuction.push_back(middlesuction);
    pointsuction.push_back(leadingedge3.calcPoint(0.0));
    pointsuction.push_back(leadingedge3.calcPoint(.01));
    SPLINE suctionside(pointsuction);
    deque<POINT> suctionside_dis = discretizeLineS(suctionside, nsuction, true, true,
    		suctionside.controlPt[1].s,suctionside.controlPt[suctionside.iEnd-1].s  );

    //pressureside
    deque<POINT> pointpressure;
    pointpressure.push_back(trailingedge3.calcPoint(0.01));
    pointpressure.push_back(trailingedge3.calcPoint(0.0));
    pointpressure.push_back(middlepressure);
    pointpressure.push_back(leadingedge3.calcPoint(1.0));
    pointpressure.push_back(leadingedge3.calcPoint(0.99));
    SPLINE pressureside(pointpressure);
    deque<POINT> pressureside_dis = discretizeLineS(pressureside, npressure, true, false,
     		pressureside.controlPt[1].s,suctionside.controlPt[pressureside.iEnd-1].s  );

    //-----------------------make the blade boundary points
    deque<POINT> tmp_blade[] = {suctionside_dis, leadingedge_dis , pressureside_dis, trailingedge_dis };
    deque< deque<POINT> > tmp_blade2 (tmp_blade, tmp_blade+sizeof(tmp_blade)/ sizeof(deque<POINT>));
    deque<POINT> meshpts = addDeques(tmp_blade2);

    for (int i=0; i<meshpts.size();i++)
    {
    	meshpts[i].rotateDegZ(-stagger_deg/2);
    }

    // ------------------------- structured mesh around the blade
    int nblade = meshpts.size();
    double distrBladeRad[nbladeBL];

    firstLengthDistr(first_l, nbladeBL, distrBladeRad);
    STRUCTMESH bladeBL(nblade, nbladeBL);
    deque<POINT> ext_BL_pts;
    for (int i=0; i<nblade; i++)
    {
      POINT tan;    //tangent vector to the surface of the blade profile
      if    ((i==0) || (i==nblade-1)) tan = meshpts[1] - meshpts[meshpts.size()-2];
      else                            tan = meshpts[i+1] - meshpts[i-1];

      POINT rad(tan.y, -tan.x);     //rad is perpendicular to tan and therefore to the blade profile
      double mag = sqrt(rad.x*rad.x + rad.y*rad.y);
      rad.x /= mag;
      rad.y /= mag;

      for (int j=0; j<nbladeBL; j++)
      {

        bladeBL.mesh[i][j] = meshpts[i] + distrBladeRad[j]*thickBL*rad;
        POINT gus(bladeBL.mesh[i][j]) ;
        if ((isnan(gus.x)!=0) || (isnan(gus.y)!=0))
        {firstLengthDistr(first_l, nbladeBL, distrBladeRad);     //only first lenght distribution
          cout<<"i is " << i <<endl;
        }
      }
      POINT tmp1(bladeBL.mesh[i][nbladeBL-1]);
      ext_BL_pts.push_back(tmp1);
    }



    //---------------------construct the boundary
    //make the corners
    POINT a(cos( angle_rad/2 )*Rin , -sin(angle_rad/2)*Rin, 0.0);
    POINT b(cos( angle_rad/2 )*Rout , -sin(angle_rad/2)*Rout, 0.0);

    //make the outer boundary
    SPLINE innerline = makeRadius(a,Rin,angle_deg);
    SPLINE outerline = makeRadius(b,Rout,angle_deg);
    deque<POINT> points3; points3.push_back(a); points3.push_back(b);
    SPLINE bottomline(points3);
    SPLINE topline = bottomline; topline.rotateDegZ(angle_deg);


    //discretize the boundary
    deque<POINT> innerpoints = discretizeLine(innerline, nboundinner, true, false);
    deque<POINT> outerpoints = discretizeLine(outerline, nboundouter, true, true);
    deque<POINT> toppoints = discretizeLine(topline, nboundtop, false, true);
    deque<POINT> bottompoints = discretizeLine(bottomline, nboundbot, false, false);
    deque<POINT> tmp_boundary[] = {innerpoints, toppoints, outerpoints, bottompoints};
    deque< deque<POINT> > tmp_boundary2 (tmp_boundary, tmp_boundary+sizeof(tmp_boundary)/ sizeof(deque<POINT>));

    deque<POINT> boundary = addDeques(tmp_boundary2);

    //make unstructured mesh
    TRIANGULATE mesh;
    mesh.extBoundary = boundary;
    mesh.triangParameters = "pq30a0.01FDY";
    HOLE hole;
    hole.holesPoints = ext_BL_pts;
    rlead.rotateDegZ(-stagger_deg/2);
    hole.insidePoints = rlead;
    mesh.holes.push_back(hole);

    mesh.unstructuredMesh2D();

    //addToDisplay(suctionside);
    //deque<POINT> testt;
    //POINT neww(-s.x, -s.y,0.0);
    //testt.push_back(neww);
    //testt.push_back(rmiddle);
    //testt.push_back(s);
    //LINE linetest(testt);
    //addToDisplay(linetest);

    addToDisplay(bladeBL);
    addToDisplay(test);
    //bezierline.rotateDegZ(-stagger_deg/2);
    //addToDisplay(bezierline);
    //addToDisplay(rmiddle);
    //addToDisplay(blade);
    //addToDisplay(suctionside_dis);
    //addToDisplay(pressureside_dis);

    //addToDisplay(pressureside);


    //addToDisplay(camberline);
    middlesuction.rotateDegZ(-stagger_deg/2);
    middlepressure.rotateDegZ(-stagger_deg/2);
    addToDisplay(middlesuction);
    addToDisplay(middlepressure);
    //addToDisplay(trailingedge_dis[1]);
    //addToDisplay(leadingedge_dis[1]);

    //addToDisplay(rtrail);
    //addToDisplay(rmiddle);
    //addToDisplay(leadingedge3);
    //addToDisplay(trailingedge3);
    //addToDisplay(bezierpointdis);

    //addToDisplay(trailingedge);


//
    addToDisplay(mesh);

//
    addToDisplay(innerline);
    addToDisplay(outerline);
    addToDisplay(topline);
    addToDisplay(bottomline);
//    addToDisplay(innerpoints);
//	addToDisplay(outerpoints);
//	addToDisplay(top);
//	addToDisplay(outerpoints[1]);
//	addToDisplay(innerpoints[1]);
//	addToDisplay(toppoints[1]);
//	addToDisplay(bottompoints[1]);
    //
//    addToDisplay(e);
//    addToDisplay(line2);
//    addToDisplay(line3);
//    addToDisplay(line4);


//POINT c(cos( angle/2 )*Rout , sin(angle/2)*Rout, 0.0);
//    deque<POINT> points;

//    for (int i=0; i<50; i++)
//    	points.pushback
//    points.push_back(a);
//	points.push_back(b);
//	points.push_back(c);

//	points.push_back(d);
//    points.push_back(a);
    //draw the lines
//    LINE ab(points);
//    addToDisplay(ab);
//    addToDisplay(a);
//    addToDisplay(b);
//    addToDisplay(c);
//    addToDisplay(d);
//--------------------------------------------------------
    bCurrTransform = bTransform = false;
    printHelp();
    flagMakeMesh = false;
  }

};




// ------------------------------------------------------------------
//
// Singleton class for opengl call back functions
//
// ------------------------------------------------------------------
class Singleton : public RADIALCOMP
{
public: // singleton accessors
  static Singleton* instanciate(const char *name)
  {
    if(instance == NULL)
      instance = new Singleton(name);//mesh);
    return instance;
  }

  static Singleton* getInstance()
  {
    if (instance == NULL)  throw(-1);
    else return instance;
  }

  static void reMakeMesh()
  {
    if (flagMakeMesh == true) getInstance()->makeMesh();
    getInstance()->coordTransformForVisualization();
    glutPostRedisplay();
  }

  static void display()
  {
    getInstance()->draw();
  }

private:
  static Singleton* instance;  // static instance of EDDY

protected:
  Singleton(const char *name) : RADIALCOMP(name) {} // overloaded constructor
};

// set singleton instance to zero
Singleton * Singleton::instance = NULL;





int main(int argc, char *argv[])
{
  initOpenGl(argc, argv);

  Singleton *radialComp = Singleton::instanciate("param.dat");

  glutDisplayFunc(radialComp->display);
  glutIdleFunc(radialComp->reMakeMesh);
  glutMainLoop();

  return 0;
}






