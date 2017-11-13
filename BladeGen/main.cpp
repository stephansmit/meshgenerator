
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
#include <fstream>
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
  RADIALCOMP(const char *name) : MESHTOOLS(name) {
	  const char* name2 = name;
  }


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


  SPLINE makeRadius(POINT startpoint, double radius, double angle) {
	  int n = 1000;
	  deque<POINT> points;
	  for (int i=0; i<=n; i++) {
		POINT newpoint(rotateDegZ(angle-angle/(double)n*i, POINT(0.0, 0.0), startpoint));
		points.push_back(newpoint);
	  }
	  SPLINE transline(points);
	  return(transline);
  }

  deque<POINT> addDeques(deque< deque<POINT> > deques) {
	  deque<POINT> points;
	  for (int i=0; i<deques.size(); i++) {
		  for (int j=0; j<deques[i].size(); j++)
		  	  {
		  		  points.push_front(deques[i][j]);
		  	  }
	  }
	  return points;
  }

  deque<POINT> discretizeSpline(SPLINE line, int npoints, bool startend, bool reverse) {
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
	  for(int i=k; i< (npoints+l); i++) {
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

  deque<POINT> discretizeLineS(SPLINE line, int npoints, bool startend, bool reverse, double sstart, double send) {
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
  	  for(int i=k; i< (npoints+l); i++){
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

  deque<POINT> discretizeBezier(BEZIER line, int npoints, bool startend, bool reverse) {
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
	  for(int i=k; i< (npoints+l); i++){
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

  deque<POINT> translateLine(SPLINE line, POINT point) {
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
    //clearParamMap();
    //addParamsFromFile("param.dat");
	//ParamMap();
    clearVisu();


    //--------------------------------------------------------

    //---------------------- read all the parameters
    string output_filename_blade_ed = getStringParam("output_filename_blade_ed");
    string output_filename_blade = getStringParam("output_filename_blade");
    string output_filename_boundary = getStringParam("output_filename_boundary");
    string output_filename_topboundary = getStringParam("output_filename_topboundary");
    string output_filename_botboundary = getStringParam("output_filename_botboundary");
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
    double posthickness = getDoubleParam("posthickness");
    double blead = getDoubleParam("btrail");
    double btrail = getDoubleParam("blead");
    double positionblade = getDoubleParam("positionblade");

    int ntrail = getIntParam("ntrail");
    int nlead = getIntParam("nlead");
    int nsuction = getIntParam("nsuction");
    int npressure = getIntParam("npressure");
    int nboundtop0=getIntParam("nboundtop0");
    int nboundtop1=getIntParam("nboundtop1");
    int nboundtop2=getIntParam("nboundtop2");
    int nboundbot0=getIntParam("nboundbot0");
    int nboundbot1=getIntParam("nboundbot1");
    int nboundbot2=getIntParam("nboundbot2");
    int nboundinner=getIntParam("nboundinner");
    int nboundouter=getIntParam("nboundouter");



    double angle_rad = (angle_deg/180.)*M_PI;
    double alphalead_rad=(alphalead/180.)*M_PI;
    double alphatrail_rad=(alphatrail/180.)*M_PI;
    double stagger_rad  = (stagger_deg/180.)*M_PI;



    // ----------make the leading and trailing edge
    //leading
    POINT rlead(Rleading,0.0,0.0);
    POINT p(0.0,-Radleading,0.0);
    SPLINE leadingedge = makeRadius(p, Radleading,180.0);
    leadingedge.rotateDegZ(-alphalead);
    deque<POINT> leadingedge2 = translateLine(leadingedge,rlead);
    SPLINE leadingedge3(leadingedge2);
    deque<POINT> leadingedge_dis = discretizeSpline(leadingedge3, nlead, false, true);
    //trailing
    POINT rtrail(Rtrailing,0.0,0.0);
    rtrail.rotateDegZ(stagger_deg);
    POINT p1(0.0,-Radtrailing,0.0);
    SPLINE trailingedge = makeRadius(p1, Radtrailing,180.0);
    trailingedge.rotateDegZ(180+alphatrail);
    deque<POINT> trailingedge2 = translateLine(trailingedge,rtrail);
    SPLINE trailingedge3(trailingedge2);
    deque<POINT> trailingedge_dis = discretizeSpline(trailingedge3, ntrail, false, true);



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
    deque<POINT> bezierpoints2;

    POINT tmppoint1 = beziertrailline.calcPoint(btrail);
    POINT tmppoint2 = bezierleadline.calcPoint(blead);
    bezierpoints2.push_back(rlead);bezierpoints2.push_back(tmppoint2); bezierpoints2.push_back(tmppoint1);bezierpoints2.push_back(rtrail);
    BEZIER beziercurve2(bezierpoints2);

    //construct the bezierline
    deque<POINT> bezierpoints;
    bezierpoints.push_back(rlead);bezierpoints.push_back(intersect);bezierpoints.push_back(rtrail);
    BEZIER beziercurve(bezierpoints);
    deque<POINT> bezierpointdis = discretizeBezier(beziercurve2,1000, true, false);
    SPLINE bezierline(bezierpointdis);
    POINT rmiddle_tmp = bezierline.calcPoint(posthickness);
    SPLINE test = bezierline;



    //------------------------calculate the thickness
    POINT s = bezierline.calcNorm2D(posthickness);
    POINT rmiddle = rmiddle_tmp; //+ thicknesp* s;
    POINT middlesuction = rmiddle + thickness* s;//(s-rmiddle);
    POINT middlepressure = rmiddle - thickness* s;//(s-rmiddle);



    //-----------------------construct the blade
    //suctionside
    deque<POINT> pointsuction;
    pointsuction.push_back(trailingedge3.calcPoint(.99));
    pointsuction.push_back(trailingedge3.calcPoint(1.0));
    pointsuction.push_back(middlesuction);
    pointsuction.push_back(leadingedge3.calcPoint(0.0));
    pointsuction.push_back(leadingedge3.calcPoint(.01));
    SPLINE suctionside(pointsuction);
    deque<POINT> suctionside_dis = discretizeLineS(suctionside, nsuction, false, true,
    		suctionside.controlPt[1].s,suctionside.controlPt[suctionside.iEnd-1].s  );
    //pressureside
    deque<POINT> pointpressure;
    pointpressure.push_back(trailingedge3.calcPoint(0.01));
    pointpressure.push_back(trailingedge3.calcPoint(0.0));
    pointpressure.push_back(middlepressure);
    pointpressure.push_back(leadingedge3.calcPoint(1.0));
    pointpressure.push_back(leadingedge3.calcPoint(0.99));
    SPLINE pressureside(pointpressure);
    deque<POINT> pressureside_dis = discretizeLineS(pressureside, npressure, false, false,
     		pressureside.controlPt[1].s,pressureside.controlPt[pressureside.iEnd-1].s  );


    //make the camberlines
    deque<POINT> camberlinepoints;
    for (int i=0; i<=100; i++){
    	double s = (double)i/100;
    	camberlinepoints.push_back((pressureside.calcPoint(s)+suctionside.calcPoint(s))/2.0);
    }
    SPLINE camberline_tmp(camberlinepoints);
    POINT point_tmp(0,positionblade,0);
    SPLINE camberline(translateLine(camberline_tmp, point_tmp));
    SPLINE camberlinebot;
    SPLINE camberlinetop;
    camberlinebot = camberline;
    camberlinetop = camberline;
    camberlinetop.rotateDegZ(-angle_deg/2);
    camberlinebot.rotateDegZ(angle_deg/2);
    deque<POINT> tmp_blade[] = {suctionside_dis, leadingedge_dis , pressureside_dis, trailingedge_dis };
    deque< deque<POINT> > tmp_blade2 (tmp_blade, tmp_blade+sizeof(tmp_blade)/ sizeof(deque<POINT>));
    deque<POINT> meshpts = addDeques(tmp_blade2);
    POINT tmp2 = meshpts[0];
    meshpts.push_back(tmp2);
    for (int i=0; i<meshpts.size(); i++){
    	meshpts[i]=meshpts[i]+point_tmp;
    }


    SPLINE bladeline(meshpts);
    deque<POINT> bladeline_equidis = discretizeSpline(bladeline, 10000, true, false);
    addToDisplay(bladeline_equidis);

    //---------------------construct the boundary

    //make the corners
    POINT bstarttop(cos( angle_rad/2)*Rin , -sin(angle_rad/2)*Rin, 0.0);
    POINT bendtop(cos( angle_rad/2 )*Rout , -sin(angle_rad/2)*Rout, 0.0);
    POINT bstartctop = camberlinetop.calcPoint(0);
    deque<POINT> btoppoints;
    btoppoints.push_back(bstarttop);btoppoints.push_back(bstartctop);
    SPLINE btopline0(btoppoints);
    SPLINE bbotline0 = btopline0;
    bbotline0.rotateDegZ(angle_deg);
    POINT bendctop = camberlinetop.calcPoint(1);
    deque<POINT> btoppoints1;
    btoppoints1.push_back(bendctop);btoppoints1.push_back(bendtop);
    SPLINE btopline2(btoppoints1);
    SPLINE bbotline2 = btopline2;
    bbotline2.rotateDegZ(angle_deg);
    //make the outer radial boundary
    SPLINE innerline = makeRadius(bstarttop,Rin,angle_deg);
    SPLINE outerline = makeRadius(bendtop,Rout,angle_deg);
    deque<POINT> points3; points3.push_back(bstarttop); points3.push_back(bendtop);
    SPLINE bottomline(points3);
    SPLINE topline = bottomline; topline.rotateDegZ(angle_deg);
    //discretize the boundary
    deque<POINT> innerpoints = discretizeSpline(innerline, nboundinner, false, true);
    deque<POINT> outerpoints = discretizeSpline(outerline, nboundouter, false, false);
    deque<POINT> toppoints0 = discretizeSpline(btopline0, nboundtop0, true, true);
    deque<POINT> toppoints1 = discretizeSpline(camberlinetop, nboundtop1, false, true);
    deque<POINT> toppoints2 = discretizeSpline(btopline2, nboundtop2, true, true);
    deque<POINT> botpoints0 = discretizeSpline(bbotline0, nboundbot0, true, false);
    deque<POINT> botpoints1 = discretizeSpline(camberlinebot, nboundbot1, false, false);
    deque<POINT> botpoints2 = discretizeSpline(bbotline2, nboundbot2, true, false);
    deque<POINT> tmp_boundary[] = {innerpoints, toppoints0, toppoints1, toppoints2, outerpoints, botpoints2, botpoints1, botpoints0};
    deque< deque<POINT> > tmp_boundary2 (tmp_boundary, tmp_boundary+sizeof(tmp_boundary)/ sizeof(deque<POINT>));
    deque<POINT> boundary = addDeques(tmp_boundary2);


    deque<POINT> tmp_periodicboundtop1[] = {toppoints0, toppoints1, toppoints2};
    deque< deque<POINT> > tmp_periodicboundtop2 (tmp_periodicboundtop1, tmp_periodicboundtop1+sizeof(tmp_periodicboundtop1)/ sizeof(deque<POINT>));
    deque<POINT> periodicboundtop = addDeques(tmp_periodicboundtop2);


    deque<POINT> tmp_periodicboundbot1[] = {botpoints2, botpoints1, botpoints0};
    deque< deque<POINT> > tmp_periodicboundbot2 (tmp_periodicboundtop1, tmp_periodicboundbot1+sizeof(tmp_periodicboundbot1)/ sizeof(deque<POINT>));
    deque<POINT> periodicboundbot = addDeques(tmp_periodicboundbot2);




    // --------------------- write results to file

    //equidistant bladeshape used for parameter analysis
    ofstream myfile;
    myfile.open (output_filename_blade_ed.c_str());
    myfile << "x,y\n";
    for (int i=0; i<bladeline_equidis.size(); i++) {
    	myfile << bladeline_equidis[i].x << ", " << bladeline_equidis[i].y << "\n";
    }
    myfile.close();

    //non-equidistant bladeshape + boundary  used for meshgen
	myfile.open (output_filename_blade.c_str());
	myfile << "x,y\n";
	for (int i=0; i<meshpts.size(); i++) {
		myfile << meshpts[i].x << ", " << meshpts[i].y << "\n";
	}
	myfile.close();

	//non-equidistant boundaries  used for meshgen
	myfile.open (output_filename_boundary.c_str());
	myfile << "x,y\n";
	for (int i=0; i<boundary.size(); i++) {
		myfile << boundary[i].x << ", " << boundary[i].y << "\n";
	}
	myfile.close();

	myfile.open (output_filename_topboundary.c_str());
	myfile << "x,y\n";
	for (int i=0; i<periodicboundtop.size(); i++) {
		myfile << periodicboundtop[i].x << ", " << periodicboundtop[i].y << "\n";
	}
	myfile.close();

	myfile.open (output_filename_botboundary.c_str());
	myfile << "x,y\n";
	for (int i=0; i<periodicboundbot.size(); i++) {
		myfile << periodicboundbot[i].x << ", " << periodicboundbot[i].y << "\n";
	}
	myfile.close();

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
//	makeMesh();
  //initOpenGl(argc, argv);
//
  Singleton *radialComp = Singleton::instanciate(argv[1]);
  radialComp->makeMesh();
//
//  glutDisplayFunc(radialComp->display);
//  glutIdleFunc(radialComp->reMakeMesh);
//  glutMainLoop();

  return 0;
}






