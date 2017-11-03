
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
  double calcheight(double a, double b, double c, POINT point)
  {
	  double height;
	  double r = sqrt(pow(point.x,2) + pow(point.y,2));
	  height = a + b*r + c*pow(r,2);
	  return height;
  }

  UNSTRUCTMESH buildMesh3D_gus(UNSTRUCTMESH &unstMesh, const double dz, const int nLayers, const double firstLength, const double lastLength, const int startLayer, const string &interpScheme)
  {
    /* Made by Gustavo J. Otero Date: 01-04-2016
     * For this function generates a 3D mesh from a 2D meshes. It first generates a second 2D mesh at a distance "dz"
     * This process can be also done for a deque of <UNSTRUCTMESH> layer2D but it is necessary that the number of elements match
     *
     *
     */
    UNSTRUCTMESH mesh3D;
    deque<int> elem_v_2D, elem_i_2D;


    deque< deque<POINT> > points2D;
    deque<UNSTRUCTMESH> layers2D;

    // Copying the 2D mesh and the nodes that construct them
    for (int i=0; i<2;i++)
    {
      UNSTRUCTMESH tmp_mesh;
      deque<POINT> tmpPts;

      tmp_mesh= unstMesh;
      //    tmp_mesh(unstMesh);
      for (int j=0;j<tmp_mesh.nodes.size();j++)
      {
        tmp_mesh.nodes[j].pt=unstMesh.nodes[j].pt+POINT(0.0,0.0,i*dz);        //making a copy of the 2D mesh at a certain dz up!!
        tmpPts.push_back(tmp_mesh.nodes[j].pt);
      }
      points2D.push_back(tmpPts);
      layers2D.push_back(tmp_mesh);
    }

    //  addToDisplay(layers2D[0]);
    //  addToDisplay(layers2D[1]);
    //  addToDisplay(points2D[0]);
    //  addToDisplay(points2D[1]);


    elem_v_2D = layers2D[0].elem_v;
    elem_i_2D = layers2D[0].elem_i;


    int nvert = (int)points2D[0].size();
    int nelem = (int)elem_i_2D.size()-1;

    // NODES OF THE 3D MESH

    // Cpying the nodes in the other direction, in this case z!
    deque< deque<POINT> > nodes3D;    //number of nodes in the 2D mesh x number of layers in z


    for (int i=0; i<nvert; i++)
    {
      deque<POINT> nodes_i;

      deque<POINT> splinepts;
      for (int s=0; s<points2D.size(); s++)
        splinepts.push_back(points2D[s][i]);


      LINE span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      firstLastLengthDistr(dz/nLayers, dz/nLayers, nLayers, distr);      //First lenght and last lenght can be used for boundary layer in the spanwise direction
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));

      //      if (interpScheme == "LINE")
      //          {
      //            LINE span(splinepts);
      //            // define distribution in the spanwise direction
      //            double distr[nLayers];
      //            firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      //            for (int l=0; l<nLayers; l++)
      //              nodes_i.push_back(span.calcPoint(distr[l]));
      //          }
      //          else if (interpScheme == "SPLINE")
      //          {
      //            SPLINE span(splinepts);
      //            // define distribution in the spanwise direction
      //            double distr[nLayers];
      //            firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      //            for (int l=0; l<nLayers; l++)
      //              nodes_i.push_back(span.calcPoint(distr[l]));
      //          }
      //          else if (interpScheme == "BEZIER")
      //          {
      //            BEZIER span(splinepts);
      //            // define distribution in the spanwise direction
      //            double distr[nLayers];
      //            firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      //            for (int l=0; l<nLayers; l++)
      //              nodes_i.push_back(span.calcPoint(distr[l]));
      //          }
      //          else
      //          {
      //            cout << "error: interpolation " << interpScheme << " not implemented" << endl;
      //            throw(-10);
      //          }
      nodes3D.push_back(nodes_i);
    }
    for (int l=0; l<nLayers; l++)
      for (int i=0; i<nvert; i++)
      {
        NODE tmp;
        tmp.pt = nodes3D[i][l];
        mesh3D.nodes.push_back(tmp);
      }
    //   mesh3D.nno_i = mesh3D.nodes.size() - nvert;
    // ELEMENTS OF THE 3D MESH
    int k=0;
    mesh3D.elem_i.push_back(k);
    for (int l=0; l<nLayers-1; l++)
      for (int c=0; c<nelem; c++)
      {
        deque<int> tmpElem;
        for (int i=elem_i_2D[c]; i<elem_i_2D[c+1]; i++)
          mesh3D.elem_v.push_back(elem_v_2D[i]+(l+1)*nvert);

        for (int i=elem_i_2D[c]; i<elem_i_2D[c+1]; i++)
          mesh3D.elem_v.push_back(elem_v_2D[i]+(l)*nvert);

        k += 2*(elem_i_2D[c+1]-elem_i_2D[c]);
        mesh3D.elem_i.push_back(k);
      }
    mesh3D.nel_i = mesh3D.elem_i.size() - 1 - nelem;

    return(mesh3D);
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
    POINT rmiddle_tmp = bezierline.calcPoint(posthickness);
    //bezierline.rotateDegZ(-stagger_deg/2);
    SPLINE test = bezierline;
    test.rotateDegZ(-stagger_deg/2);

    //------------------------calculate the thickness
    POINT s = bezierline.calcNorm2D(posthickness);
    POINT rmiddle = rmiddle_tmp + thicknesp* s;
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

    //-----------------------make the blade boundary points
    deque<POINT> tmp_blade[] = {suctionside_dis, leadingedge_dis , pressureside_dis, trailingedge_dis };
    deque< deque<POINT> > tmp_blade2 (tmp_blade, tmp_blade+sizeof(tmp_blade)/ sizeof(deque<POINT>));
    deque<POINT> meshpts = addDeques(tmp_blade2);
    POINT tmp2 = meshpts[0];
    meshpts.push_back(tmp2);
    for (int i=0; i<meshpts.size();i++)
    {
    	meshpts[i].rotateDegZ(-stagger_deg/2);
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

    //make height distribution


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

    UNSTRUCTMESH unstructbladeBL=(bladeBL);

    //make unstructured mesh
    TRIANGULATE mesh;
    mesh.extBoundary = boundary;
    HOLE hole;
    hole.holesPoints = ext_BL_pts;
    rlead.rotateDegZ(-stagger_deg/2);
    hole.insidePoints = rlead;
    mesh.holes.push_back(hole);

    //mesh.triangParameters = "pq30a0.01FDY";
    double triangSize = getDoubleParam("TRIAGSIZE");

    char param[200];
    sprintf(param, "pq30a%fFDY", triangSize);
    printf("size = %s", param);

    mesh.triangParameters = param;
    int smooth = getIntParam("smooth");

    mesh.unstructuredMesh2D();

    //UNSTRUCTMESH combinedmesh =  mesh;
    //mesh.smoothMesh(smooth, 1e-12);
    UNSTRUCTMESH mesh3d = mesh2Dto3D(mesh,0.1,2 );


    addToDisplay(mesh3d);



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






