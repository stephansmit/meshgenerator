
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

  void makeMesh()
  {
    // make sure to re-read the file and to delete copies in visu
    clearParamMap();
    addParamsFromFile("param.dat");
    clearVisu();

//    POINT a(0.0, 0.0);
//    POINT b(1.0, 0.0);
//    POINT c(1.0, 2.0);
//    POINT d(0.0, 2.0);
//
//
////
////    addToDisplay(a);
////    addToDisplay(b);
////    addToDisplay(c);
////    addToDisplay(d);
//
//
//    deque<POINT> tmp;
//    tmp.push_back(a);
//    tmp.push_back(b);
//    tmp.push_back(c);
//    tmp.push_back(d);
//
//    addToDisplay(tmp);
//    LINE l1(tmp);
//
//
//
//    int imax = 10;
//    int jmax = 20;
//    STRUCTMESH strMesh(imax,jmax);
//
//
//    for(int i=0; i<imax; i++)
//      for(int j=0; j<jmax; j++)
//      {
//        POINT pA = a + (double)i/(double)(imax-1)*(b-a);
//        POINT pB = d + (double)i/(double)(imax-1)*(c-d);
//        strMesh.mesh[i][j] = pA + (double)j/(double)(jmax-1)*(pB-pA);
//      }
//
////    addToDisplay(strMesh);
//
//    UNSTRUCTMESH unstMesh(strMesh);
//
//
//    addToDisplay(unstMesh);

    // ==========================================================================================
    // READ HUB AND SHROUD POINTS
    // ==========================================================================================
    LINE hubLine("BladeGen_hub2.curve");
    LINE shrLine("BladeGen_shroud2.curve");

    // define hub and shr lines on the r-z plane
    deque<POINT> tmp_pts;
    for (int i=0; i<hubLine.controlPt.size(); i++) tmp_pts.push_back(POINT(hubLine[i].radZ(), 0.0));
    LINE tmp_hub(tmp_pts);
    tmp_pts.clear();

    for (int i=0; i<shrLine.controlPt.size(); i++) tmp_pts.push_back(POINT(shrLine[i].radZ(), shrLine[i].z));
    LINE tmp_shr(tmp_pts);
    tmp_pts.clear();

    // recompute hub and shroud lines
    int nnew = 350;
    deque<POINT> new_hub, new_shr;
    for (int i=0; i<nnew; i++)
    {
      new_hub.push_back(tmp_hub.calcPoint((double)i/(double)(nnew-1)));
      new_shr.push_back(tmp_shr.calcPoint((double)i/(double)(nnew-1)));
    }

    hubLine = LINE(new_hub);
    shrLine = LINE(new_shr);


    hubLine.writeCtrPts("hubLine.dat");
    shrLine.writeCtrPts("shrLine.dat");

//    initTransform(hubLine, shrLine);

//    throw(-1);


    deque<POINT> halfSpan_pts;
    for (int i=0; i<nnew; i++)
      halfSpan_pts.push_back(0.5*(hubLine.controlPt[i]+shrLine.controlPt[i]));
    LINE midspan(halfSpan_pts);


    double radInt = shrLine.controlPt[nnew-1].x; //0.0873*1000.;
    double radOut = shrLine.controlPt[0].x;


    // ==========================================================================================
    // READ PROFILE POINTS
    // ==========================================================================================
    deque<BLADE> blades;
    if (!readAnsysBladeGenGeometry(getStringParam("fileName").c_str(), blades))     throw(-1);

    // ==========================================================================================
    // ADJUST POINTS
    // ==========================================================================================


    blades[0].profile[0].writeCtrPts("blade0.dat");
    blades[0].profile[1].writeCtrPts("blade0.25.dat");
    blades[0].profile[2].writeCtrPts("blade0.5.dat");
    blades[0].profile[3].writeCtrPts("blade0.75.dat");
    blades[0].profile[4].writeCtrPts("blade1.dat");


    double fact = 119.0;

//    double phiRot = 360./36.0; // 1:2
    double phiRot = 360./43.; // 2:5
    double phiSta = 360./18.;


    SPLINE blade_sp( blades[0].profile[2].controlPt);
    LINE blade_li( blades[0].profile[2].controlPt);

    deque<POINT> pts_le;
    pts_le.push_back(blades[0].profile[2].controlPt[42]);
    pts_le.push_back(blades[0].profile[2].controlPt[43]);
    pts_le.push_back(blades[0].profile[2].controlPt[44]);

    pts_le.push_back(blades[0].profile[2].controlPt[45]);
    pts_le.push_back(blades[0].profile[2].controlPt[46]);
    pts_le.push_back(blades[0].profile[2].controlPt[47]);
    pts_le.push_back(blades[0].profile[2].controlPt[54]);
    pts_le.push_back(blades[0].profile[2].controlPt[55]);
    pts_le.push_back(blades[0].profile[2].controlPt[56]);
    pts_le.push_back(blades[0].profile[2].controlPt[57]);
    pts_le.push_back(blades[0].profile[2].controlPt[58]);
    pts_le.push_back(blades[0].profile[2].controlPt[59]);
    pts_le.push_back(blades[0].profile[2].controlPt[60]);
    pts_le.push_back(blades[0].profile[2].controlPt[61]);
    pts_le.push_back(blades[0].profile[2].controlPt[62]);
    pts_le.push_back(blades[0].profile[2].controlPt[63]);
    pts_le.push_back(blades[0].profile[2].controlPt[64]);
    pts_le.push_back(blades[0].profile[2].controlPt[65]);
    pts_le.push_back(blades[0].profile[2].controlPt[66]);
    pts_le.push_back(blades[0].profile[2].controlPt[73]);
    pts_le.push_back(blades[0].profile[2].controlPt[74]);
    pts_le.push_back(blades[0].profile[2].controlPt[75]);
    pts_le.push_back(blades[0].profile[2].controlPt[76]);
    pts_le.push_back(blades[0].profile[2].controlPt[77]);

    pts_le.push_back(blades[0].profile[2].controlPt[78]);
    pts_le.push_back(blades[0].profile[2].controlPt[79]);
    pts_le.push_back(blades[0].profile[2].controlPt[80]);

    BEZIER le_bez(pts_le);

    int ind1=42, ind2=80;

    deque<POINT> mod_profile_pts;
    for (int i=0; i<ind1; i++)  mod_profile_pts.push_back(blades[0].profile[2].controlPt[i]);
    for (int i=0; i<60; i++)    mod_profile_pts.push_back(le_bez.calcPoint((double)i/59.0));
    for (int i=ind2+1; i<blades[0].profile[2].controlPt.size(); i++)  mod_profile_pts.push_back(blades[0].profile[2].controlPt[i]);


    SPLINE tmp(mod_profile_pts);
//    addToDisplay(tmp);
//    addToDisplay(blades[0].profile[2].controlPt);


    // ==========================================================================================
    // ROTOR MESH
    // ==========================================================================================
#if 1

//    LINE perMinLine(halfSpan_pts);
//    for (int i=0; i<perMinLine.controlPt.size(); i++) perMinLine.controlPt[i].rotateDegZ(0.5);
//    LINE perPluLine = perMinLine;
//    for (int i=0; i<perPluLine.controlPt.size(); i++) perPluLine.controlPt[i].rotateDegZ(phiRot);
//    for (int i=0; i<perMinLine.controlPt.size(); i++) perMinLine.controlPt[i].z = perPluLine.controlPt[i].z = 0.0;
//
//
//    addToDisplay(perPluLine);
//    addToDisplay(perMinLine);


    // ==============================
    // BOUNDARY LAYER MESH

    int rotBlPts = getIntParam("ROT_BLPTS");
    int rotSSPts = getIntParam("ROT_SSPTS");
    int rotPSPts = getIntParam("ROT_PSPTS");

    double rotBlFl = getDoubleParam("ROT_BLFL");
    double rotBlTh = getDoubleParam("ROT_BLTH");

//    deque<POINT> profile(mainBl.profile[2].controlPt);
    deque<POINT> profile(mod_profile_pts);
    profile.push_back(profile[0]);
    profile.push_back(profile[1]);
    profile.push_back(profile[2]);
    for (int i=0; i<profile.size(); i++)   profile[i].z = 0.0;

    deque<POINT> ss, ps;
    int indSplit = 69; // 61
    for (int i=0; i<=indSplit; i++)             ss.push_back(profile[i]);
    for (int i=indSplit-2; i<profile.size(); i++) ps.push_back(profile[i]);

    SPLINE ssSpline(ss, 1, ss.size()-2);
    SPLINE psSpline(ps, 1, ps.size()-2);

    STRUCTMESH omeshSS = makeMeshEdge(ssSpline, rotSSPts, rotBlPts, -rotBlTh, rotBlFl, 1, 1.5, 0.4);//, 0.001, 0.0008);
    STRUCTMESH omeshPS = makeMeshEdge(psSpline, rotPSPts, rotBlPts, -rotBlTh, rotBlFl, 1, 1.5, 0.6);//, 0.001, 0.0008);
    for (int j=0; j<omeshSS.jmax; j++) omeshSS.mesh[0][j] = omeshPS.mesh[omeshPS.imax-1][j];
    for (int j=0; j<omeshSS.jmax; j++) omeshSS.mesh[omeshSS.imax-1][j] = omeshPS.mesh[0][j];

    // write blade points
    FILE *file=fopen("blade_rot.dat","wt");
    for (int i=0; i<omeshSS.imax; i++)  fprintf(file, "%1.6le\t%1.6le\n", omeshSS.mesh[i][0].x/1000, -omeshSS.mesh[i][0].y/1000, 1.0);
    for (int i=0; i<omeshPS.imax; i++)  fprintf(file, "%1.6le\t%1.6le\n", omeshPS.mesh[i][0].x/1000, -omeshPS.mesh[i][0].y/1000, 1.0);
    fclose(file);
//    throw(-1);

//    addToDisplay(omeshSS);
//    addToDisplay(omeshPS);

#if 0
    //========================================================================//
    /* calculate area distribution */
    SPLINE ssSpl = ssSpline;
    SPLINE psSpl = psSpline;
    for (int i=0; i<ssSpl.controlPt.size(); i++)
      ssSpl.controlPt[i].rotateDegZ(-phiRot);
    ssSpl.init();

//    addToDisplay(ssSpl);
//    addToDisplay(psSpl);


    int iStart = 8;
    int iEnd   = psSpl.controlPt.size()-8;
    double sStart = psSpl.controlPt[iStart].s;
    double sEnd   = psSpl.controlPt[iEnd].s;
    int npts = 200;
    int nsplit = 6000;
    int nsplit_area = 500;

    double MIN = midspan[0].radZ();
    double MAX = midspan[nnew-1].radZ();

    double area[nnew];
    for (int i=0; i<nnew; i++)
    {
      POINT tmp = hubLine[i] - shrLine[i];
      area[i] = sqrt(tmp.x*tmp.x + tmp.y*tmp.y + tmp.z*tmp.z);
    }

    deque<POINT> mid_points;
    deque<double> area_distr;

    for (int i=0; i<npts; i++)
    {
      double s = sStart + (sEnd-sStart)*(double)i/(double)(npts-1);

      POINT ps_pt = psSpl.calcPoint(s);

      double dist = 1.e20;
      double s_min = 0.0;

      for (int j=0; j<nsplit; j++)
      {
        double s_ss = (double)j/(double)(nsplit-1);
        POINT tmp = ssSpl.calcPoint(s_ss);
        double dist_ = mag(tmp-ps_pt);
        if (dist_<dist)
        {
          s_min = s_ss;
          dist = dist_;
        }
      }

      POINT ss_pt = ssSpl.calcPoint(s_min);
      LINE tmp(ss_pt,ps_pt);

      POINT mid_pt = tmp.calcPoint(0.5);
      mid_points.push_back(mid_pt);

      // calculate area
      double areai = 0.0;

      for (int k=0; k<nsplit_area-1; k++)
      {
        POINT pA = tmp.calcPoint((double)k/(double)(nsplit_area-1));
        POINT pB = tmp.calcPoint((double)(k+1)/(double)(nsplit_area-1));
        POINT pC = 0.5*(pA+pB);
        double rad_i = pC.radZ();
        int ind=1;
        while (midspan[ind].radZ()<rad_i) ind++;
        double fact = (rad_i - midspan[ind-1].radZ()) / (midspan[ind].radZ() - midspan[ind-1].radZ() + 1.0e-20);
        areai += (area[ind-1] + fact*(area[ind]-area[ind-1]))*mag(pA-pB);
      }

      area_distr.push_back(areai);

//      addToDisplay(tmp);

    }

    SPLINE mid_pts_spl(mid_points);

    FILE *fp_ = fopen("passage_area_rotor_43.dat", "wt");
    for (int i=0; i<mid_points.size(); i++)
      fprintf(fp_,"%1.6le\t%1.6le\n", mid_pts_spl[i].s, area_distr[i]);
    fclose(fp_);

    /* end */
    //========================================================================//
#endif


//
//    // write face centers
//    FILE *fp = fopen("xfa_ss.dat", "wt");
//    for (int i=0; i<omeshSS.imax-1; i++)
//    {
//      POINT xfa = 0.5*(omeshSS.mesh[i][0]+omeshSS.mesh[i+1][0]);
//      fprintf(fp, "%1.10le\t%1.10le\n", xfa.x/1000.0, -xfa.y/1000.0);
//    }
//    fclose(fp);
//
//    fp = fopen("xfa_ps.dat", "wt");
//    for (int i=0; i<omeshPS.imax-1; i++)
//    {
//      POINT xfa = 0.5*(omeshPS.mesh[i][0]+omeshPS.mesh[i+1][0]);
//      fprintf(fp, "%1.10le\t%1.10le\n", xfa.x/1000.0, -xfa.y/1000.0);
//    }
//    fclose(fp);

//    addToDisplay(ssSpline.controlPt[0]);
//    addToDisplay(psSpline.controlPt[0]);




    // ===============================
    // Define the periodic boundary
#if 0
    int rotPer = getIntParam("ROT_PER");
    int rotIn  = getIntParam("ROT_IN");
    int rotOut = getIntParam("ROT_OUT");


    int ntmp = 50;
    deque<POINT> camber;
    for (int i=0; i<ntmp; i++)
    {
      double s = (double)i/(double)(ntmp-1);
      POINT pA = psSpline.calcPoint(s);
      POINT pB = ssSpline.calcPoint(1.-s);
      camber.push_back(0.5*(pA+pB));
    }

    camber.push_front(POINT(sqrt(radInt*radInt - camber[0].y*camber[0].y), camber[0].y));
    camber.push_back(POINT(sqrt(radOut*radOut - camber[camber.size()-1].y*camber[camber.size()-1].y), camber[camber.size()-1].y));

    // rotate half pitch
    for (int i=0; i<camber.size(); i++)
      camber[i].rotateDegZ(phiRot/2.);

    SPLINE perSpline(camber);


    deque<POINT> rotExtBound;
    for (int i=0; i<rotPer; i++)  rotExtBound.push_back(perSpline.calcPoint((double)i/(double)(rotPer-1.)));
    for (int i=1; i<rotOut; i++)
    {
      POINT tmp = perSpline.controlPt[perSpline.controlPt.size()-1];
      tmp.rotateDegZ(-phiRot*(double)i/(double)(rotOut-1.));
      rotExtBound.push_back(tmp);
    }

    for (int i=0; i<camber.size(); i++) camber[i].rotateDegZ(-phiRot);
    SPLINE perMinSpline(camber);

    for (int i=1; i<rotPer; i++)  rotExtBound.push_back(perMinSpline.calcPoint(1.-(double)i/(double)(rotPer-1.)));
    for (int i=1; i<rotOut-1; i++)
    {
      POINT tmp = perMinSpline.controlPt[0];
      tmp.rotateDegZ(phiRot*(double)i/(double)(rotOut-1.));
      rotExtBound.push_back(tmp);
    }

    LINE rotExtBoundLine(rotExtBound);
    rotExtBoundLine.writeCtrPts("rot_ext.dat");

    deque<POINT> bladeBound;
    for (int i=0; i<omeshPS.imax; i++) bladeBound.push_back(omeshPS.mesh[i][omeshPS.jmax-1]);
    for (int i=1; i<omeshSS.imax-1; i++) bladeBound.push_back(omeshSS.mesh[i][omeshSS.jmax-1]);
    LINE bladeBoundLine(bladeBound);
    bladeBoundLine.writeCtrPts("rot_int.dat");

#endif


    // ==============================
    // Complete rotor mesh
#if 1
    {
    UNSTRUCTMESH rotorUnst;
    rotorUnst.readGambitNeu2D("rotorUnst.neu");
//    rotorUnst.readGambitNeu2D("rotorUnst1-2.neu");

    UNSTRUCTMESH rotor2D = rotorUnst + (UNSTRUCTMESH)omeshPS + (UNSTRUCTMESH)omeshSS;
    UNSTRUCTMESH rotFinal = rotor2D;
    UNSTRUCTMESH MESH_ROTOR3D;

//    int nBlades = 5;
//    for (int k=1; k<nBlades; k++)
//    {
//      UNSTRUCTMESH rotor2D_copy = rotor2D;
//      for (int i=0; i<rotor2D_copy.nodes.size(); i++)
//        rotor2D_copy.nodes[i].pt.rotateDegZ(phiRot*k);
//
//      rotFinal += rotor2D_copy;
//
//    }
    addToDisplay(rotFinal);

    rotFinal.writeGambitNeu("rotMesh2D.neu", 2);



    // ==============================
    // Make the mesh 3D
    double area[nnew];
    for (int i=0; i<nnew; i++) area[i] = shrLine[i].y;
    // normalize
//    double area_ = 0.0081*1000.0;
    double area_ = 0.05;
    for (int i=0; i<nnew; i++)  area[i] /= area[nnew-1]/(area_);



    for (int i=0; i<nnew; i++) midspan[i].z = midspan[i].y;
    midspan.writeCtrPts("midspan.dat");

//    double MIN = midspan[0].x;
//    double MAX = midspan[nnew-1].x;


    deque<POINT> topNodes, botNodes;

    // variable area
    for (int i=0; i<rotFinal.nodes.size(); i++)
    {
      double rad_i = rotFinal.nodes[i].pt.radZ();

      int ind = 1;
      while (shrLine[ind].x<rad_i && ind<nnew-1) ind++;
      double areai = shrLine[ind-1].y + (rad_i-shrLine[ind-1].x)/(shrLine[ind].x-shrLine[ind-1].x)*(shrLine[ind].y-shrLine[ind-1].y);
//      double areai = area[ind-1] + (rad_i-shrLine[ind-1].x)/(shrLine[ind].x-shrLine[ind-1].x)*(area[ind]-area[ind-1]);

      topNodes.push_back(rotFinal.nodes[i].pt+POINT(0.0,0.0,areai/2.));
      botNodes.push_back(rotFinal.nodes[i].pt+POINT(0.0,0.0,-areai/2.));
    }


    for (int i=0; i<botNodes.size(); i++)
      MESH_ROTOR3D.push_back_node(botNodes[i]);
    for (int i=0; i<topNodes.size(); i++)
      MESH_ROTOR3D.push_back_node(topNodes[i]);


    MESH_ROTOR3D.elem_i.push_back(0);

    for (int i=1; i<rotFinal.elem_i.size(); i++)
    {
      int startJ = rotFinal.elem_i[i-1];
      int endJ = rotFinal.elem_i[i];
      for (int j=startJ; j<endJ; j++)
        MESH_ROTOR3D.elem_v.push_back(rotFinal.elem_v[j]);
      for (int j=startJ; j<endJ; j++)
        MESH_ROTOR3D.elem_v.push_back(rotFinal.elem_v[j]+rotFinal.nodes.size());

      MESH_ROTOR3D.elem_i.push_back(2*(endJ-startJ) + MESH_ROTOR3D.elem_i[MESH_ROTOR3D.elem_i.size()-1]);
    }

    // ===================================== //
    // put it back to its place
    for (int i=0; i<MESH_ROTOR3D.nodes.size(); i++)
    {
      double rad_i = MESH_ROTOR3D.nodes[i].pt.radZ();
//      if (rad_i>MAX) rad_i = MAX;
//      if (rad_i<MIN) rad_i = MIN;

      int ind=1;
      while (midspan[ind].x<rad_i && ind<nnew-1) ind++;
      double fact = (rad_i - midspan[ind-1].x) / (midspan[ind].x - midspan[ind-1].x);
      double dz = midspan.calcPoint(midspan[ind-1].s + fact*(midspan[ind].s-midspan[ind-1].s)).y;

      MESH_ROTOR3D.nodes[i].pt.z += dz;
    }

    for (int i=0; i<rotFinal.nodes.size(); i++)
    {
//      rotFinal.nodes[i].pt =  rotFinal.nodes[i].pt * 1.0/fact;
      rotFinal.nodes[i].pt.y = -rotFinal.nodes[i].pt.y;// * 1.0/fact;
    }
    for (int i=0; i<MESH_ROTOR3D.nodes.size(); i++)
    {
//      MESH_ROTOR3D.nodes[i].pt =  MESH_ROTOR3D.nodes[i].pt * 1.0/fact;
      MESH_ROTOR3D.nodes[i].pt.y = -MESH_ROTOR3D.nodes[i].pt.y;// * 1.0/fact;
    }

    rotFinal.writeGambitNeu("rotMesh2D_43.neu", 2);
    MESH_ROTOR3D.writeTecplot("rotMesh3D_43.plt",3);
    MESH_ROTOR3D.writeGambitNeu("rotMesh3D_43.neu", 3);
//    rotFinal.writeGambitNeu("rotMesh2D_2-5.neu", 2);
//    MESH_ROTOR3D.writeTecplot("rotMesh3D_2-5.plt",3);
//    MESH_ROTOR3D.writeGambitNeu("rotMesh3D_2-5.neu", 3);
//    rotFinal.writeGambitNeu("rotMesh2D_2-5.neu", 2);
//    MESH_ROTOR3D.writeTecplot("rotMesh3D_2-5_flat.plt",3);
//    MESH_ROTOR3D.writeGambitNeu("rotMesh3D_2-5_flat.neu", 3);

    }
#endif



#endif



    // ==========================================================================================
    // STATOR MESH
    // ==========================================================================================
#if 0

    // ===================
    // Define lines
    int isplit1 = 250;
    int isplit2 = 20;
//    LINE newstatleft;
//    for (int i=isplit1; i>=0; i--)
//      newstatleft.controlPt.push_back(statorRight.controlPt[i]);
//    for (int i=1; i<=statorLeft.size()-isplit2+1; i++)
//      newstatleft.controlPt.push_back(statorLeft.controlPt[i]);
//    newstatleft.writeCtrPts("statleft");
    LINE newstatleft("statleft");

//    LINE newstatright;
//    for (int i=isplit1-2; i<statorRight.size(); i++)
//      newstatright.controlPt.push_back(statorRight.controlPt[i]);
//    for (int i=statorLeft.size()-2; i>=statorLeft.size()-isplit2-1; i--)
//      newstatright.controlPt.push_back(statorLeft.controlPt[i]);
//    newstatright.rotateDegZ(rotateright);
//    newstatright.writeCtrPts("statright");
    LINE newstatright("statright");

//    LINE newstatleft1, newstatleft2;
//    for (int i=0; i<newstatleft.size()-1229+1; i++)                   newstatleft1.controlPt.push_back(newstatleft[i]);
//    for (int i=newstatleft.size()-1229-2; i<newstatleft.size(); i++)  newstatleft2.controlPt.push_back(newstatleft[i]);
//    newstatleft1.writeCtrPts("statleft1");
//    newstatleft2.writeCtrPts("statleft2");
    LINE newstatleft1("statleft1");
    LINE newstatleft2("statleft2");

//    LINE newstatright1, newstatright2;
//    for (int i=newstatright.size()-463; i>=0; i--)                        newstatright2.controlPt.push_back(newstatright[i]);
//    for (int i=newstatright.size()-1; i>=newstatright.size()-463-2; i--)  newstatright1.controlPt.push_back(newstatright[i]);
//    newstatright1.writeCtrPts("statright1");
//    newstatright2.writeCtrPts("statright2");
    LINE newstatright1("statright1");
    LINE newstatright2("statright2");


//    addToDisplay(newstatright1);
//    addToDisplay(newstatright2);
//    addToDisplay(newstatleft1);
//    addToDisplay(newstatleft2);


    // Read mesh parameters
    int staLeft1 = getIntParam("STA_LEFT1");
    int staLeft2 = getIntParam("STA_LEFT2");
    int staRight1 = getIntParam("STA_RIGHT1");
    int staRight2 = getIntParam("STA_RIGHT2");
    int staBlPts = getIntParam("STA_BLPTS");
    double staBlTh = getDoubleParam("STA_BLTH");
    double staBlFl = getDoubleParam("STA_BLFL");



    double par1 = getDoubleParam("MESH1_PAR1");
    double par2 = getDoubleParam("MESH1_PAR2");
    SPLINE newstatleft1spl(newstatleft1.controlPt, 1, newstatleft1.size()-2);
    STRUCTMESH mesh1 = makeMeshEdgeSpecial(newstatleft1spl, staLeft1, staBlPts, staBlTh, staBlFl, 3,  par1, par2, 6.0, 1.0);

    par1 = getDoubleParam("MESH2_PAR1");
    par2 = getDoubleParam("MESH2_PAR2");
    SPLINE newstatleft2spl(newstatleft2.controlPt, 1, newstatleft2.size()-2);
    STRUCTMESH mesh2 = makeMeshEdge(newstatleft2spl, staLeft2, staBlPts, staBlTh, staBlFl, 1, par1, par2);


    par1 = getDoubleParam("MESH3_PAR1");
    par2 = getDoubleParam("MESH3_PAR2");
    SPLINE newstatright1spl(newstatright1.controlPt, 1, newstatright1.size()-2);
    STRUCTMESH mesh3 = makeMeshEdge(newstatright1spl, staRight1, staBlPts, staBlTh, staBlFl, 1, par1, par2);


    par1 = getDoubleParam("MESH4_PAR1");
    par2 = getDoubleParam("MESH4_PAR2");
    SPLINE newstatright2spl(newstatright2.controlPt, 1, newstatright2.size()-2);
    STRUCTMESH mesh4 = makeMeshEdgeSpecial(newstatright2spl, staRight2, staBlPts, staBlTh, staBlFl, 3, par1, par2, 1.0, 6.0);


    // make sure the structured meshes have no gap at their boundaries
    for (int j=0; j<mesh1.jmax; j++)  mesh1.mesh[mesh1.imax-1][j] = mesh2.mesh[0][j];
    for (int j=0; j<mesh4.jmax; j++)  mesh4.mesh[0][j] = mesh3.mesh[mesh3.imax-1][j];
    for (int j=0; j<mesh1.jmax; j++)
    {
      POINT pttmp(mesh1.mesh[0][j]);
      pttmp.rotateDegZ(phiSta);
      mesh4.mesh[mesh4.imax-1][j] = pttmp;
    }
    for (int j=0; j<mesh2.jmax; j++)
    {
      POINT pttmp(mesh2.mesh[mesh2.imax-1][j]);
      pttmp.rotateDegZ(phiSta);
      mesh3.mesh[0][j] = pttmp;
    }

//    addToDisplay(mesh1);
//    addToDisplay(mesh2);
//    addToDisplay(mesh3);
//    addToDisplay(mesh4);

    STRUCTMESH mesh1_copy = mesh1;
    STRUCTMESH mesh2_copy = mesh2;

    for (int i=0; i<mesh1_copy.imax; i++)
      for (int j=0; j<mesh1_copy.jmax; j++)
        mesh1_copy.mesh[i][j].rotateDegZ(phiSta);

    for (int i=0; i<mesh2_copy.imax; i++)
      for (int j=0; j<mesh2_copy.jmax; j++)
        mesh2_copy.mesh[i][j].rotateDegZ(phiSta);

//    addToDisplay(mesh1_copy);
//    addToDisplay(mesh2_copy);




    // =====================================================
    // Define external boundary for unstructured mesh
#if 0
    int staPer1 = getIntParam("STA_PER1");
    int staPer2 = getIntParam("STA_PER2");
    int staThroat = getIntParam("STA_THROAT");
    int staIn = getIntParam("STA_IN");
    int staOut = getIntParam("STA_OUT");

    // inlet section
    POINT corner(mesh1.mesh[0][mesh1.jmax-1]);
    deque<POINT> perInMin, perInPlu;
    for (int i=0; i<staPer1; i++)
    {
      POINT tmp = corner + (double)i/(double)(staPer1-1)*corner*0.11/mag(corner);
//      POINT tmp = corner + (double)i/(double)(staPer1-1)*corner*radInt/119/mag(corner);
      perInMin.push_back(tmp);
      tmp.rotateDegZ(phiSta);
      perInPlu.push_back(tmp);
    }

    addToDisplay(perInMin);
    addToDisplay(perInPlu);

    deque<POINT> inletSta;
    for (int i=0; i<staIn; i++)
    {
      POINT tmp = perInPlu[perInPlu.size()-1];
      tmp.rotateDegZ(-(double)i/(double)(staIn-1)*phiSta);
      inletSta.push_back(tmp);
    }
    addToDisplay(inletSta);

    // throat line
    deque<POINT> throatLine;
    for (int i=0; i<staThroat; i++)
    {
      POINT pt1 = mesh1.mesh[mesh1.imax-1][mesh1.jmax-1];
      POINT pt2 = mesh4.mesh[0][mesh1.jmax-1];
      POINT pt = pt1 + (double)i/(double)(staThroat-1)*(pt2-pt1);
      throatLine.push_back(pt);
    }
    addToDisplay(throatLine);

    //unstructured mesh
    deque<POINT> unstrInBoundary;
    for (int i=1; i<inletSta.size(); i++)      unstrInBoundary.push_back(inletSta[i]);
    for (int i=perInMin.size()-2; i>=0; i--)   unstrInBoundary.push_back(perInMin[i]);
    for (int i=1; i<mesh1.imax; i++)           unstrInBoundary.push_back(mesh1.mesh[i][mesh1.jmax-1]);
    for (int i=1; i<throatLine.size(); i++)    unstrInBoundary.push_back(throatLine[i]);
    for (int i=1; i<mesh4.imax; i++)           unstrInBoundary.push_back(mesh4.mesh[i][mesh4.jmax-1]);
    for (int i=1; i<perInPlu.size(); i++)      unstrInBoundary.push_back(perInPlu[i]);

    LINE tmpLine1(unstrInBoundary);
    tmpLine1.writeCtrPts("sta_unstr_in.dat");

    addToDisplay(unstrInBoundary);


    // outlet
    deque<POINT> PerOutMin, PerOutPlu;
    corner = mesh2.mesh[mesh2.imax-1][mesh2.jmax-1];

    double delta = fabs(corner.radZ()-radInt/fact);

    for (int i=0; i<staPer2; i++)
    {
      POINT tmp = corner - (double)i/(double)(staPer2-1)*corner/mag(corner)*delta;
      PerOutMin.push_back(tmp);
      tmp.rotateDegZ(phiSta);
      PerOutPlu.push_back(tmp);
    }
    addToDisplay(PerOutMin);
    addToDisplay(PerOutPlu);

    deque<POINT> outletSta;
    for (int i=0; i<staOut; i++)
    {
      POINT tmp = PerOutMin[PerOutMin.size()-1];
      tmp.rotateDegZ((double)i/(double)(staOut-1)*phiSta);
      outletSta.push_back(tmp);
    }

    deque<POINT> unstrOutBoundary;
    for (int i=throatLine.size()-2; i>=0; i--)    unstrOutBoundary.push_back(throatLine[i]);
    for (int i=1; i<mesh2.imax; i++)              unstrOutBoundary.push_back(mesh2.mesh[i][mesh2.jmax-1]);
    for (int i=1; i<PerOutMin.size(); i++)        unstrOutBoundary.push_back(PerOutMin[i]);
    for (int i=1; i<outletSta.size(); i++)        unstrOutBoundary.push_back(outletSta[i]);
    for (int i=PerOutPlu.size()-2; i>=0; i--)     unstrOutBoundary.push_back(PerOutPlu[i]);
    for (int i=1; i<mesh3.imax; i++)              unstrOutBoundary.push_back(mesh3.mesh[i][mesh3.jmax-1]);

    LINE tmpLine2(unstrOutBoundary);
    tmpLine2.writeCtrPts("sta_unstr_out.dat");

    addToDisplay(unstrOutBoundary);
#endif


    // ===================================
    // Read unstructured mesh
#if 0
    UNSTRUCTMESH staUnstIn;
    staUnstIn.readGambitNeu2D("sta_unst_in.neu");

    UNSTRUCTMESH staUnstOut;
    staUnstOut.readGambitNeu2D("sta_unst_out.neu");


    UNSTRUCTMESH sta2D = staUnstIn + staUnstOut + (UNSTRUCTMESH)mesh1 + (UNSTRUCTMESH)mesh2 + (UNSTRUCTMESH)mesh3 + (UNSTRUCTMESH)mesh4;
    UNSTRUCTMESH sta2D_copy = sta2D;
    UNSTRUCTMESH staFinal = sta2D;

//    for (int i=0; i<sta2D_copy.nodes.size(); i++)
//      sta2D_copy.nodes[i].pt.rotateDegZ(phiSta);
//    staFinal += sta2D_copy;
//    addToDisplay(staFinal);

    staFinal.writeGambitNeu("staFinal2D_2-5.neu",2);
//    staFinal.writeGambitNeu("staFinal2D_1-2.neu",2);

    cout << "number of elements: " << staFinal.elem_i.size()-1 << endl;


    // ==============================
    // Make the mesh 3D

    double dz = 0.5*(hubLine[hubLine.size()-1] + shrLine[shrLine.size()-1]).y;
//    double area = 0.0081*1000.0;
    double area = 0.05;

    UNSTRUCTMESH MESH_STATOR3D;

    deque<POINT> topNodes, botNodes;

    // variable area
    for (int i=0; i<staFinal.nodes.size(); i++)
    {
      topNodes.push_back(staFinal.nodes[i].pt+POINT(0.0,0.0, dz + area/2.0));
      botNodes.push_back(staFinal.nodes[i].pt+POINT(0.0,0.0, dz - area/2.0));
    }

    for (int i=0; i<botNodes.size(); i++)
      MESH_STATOR3D.push_back_node(botNodes[i]);
    for (int i=0; i<topNodes.size(); i++)
      MESH_STATOR3D.push_back_node(topNodes[i]);


    MESH_STATOR3D.elem_i.push_back(0);

    for (int i=1; i<staFinal.elem_i.size(); i++)
    {
      int startJ = staFinal.elem_i[i-1];
      int endJ = staFinal.elem_i[i];
      for (int j=startJ; j<endJ; j++)
        MESH_STATOR3D.elem_v.push_back(staFinal.elem_v[j]);
      for (int j=startJ; j<endJ; j++)
        MESH_STATOR3D.elem_v.push_back(staFinal.elem_v[j]+staFinal.nodes.size());

      MESH_STATOR3D.elem_i.push_back(2*(endJ-startJ) + MESH_STATOR3D.elem_i[MESH_STATOR3D.elem_i.size()-1]);
    }


    for (int i=0; i<MESH_STATOR3D.nodes.size(); i++)
    {
      MESH_STATOR3D.nodes[i].pt.x = MESH_STATOR3D.nodes[i].pt.x*fact;
      MESH_STATOR3D.nodes[i].pt.y = MESH_STATOR3D.nodes[i].pt.y*fact;
    }

    // ===============================================================================
    double alpha1 = atan2(-perMinSpline[0].y, perMinSpline[0].x)*180.0/4.0/atan(1.0); // I change the sign of y
    double alpha2 = atan2(outletSta[staOut-1].y, outletSta[staOut-1].x)*180.0/4.0/atan(1.0);

//    printf(" > RADIUS: %1.10le\n", sqrt(inlet[0].x*inlet[0].x + inlet[0].y*inlet[0].y));

    double delta_rotor = alpha1 - alpha2;

    for (int i=0; i<MESH_STATOR3D.nodes.size(); i++)
//      MESH_STATOR3D.nodes[i].pt.rotateDegZ(delta_rotor);
      MESH_STATOR3D.nodes[i].pt.rotateDegZ(delta_rotor-phiSta);

    MESH_STATOR3D.writeTecplot("staMesh3D_2-5.plt",3);
    MESH_STATOR3D.writeGambitNeu("staMesh3D_2-5.neu", 3);
//    MESH_STATOR3D.writeTecplot("staMesh3D_1-2.plt",3);
//    MESH_STATOR3D.writeGambitNeu("staMesh3D_1-2.neu", 3);
//    addToDisplay(MESH_STATOR3D);



    // write faces
//    FILE *ff=fopen("xfa_ps_sta.dat", "wt");
//    for (int i=0; i<mesh3.imax-1; i++)
//    {
//      POINT fa = 0.5*(mesh3.mesh[i][0] + mesh3.mesh[i+1][0]);
//      fa.rotateDegZ(delta_rotor-1.0*phiSta);
//      fprintf(ff, "%1.10le\t%1.10le\n", fa.x*fact/1000, fa.y*fact/1000);
//    }
//    for (int i=0; i<mesh4.imax-1; i++)
//    {
//      POINT fa = 0.5*(mesh4.mesh[i][0] + mesh4.mesh[i+1][0]);
//      fa.rotateDegZ(delta_rotor-1.0*phiSta);
//      fprintf(ff, "%1.10le\t%1.10le\n", fa.x*fact/1000, fa.y*fact/1000);
//    }
//    fclose(ff);
//
//
//    ff=fopen("xfa_ss_sta.dat", "wt");
//    for (int i=0; i<mesh1_copy.imax-1; i++)
//    {
//      POINT fa = 0.5*(mesh1_copy.mesh[i][0] + mesh1_copy.mesh[i+1][0]);
//      fa.rotateDegZ(delta_rotor-1.*phiSta);
//      fprintf(ff, "%1.10le\t%1.10le\n", fa.x*fact/1000, fa.y*fact/1000);
//    }
//    for (int i=0; i<mesh2_copy.imax-1; i++)
//    {
//      POINT fa = 0.5*(mesh2_copy.mesh[i][0] + mesh2_copy.mesh[i+1][0]);
//      fa.rotateDegZ(delta_rotor-1.*phiSta);
//      fprintf(ff, "%1.10le\t%1.10le\n", fa.x*fact/1000, fa.y*fact/1000);
//    }
//    fclose(ff);


#endif



#endif


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






