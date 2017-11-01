
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
//--------------------------------------------------------



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






