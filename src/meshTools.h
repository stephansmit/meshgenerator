/*
 * meshTools.h
 *
 *  Created on: Mar 1, 2012
 *      Author: enrico
 */

#ifndef MESHTOOLS_H_
#define MESHTOOLS_H_

#include "transform.h"




typedef struct
{
  deque<LINE> profile;
  deque<double> span;
  int iLE, iTE1, iTE2;
}BLADE;


typedef struct
{
  UNSTRUCTMESH faceMesh;
  STRUCTMESH mainBL, spliBL, inlet, outlet;
  UNSTRUCTMESH mainTip, spliTip;
  TRIANGULATE unstr;
  double span;
}SPANLEVEL;



class MESHTOOLS : public TRANSFORM, public ParamMap
{
public:
	MESHTOOLS(const char *name) : ParamMap(name) {}

public:

	/*!
	 * \brief read ANSYSBladeGen turbomachinery geometry
	 */
	int readAnsysBladeGenGeometry(string fname, deque<BLADE> &blade);

	/*!
	 * \brief read ANSYSBladeGen blade
	 */
	int readAnsysBladeGenBlade(ifstream &ifile, BLADE &blade, string &lastLineRead);

	/*!
	 * \brief read ANSYSBladeGen profile
	 */
	void readAnsysBladeGenProfile(ifstream &ifile, LINE &prof, int &iLE, int &iTE1, int &iTE2);


	SPLINE arbitraryClustering(vector<Param> *distrPtsFile);


  SPLINE splineClustering(vector<Param> *distrPtsFile);


	/*! \brief Define the tip profile of a blade
	 *
	 *  \param hubLine line which defines the hub profile
	 *  \param shrLine line which defines the shroud profile
	 *  \param tipClearance value of the tip clearance [m]
	 *  \param bladeTip control points of the blade's tip profile
	 */
	void tipProfile(LINE &hubLine, LINE &shrLine, const double tipClearance, BLADE &blade, LINE &bladeTip);


	/*! \brief split the blade into 4 splines: trailing edge, leading edge, bottom and top
	 *
	 * 	\param filePts points defining the blade's profile
	 * 	\param ind1 index identifying the leading edge
	 * 	\param ind2 index identifying the beginning of the trailing edge
	 * 	\param ind3 index identifying the end of the trailing edge
	 * 	\param cLe percentage of the chord defining the leading edge region
	 *
	 */
	void splitBladeSplines(const deque<POINT> &filePts, const int ind1, const int ind2, const int ind3, const double cLE,
			SPLINE &leadEdge, SPLINE &trailEdge, SPLINE &bottom, SPLINE &top);

	/*! \brief rotates nodes of a structured mesh within a range of i
	 *
	 * 	\param mesh structured mesh
	 * 	\param iStart starting index i
	 * 	\param iEnd ending index i
	 * 	\param phiFact multiplying factor of the angle of rotation (<1)
	 *
	 */
	void rotateStructuredNodes(STRUCTMESH &mesh, const int iStart, const int iEnd, const double phiFact);

	/*! \brief makes 2D boundary layer structured mesh
	 *
	 * 	\param spline spline defining the edge
	 * 	\param thickness thickness of the boundary layer
	 * 	\param firstLength dimension of the first cell [dimensional]
	 * 	\param str stretching of the distribution (tanh) of nodes along the edge
	 * 	\param center center of the mesh distribution (tanh) of nodes along the edge
	 * 	\param imax number of nodes along the edge
	 * 	\param jmax number of nodes in radial direction
	 *
	 */
	STRUCTMESH makeMeshEdge(const SPLINE &line, const int imax, const int jmax, const double thickness, const double firstLength, const int interp, const double par1, const double par2);

	STRUCTMESH makeMeshEdgeSpecial(const SPLINE &line, //const SPLINE &line1distr,
	    const int imax, const int jmax, const double thickness,
	    const double firstLength, const int interp, const double par1, const double par2, const double h1, const double h2);


	/*! \brief makes 2D block mesh between two splines
	 *
	 * 	\param bottom lower spline
	 * 	\param top higher spline
	 * 	\param imax number of nodes along the spline
	 * 	\param jmax number of nodes in radial direction
	 *
	 */
	STRUCTMESH structBlockMesh2D(SPLINE &bottom, SPLINE &top, const int imax, const int jmax);


	/*! \brief makes the 2D boundary layer mesh around a blade with a sharp trailing edge (e.g. radial compressor)
	 *
	 * 	\param filePts points defining the blade
	 * 	\param ind1 index identifying the leading edge
	 * 	\param ind2 index identifying the beginning of the trailing edge
	 * 	\param ind3 index identifying the end of the trailing edge
	 * 	\param thickness thickness of the boundary layer
	 * 	\param firstLength dimension of the first cell [dimensional]
	 * 	\param cLe percentage of the chord defining the leading edge region
	 * 	\param nLE number of nodes in the leading edge region
	 * 	\param nBOT number of nodes in the bottom region
	 * 	\param nTE number of nodes in the trailing edge region
	 * 	\param nTOP number of nodes in the top region
	 * 	\param jmax number of nodes in radial direction
	 * 	\param strLE stretching of the distribution (tanh) of nodes along the leading edge
	 * 	\param centerLE center of the distribution (tanh) of nodes along the leading edge
	 * 	\param strBOT stretching of the distribution (tanh) of nodes along the bottom
	 * 	\param centerBOT center of the distribution (tanh) of nodes along the bottom
	 *  \param strTE stretching of the distribution (tanh) of nodes along the trailing edge
	 * 	\param centerTE center of the distribution (tanh) of nodes along the trailing edge
	 * 	\param strTOP stretching of the distribution (tanh) of nodes along the top
	 * 	\param centerTOP center of the distribution (tanh) of nodes along the top
	 * 	\param nRotTE number of "columns" j of nodes to rotate near the trailing edge
	 * 	\param phiFactTOP multiplying factor of the angle of rotation (<1) on the top
	 * 	\param phiFactBOT multiplying factor of the angle of rotation (<1) on the bottom
	 *
	 */
	STRUCTMESH makeSharpTeBladeMesh(const deque<POINT> &filePts, const int ind1, const int ind2, const int ind3, const double thickness, const double firstLength, const double cLE,
      														const int nLE, const int nBOT, const int nTE, const int nTOP, const int jmax,
      														const double strLE, const double centerLE, const double strBOT, const double centerBOT, const double strTE, const double centerTE,
      														const double strTOP, const double centerTOP, const int nRotTE, const double phiFactTOP, const double phiFactBOT);



  /*! \brief makes the 2D boundary layer mesh around a blade with a sharp trailing edge (e.g. radial compressor) and wake
   *
   *  \param filePts points defining the blade
   *  \param ind1 index identifying the leading edge
   *  \param ind2 index identifying the beginning of the trailing edge
   *  \param ind3 index identifying the end of the trailing edge
   *  \param thickness thickness of the boundary layer
   *  \param firstLength dimension of the first cell [dimensional]
   *  \param cLe percentage of the chord defining the leading edge region
   *  \param nLE number of nodes in the leading edge region
   *  \param nBOT number of nodes in the bottom region
   *  \param nTE number of nodes in the trailing edge region
   *  \param nTOP number of nodes in the top region
   *  \param jmax number of nodes in radial direction
   *  \param strLE stretching of the distribution (tanh) of nodes along the leading edge
   *  \param centerLE center of the distribution (tanh) of nodes along the leading edge
   *  \param strBOT stretching of the distribution (tanh) of nodes along the bottom
   *  \param centerBOT center of the distribution (tanh) of nodes along the bottom
   *  \param strTE stretching of the distribution (tanh) of nodes along the trailing edge
   *  \param centerTE center of the distribution (tanh) of nodes along the trailing edge
   *  \param strTOP stretching of the distribution (tanh) of nodes along the top
   *  \param centerTOP center of the distribution (tanh) of nodes along the top
   *  \param nRotTE number of "columns" j of nodes to rotate near the trailing edge
   *  \param wakeLength maximum m coordinate of the wake
   *
   */
  UNSTRUCTMESH makeSharpTeBladeWakeMesh(const deque<POINT> &filePts, const int ind1, const int ind2, const int ind3, const double thickness, const double firstLength, const double cLE,
                                  const int nLE, const int nBOT, const int nTE, const int nTOP, const int jmax,
                                  const double strLE, const double centerLE, const double strBOT, const double centerBOT, const double strTE, const double centerTE,
                                  const double strTOP, const double centerTOP, const int nPTS, const int nWake, const double wakeLength);

  STRUCTMESH makeWakeMesh(STRUCTMESH &meshTOP, STRUCTMESH &meshBOT, STRUCTMESH &meshTE, const double firstLength, const int nWake, const int nPTS, const double wakeLength);


	/*! \brief makes 2D mesh of the tip of a blade
	 *
	 * 	\param nLE number of nodes in the leading edge region
	 * 	\param nBOT number of nodes in the bottom region
	 * 	\param nTE number of nodes in the trailing edge region
	 * 	\param nTOP number of nodes in the top region
	 * 	\param jmaxTC2 number of nodes in the radial direction of the structured mesh in the trailing edge mesh of the tip
	 * 	\param meshBlade boundary layer mesh of the blade at the same span of the tip
	 * 	\param tcSpan span (mean value, adimensional) which defines the location of the tip of the blade
	 *
	 */
	UNSTRUCTMESH makeTipMesh(const int nLE, const int nBOT, const int nTE, const int nTOP, const int jmaxTC2, const STRUCTMESH meshBlade, const double tcSpan);


  /*! \brief Defines the external boundary of the unstructured mesh, this can be different from the periodic boundary if a sructured mesh is placed at inlet and/or outlet
   *
   * \param DistrPerBC point's distribution along the periodic boundary
   * \param PerBCp Upper periodic boundary
   * \param PerBCm Lower periodic boundary
   * \param nPerBC number of nodes on the periodic boundary
   * \param cutIn index defining the start of the boundary for the unstructured mesh on the periodic boundary. should be set to 0 if no structured mesh is present at inlet
   * \param imaxOut number of nodes of the structured mesh at outlet on the periodic boundary. should be set to 0 if no structured mesh is present at outlet
   * \param jmaxOut number of nodes on the vertical edge at outlet
   * \param jmaxIn number of nodes on the vertical edge at inlet
   * \param list of all the nodes defining the external boundary of the unstructured mesh
   *
   */
  deque<POINT> defineExternalBoundaryNodes(SPLINE &DistrPerBC, SPLINE &PerBC, const double mPhi, const double pPhi, const int nPerBC, const int jmaxOut, const int jmaxIn);





	/*! \brief makes 2D mesh for a level defined by a fixed adimensional span
	 *
	 */
  SPANLEVEL meshSpanLevel(bool tipMesh, const deque<POINT> &mainBlade, const deque<POINT> &spliBlade, const double span,
                     const int iLEM, const int iTE1M, const int iTE2M, const int iLES, const int iTE1S, const int iTE2S,
                     const int nLeMain, const int nTopMain, const int nBotMain, const int nTeMain, const int jmaxMain, const double thicknessMain, const double cLeMain,
                     const int nLeSpli, const int nTopSpli, const int nBotSpli, const int nTeSpli, const int jmaxSpli, const double thicknessSpli, const double cLeSpli,
                     const double firstLengthMain, const double firstLengthSpli, const int nPerBC, const int imaxIn, const int jmaxIn, const int imaxOut, const int jmaxOut,
                     deque<SPLINE> &spliHubPerBC, deque<SPLINE> &spliShrPerBC, const double mPhi, const double pPhi, deque<POINT> &extpointsHub, deque<POINT> &extpointsShr,
                     TRIANGULATE hubMesh = TRIANGULATE());


	/*! \brief move a mesh from constant adimensional span to constant dimensional span (needed for tip mesh)
	 *
	 *  \param tcMPS line on which the mesh must be projected (not constant adimensional span S)
	 *  \param spanLevel span level mesh that has to be moved
	 *
	 */
	void moveMeshSpan(deque<POINT> &tcMPS, SPANLEVEL &spanLevel);



  /*! \brief interpolates between fixed layers to create a 3D mesh in the span direction
   *
   *  \param layers2D fixed layers to use as "control points" for the interpolation
   *  \param nLayers numbers of points in spanwise direction
   *  \param hubLength dimension of the first cell at the hub
   *  \param shrLength dimension of the first cell at the shroud
   *  \param markIndex
   *  \param interpScheme interpolation scheme
   *
   */
  UNSTRUCTMESH buildMesh3D(deque<SPANLEVEL> &layers2D, const int nLayers, const double firstLength, const double lastLength, const int startLayer, const string &interpScheme);

  UNSTRUCTMESH mesh2Dto3D(const UNSTRUCTMESH mesh2D, const double a,const double b,const double c, const double Rmax, const int kmax);

  UNSTRUCTMESH buildMesh3DNew(deque<UNSTRUCTMESH> &layers2D, const int nLayers, const double firstLength, const double lastLength, const string &interpScheme);  //GUSTAVO 15/06/16

  UNSTRUCTMESH buildMesh3DNewwithFace(deque<UNSTRUCTMESH> &layers2D, const int nLayers, const double firstLength, const double lastLength, const string &interpScheme) ;  //GUSTAVO 15/06/16

  UNSTRUCTMESH buildMesh3DNewwithFace_obj(deque<UNSTRUCTMESH> &layers2D, const int nLayers, const double firstLength, const double lastLength, const string &interpScheme) ;  //GUSTAVO 22/08/16


  void writeFluentMsh_Interface(const char *filename, const UNSTRUCTMESH meshL, const UNSTRUCTMESH meshR, deque<string> bc_names );  //GUSTAVO 15/06/16

  void writeFluentMsh_Interface_obj(const char *filename, const UNSTRUCTMESH meshL, const UNSTRUCTMESH meshR, deque<string> bc_names );  //GUSTAVO 22/08/16

  void readPts2D(deque<POINT> &pts, const string &name)
  {
    FILE *fp = fopen(name.c_str(), "rt");
    if (fp == NULL)
    {
      printf("couldn't open %s\n", name.c_str());
      throw(-1);
    }

    POINT pt;
    while(fscanf(fp, "%lf%lf", &pt.x, &pt.y) != EOF)  { pts.push_back(pt); }
    fclose(fp);
  }


};

#endif /* MESHTOOLS_H_ */
