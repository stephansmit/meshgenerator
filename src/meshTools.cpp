/*
 * meshTools.cpp
 *
 *  Created on: Mar 19, 2012
 *      Author: enrico
 */
#include "Param.h"
#include "transform.h"
#include "struct.h"
#include "unstruct.h"
#include "meshTools.h"

#include "fstream"
using namespace std;


int MESHTOOLS::readAnsysBladeGenGeometry(string fname, deque<BLADE> &blade)
{
  ifstream ifile(fname.c_str());
  if (!ifile.is_open())
  {
    cerr << "couldn't find file: " << fname << endl;
    return 0;
  }

  string lastLineRead;
  getline(ifile, lastLineRead);
  while(lastLineRead.find("Blade") != string::npos)
  {
    BLADE bladeRead;
    readAnsysBladeGenBlade(ifile, bladeRead, lastLineRead);
    blade.push_back(bladeRead);
  }
  ifile.close();

  cout << "Summary of BladeGen file read:\nFound "<< blade.size() << " blades\n";
  for (int b=0; b<blade.size(); b++)
  {
    cout << "  Blade " << b << " with " << blade[b].profile.size() << " profiles\n";
    for (int p=0; p<blade[b].profile.size(); p++)
      cout << "    Profile " << p << " has " << blade[b].profile[p].controlPt.size() << " points\n";
    cout << "    iLE/iTE1/iTE2 are " << blade[b].iLE << ", " << blade[b].iTE1 << ", " << blade[b].iTE2 << endl;
  }

  return 1;
}

int MESHTOOLS::readAnsysBladeGenBlade(ifstream &ifile, BLADE &blade, string &lastLineRead)
{
  string line;
  getline (ifile, line);
  do
  {
    vector<string> token;
    tokenizeString(token, line, " %#\n\r");

    // get span
    double value;
    if (from_string<double>(value, token[3], std::dec)) blade.span.push_back(value/100.0);
    else
    {
      cerr << "Error: value at pos: 4 of /" <<  line << "/ is not double." << endl;
      return 0;
    }

    LINE profile;
    int iLE, iTE, iTE2;
    readAnsysBladeGenProfile(ifile, profile, blade.iLE, blade.iTE1, blade.iTE2);
    blade.profile.push_back(profile);

    line.clear();
    getline (ifile, line);
  }while(line.find("Profile") != string::npos);
  lastLineRead = line;
  return 1;
}

void MESHTOOLS::readAnsysBladeGenProfile(ifstream &ifile, LINE &prof, int &iLE, int &iTE1, int &iTE2)
{
  string line;
  getline(ifile, line);
  do
  {
    vector<string> token;
    tokenizeString(token, line, " %#\n\r");

    if (token.size() == 4)
    {
      if ((token[3] == "LE")  || (token[3] == "%LE"))  iLE  = prof.size();
      if ((token[3] == "TE")  || (token[3] == "%TE"))  iTE1 = prof.size();
      if ((token[3] == "TE2") || (token[3] == "%TE2")) iTE2 = prof.size();
    }

    double x, y, z;
    from_string<double>(x, token[0], std::dec);
    from_string<double>(y, token[1], std::dec);
    from_string<double>(z, token[2], std::dec);
    prof.controlPt.push_back(POINT(x, y, z));

    line.clear();
    getline (ifile, line);
  }while(line.size()>1);
}

SPLINE MESHTOOLS::arbitraryClustering(vector<Param> *distrPtsFile)
{
  SPLINE distr;
  for (int i=0; i<(*distrPtsFile).size(); i++)
    distr.controlPt.push_back(POINT((*distrPtsFile)[i].getDouble(1), (*distrPtsFile)[i].getDouble(2), 0.0));

  // set the spline interpolant to x and initialize (calc derivative)
  for (int i=0; i<distr.controlPt.size(); i++)
    distr.controlPt[i].s = distr.controlPt[i].x;
  distr.calcDerivative();
  return(distr);
}


SPLINE MESHTOOLS::splineClustering(vector<Param> *distrPtsFile)
{
  SPLINE distr;

  // if you make this dx larger you increase the region in which points are clustered
  // (dy is the one specified on a wider range of x)
  double dx = 0.01;

  distr.controlPt.push_back(POINT(0.0,0.0,0.0));

  for (int i=0; i<(*distrPtsFile).size(); i++)
  {
    double x = (*distrPtsFile)[i].getDouble(1);
    double dy = (*distrPtsFile)[i].getDouble(2);
    double y = x;
    distr.controlPt.push_back(POINT(x-dx,y-dy*dx,0.0));
    distr.controlPt.push_back(POINT(x,y,0.0));
    distr.controlPt.push_back(POINT(x+dx,y+dy*dx,0.0));
  }
  distr.controlPt.push_back(POINT(1.0,1.0,0.0));

  distr.iStart = 0;
  distr.iEnd = distr.controlPt.size()-1;

  for (int i=0; i<distr.controlPt.size(); i++)
    distr.controlPt[i].s = distr.controlPt[i].x;
  distr.calcDerivative();

  return(distr);
}


void MESHTOOLS::tipProfile(LINE &hubLine, LINE &shrLine, const double tipClearance, BLADE &blade, LINE &bladeTip)
{
  // transform shroud points to RZ coordinates, set third coord to zero to perform intersection in RZ plane
  deque<POINT> shrRZ;
  for (int i=0; i<shrLine.size(); i++)
    shrRZ.push_back(POINT(shrLine[i].radZ(), shrLine[i].z, 0.0));
  // define spline in order to get the normal vec at control points
  SPLINE shroudSplineRZ(shrRZ);

  // parallel offset of shroud to tip clearance
  deque<POINT> tipclPtsRZ;
  for (int i=0; i<shroudSplineRZ.size(); i++)
  {
    VEC norm = shroudSplineRZ.calcNorm2D(shroudSplineRZ[i].s);
    // push_back tip clearance offset point
    tipclPtsRZ.push_back(shroudSplineRZ[i] - tipClearance*norm);
  }
  // define spline for intersection used below
  SPLINE tipclSplineRZ(tipclPtsRZ);


  // tip clearance profile for main and splitter blade
  for (int i=0; i<blade.profile[0].size(); i++)
  {
    // transform blade points to RZPhi coordinates
    deque<POINT> tmpPts, tmpPts2;
    for (int s=0; s<blade.profile.size(); s++)
      tmpPts.push_back(POINT(blade.profile[s].controlPt[i].radZ(), blade.profile[s].controlPt[i].z, blade.profile[s].controlPt[i].phiRadZ()/1000.));
    SPLINE tmpSpl(tmpPts);


    // intersect blade with tip clearance spline in RZ
    POINT tmpPoint = tmpSpl.intersectXY(tipclSplineRZ, 1, getDoubleParam("INTERSECT_ACCURACY","1.e-8"));
    bladeTip.controlPt.push_back(POINT(tmpPoint.x,tmpPoint.y,tmpPoint.z*1000.));
  }

  // transform new tip clearance profile back to xyz coordinates
  transformRZPHI_to_XYZ(bladeTip.controlPt);
}



void MESHTOOLS::splitBladeSplines(const deque<POINT> &filePts, const int ind1, const int ind2, const int ind3, const double cLE,
																	SPLINE &leadEdge, SPLINE &trailEdge, SPLINE &bottom, SPLINE &top)
{
	// ==========================================================================================
	// SPLIT IN 3 SPLINES
	// ==========================================================================================
	deque<POINT> botPts, tePts, topPts;

	for (int i=0; i<=ind2; i++)
		botPts.push_back(filePts[i]);

	for (int i=filePts.size()-1; i>=ind3; i--)
			topPts.push_back(filePts[i]);

	for (int i=ind2; i<=ind3; i++)
			tePts.push_back(filePts[i]);

	// ==========================================================================================
	// SPLIT IN LEADING EDGE, BOTTOM, TRAILING EDGE AND TOP
	// ==========================================================================================
	SPLINE bot_(botPts), top_(topPts);
	SPLINE trailEdgeTmp(tePts, 0, tePts.size()-1);


	// Cut the leading edge projecting points on the chord
	int np = 1000;
	POINT le = filePts[0], te = trailEdgeTmp.calcPoint(0.5);
	double chord = sqrt((le.x-te.x)*(le.x-te.x) + (le.y-te.y)*(le.y-te.y));

	double m = (te.y-le.y)/(te.x-le.x) + 1e-20;
	double c = 0.0;


	// BOTTOM
	int j_ = 1;
	while ((c<=cLE) && (j_<np))
	{
		POINT tmp = bot_.calcPoint((double)j_/(double)(np-1));
		double xbar = (tmp.y-le.y + tmp.x/m + m*le.x) / (m + 1/m);
		double ybar = m*(xbar-le.x) + le.y;
		c = sqrt((xbar-le.x)*(xbar-le.x) + (ybar-le.y)*(ybar-le.y)) / chord;
		j_++;
	}
	POINT pB = bot_.calcPoint((double)j_/(double)(np-1));
	double spB = (double)j_/(double)(np-1);

	// TOP
	j_ = 1;
	c = 0.0;
	while ((c<=cLE) && (j_<np))
	{
		POINT tmp = top_.calcPoint((double)j_/(double)(np-1));
		double xbar = (tmp.y-le.y + tmp.x/m + m*le.x) / (m + 1/m);
		double ybar = m*(xbar-le.x) + le.y;
		c = sqrt((xbar-le.x)*(xbar-le.x) + (ybar-le.y)*(ybar-le.y)) / chord;
		j_++;
	}
	POINT pT = top_.calcPoint((double)j_/(double)(np-1));
	double spT = (double)j_/(double)(np-1);

	// Cut the two splines in spB and spT
	deque<POINT> LEPtsT, LEPtsB, LEPts, BPts, TPts;
	top_.split(LEPtsT, TPts, spT);
	bot_.split(LEPtsB, BPts, spB);

	SPLINE bottomTmp(BPts);

	deque<POINT> newTPts;
	for (int i=TPts.size()-1; i>=0; i--)
		newTPts.push_back(TPts[i]);
	SPLINE topTmp(newTPts);


	LEPtsT.push_back(TPts[1]);
	for (int i=LEPtsT.size()-1; i>=0; i--)
			LEPts.push_back(LEPtsT[i]);
	for (int i=0; i<LEPtsB.size(); i++)
			LEPts.push_back(LEPtsB[i]);
	LEPts.push_back(BPts[1]);
	SPLINE leadEdgeTmp(LEPts, 1, LEPts.size()-2);

	// assign new splines
	leadEdge = leadEdgeTmp;
	trailEdge = trailEdgeTmp;
	bottom = bottomTmp;
	top = topTmp;
}


void MESHTOOLS::rotateStructuredNodes(STRUCTMESH &mesh, const int iStart, const int iEnd, const double phiFact)
{
	int iN = abs(iStart-iEnd);
	if (iStart<iEnd)
	{
		for (int i=iStart; i<=iEnd; i++)
		{
			POINT pt2 = mesh.mesh[i][ mesh.jmax-1]- mesh.mesh[i][0];
			POINT pt1 = mesh.mesh[i-1][ mesh.jmax-1]- mesh.mesh[i-1][0];
			double phi = acos((pt1.x*pt2.x + pt1.y*pt2.y) / sqrt(pt1.x*pt1.x + pt1.y*pt1.y) / sqrt(pt2.x*pt2.x + pt2.y*pt2.y))*180/4/atan(1.0);
			for (int j=0; j<mesh.jmax; j++)
				mesh.mesh[i][j] = rotateDegZ(-phiFact*phi*(1.-(double)i/(double)(iN)/3.), mesh.mesh[i][0], mesh.mesh[i][j]);
		}
	}
	else
	{
		for (int i=iStart; i>=iEnd; i--)
		{
			POINT pt1 = mesh.mesh[i][mesh.jmax-1]- mesh.mesh[i][0];
			POINT pt2 = mesh.mesh[i+1][mesh.jmax-1]- mesh.mesh[i+1][0];
			double phi = acos((pt1.x*pt2.x + pt1.y*pt2.y) / sqrt(pt1.x*pt1.x + pt1.y*pt1.y) / sqrt(pt2.x*pt2.x + pt2.y*pt2.y))*180/4/atan(1.0);
			for (int j=0; j<mesh.jmax; j++)
				mesh.mesh[i][j] = rotateDegZ(phiFact*phi*(1-(double)(iStart-i)/(double)(iN)/3.), mesh.mesh[i][0], mesh.mesh[i][j]);
		}
	}
}



STRUCTMESH MESHTOOLS::makeMeshEdge(const SPLINE &line, const int imax, const int jmax, const double thickness, const double firstLength,
                                   const int interp, const double par1, const double par2)
//STRUCTMESH MESHTOOLS::makeMeshEdge(const SPLINE &line, const double thickness, const double firstLength, const double str, const double center, const int imax, const int jmax)
{
	// nodes evenly spaced
	double t[imax];

	for (int i=0; i<imax; i++)
		t[i] = (double)i/(double)(imax-1);


	if (interp==1) // tanh
	{
    for (int i=0; i<imax; i++)
      t[i] = stretchingTanh(t[i], par1, par2);
	}
	else if (interp==2) // atanh
	{
    for (int i=0; i<imax; i++)
      t[i] = stretchingAtanh(t[i], par2, par1);
	}
	else if (interp==3) // first-last length
	{
	  firstLastLengthDistr(par1/line.length, par2/line.length, imax, t);
	}


	// thickness of the mesh
	LINE line1(imax); // this is the line along the blade geometry
	LINE line2(imax); // this is the line with the offset=thickness
	deque <POINT> pts_line2;

	for (int i=0; i<line1.size(); i++)
	{
		POINT pt = line.calcPoint(t[i]);
		POINT norm = line.calcNorm2D(t[i]);

		line1[i] = pt;
		line2[i] = pt + thickness*norm;  // point on the offset line
//		pts_line2.push_back(pt + thickness*norm);  // point on the offset line
	}

//	LINE tmp_line2(pts_line2); // this is the line with the offset=thickness
//	for (int i=0; i<line1.size(); i++)
//	  line2[i] =  tmp_line2.calcPoint(t[i]);

	STRUCTMESH block(imax,jmax);

	for (int i=0; i<line1.size(); i++)
	{
		POINT p1 = line1[i];
		POINT p2 = line2[i];

		double radDist[jmax];
		firstLengthDistr(firstLength/mag(p2-p1), jmax, radDist);

		for (int j=0; j<block.jmax; j++)
		{
			block.mesh[i][j] = p1 + radDist[j]*(p2-p1);
			if(isnan(block.mesh[i][j].x) || isnan(block.mesh[i][j].y))
			  cout<< "Not a Number"<< i <<" <-i, j-> "<< j <<" radDist: " << radDist[j]<<endl;
		}
	}

	return(block);
}

STRUCTMESH MESHTOOLS::makeMeshEdgeSpecial(const SPLINE &line, //const SPLINE &line1distr,
    const int imax, const int jmax, const double thickness, const double firstLength,
    const int interp, const double par1, const double par2, const double h1, const double h2)
//STRUCTMESH MESHTOOLS::makeMeshEdge(const SPLINE &line, const double thickness, const double firstLength, const double str, const double center, const int imax, const int jmax)
{
  // nodes evenly spaced
  double t[imax];

  for (int i=0; i<imax; i++)
    t[i] = (double)i/(double)(imax-1);


  if (interp==1) // tanh
  {
    for (int i=0; i<imax; i++)
      t[i] = stretchingTanh(t[i], par1, par2);
  }
  else if (interp==2) // atanh
  {
    for (int i=0; i<imax; i++)
      t[i] = stretchingAtanh(t[i], par2, par1);
  }
  else if (interp==3) // first-last length
  {
    firstLastLengthDistr(par1/line.length, par2/line.length, imax, t);
  }


  // thickness of the mesh
  LINE line1(imax); // this is the line along the blade geometry
  LINE line2(imax); // this is the line with the offset=thickness
  for (int i=0; i<line1.size(); i++)
  {
    POINT pt = line.calcPoint(t[i]);
    POINT norm = line.calcNorm2D(t[i]);

    double fact = h1 + t[i]*(h2-h1);
//    double fact = t[min(i+1, line1.size()-1)]/t[1];
//    double fact = 1.0;
//    if (i<line1.size()-2) fact = (t[i+1]-t[i])/(t[line1.size()-1]-t[line1.size()-2]);

    line1[i] = pt;
    line2[i] = pt + thickness*norm*fact;  // point on the offset line
  }

  STRUCTMESH block(imax,jmax);

  for (int i=0; i<line1.size(); i++)
  {
    POINT p1 = line1[i];
    POINT p2 = line2[i];

    double radDist[jmax];
    firstLengthDistr(firstLength/mag(p2-p1), jmax, radDist);

    for (int j=0; j<block.jmax; j++)
      block.mesh[i][j] = p1 + radDist[j]*(p2-p1);
  }

  return(block);
}


STRUCTMESH MESHTOOLS::structBlockMesh2D(SPLINE &bottom, SPLINE &top, const int imax, const int jmax)
{
	STRUCTMESH block2D(imax,jmax);

	for (int i=0; i<imax; i++)
	{
		double ti = (double)i/(double)(imax-1.);

		POINT p1, p2;
		p1 = bottom.calcPoint(ti);
		p2 = top.calcPoint(ti);

		for (int j=0; j<jmax; j++)
		{
			double tj = (double)j/(double)(jmax-1.);
			POINT pt = p1 + tj*(p2 - p1);
			block2D.mesh[i][j] = pt;
		}
	}

	return (block2D);
}



STRUCTMESH MESHTOOLS::makeSharpTeBladeMesh(const deque<POINT> &filePts, const int ind1, const int ind2, const int ind3, const double thickness, const double firstLength, const double cLE,
		                                const int nLE, const int nBOT, const int nTE, const int nTOP, const int jmax,
		                                const double strLE, const double centerLE, const double strBOT, const double centerBOT, const double strTE, const double centerTE,
		                                const double strTOP, const double centerTOP, const int nRotTE, const double phiFactTOP, const double phiFactBOT)
{
  // ==========================================================================================
  // SPLIT THE PROFILE INTO 4 SPLINES: LEADING EDGE, TRAILING EDGE, BOTTOM AND TOP
  // ==========================================================================================
	SPLINE leadEdge, trailEdge, bottom, top;
//	deque<POINT> tmpPts = filePts;
	splitBladeSplines(filePts, ind1, ind2, ind3, cLE, leadEdge, trailEdge, bottom, top);

	// ==========================================================================================
  // MESH EDGES
  // ==========================================================================================
  STRUCTMESH meshLE = makeMeshEdge(leadEdge, nLE, jmax, thickness, firstLength, 2, strLE, centerLE);
  STRUCTMESH meshTE = makeMeshEdge(trailEdge, nTE, jmax, thickness, firstLength, 1, strTE, centerTE);

  STRUCTMESH meshTOP = makeMeshEdge(top, nTOP, jmax, thickness, firstLength, 3, 4*mag(meshTE.mesh[nTE-1][0]-meshTE.mesh[nTE-2][0]), mag(meshLE.mesh[0][0]-meshLE.mesh[1][0]));
  STRUCTMESH meshBOT = makeMeshEdge(bottom, nBOT, jmax, thickness, firstLength, 3, mag(meshLE.mesh[nLE-1][0]-meshLE.mesh[nLE-2][0]), 4*mag(meshTE.mesh[0][0]-meshTE.mesh[1][0]));

	// ==========================================================================================
  // ROTATE NODES AT THE TRAILING EDGE TO HAVE A SMOOTH TRANSITION
  // ==========================================================================================
	POINT normTe = trailEdge.calcNorm2D(0.5);

//	 define the angle between the two vectors
	POINT pt1, pt2;
	double phi;
	pt1 = meshTE.mesh[nTE-1][jmax-1]- meshTE.mesh[nTE-1][0];
	pt2 = meshTOP.mesh[0][jmax-1] - meshTOP.mesh[0][0];
	phi = acos((pt1.x*pt2.x + pt1.y*pt2.y) / sqrt(pt1.x*pt1.x + pt1.y*pt1.y) / sqrt(pt2.x*pt2.x + pt2.y*pt2.y))*180/4/atan(1.0);
	for (int j=0; j<jmax; j++)
		meshTOP.mesh[0][j] = rotateDegZ(-0.5*phi, meshTOP.mesh[0][0], meshTOP.mesh[0][j]);

	rotateStructuredNodes(meshTOP, 1, nRotTE, phiFactTOP);

	pt1 = meshBOT.mesh[nBOT-1][jmax-1]- meshBOT.mesh[nBOT-1][0];
	pt2 = meshTE.mesh[0][jmax-1] - meshTE.mesh[0][0];
	phi = acos((pt1.x*pt2.x + pt1.y*pt2.y) / sqrt(pt1.x*pt1.x + pt1.y*pt1.y) / sqrt(pt2.x*pt2.x + pt2.y*pt2.y))*180/4/atan(1.0);
	for (int j=0; j<jmax; j++)
	meshBOT.mesh[nBOT-1][j] = rotateDegZ(0.5*phi, meshBOT.mesh[nBOT-1][0], meshBOT.mesh[nBOT-1][j]);

	rotateStructuredNodes(meshBOT, nBOT-2, nBOT-nRotTE-1, phiFactBOT);

	for (int j=0; j<jmax; j++)
		meshTE.mesh[0][j] = meshBOT.mesh[nBOT-1][j];
	for (int j=0; j<jmax; j++)
		meshTE.mesh[nTE-1][j] = meshTOP.mesh[0][j];

	// Rotate also the trailing edge
	if (nTE%2==0)
	{
		rotateStructuredNodes(meshTE, nTE-2, nTE/2, phiFactTOP);
		rotateStructuredNodes(meshTE, 1, (nTE-2)/2, phiFactBOT);
	}
	else
  {
    rotateStructuredNodes(meshTE, nTE-2, (nTE-1)/2, phiFactTOP);
    rotateStructuredNodes(meshTE, 1, (nTE-3)/2, phiFactBOT);
  }


	// ==========================================================================================
  // MESH TRAILING EDGE (version with no rotation of the trailing edge nodes)
  // ==========================================================================================
//	for (int i=0; i<nTE; i++)
//	{
//		double t = (double)i/(double)(nTE-1);
//			t = stretching(t, strTE, centerTE);
//		POINT pt = trailEdge.calcPoint(t);
//		meshTE.mesh[i][0] = pt;
//	}
//
//	for (int i=0; i<nTE; i++)
//		for (int j=1; j<jmax; j++)
//		{
//			double t = (double)i/(double)(nTE-1);
//			t = stretching(t, strTE, centerTE);
//			POINT pt = meshBOT.mesh[nBOT-1][j] + t*(meshTOP.mesh[0][j] - meshBOT.mesh[nBOT-1][j]);
//			meshTE.mesh[i][j] = pt;
//		}



	// ==========================================================================================
  // ASSEMBLE THE COMPLETE MESH
  // ==========================================================================================
	int imax = nTOP+nLE+nBOT-2+nTE-2+1;	// +1 just to close the mesh, when go unstruct double points are removed
	STRUCTMESH meshBlade(imax, jmax);

	for (int i=0; i<nLE; i++)
		for (int j=0; j<jmax; j++)
			meshBlade.mesh[i][j] = meshLE.mesh[i][j];

	for (int i=nLE-1; i<nBOT+nLE-1; i++)
		for (int j=0; j<jmax; j++)
			meshBlade.mesh[i][j] = meshBOT.mesh[i-nLE+1][j];

	for (int i=nBOT+nLE-1; i<nBOT+nLE+nTE-3; i++)
		for (int j=0; j<jmax; j++)
			meshBlade.mesh[i][j] = meshTE.mesh[i-nBOT-nLE+2][j];

	for (int i=nTE+nLE+nBOT-3; i<nTOP+nLE+nBOT+nTE-3; i++)
		for (int j=0; j<jmax; j++)
			meshBlade.mesh[i][j] = meshTOP.mesh[i-nTE-nLE-nBOT+3][j];

	for (int j=0; j<jmax; j++)
		meshBlade.mesh[imax-1][j] = meshBlade.mesh[0][j];


	for (int i=nLE+nBOT-nRotTE-2; i<nLE+nBOT+nTE+nRotTE-2; i++)
	{
	  LINE tmpline(jmax);
	  for (int j=0; j<jmax; j++)
	    tmpline[j] = meshBlade.mesh[i][j];

	  for (int j=jmax/2; j<jmax; j++)
	  {
	    POINT vec = tmpline[j]-tmpline[j-1];
	    meshBlade.mesh[i][j] = meshBlade.mesh[i][j-1] + 0.4*vec;
//	    meshBlade.mesh[i][j] = meshBlade.mesh[i][j-1]
//              + 0.6*mag(meshBlade.mesh[i+1][j]-meshBlade.mesh[i-1][j])*vec/mag(vec);
	  }
	}



//  LINE tmpline(imax);
//  for (int iter=0; iter<50; iter++)
//  {
//    for (int j=jmax-3; j<jmax; j++)
//    {
//      for (int i=0; i<imax; i++)
//        tmpline[i] = meshBlade.mesh[i][j];
//
//      for (int i=nLE+nBOT-30; i<nLE+nBOT+nTE+20; i++)
//      {
//        POINT pnew = 0.5*(meshBlade.mesh[i+1][j]+meshBlade.mesh[i-1][j]);
////        tmpline[i] = tmpline[i] + 0.4*pow(j/(double)(jmax-1),5.)*(pnew-tmpline[i]);
//        tmpline[i] = tmpline[i] + 0.4*(pnew-tmpline[i]);
//        POINT nvec = tmpline[i] - meshBlade.mesh[i][j-1];
//        double t = mag(meshBlade.mesh[i][j]-meshBlade.mesh[i][j-1]);
//        nvec = nvec/mag(nvec);
//        tmpline.controlPt[i] = meshBlade.mesh[i][j-1] + t*nvec;
//      }
//
//      for (int i=0; i<imax; i++)
//        meshBlade.mesh[i][j] = tmpline[i];
//    }
//  }


	return (meshBlade);
}


STRUCTMESH MESHTOOLS::makeWakeMesh(STRUCTMESH &meshTOP, STRUCTMESH &meshBOT, STRUCTMESH &meshTE, const double firstLength, const int nWake, const int nPTS, const double wakeLength)
{
  int nTOP = meshTOP.imax;
  int nTE  = meshTE.imax;
  int nBOT = meshBOT.imax;
  int jmax = meshBOT.jmax;


  double m = (meshTOP.mesh[0][0].y-meshTOP.mesh[1][0].y)/(meshTOP.mesh[0][0].x-meshTOP.mesh[1][0].x);
  double y1 = meshTOP.mesh[0][0].y + m*(wakeLength-meshTOP.mesh[0][0].x);

  m = (meshTOP.mesh[0][jmax-1].y-meshTOP.mesh[1][jmax-1].y)/(meshTOP.mesh[0][jmax-1].x-meshTOP.mesh[1][jmax-1].x);
  double y2 = meshTOP.mesh[0][jmax-1].y + m*(wakeLength-meshTOP.mesh[0][jmax-1].x);


  m = (meshBOT.mesh[nBOT-1][0].y-meshBOT.mesh[nBOT-2][0].y)/(meshBOT.mesh[nBOT-1][0].x-meshBOT.mesh[nBOT-2][0].x);
  double y3 = meshBOT.mesh[nBOT-1][0].y + m*(wakeLength-meshBOT.mesh[nBOT-1][0].x);

  m = (meshBOT.mesh[nBOT-1][jmax-1].y-meshBOT.mesh[nBOT-2][jmax-1].y)/(meshBOT.mesh[nBOT-1][jmax-1].x-meshBOT.mesh[nBOT-2][jmax-1].x);
  double y4 = meshBOT.mesh[nBOT-1][jmax-1].y + m*(wakeLength-meshBOT.mesh[nBOT-1][jmax-1].x);


  double z = meshTE.mesh[0][0].z;

  //------------------------------------------------------
  // central part of the wake
  //------------------------------------------------------
  STRUCTMESH centWake(nTE,nWake);

  deque<POINT> line1, line2;
  for (int i=0; i<nTE; i++)
    line1.push_back(meshTE.mesh[i][0]);

  LINE wakeEnd(POINT(wakeLength,y3,z),POINT(wakeLength,y1,z));
  for (int i=0; i<nTE; i++)
  {
    double ti = (double)i/(double)(nTE-1);
    ti = stretchingTanh(ti,0.5,2.2);
    line2.push_back(wakeEnd.calcPoint(ti));
  }

  for (int i=0; i<nTE; i++)
  {
    double distr[nWake];
    firstLengthDistr(3e-5/mag(line1[i]-line2[i]), nWake, distr);

    for (int j=0; j<nWake; j++)
      centWake.mesh[i][j] = line1[i] + distr[j]*(line2[i]-line1[i]);
  }


  //------------------------------------------------------
  // upper part of the wake
  //------------------------------------------------------
  deque<POINT> left, right, up, down;

  for(int i=0; i<nWake; i++)
    down.push_back(centWake.mesh[nTE-1][i]);

  double distr2[jmax];
  firstLengthDistr((centWake.mesh[nTE-1][nWake-1].y-centWake.mesh[nTE-2][nWake-1].y)/(y2-y1), jmax, distr2);
  for(int i=0; i<jmax; i++)
    right.push_back(POINT(wakeLength,y1+distr2[i]*(y2-y1),z));

  double distr3[nWake];
  firstLengthDistr(10*mag(centWake.mesh[nTE-1][0]-centWake.mesh[nTE-1][1])/mag(meshTOP.mesh[0][jmax-1]-POINT(wakeLength,y2,z)), nWake, distr3);
  for(int i=0; i<nWake; i++)
    up.push_back(meshTOP.mesh[0][jmax-1]+distr3[i]*(POINT(wakeLength,y2,z) - meshTOP.mesh[0][jmax-1]));


  // define external upper spline
  deque<POINT> pts;
  for (int i=nTOP-1; i>=0; i--)
    pts.push_back(meshTOP.mesh[i][jmax-1]);
  for (int i=0; i<nWake; i++)
    pts.push_back(up[i]);

  int iStart = nTOP-nPTS;

  SPLINE extSpl(pts,iStart,pts.size()-1);

  pts.clear();
  // redistribute points on the spline
  for (int i=0; i<nPTS+nWake; i++)
    pts.push_back(extSpl.calcPoint((double)i/(double)(nPTS+nWake-2)));

  for (int i=0; i<nPTS+nWake; i++)
    extSpl[i+iStart] = pts[i];
  extSpl.init();

  // recompute nodes of meshTOP
  for (int i=0; i<nPTS; i++)
  {
    double distr[jmax];
    firstLengthDistr(firstLength/mag(meshTOP.mesh[nPTS-i-1][0]-extSpl[iStart+i]), jmax, distr);
    for (int j=0; j<jmax; j++)
      meshTOP.mesh[nPTS-i-1][j] = meshTOP.mesh[nPTS-i-1][0] + distr[j]*(extSpl[iStart+i]-meshTOP.mesh[nPTS-i-1][0]);
  }

  up.clear();

  // new left line
  for (int i=0; i<jmax; i++)
    left.push_back(meshTOP.mesh[0][i]);

  // new up line
  for (int i=0; i<nWake; i++)
    up.push_back(extSpl[iStart+nPTS+i-1]);



  STRUCTMESH upWake(nWake,jmax);
  for (int i=0; i<nWake; i++)
    for (int j=0; j<jmax; j++)
    {

      double K = mag(left[j]-left[0])/mag(left[jmax-1]-left[0]) +
                  pow((double)i/(double)(nWake-1),3.)*(mag(right[j]-right[0])/mag(right[jmax-1]-right[0])- mag(left[j]-left[0])/mag(left[jmax-1]-left[0]));
      upWake.mesh[i][j] = down[i] + K*(up[i]-down[i]);
    }


  //------------------------------------------------------
  // lower part of the wake
  //------------------------------------------------------
  left.clear();
  right.clear();
  up.clear();
  down.clear();

  for(int i=0; i<nWake; i++)
    up.push_back(centWake.mesh[0][i]);

  double distr4[jmax];
  firstLengthDistr((centWake.mesh[1][nWake-1].y-centWake.mesh[0][nWake-1].y)/(y3-y4), jmax, distr4);
  for(int i=0; i<jmax; i++)
    right.push_back(POINT(wakeLength,y3+distr4[i]*(y4-y3),z));

  double distr5[nWake];
  firstLengthDistr(10*mag(centWake.mesh[0][0]-centWake.mesh[0][1])/mag(meshBOT.mesh[nBOT-1][jmax-1]-POINT(wakeLength,y4,z)), nWake, distr5);
  for(int i=0; i<nWake; i++)
    down.push_back(meshBOT.mesh[nBOT-1][jmax-1]+distr5[i]*(POINT(wakeLength,y4,z) - meshBOT.mesh[nBOT-1][jmax-1]));


  // define external lower spline
  pts.clear();
  for (int i=0; i<nBOT; i++)
    pts.push_back(meshBOT.mesh[i][jmax-1]);
  for (int i=0; i<nWake; i++)
    pts.push_back(down[i]);

  iStart = nBOT-nPTS;

  SPLINE extSpl2(pts,iStart,pts.size()-1);

  // redistribute points on the spline
  pts.clear();
  for (int i=0; i<nPTS+nWake; i++)
    pts.push_back(extSpl2.calcPoint((double)i/(double)(nPTS+nWake-2)));
  for (int i=0; i<nPTS+nWake; i++)
    extSpl2[i+iStart] = pts[i];
  extSpl2.init();

  // recompute nodes of meshBOT
  for (int i=nBOT-nPTS; i<nBOT; i++)
  {
    double distr[jmax];
    firstLengthDistr(firstLength/mag(meshBOT.mesh[i][0]-extSpl2[i]), jmax, distr);
    for (int j=0; j<jmax; j++)
      meshBOT.mesh[i][j] = meshBOT.mesh[i][0] + distr[j]*(extSpl2[i]-meshBOT.mesh[i][0]);
  }

  down.clear();

  // new left line
  for (int i=0; i<jmax; i++)
    left.push_back(meshBOT.mesh[nBOT-1][i]);

  // new up line
  for (int i=0; i<nWake; i++)
    down.push_back(extSpl2[nBOT-1+i]);



  STRUCTMESH downWake(nWake,jmax);
  for (int i=0; i<nWake; i++)
    for (int j=0; j<jmax; j++)
    {

      double K = mag(left[j]-left[0])/mag(left[jmax-1]-left[0]) +
                  pow((double)i/(double)(nWake-1),3.)*(mag(right[j]-right[0])/mag(right[jmax-1]-right[0])- mag(left[j]-left[0])/mag(left[jmax-1]-left[0]));
      downWake.mesh[i][j] = up[i] + K*(down[i]-up[i]);
    }

  // ==========================================================================================
  // ASSEMBLE THE MESH
  // ==========================================================================================
  STRUCTMESH meshWake(2*jmax+nTE-2,nWake);
  for (int i=0; i<jmax; i++)
    for (int j=0; j<nWake; j++)
      meshWake.mesh[i][j] = downWake.mesh[j][jmax-1-i];

  for (int i=1; i<nTE-1; i++)
    for (int j=0; j<nWake; j++)
      meshWake.mesh[jmax-1+i][j] = centWake.mesh[i][j];

  for (int i=0; i<jmax; i++)
    for (int j=0; j<nWake; j++)
      meshWake.mesh[jmax-1+nTE-1+i][j] = upWake.mesh[j][i];

  return(meshWake);
}



UNSTRUCTMESH MESHTOOLS::makeSharpTeBladeWakeMesh(const deque<POINT> &filePts, const int ind1, const int ind2, const int ind3, const double thickness, const double firstLength,
                                const double cLE, const int nLE, const int nBOT, const int nTE, const int nTOP, const int jmax,
                                const double strLE, const double centerLE, const double strBOT, const double centerBOT, const double strTE, const double centerTE,
                                const double strTOP, const double centerTOP, const int nPTS, const int nWake, const double wakeLength)
{
  // ==========================================================================================
  // SPLIT THE PROFILE INTO 4 SPLINES: LEADING EDGE, TRAILING EDGE, BOTTOM AND TOP
  // ==========================================================================================
  SPLINE leadEdge, trailEdge, bottom, top;
  splitBladeSplines(filePts, ind1, ind2, ind3, cLE, leadEdge, trailEdge, bottom, top);

  // ==========================================================================================
  // MESH EDGES
  // ==========================================================================================
  STRUCTMESH meshLE = makeMeshEdge(leadEdge, nLE, jmax, thickness, firstLength, 2, strLE, centerLE);
  STRUCTMESH meshTE = makeMeshEdge(trailEdge, nTE, jmax, thickness, firstLength, 1, strTE, centerTE);

  STRUCTMESH meshTOP = makeMeshEdge(top, nTOP, jmax, thickness, firstLength, 3, 10*mag(meshTE.mesh[nTE-1][0]-meshTE.mesh[nTE-2][0]), mag(meshLE.mesh[0][0]-meshLE.mesh[1][0]));
  STRUCTMESH meshBOT = makeMeshEdge(bottom, nBOT, jmax, thickness, firstLength, 3, mag(meshLE.mesh[nLE-1][0]-meshLE.mesh[nLE-2][0]), 10*mag(meshTE.mesh[0][0]-meshTE.mesh[1][0]));

  // ==========================================================================================
  // MESH WAKE
  // ==========================================================================================
  STRUCTMESH meshWake = makeWakeMesh(meshTOP, meshBOT, meshTE, firstLength, nWake, nPTS, wakeLength);


  // ==========================================================================================
  // ASSEMBLE THE COMPLETE MESH
  // ==========================================================================================
  // define new external boundary
  deque<POINT> newExt;
  for (int i=0; i<nWake; i++)
    newExt.push_back(meshWake.mesh[meshWake.imax-1][nWake-1-i]);
  for (int i=1; i<nTOP; i++)
    newExt.push_back(meshTOP.mesh[i][jmax-1]);
  for (int i=1; i<nLE; i++)
    newExt.push_back(meshLE.mesh[i][jmax-1]);
  for (int i=1; i<nBOT; i++)
    newExt.push_back(meshBOT.mesh[i][jmax-1]);
  for (int i=1; i<nWake; i++)
    newExt.push_back(meshWake.mesh[0][i]);
  for (int i=1; i<meshWake.imax-1; i++)
    newExt.push_back(meshWake.mesh[i][nWake-1]);


  UNSTRUCTMESH meshBlade = (UNSTRUCTMESH)meshTOP + (UNSTRUCTMESH)meshLE + (UNSTRUCTMESH)meshBOT + (UNSTRUCTMESH)meshWake;

  meshBlade.moveBoundaryNodesEnd(newExt);
  meshBlade.moveBoundaryElementsEnd(newExt);

  return(meshBlade);
}






UNSTRUCTMESH MESHTOOLS::makeTipMesh(const int nLE, const int nBOT, const int nTE, const int nTOP, const int jmaxTC2, const STRUCTMESH meshBlade, const double tcSpan)
{
	int imaxTC1 = nTOP, jmaxTC1 = nTE;

	STRUCTMESH tipCle1(imaxTC1, jmaxTC1);

	POINT p1TE = meshBlade.mesh[nLE+nBOT-2][0], p2TE = meshBlade.mesh[nLE+nBOT+nTE-3][0];
	POINT dpTE = p2TE - p1TE;


	for (int i=0; i<imaxTC1; i++)
		tipCle1.mesh[i][0] = meshBlade.mesh[nLE+nBOT-2-i][0];

	deque<POINT> tmpPt;
	for (int i=nLE+nBOT-2; i<nLE+nBOT-2+nTE; i++ )
		tmpPt.push_back(meshBlade.mesh[i][0]);
	SPLINE te(tmpPt);

	for (int i=0; i<imaxTC1-1; i++)
	{
		POINT p1 = meshBlade.mesh[nLE+nBOT-2-i][0];
		POINT p2 = meshBlade.mesh[nLE+nBOT+nTE-3+i][0];

		for (int j=1; j<jmaxTC1; j++) // add manually the top line
		{
			POINT p4 = meshBlade.mesh[nLE+nBOT+j-2][0];
			POINT dpTEj = p4 - p1TE;
			double t = mag(dpTEj)/mag(dpTE);

			POINT pt = p1 + t *(p2-p1);
			if (i==0)
				pt = te.calcPoint(t);
			tipCle1.mesh[i][j] = pt;
		}
	}

	POINT p1 = meshBlade.mesh[nLE-1][0];
	POINT p2 = meshBlade.mesh[0][0];
	for (int j=1; j<jmaxTC1; j++) // add manually the top line
	{
		POINT p4 = meshBlade.mesh[nLE+nBOT+j-2][0];
		POINT dpTEj = p4 - p1TE;
		double t = mag(dpTEj)/mag(dpTE);

		POINT pt = p1 + t *(p2-p1);
		tipCle1.mesh[imaxTC1-1][j] = pt;
	}



	// PART 2: LE
	int imaxTC2 = nLE;
	STRUCTMESH tipCle2(imaxTC2, jmaxTC2);

	for (int i=0; i<imaxTC2; i++)
		tipCle2.mesh[i][0] = meshBlade.mesh[i][0];

	for (int j=0; j<jmaxTC2; j++)
		tipCle2.mesh[0][j] = tipCle1.mesh[imaxTC1-1][jmaxTC1-1-j];

	for (int i=1; i<imaxTC2-1; i++)
	{
		POINT p1 = meshBlade.mesh[i][0];
		POINT p2 = meshBlade.mesh[i][1];
		POINT dp = p2-p1;

		for (int j=1; j<jmaxTC2; j++) // add manually the top line
		{
			double t = mag(tipCle1.mesh[imaxTC1-1][j]-tipCle1.mesh[imaxTC1-1][0]);

			POINT pt = p1 - t *dp/mag(dp);
			tipCle2.mesh[i][j] = pt;
		}
	}

	for (int j=0; j<jmaxTC2; j++)
		tipCle2.mesh[imaxTC2-1][j] = tipCle1.mesh[imaxTC1-1][j];


	// smooth structured tip clearance mesh at the leading edge
//	LINE tmpline(imaxTC2);
//	for (int iter=0; iter<50; iter++)
//	{
//		for (int j=1; j<jmaxTC2; j++)
//		{
//			for (int i=0; i<imaxTC2; i++)
//				tmpline[i] = tipCle2.mesh[i][j];
//
//			for (int i=1; i<imaxTC2-1; i++)
//			{
//				POINT pnew = 0.5*(tipCle2.mesh[i+1][j]+tipCle2.mesh[i-1][j]);
//				tmpline[i] = tmpline[i] + 0.5*(pnew-tmpline[i]);
////					tmpline[i] = tmpline[i] + (double)(jmaxTC2-j)/(double)(jmaxTC2-1)*(pnew-tmpline[i]);
//				POINT nvec = tmpline[i] - tipCle2.mesh[i][j-1];
//				double t = mag(tipCle2.mesh[i][j]-tipCle2.mesh[i][j-1]);
//				nvec = nvec/mag(nvec);
//				tmpline.controlPt[i] = tipCle2.mesh[i][j-1] + t*nvec;
//			}
//
//			for (int i=0; i<imaxTC2; i++)
//				tipCle2.mesh[i][j] = tmpline[i];
//		}
//	}

	//TRIANGLES
	TRIANGULATE tipCle3;

  deque<POINT> tmp;
  for (int i=0; i<tipCle2.imax; i++)
    tmp.push_back(tipCle2.mesh[i][tipCle2.jmax-1]);


  for (int j=jmaxTC2; j<jmaxTC1-jmaxTC2; j++)
    tmp.push_back(tipCle1.mesh[imaxTC1-1][j]);


  tipCle3.triangParameters = "pq30a0.000002FDY";
  tipCle3.extBoundary = tmp;
  tipCle3.unstructuredMesh2D();
  tipCle3.smoothMesh(getIntParam("nIterSmooth"), getDoubleParam("tollSmooth"));
  for (int j=0; j<tipCle3.nodes.size(); j++)
  	tipCle3.nodes[j].pt.z = tcSpan;

	UNSTRUCTMESH meshTip = (UNSTRUCTMESH)tipCle1 + (UNSTRUCTMESH)tipCle2 + tipCle3;

  deque<POINT> extbnd;
  for (int i=0; i<meshBlade.imax; i++)
  	extbnd.push_back(meshBlade.mesh[i][0]);

//  meshTip.moveBoundaryNodesEnd(extbnd);
//  meshTip.moveBoundaryElementsEnd(extbnd);

	return (meshTip);
}




deque<POINT> MESHTOOLS::defineExternalBoundaryNodes(SPLINE &DistrPerBC, SPLINE &PerBC, const double mPhi, const double pPhi, const int nPerBC, const int jmaxOut, const int jmaxIn)
{
  SPLINE PerBCm = PerBC + POINT(0.0, mPhi*4.*atan(1.0)/180.0, 0.0);
  SPLINE PerBCp = PerBC + POINT(0.0, pPhi*4.*atan(1.0)/180.0, 0.0);
  deque<POINT> extpoints;

  for (int i=0; i<nPerBC; i++)
  {
    double ti = (double)i/(double)(nPerBC-1);
    ti = DistrPerBC.calcPoint(ti).y;
    POINT pt = PerBCm.calcPoint(ti);
    extpoints.push_back(pt);
  }

  // second edge (vertical out)
  for (int j=1; j<jmaxOut-1; j++)
  {
    double tj = (double)j/(double)(jmaxOut-1);
    POINT pt = PerBCm.calcPoint(1.0) + tj*(PerBCp.calcPoint(1.0) - PerBCm.calcPoint(1.0));
    extpoints.push_back(pt);
  }

  // third edge (upper periodic)
  for (int i=nPerBC-1; i>=0; i--)
  {
    double ti = (double)i/(double)(nPerBC-1);
    ti = DistrPerBC.calcPoint(ti).y;
    POINT pt = PerBCp.calcPoint(ti);
    extpoints.push_back(pt);
  }

  // forth edge (vertical inlet)
  for (int j=jmaxIn-2; j>=1; j--)
  {
    double tj = (double)j/(double)(jmaxIn-1);
    POINT pt = PerBCm.calcPoint(0.0) + tj*(PerBCp.calcPoint(0.0) - PerBCm.calcPoint(0.0));
    extpoints.push_back(pt);
  }
  return(extpoints);
}






SPANLEVEL MESHTOOLS::meshSpanLevel(bool tipMesh, const deque<POINT> &mainBlade, const deque<POINT> &spliBlade, const double span,
                   const int iLEM, const int iTE1M, const int iTE2M, const int iLES, const int iTE1S, const int iTE2S,
                   const int nLeMain, const int nTopMain, const int nBotMain, const int nTeMain, const int jmaxMain, const double thicknessMain, const double cLeMain,
                   const int nLeSpli, const int nTopSpli, const int nBotSpli, const int nTeSpli, const int jmaxSpli, const double thicknessSpli, const double cLeSpli,
                   const double firstLengthMain, const double firstLengthSpli, const int nPerBC, const int imaxIn, const int jmaxIn, const int imaxOut, const int jmaxOut,
                   deque<SPLINE> &spliHubPerBC, deque<SPLINE> &spliShrPerBC, const double mPhi, const double pPhi, deque<POINT> &extpointsHub, deque<POINT> &extpointsShr,
                   TRIANGULATE hubMesh)
{
  // ==========================================================================================
  // READ DISTRIBUTIONS PARAMETERS FOR MAIN AND SPLITTER BLADES
  // ==========================================================================================
  double strLeMain = getDoubleParam("strLeMain");
  double centerLeMain = getDoubleParam("centerLeMain");
  double strTeMain = getDoubleParam("strTeMain");
  double centerTeMain = getDoubleParam("centerTeMain");
  double strTopMain = getDoubleParam("strTopMain");
  double centerTopMain = getDoubleParam("centerTopMain");
  double strBotMain = getDoubleParam("strBotMain");
  double centerBotMain = getDoubleParam("centerBotMain");
  int nRotTeMain = getIntParam("nRotTeMain");
  double phiFactTopMain = getDoubleParam("phiFactTopMain");
  double phiFactBotMain = getDoubleParam("phiFactBotMain");
	vector<Param> (*pointInsideMain) = getParamVector("pointInsideMain");


  double strLeSpli = getDoubleParam("strLeSpli");
  double centerLeSpli = getDoubleParam("centerLeSpli");
  double strTeSpli = getDoubleParam("strTeSpli");
  double centerTeSpli = getDoubleParam("centerTeSpli");
  double strTopSpli = getDoubleParam("strTopSpli");
  double centerTopSpli = getDoubleParam("centerTopSpli");
  double strBotSpli = getDoubleParam("strBotSpli");
  double centerBotSpli = getDoubleParam("centerBotSpli");
  int nRotTeSpli = getIntParam("nRotTeSpli");
  double phiFactTopSpli = getDoubleParam("phiFactTopSpli");
  double phiFactBotSpli = getDoubleParam("phiFactBotSpli");
	vector<Param> (*pointInsideSpli) = getParamVector("pointInsideSpli");

	SPANLEVEL spanLevel;

  // ==========================================================================================
  // BOUNDARY LAYER BLADE MESH
  // ==========================================================================================
	spanLevel.span = span;

	spanLevel.mainBL = makeSharpTeBladeMesh(mainBlade, iLEM, iTE1M, iTE2M, thicknessMain, firstLengthMain, cLeMain, nLeMain, nBotMain, nTeMain, nTopMain, jmaxMain,
      	strLeMain, centerLeMain, strBotMain, centerBotMain, strTeMain, centerTeMain, strTopMain, centerTopMain, nRotTeMain, phiFactTopMain, phiFactBotMain);
	spanLevel.spliBL = makeSharpTeBladeMesh(spliBlade, iLEM, iTE1M, iTE2M, thicknessSpli, firstLengthSpli, cLeSpli, nLeSpli, nBotSpli, nTeSpli, nTopSpli, jmaxSpli,
      	strLeSpli, centerLeSpli, strBotSpli, centerBotSpli, strTeSpli, centerTeSpli, strTopSpli, centerTopSpli, nRotTeSpli, phiFactTopSpli, phiFactBotSpli);

	int imaxMain = spanLevel.mainBL.imax;
	int imaxSpli = spanLevel.spliBL.imax;

	// ==========================================================================================
	// STRUCTURED MESH AT INLET AN OUTLET
	// ==========================================================================================
	SPLINE inHubPerBCm = spliHubPerBC[0]+POINT(0.0,mPhi*4.*atan(1.0)/180.0,0.0);
	SPLINE inHubPerBCp = spliHubPerBC[0]+POINT(0.0,pPhi*4.*atan(1.0)/180.0,0.0);
	STRUCTMESH inHub = structBlockMesh2D(inHubPerBCm, inHubPerBCp, imaxIn, jmaxIn);

	if (span==0.0)
		spanLevel.inlet = inHub;
	else
	{
	  SPLINE inShrPerBCm = spliShrPerBC[0]+POINT(0.0,mPhi*4.*atan(1.0)/180.0,0.0);
	  SPLINE inShrPerBCp = spliShrPerBC[0]+POINT(0.0,pPhi*4.*atan(1.0)/180.0,0.0);
	  STRUCTMESH inShr = structBlockMesh2D(inShrPerBCm, inShrPerBCp, imaxIn, jmaxIn);

		STRUCTMESH inMsh(imaxIn, jmaxIn);
		for (int i=0; i<imaxIn; i++)
			for (int j=0; j<jmaxIn; j++)
				inMsh.mesh[i][j] = inHub.mesh[i][j] + span*(inShr.mesh[i][j]-inHub.mesh[i][j]);
		spanLevel.inlet = inMsh;
	}


  SPLINE outHubPerBCm = spliHubPerBC[2]+POINT(0.0,mPhi*4.*atan(1.0)/180.0,0.0);
  SPLINE outHubPerBCp = spliHubPerBC[2]+POINT(0.0,pPhi*4.*atan(1.0)/180.0,0.0);
  STRUCTMESH outHub = structBlockMesh2D(outHubPerBCm, outHubPerBCp, imaxOut, jmaxOut);

	if (span==0.0)
		spanLevel.outlet = outHub;
	else
	{
	  SPLINE outShrPerBCm = spliShrPerBC[2]+POINT(0.0,mPhi*4.*atan(1.0)/180.0,0.0);
	  SPLINE outShrPerBCp = spliShrPerBC[2]+POINT(0.0,pPhi*4.*atan(1.0)/180.0,0.0);
	  STRUCTMESH outShr = structBlockMesh2D(outShrPerBCm, outShrPerBCp, imaxOut, jmaxOut);

	  STRUCTMESH outMsh(imaxOut, jmaxOut);
		for (int i=0; i<imaxOut; i++)
			for (int j=0; j<jmaxOut; j++)
				outMsh.mesh[i][j] = outHub.mesh[i][j] + span*(outShr.mesh[i][j]-outHub.mesh[i][j]);
		spanLevel.outlet = outMsh;
	}


	// ==========================================================================================
	// UNSTRUCTURED MESH
	// ==========================================================================================
	deque<POINT> extpoints;
	for (int j=0; j<extpointsHub.size(); j++)
		extpoints.push_back(extpointsHub[j] + span*(extpointsShr[j]-extpointsHub[j]));

	spanLevel.unstr = TRIANGULATE();

	if (span==0.0)
	{
		// assign holes points
		HOLE holeMainBlade;
		holeMainBlade.insidePoints = POINT((*pointInsideMain)[0].getDouble(1), (*pointInsideMain)[0].getDouble(2));
		for (int i=0; i<imaxMain-1; i++)
			holeMainBlade.holesPoints.push_back(spanLevel.mainBL.mesh[i][jmaxMain-1] );
		spanLevel.unstr.holes.push_back(holeMainBlade);

		HOLE holeSpliBlade;
		holeSpliBlade.insidePoints = POINT((*pointInsideSpli)[0].getDouble(1), (*pointInsideSpli)[0].getDouble(2));
		for (int i=0; i<imaxSpli-1; i++)
			holeSpliBlade.holesPoints.push_back(spanLevel.spliBL.mesh[i][jmaxSpli-1] );
		spanLevel.unstr.holes.push_back(holeSpliBlade);


//				char triPar[100] = getStringParam("triangleParameters").;
		spanLevel.unstr.triangParameters = "pq30a0.0003FDY";
		spanLevel.unstr.extBoundary = extpointsHub;
		spanLevel.unstr.unstructuredMesh2D();
		spanLevel.unstr.smoothMesh(getIntParam("nIterSmooth"), getDoubleParam("tollSmooth"));
	}
	else
	{
		spanLevel.unstr = hubMesh;

		spanLevel.unstr.replacePoints(extpointsHub, extpoints);

		deque<POINT> tmp2;
		for (int i=0; i<imaxMain-1; i++)
			tmp2.push_back(spanLevel.mainBL.mesh[i][jmaxMain-1]);
		spanLevel.unstr.replacePoints(spanLevel.unstr.holes[0].holesPoints, tmp2);

		tmp2.clear();
		for (int i=0; i<imaxSpli-1; i++)
			tmp2.push_back(spanLevel.spliBL.mesh[i][jmaxSpli-1]);
		spanLevel.unstr.replacePoints(spanLevel.unstr.holes[1].holesPoints, tmp2);

		for (int j=0; j<spanLevel.unstr.nodes.size(); j++)
			spanLevel.unstr.nodes[j].pt.z = span;

		spanLevel.unstr.smoothMesh(getIntParam("nIterSmooth"), getDoubleParam("tollSmooth"));
	}

	spanLevel.faceMesh = spanLevel.unstr + (UNSTRUCTMESH)spanLevel.inlet + (UNSTRUCTMESH)spanLevel.outlet + (UNSTRUCTMESH)spanLevel.mainBL + (UNSTRUCTMESH)spanLevel.spliBL;

	if (tipMesh == true)
	{
		spanLevel.mainTip = makeTipMesh(nLeMain, nBotMain, nTeMain, nTopMain, getDoubleParam("jmaxTCMain"), spanLevel.mainBL, span);
		spanLevel.spliTip = makeTipMesh(nLeSpli, nBotSpli, nTeSpli, nTopSpli, getDoubleParam("jmaxTCSpli"), spanLevel.spliBL, span);

		spanLevel.faceMesh += spanLevel.mainTip + spanLevel.spliTip;
	}

	return(spanLevel);
}




void MESHTOOLS::moveMeshSpan(deque<POINT> &tcMPS, SPANLEVEL &spanLevel)
{
	POINT pmin = spanLevel.inlet.mesh[0][0].x;
	POINT pmax = spanLevel.outlet.mesh[spanLevel.outlet.imax-1][0];

	deque<POINT> splPts;
	for (int i=1; i<tcMPS.size()-1; i++)
		splPts.push_back(POINT(tcMPS[i].x, tcMPS[i].z, 0));
		SPLINE splineTcMPS(splPts);
// 		splineTcMPS.writeCtrPts("check.dat");


	double fact1 = pmin.x/tcMPS[0].x;
	double fact2 = (pmax.x-1)/(tcMPS[tcMPS.size()-1].x-1);

	// move complete mesh
	for (int i=0; i<spanLevel.faceMesh.nodes.size(); i++)
	{
		if (spanLevel.faceMesh.nodes[i].pt.x<0)
			{spanLevel.faceMesh.nodes[i].pt.z = tcMPS[0].z; spanLevel.faceMesh.nodes[i].pt.x = spanLevel.faceMesh.nodes[i].pt.x/fact1;}
		else if (spanLevel.faceMesh.nodes[i].pt.x>1)
			{spanLevel.faceMesh.nodes[i].pt.z = tcMPS[tcMPS.size()-1].z; spanLevel.faceMesh.nodes[i].pt.x = 1 + (spanLevel.faceMesh.nodes[i].pt.x-1)/fact2;}
		else
		{
			int ind=1;
//				while ((spanLevel.faceMesh.nodes[i].pt.x>tcMPS[ind].x) && (ind<tcMPS.size()-1)) ind++;
			while ((spanLevel.faceMesh.nodes[i].pt.x>splineTcMPS.controlPt[ind].x) && (ind<splineTcMPS.controlPt.size()-1)) ind++;

			double fact = (spanLevel.faceMesh.nodes[i].pt.x - splineTcMPS.controlPt[ind-1].x) / (splineTcMPS.controlPt[ind].x - splineTcMPS.controlPt[ind-1].x);
			if (fact>1) throw(-1);

			POINT tmpPt = splineTcMPS.calcPoint(splineTcMPS.controlPt[ind-1].s+fact*(splineTcMPS.controlPt[ind].s-splineTcMPS.controlPt[ind-1].s));
			spanLevel.faceMesh.nodes[i].pt.z = tmpPt.y;
//				spanLevel.faceMesh.nodes[i].pt.z = tcMPS[ind-1].z + fact*(tcMPS[ind].z-tcMPS[ind-1].z);
		}
	}
}


UNSTRUCTMESH MESHTOOLS::buildMesh3D(deque<SPANLEVEL> &layers2D, const int nLayers, const double firstLength, const double lastLength, const int startLayer, const string &interpScheme)
{


//	markIndex *= nLayers;

	UNSTRUCTMESH mesh3D;

	deque<int> elem_v_2D, elem_i_2D;
	deque< deque<POINT> > points2D;

	elem_v_2D = layers2D[0].faceMesh.elem_v;
	elem_i_2D = layers2D[0].faceMesh.elem_i;

	if (startLayer==0)
	{
	  for (int s=0; s<layers2D.size(); s++)
    {
      deque<POINT> tmpPts;
      for (int i=0; i<layers2D[s].faceMesh.nodes.size(); i++)
        tmpPts.push_back(layers2D[s].faceMesh.nodes[i].pt);
      points2D.push_back(tmpPts);
    }
	}
	else if (startLayer==1)
	{
    for (int s=layers2D.size()-1; s>=0; s--)
    {
      deque<POINT> tmpPts;
      for (int i=0; i<layers2D[s].faceMesh.nodes.size(); i++)
        tmpPts.push_back(layers2D[s].faceMesh.nodes[i].pt);
      points2D.push_back(tmpPts);
    }
	}
	else {printf("error in startLayer definition: must be either 0 or 1"); throw(-111);}


  int nvert = (int)points2D[0].size();
  int nelem = (int)elem_i_2D.size()-1;


  // NODES OF THE 3D MESH
  deque< deque<POINT> > nodes3D;

  for (int i=0; i<nvert; i++)
  {
    deque<POINT> nodes_i;

    deque<POINT> splinepts;
    for (int s=0; s<points2D.size(); s++)
      splinepts.push_back(points2D[s][i]);

    if (interpScheme == "LINE")
    {
      LINE span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));
    }
    else if (interpScheme == "SPLINE")
    {
      SPLINE span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));
    }
    else if (interpScheme == "BEZIER")
    {
      BEZIER span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));
    }
    else
    {
      cout << "error: interpolation " << interpScheme << " not implemented" << endl;
      throw(-10);
    }

    nodes3D.push_back(nodes_i);
  }

  for (int l=0; l<nLayers; l++)
    for (int i=0; i<nvert; i++)
    {
      NODE tmp;
      tmp.pt = nodes3D[i][l];
      mesh3D.nodes.push_back(tmp);
    }

//  for (int l=0; l<nLayers; l++)
//  	for (int i=0; i<nvert; i++)
//		{
//
//			deque<POINT> splinepts;
//			for (int s=0; s<points2D.size(); s++)
//				splinepts.push_back(points2D[s][i]);
//
//			POINT pt;
//
//
//
//
//
//			if (interpScheme == "LINE")
//			{
//				LINE span(splinepts);
//	      // define distribution in the spanwise direction
//	      double distr[nLayers];
//	      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
//        pt = span.calcPoint(distr[l]);
//
//			}
//			else if (interpScheme == "SPLINE")
//			{
//				SPLINE span(splinepts);
//	      // define distribution in the spanwise direction
//	      double distr[nLayers];
//	      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
//				pt = span.calcPoint(distr[l]);
//			}
//			else if (interpScheme == "BEZIER")
//			{
//				BEZIER span(splinepts);
//        // define distribution in the spanwise direction
//        double distr[nLayers];
//        firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
//				pt = span.calcPoint(distr[l]);
//			}
//			else
//			{
//				cout << "error: interpolation " << interpScheme << " not implemented" << endl;
//				throw(-10);
//			}
//
//			NODE tmpNod;
//			tmpNod.pt = pt;
//
//
//			mesh3D.nodes.push_back(tmpNod);
//		}
  mesh3D.nno_i = mesh3D.nodes.size() - nvert;




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



//
//UNSTRUCTMESH MESHTOOLS::mesh2Dto3D(const UNSTRUCTMESH mesh2D, const double dz)
//{
//  UNSTRUCTMESH mesh3D;
//
//  for (int i=0; i<mesh2D.nodes.size(); i++)
//    mesh3D.nodes.push_back(mesh2D.nodes[i]);
//  for (int i=0; i<mesh2D.nodes.size(); i++)
//    mesh3D.push_back_node(mesh2D.nodes[i].pt+POINT(0.0,0.0,dz));
//
//  mesh3D.elem_i.push_back(0);
//
//  for (int i=1; i<mesh2D.elem_i.size(); i++)
//  {
//    int startJ = mesh2D.elem_i[i-1];
//    int endJ = mesh2D.elem_i[i];
//    for (int j=startJ; j<endJ; j++)
//      mesh3D.elem_v.push_back(mesh2D.elem_v[j]);
//    for (int j=startJ; j<endJ; j++)
//      mesh3D.elem_v.push_back(mesh2D.elem_v[j]+mesh2D.nodes.size());
//
//    mesh3D.elem_i.push_back(2*(endJ-startJ) + mesh3D.elem_i[mesh3D.elem_i.size()-1]);
//  }
//  return(mesh3D);
//}
UNSTRUCTMESH MESHTOOLS::mesh2Dto3D(const UNSTRUCTMESH mesh2D, const double a, const double b, const double c, const double Rmax, const int kmax=2)
{
   UNSTRUCTMESH mesh3D;

   //--------------------------------
   // Copying the elements
   deque<int> elem_v_2D, elem_i_2D;
   elem_v_2D = mesh2D.elem_v;
   elem_i_2D = mesh2D.elem_i;
   //Node count: all nodes and internal nodes
   int nvert =     (int)mesh2D.nodes.size();
   int nno_i =     (int)mesh2D.nno_i;
   int nno_bc =    nvert-nno_i;
   //Element count: all elements and internal elements
   int nelem =     (int)elem_i_2D.size()-1;
   int nel_i =     (int)mesh2D.nel_i;
   int nel_bc =    nelem-nel_i;
   //--------------------------------


//****************************************************************
   //      3D NODES
   // 1. Internal nodes of the 2D mesh
   //    1.1 Bottom internal nodes
   //    1.2 Top internal nodes
   // 2. Boundary nodes of the 2D mesh
   //    2.1 Bottom boundary nodes
   //    2.2 Top boundary nodes
   // They are divided like this to make the face generation easier...


   //=============================================================
   //-----------------------------------------
   //    1. internal nodes of the 2D mesh
   //         1.1 From bottom to top  nodes
   for (int k=0; k<kmax; k++)
     for (int i=0; i<nno_i; i++)
     {
       NODE tmp = mesh2D.nodes[i];
       POINT tmppoint = tmp.pt;

       tmp.pt += POINT(0.0,0.0,(double)k/((double)kmax-1)*(a+b*(Rmax-tmppoint.radZ())+ c*(Rmax-tmppoint.radZ())*(Rmax-tmppoint.radZ())));
       mesh3D.nodes.push_back(tmp);
     }
   //-----------------------------------------
   //    2. boundary nodes of the 2D mesh
   //         2.1 Bottom  nodes
   for (int k=0; k<kmax; k++)
     for (int i=nno_i; i<nvert; i++)
     {
       NODE tmp = mesh2D.nodes[i];
       POINT tmppoint = tmp.pt;
       tmp.pt += POINT(0.0,0.0,(double)k/((double)kmax-1)*(a+b*(Rmax-tmppoint.radZ())+ c*(Rmax-tmppoint.radZ())*(Rmax-tmppoint.radZ())));
       mesh3D.nodes.push_back(tmp);
     }
   int nno_i_3D = kmax*nno_i;
   int nvert_3D = kmax*nvert;

   if((int)mesh3D.nodes.size()- nvert_3D != 0)
   {
     cout<<"Error, the number of vertices of each mesh is not the same"<<endl;
     //      cout<<"Nvert: "<<nvert_3D<<" mesh3D:
     //"<<mesh3D.nodes.size()<<endl;
     throw(-1);
   }


//****************************************************************
   // ELEMENTS OF THE 3D MESH,
   int j=0;        //element count

   //-------------------
   // 1.1 Internal elements of the original 2D mesh
   mesh3D.elem_i.push_back(j);
   for (int k=0; k<kmax-1; k++)
     for (int c=0; c<nelem; c++) //nel_i //nelem
     {
       for (int i=elem_i_2D[c]; i<elem_i_2D[c+1]; i++)
         if(elem_v_2D[i]>=nno_i)
mesh3D.elem_v.push_back(elem_v_2D[i]+(k)  *nno_bc +(kmax-1)*nno_i);
         else mesh3D.elem_v.push_back(elem_v_2D[i]+(k)  *nno_i);

       for (int i=elem_i_2D[c]; i<elem_i_2D[c+1]; i++)
         if(elem_v_2D[i]>=nno_i)
mesh3D.elem_v.push_back(elem_v_2D[i]+(k+1)*nno_bc +(kmax-1)*nno_i);
         else mesh3D.elem_v.push_back(elem_v_2D[i]+(k+1)*nno_i);


       j += 2*(elem_i_2D[c+1]-elem_i_2D[c]);
       mesh3D.elem_i.push_back(j);
     }



   return(mesh3D);

}
UNSTRUCTMESH MESHTOOLS::buildMesh3DNew(deque<UNSTRUCTMESH> &layers2D, const int nLayers, const double firstLength, const double lastLength, const string &interpScheme) //, const int startLayer, const string &interpScheme)   //GUSTAVO 15/06/16
{

  /* Made by Gustavo J. Otero Date: 01-04-2016
   * For this function generates a 3D mesh from two or more meshes (they do not need to be 2D_. That is why the input is a deque of meshes.
   *
   *
   *
   */
  UNSTRUCTMESH mesh3D;
  deque<int> elem_v_2D, elem_i_2D;

  for(int i=0; i<layers2D.size(); i++)
  {
    cout<<"Number of vertices in original mesh layer "<< i << ": "<< layers2D[0].nodes.size() <<endl;
    if((i>0) && i<layers2D.size()-1)
      if(layers2D[i-1].nodes.size()- layers2D[i].nodes.size() != 0)
      {
        cout<<"Error, the number of vertices of each mesh is not the same"<<endl;
        throw(-1);
      }
  }

//  cout<<"Number of vertices in original mesh layer 1  "<< layers2D[0].nodes.size() <<endl;
//  cout<<"Number of vertices in original mesh layer 2  "<< layers2D[1].nodes.size() <<endl;
//  if(layers2D[0].nodes.size()- layers2D[1].nodes.size() != 0)
//  {
//    cout<<"Error, the number of vertices of each mesh is not the same"<<endl;
//    throw(-1);
//  }


  deque< deque<POINT> > points2D;

  // Copying the 2D mesh and the nodes that construct them
  for (int i=0; i<layers2D.size();i++)
  {
    deque<POINT> tmpPts;
    for (int j=0;j<layers2D[0].nodes.size();j++)
      tmpPts.push_back(layers2D[i].nodes[j].pt);
    points2D.push_back(tmpPts);
  }

  elem_v_2D = layers2D[0].elem_v;
  elem_i_2D = layers2D[0].elem_i;



  int nvert = (int)points2D[0].size();
  int nelem = (int)elem_i_2D.size()-1;

  // NODES OF THE 3D MESH

  // Cpying the nodes in the other direction, in this case z!
  deque< deque<POINT> > nodes3D;    //number of nodes in the 2D mesh x number of layers in z

//  string interpScheme="SPLINE";

  for (int i=0; i<nvert; i++)
  {
    deque<POINT> nodes_i;

    deque<POINT> splinepts;
    for (int s=0; s<points2D.size(); s++)
      splinepts.push_back(points2D[s][i]);


    if (interpScheme == "LINE")
    {
      LINE span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      for (int i=0; i<nLayers; i++)
        distr[i] = (double)i/(double)(nLayers-1);
//      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));
    }
    else if (interpScheme == "SPLINE")
    {
      SPLINE span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      for (int i=0; i<nLayers; i++)
        distr[i] = (double)i/(double)(nLayers-1);
//      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));
    }
    else if (interpScheme == "BEZIER")
    {
      BEZIER span(splinepts);
      // define distribution in the spanwise direction
      double distr[nLayers];
      for (int i=0; i<nLayers; i++)
        distr[i] = (double)i/(double)(nLayers-1);
//      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
      for (int l=0; l<nLayers; l++)
        nodes_i.push_back(span.calcPoint(distr[l]));
    }
    else
    {
      cout << "error: interpolation " << interpScheme << " not implemented" << endl;
      throw(-10);
    }

    nodes3D.push_back(nodes_i);








    /*
    deque<POINT> nodes_i;

    deque<POINT> splinepts;
    for (int s=0; s<points2D.size(); s++)
      splinepts.push_back(points2D[s][i]);

    LINE span(splinepts);
    // define distribution in the spanwise direction
    double distr[nLayers];
    //firstLastLengthDistr(dz/nLayers, dz/nLayers, nLayers, distr);      //First lenght and last lenght can be used for boundary layer in the spanwise direction

    //-------------------------------------------------------
    double z1=layers2D[1].nodes[i].pt.z;
    double z0=layers2D[0].nodes[i].pt.z;
    double span_distance=z1-z0; //layers2D[i].nodes[j].pt


    for (int i=0; i<nLayers; i++)
      distr[i] = (double)i/(double)(nLayers-1);
    //firstLastLengthDistr(span_distance/nLayers, span_distance/nLayers, nLayers, distr);
    //-------------------------------------------------------

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
    */
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

UNSTRUCTMESH MESHTOOLS::buildMesh3DNewwithFace(deque<UNSTRUCTMESH> &layers2D, const int nLayers, const double firstLength, const double lastLength, const string &interpScheme)    //GUSTAVO 15/06/16
 {
  /* Made by Gustavo J. Otero Date: 01-04-2016
     * For this function generates a 3D mesh from two or more meshes (they do not need to be 2D_. That is why the input is a deque of meshes.
     *  It includes the generatation of the faces, later use for the creating of the Fluent *msh file
     *
     *
     */
   UNSTRUCTMESH mesh3D;
   deque<int> elem_v_2D, elem_i_2D;
   deque<int> face_elem_2D, face_i_2D;

//   cout<<"Number of vertices in original mesh   "<< layers2D[0].nodes.size() <<endl;
//   cout<<"Number of elements in original mesh   "<< layers2D[0].elem_i.size()-1 <<endl;
//   cout<<"Number of faces in original mesh   "<< layers2D[0].face_i.size() <<endl;

    for(int i=0; i<layers2D.size(); i++)
    {
      cout<<"Number of vertices in original mesh layer "    << i << ": "<< layers2D[i].nodes.size() <<endl;
      cout<<"Number of elements in original mesh layer "    << i << ": "<< layers2D[i].elem_i.size()-1 <<endl;
      cout<<"Number of faces in original mesh layer "       << i << ": "<< layers2D[i].face_i.size() <<endl;

      if((i>0) && i<layers2D.size()-1)
        if(layers2D[i-1].nodes.size()- layers2D[i].nodes.size() != 0)
        {
          cout<<"Error, the number of vertices of each mesh is not the same"<<endl;
          throw(-1);
        }
    }




   deque< deque<POINT> > points2D;
   //deque<UNSTRUCTMESH> layers2D;

   // Copying the 2D mesh and the nodes that construct them
   for (int i=0; i<layers2D.size();i++)
   {
     //UNSTRUCTMESH tmp_mesh;
     deque<POINT> tmpPts;

     //tmp_mesh= unstMesh;
     //    tmp_mesh(unstMesh);
     for (int j=0;j<layers2D[0].nodes.size();j++)
     {
       //tmp_mesh.nodes[j].pt=unstMesh.nodes[j].pt+POINT(0.0,0.0,i*dz);        //making a copy of the 2D mesh at a certain dz up!!
       tmpPts.push_back(layers2D[i].nodes[j].pt);
     }
     points2D.push_back(tmpPts);
     //layers2D.push_back(tmp_mesh);


   }

   elem_v_2D = layers2D[0].elem_v;
   elem_i_2D = layers2D[0].elem_i;
   face_elem_2D = layers2D[0].face_elem;
   face_i_2D = layers2D[0].face_i;



   int nvert = (int)points2D[0].size();
   int nelem = (int)elem_i_2D.size()-1;
   int nface = (int)face_i_2D.size();

   // NODES OF THE 3D MESH

   // Cpying the nodes in the other direction, in this case z!
   deque< deque<POINT> > nodes3D;    //number of nodes in the 2D mesh x number of layers in z


//   for (int i=0; i<nvert; i++)
//   {
//     deque<POINT> nodes_i;
//
//     deque<POINT> splinepts;
//     for (int s=0; s<points2D.size(); s++)
//       splinepts.push_back(points2D[s][i]);
//
//     LINE span(splinepts);
//     // define distribution in the spanwise direction
//     double distr[nLayers];
//     //firstLastLengthDistr(dz/nLayers, dz/nLayers, nLayers, distr);      //First lenght and last lenght can be used for boundary layer in the spanwise direction
//
//     //-------------------------------------------------------
//     double z1=layers2D[1].nodes[i].pt.z;
//     double z0=layers2D[0].nodes[i].pt.z;
//     double span_distance=z1-z0; //layers2D[i].nodes[j].pt
//
//
//     for (int i=0; i<nLayers; i++)
//       distr[i] = (double)i/(double)(nLayers-1);
//     //firstLastLengthDistr(span_distance/nLayers, span_distance/nLayers, nLayers, distr);
//     //-------------------------------------------------------
//
//     for (int l=0; l<nLayers; l++)
//       nodes_i.push_back(span.calcPoint(distr[l]));
//
//     nodes3D.push_back(nodes_i);
//   }
   for (int i=0; i<nvert; i++)
    {
      deque<POINT> nodes_i;

      deque<POINT> splinepts;
      for (int s=0; s<points2D.size(); s++)
        splinepts.push_back(points2D[s][i]);


      if (interpScheme == "LINE")
      {
        LINE span(splinepts);
        // define distribution in the spanwise direction
        double distr[nLayers];
        for (int i=0; i<nLayers; i++)
          distr[i] = (double)i/(double)(nLayers-1);
  //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
        for (int l=0; l<nLayers; l++)
          nodes_i.push_back(span.calcPoint(distr[l]));
      }
      else if (interpScheme == "SPLINE")
      {
        SPLINE span(splinepts);
        // define distribution in the spanwise direction
        double distr[nLayers];
        for (int i=0; i<nLayers; i++)
          distr[i] = (double)i/(double)(nLayers-1);
  //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
        for (int l=0; l<nLayers; l++)
          nodes_i.push_back(span.calcPoint(distr[l]));
      }
      else if (interpScheme == "BEZIER")
      {
        BEZIER span(splinepts);
        // define distribution in the spanwise direction
        double distr[nLayers];
        for (int i=0; i<nLayers; i++)
          distr[i] = (double)i/(double)(nLayers-1);
  //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
        for (int l=0; l<nLayers; l++)
          nodes_i.push_back(span.calcPoint(distr[l]));
      }
      else
      {
        cout << "error: interpolation " << interpScheme << " not implemented" << endl;
        throw(-10);
      }

      nodes3D.push_back(nodes_i);


    }


//   for (int i=0; i<nvert; i++)
//   {
//     deque<POINT> nodes_i;
//
//     deque<POINT> splinepts;
//     for (int s=0; s<points2D.size(); s++)
//       splinepts.push_back(points2D[s][i]);
//
//
//     if (interpScheme == "LINE")
//     {
//       LINE span(splinepts);
//       // define distribution in the spanwise direction
//       double distr[nLayers];
//       for (int i=0; i<nLayers; i++)
//         distr[i] = (double)i/(double)(nLayers-1);
//       //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
//       for (int l=0; l<nLayers; l++)
//         nodes_i.push_back(span.calcPoint(distr[l]));
//     }
//     else if (interpScheme == "SPLINE")
//     {
//       SPLINE span(splinepts);
//       // define distribution in the spanwise direction
//       double distr[nLayers];
//       for (int i=0; i<nLayers; i++)
//         distr[i] = (double)i/(double)(nLayers-1);
//       //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
//       for (int l=0; l<nLayers; l++)
//         nodes_i.push_back(span.calcPoint(distr[l]));
//     }
//     else if (interpScheme == "BEZIER")
//     {
//       BEZIER span(splinepts);
//       // define distribution in the spanwise direction
//       double distr[nLayers];
//       for (int i=0; i<nLayers; i++)
//         distr[i] = (double)i/(double)(nLayers-1);
//       //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
//       for (int l=0; l<nLayers; l++)
//         nodes_i.push_back(span.calcPoint(distr[l]));
//     }
//     else
//     {
//       cout << "error: interpolation " << interpScheme << " not implemented" << endl;
//       throw(-10);
//     }
//
//     nodes3D.push_back(nodes_i);
//   }


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


   // FACES OF THE 3D MESH      With its connectivity as the *.msh files needs
/*          ORIGINAl WORKING LOOP, BUT SLOOOOOW!!!!!!!
   for (int l=0; l<nLayers-1; l++)
     for(int i=0; i<elem_i_2D.size()-1; i++)
     {
       FACE tmpface_x1;
       for (int j=elem_i_2D[i]; j<elem_i_2D[i+1]; j++)
         tmpface_x1.node.push_back(elem_v_2D[j]+(l*nvert));

       tmpface_x1.elem_1= (i+1)+(l*nelem);
       tmpface_x1.elem_2= 0;//l*(i+1+nelem);
       //Naming the boundary at z = zmin
       if(l==0)   strcpy (tmpface_x1.name,"bottom");
//         cout<<endl;
       mesh3D.faces.push_back(tmpface_x1);

       FACE tmpface_x2;
       for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)    tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));

       tmpface_x2.elem_1= (i+1)+(l*nelem);
//         if (l!=nLayers-2) tmpface_x2.elem_2= l*(i+1+nelem);
//         else tmpface_x2.elem_2= 0;
       tmpface_x2.elem_2= 0;
       //Naming the boundary at z = zmax
       if(l==nLayers-2)   strcpy (tmpface_x2.name,"top");
       mesh3D.faces.push_back(tmpface_x2);

       for(int j=elem_i_2D[i]; j<elem_i_2D[i+1]; j++)
       {
         FACE tmpface_y;
         if(elem_i_2D[i+1]-j>1)
         {
           tmpface_y.node.push_back(elem_v_2D[j+1]+(l*nvert));
           tmpface_y.node.push_back(elem_v_2D[j]+(l*nvert));
           tmpface_y.node.push_back(elem_v_2D[j]+((l+1)*nvert));
           tmpface_y.node.push_back(elem_v_2D[j+1]+((l+1)*nvert));
         }
         else
         {
           tmpface_y.node.push_back(elem_v_2D[elem_i_2D[i]]+(l*nvert));
           tmpface_y.node.push_back(elem_v_2D[elem_i_2D[i+1]-1]+(l*nvert));
           tmpface_y.node.push_back(elem_v_2D[elem_i_2D[i+1]-1]+((l+1)*nvert));
           tmpface_y.node.push_back(elem_v_2D[elem_i_2D[i]]+((l+1)*nvert));
         }
         tmpface_y.elem_1= (i+1)+(l*nelem);
         tmpface_y.elem_2= 0;

         mesh3D.faces.push_back(tmpface_y);
       }
//         cout<<"In side face calculation" <<endl;

     }




     */

   ///WIth the with list
//     if(layers2D[0].face_i.size()!=0)
//     {
//       int count=0;
//       for (int l=0; l<nLayers-1; l++)
//       {
//         //--------------------------------------------------------------------------
//         // Faces perpendicular to the original 2D mesh
//         for(int i=0; i<layers2D[0].face_i.size(); i++)
//         {
////           FACE tmpface_y;
//
//           //Giving the nodes
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]+1]+(l*nvert));
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]]+(l*nvert));
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]]+((l+1)*nvert));
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]+1]+((l+1)*nvert));
//
//           //Index
//           mesh3D.face_i.push_back(4*count);//(4*i+l*nvert);
//           count++;
//
//           //Giving the elements
//           mesh3D.face_elem.push_back(face_elem_2D[2*i]+(l*nelem));
//
//           if(layers2D[0].face_elem[2*i+1]==0)    mesh3D.face_elem.push_back(0);
//           else                                   mesh3D.face_elem.push_back(face_elem_2D[2*i+1]+(l*nelem));
//
//           //Giving the name of the face as the one in the original 2D mesh
//           mesh3D.face_name.push_back(layers2D[0].face_name[i]);
//
//         }
//         //--------------------------------------------------------------------------
//         // Faces parallel to the original 2D mesh
//         for(int i=0; i<elem_i_2D.size()-1; i++)
//         {
//           //Giving the nodes
//           for (int j=elem_i_2D[i]; j<elem_i_2D[i+1]; j++)
//             mesh3D.face_v.push_back(elem_v_2D[j]+(l*nvert));//tmpface_x1.node.push_back(elem_v_2D[j]+(l*nvert));
//
//           //Index
//           mesh3D.face_i.push_back(4*count);//(4*i+l*nvert);
//           count++;
//
//           //Giving the elements and the names!
//           if(l==0)                                 //Faces of the bottom layer
//           {
//             mesh3D.face_elem.push_back((i+1)+(l*nelem));       //tmpface_x1.elem_1= (i+1)+(l*nelem);
//             mesh3D.face_elem.push_back(0);                     //tmpface_x1.elem_2= 0;//l*(i+1+nelem);
//             mesh3D.face_name.push_back("bottom");              //strcpy (tmpface_x1.name,"bottom");
//           }
//           else                                     //Faces In-between two layers
//           {
//             mesh3D.face_elem.push_back((i+1)+(l*nelem)); //tmpface_x1.elem_1= (i+1)+((l)*nelem); //l*(i+1+nelem);
//             mesh3D.face_elem.push_back((i+1)+((l-1)*nelem)); //tmpface_x1.elem_2= (i+1)+((l-1)*nelem);
//             mesh3D.face_name.push_back("fluid"); //strcpy (tmpface_x1.name,"fluid");
//           }
////           mesh3D.faces.push_back(tmpface_x1);
//           //--------------------------------------------------------------------------
//           if(l==nLayers-2)         //Faces of the top layer
//           {
////             FACE tmpface_x2;
//             //Giving the nodes
//             for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)
//               mesh3D.face_v.push_back(elem_v_2D[j]+((l+1)*nvert)); //tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));
//
//             //Index
//             mesh3D.face_i.push_back(4*count);//(4*i+l*nvert);
//             count++;
//
//             //Giving the elements and the names!
//             mesh3D.face_elem.push_back((i+1)+(l*nelem));       //tmpface_x2.elem_1= (i+1)+(l*nelem);
//             mesh3D.face_elem.push_back(0);                     //tmpface_x2.elem_2= 0;
//             mesh3D.face_name.push_back("top");              //strcpy (tmpface_x2.name,"top");
//
////             mesh3D.faces.push_back(tmpface_x2);
//           }
//         }
//       }
//
//
//
//       //     mesh3D.removeDoubleFaces();
//
//     }
   //Faces with Objects!
   if(layers2D[0].faces.size()!=0)
   {
     for (int l=0; l<nLayers-1; l++)
     {
       //--------------------------------------------------------------------------
       // Faces parallel to the original 2D mesh
       for(int i=0; i<layers2D[0].faces.size(); i++)
       {
         FACE tmpface_y;

         //Giving the elements
         tmpface_y.node.push_back(layers2D[0].faces[i].node[0]+(l*nvert));
         tmpface_y.node.push_back(layers2D[0].faces[i].node[1]+(l*nvert));
         tmpface_y.node.push_back(layers2D[0].faces[i].node[1]+((l+1)*nvert));
         tmpface_y.node.push_back(layers2D[0].faces[i].node[0]+((l+1)*nvert));

         //Giving the elements
         tmpface_y.elem_1= layers2D[0].faces[i].elem_1+(l*nelem);

         if(layers2D[0].faces[i].elem_2==0)     tmpface_y.elem_2= 0;
         else                                   tmpface_y.elem_2= layers2D[0].faces[i].elem_2+(l*nelem);

         //Giving the name of the face as the one in the original 2D mesh
         strcpy (tmpface_y.name,layers2D[0].faces[i].name);


         mesh3D.faces.push_back(tmpface_y);

       }
       //--------------------------------------------------------------------------
       // Faces parallel to the original 2D mesh
       for(int i=0; i<elem_i_2D.size()-1; i++)
       {
         FACE tmpface_x1;
         for (int j=elem_i_2D[i]; j<elem_i_2D[i+1]; j++)
           tmpface_x1.node.push_back(elem_v_2D[j]+(l*nvert));

         if(l==0)                                 //Faces of the bottom layer
         {
           tmpface_x1.elem_1= (i+1)+(l*nelem);
           tmpface_x1.elem_2= 0;//l*(i+1+nelem);
           strcpy (tmpface_x1.name,"bottom");
         }
         else                                     //Faces In-between two layers
         {
           tmpface_x1.elem_1= (i+1)+((l)*nelem); //l*(i+1+nelem);
           tmpface_x1.elem_2= (i+1)+((l-1)*nelem);
           strcpy (tmpface_x1.name,"fluid");
         }
         mesh3D.faces.push_back(tmpface_x1);
         //--------------------------------------------------------------------------
         if(l==nLayers-2)         //Faces of the top layer
         {
           FACE tmpface_x2;
           for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)    tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));

           tmpface_x2.elem_1= (i+1)+(l*nelem);
           tmpface_x2.elem_2= 0;
           strcpy (tmpface_x2.name,"top");

           mesh3D.faces.push_back(tmpface_x2);
         }

         // THe following piece of code will create face which are doubled!!!
         //       FACE tmpface_x2;
         //       for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)    tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));
         //       tmpface_x2.elem_1= (i+1)+(l*nelem);
         //       if(l==nLayers-2)
         //       {
         //         tmpface_x2.elem_2= 0;
         //         strcpy (tmpface_x2.name,"top");
         //       }
         //       else
         //       {
         //         tmpface_x2.elem_2= (i+1)+((l-1)*nelem); //l*(i+1+nelem);
         //         strcpy (tmpface_x2.name,"fluid");
         //       }
         //       mesh3D.faces.push_back(tmpface_x2);

         //--------------------------------------------------------------------------
         //         cout<<"In side face calculation" <<endl;
       }
     }
   }
   else
     cout<<"Inside: buildMesh3DNew(): The 2D faces objects of the mesh were not already created! Therefore the face connectivity for the *.msh was not created" <<endl;


   return(mesh3D);
 }


UNSTRUCTMESH MESHTOOLS::buildMesh3DNewwithFace_obj(deque<UNSTRUCTMESH> &layers2D, const int nLayers, const double firstLength, const double lastLength, const string &interpScheme)    //GUSTAVO 15/06/16
 {
  /* Made by Gustavo J. Otero Date: 01-04-2016
     * For this function generates a 3D mesh from two or more meshes (they do not need to be 2D_. That is why the input is a deque of meshes.
     *  It includes the generatation of the faces, later use for the creating of the Fluent *msh file
     *
     *
     */
   UNSTRUCTMESH mesh3D;
   deque<int> elem_v_2D, elem_i_2D;
   deque<int> face_elem_2D, face_i_2D;

//   cout<<"Number of vertices in original mesh   "<< layers2D[0].nodes.size() <<endl;
//   cout<<"Number of elements in original mesh   "<< layers2D[0].elem_i.size()-1 <<endl;
//   cout<<"Number of faces in original mesh   "<< layers2D[0].face_i.size() <<endl;

    for(int i=0; i<layers2D.size(); i++)
    {
      cout<<"Number of vertices in original mesh layer "    << i << ": "<< layers2D[i].nodes.size() <<endl;
      cout<<"Number of elements in original mesh layer "    << i << ": "<< layers2D[i].elem_i.size()-1 <<endl;
      cout<<"Number of faces in original mesh layer "       << i << ": "<< layers2D[i].face_i.size() <<endl;

      if((i>0) && i<layers2D.size()-1)
        if(layers2D[i-1].nodes.size()- layers2D[i].nodes.size() != 0)
        {
          cout<<"Error, the number of vertices of each mesh is not the same"<<endl;
          throw(-1);
        }
    }




   deque< deque<POINT> > points2D;
   //deque<UNSTRUCTMESH> layers2D;

   // Copying the 2D mesh and the nodes that construct them
   for (int i=0; i<layers2D.size();i++)
   {
     //UNSTRUCTMESH tmp_mesh;
     deque<POINT> tmpPts;

     //tmp_mesh= unstMesh;
     //    tmp_mesh(unstMesh);
     for (int j=0;j<layers2D[0].nodes.size();j++)
     {
       //tmp_mesh.nodes[j].pt=unstMesh.nodes[j].pt+POINT(0.0,0.0,i*dz);        //making a copy of the 2D mesh at a certain dz up!!
       tmpPts.push_back(layers2D[i].nodes[j].pt);
     }
     points2D.push_back(tmpPts);
     //layers2D.push_back(tmp_mesh);


   }

   elem_v_2D = layers2D[0].elem_v;
   elem_i_2D = layers2D[0].elem_i;
   face_elem_2D = layers2D[0].face_elem;
   face_i_2D = layers2D[0].face_i;



   int nvert = (int)points2D[0].size();
   int nelem = (int)elem_i_2D.size()-1;
   int nface = (int)face_i_2D.size();

   // NODES OF THE 3D MESH

   // Cpying the nodes in the other direction, in this case z!
   deque< deque<POINT> > nodes3D;    //number of nodes in the 2D mesh x number of layers in z

   //Distribution along the span
   double distr[nLayers];
   for (int i=0; i<nLayers; i++)
     distr[i] = (double)i/(double)(nLayers-1);
   //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);


   for (int i=0; i<nvert; i++)
   {
     deque<POINT> nodes_i;

     deque<POINT> splinepts;
     for (int s=0; s<points2D.size(); s++)
       splinepts.push_back(points2D[s][i]);


     if (interpScheme == "LINE")
     {
       LINE span(splinepts);
       // define distribution in the spanwise direction
       double distr[nLayers];
       for (int i=0; i<nLayers; i++)
         distr[i] = (double)i/(double)(nLayers-1);
       //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
       for (int l=0; l<nLayers; l++)
         nodes_i.push_back(span.calcPoint(distr[l]));
     }
     else if (interpScheme == "SPLINE")
     {
       SPLINE span(splinepts);
       // define distribution in the spanwise direction
       double distr[nLayers];
       for (int i=0; i<nLayers; i++)
         distr[i] = (double)i/(double)(nLayers-1);
       //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
       for (int l=0; l<nLayers; l++)
         nodes_i.push_back(span.calcPoint(distr[l]));
     }
     else if (interpScheme == "BEZIER")
     {
       BEZIER span(splinepts);
       // define distribution in the spanwise direction
       double distr[nLayers];
       for (int i=0; i<nLayers; i++)
         distr[i] = (double)i/(double)(nLayers-1);
       //      firstLastLengthDistr(firstLength/span.length, lastLength/span.length, nLayers, distr);
       for (int l=0; l<nLayers; l++)
         nodes_i.push_back(span.calcPoint(distr[l]));
     }
     else
     {
       cout << "error: interpolation " << interpScheme << " not implemented" << endl;
       throw(-10);
     }

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


   // FACES OF THE 3D MESH      With its connectivity as the *.msh files needs


   if(layers2D[0].faces.size()!=0)
     for (int l=0; l<nLayers-1; l++)
     {
       //--------------------------------------------------------------------------
       // Faces parallel to the original 2D mesh
       for(int i=0; i<layers2D[0].faces.size(); i++)
       {
         FACE tmpface_y;

         //Giving the elements
         tmpface_y.node.push_back(layers2D[0].faces[i].node[0]+(l*nvert));
         tmpface_y.node.push_back(layers2D[0].faces[i].node[1]+(l*nvert));
         tmpface_y.node.push_back(layers2D[0].faces[i].node[1]+((l+1)*nvert));
         tmpface_y.node.push_back(layers2D[0].faces[i].node[0]+((l+1)*nvert));

         //Giving the elements
         tmpface_y.elem_1= layers2D[0].faces[i].elem_1+(l*nelem);

         if(layers2D[0].faces[i].elem_2==0)     tmpface_y.elem_2= 0;
         else                                   tmpface_y.elem_2= layers2D[0].faces[i].elem_2+(l*nelem);

         //Giving the name of the face as the one in the original 2D mesh
         strcpy (tmpface_y.name,layers2D[0].faces[i].name);


         mesh3D.faces.push_back(tmpface_y);

       }
       //--------------------------------------------------------------------------
       // Faces parallel to the original 2D mesh
       for(int i=0; i<elem_i_2D.size()-1; i++)
       {
         FACE tmpface_x1;
         //Giving the nodes
         for (int j=elem_i_2D[i]; j<elem_i_2D[i+1]; j++)
           tmpface_x1.node.push_back(elem_v_2D[j]+(l*nvert));

         //Giving the elements and the names!
         if(l==0)                                 //Faces of the bottom layer
         {
           tmpface_x1.elem_1= (i+1)+(l*nelem);
           tmpface_x1.elem_2= 0;//l*(i+1+nelem);
           strcpy (tmpface_x1.name,"bottom");
         }
         else                                     //Faces In-between two layers
         {
           tmpface_x1.elem_1= (i+1)+((l)*nelem); //l*(i+1+nelem);
           tmpface_x1.elem_2= (i+1)+((l-1)*nelem);
           strcpy (tmpface_x1.name,"fluid");
         }
         mesh3D.faces.push_back(tmpface_x1);
         //--------------------------------------------------------------------------
         if(l==nLayers-2)         //Faces of the top layer
         {
           FACE tmpface_x2;
           //Giving the nodes
           for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)
             tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));

           tmpface_x2.elem_1= (i+1)+(l*nelem);
           tmpface_x2.elem_2= 0;
           strcpy (tmpface_x2.name,"top");

           mesh3D.faces.push_back(tmpface_x2);
         }

         // THe following piece of code will create face which are doubled!!!
         //       FACE tmpface_x2;
         //       for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)    tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));
         //       tmpface_x2.elem_1= (i+1)+(l*nelem);
         //       if(l==nLayers-2)
         //       {
         //         tmpface_x2.elem_2= 0;
         //         strcpy (tmpface_x2.name,"top");
         //       }
         //       else
         //       {
         //         tmpface_x2.elem_2= (i+1)+((l-1)*nelem); //l*(i+1+nelem);
         //         strcpy (tmpface_x2.name,"fluid");
         //       }
         //       mesh3D.faces.push_back(tmpface_x2);

         //--------------------------------------------------------------------------
         //         cout<<"In side face calculation" <<endl;
       }
     }
   else
     cout<<"Inside: buildMesh3DNew(): The 2D faces objects of the mesh were not already created! Therefore the face connectivity for the *.msh was not created" <<endl;



//     ///WIth the newwwwwww
//     if(layers2D[0].face_i.size()!=0)
//     {
//       int count=0;
//       for (int l=0; l<nLayers-1; l++)
//       {
//         //--------------------------------------------------------------------------
//         // Faces perpendicular to the original 2D mesh
//         for(int i=0; i<layers2D[0].face_i.size(); i++)
//         {
//           //           FACE tmpface_y;
//
//           //Giving the nodes
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]+1]+(l*nvert));
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]]+(l*nvert));
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]]+((l+1)*nvert));
//           mesh3D.face_v.push_back(layers2D[0].face_v[face_i_2D[i]+1]+((l+1)*nvert));
//
//           //Index
//           mesh3D.face_i.push_back(4*count);//(4*i+l*nvert);
//           count++;
//
//           //Giving the elements
//           mesh3D.face_elem.push_back(face_elem_2D[2*i]+(l*nelem));
//
//           if(layers2D[0].face_elem[2*i+1]==0)    mesh3D.face_elem.push_back(0);
//           else                                   mesh3D.face_elem.push_back(face_elem_2D[2*i+1]+(l*nelem));
//
//           //Giving the name of the face as the one in the original 2D mesh
//           mesh3D.face_name.push_back(layers2D[0].face_name[i]);
//
//         }
//         //--------------------------------------------------------------------------
//         // Faces parallel to the original 2D mesh
//         for(int i=0; i<elem_i_2D.size()-1; i++)
//         {
//           //Giving the nodes
//           for (int j=elem_i_2D[i]; j<elem_i_2D[i+1]; j++)
//             mesh3D.face_v.push_back(elem_v_2D[j]+(l*nvert));//tmpface_x1.node.push_back(elem_v_2D[j]+(l*nvert));
//
//           //Index
//           mesh3D.face_i.push_back(4*count);//(4*i+l*nvert);
//           count++;
//
//           //Giving the elements and the names!
//           if(l==0)                                 //Faces of the bottom layer
//           {
//             mesh3D.face_elem.push_back((i+1)+(l*nelem));       //tmpface_x1.elem_1= (i+1)+(l*nelem);
//             mesh3D.face_elem.push_back(0);                     //tmpface_x1.elem_2= 0;//l*(i+1+nelem);
//             mesh3D.face_name.push_back("bottom");              //strcpy (tmpface_x1.name,"bottom");
//           }
//           else                                     //Faces In-between two layers
//           {
//             mesh3D.face_elem.push_back((i+1)+(l*nelem)); //tmpface_x1.elem_1= (i+1)+((l)*nelem); //l*(i+1+nelem);
//             mesh3D.face_elem.push_back((i+1)+((l-1)*nelem)); //tmpface_x1.elem_2= (i+1)+((l-1)*nelem);
//             mesh3D.face_name.push_back("fluid"); //strcpy (tmpface_x1.name,"fluid");
//           }
//           //           mesh3D.faces.push_back(tmpface_x1);
//           //--------------------------------------------------------------------------
//           if(l==nLayers-2)         //Faces of the top layer
//           {
//             //             FACE tmpface_x2;
//             //Giving the nodes
//             for (int j=elem_i_2D[i+1]-1; j>elem_i_2D[i]-1; j--)
//               mesh3D.face_v.push_back(elem_v_2D[j]+((l+1)*nvert)); //tmpface_x2.node.push_back(elem_v_2D[j]+((l+1)*nvert));
//
//             //Index
//             mesh3D.face_i.push_back(4*count);//(4*i+l*nvert);
//             count++;
//
//             //Giving the elements and the names!
//             mesh3D.face_elem.push_back((i+1)+(l*nelem));       //tmpface_x2.elem_1= (i+1)+(l*nelem);
//             mesh3D.face_elem.push_back(0);                     //tmpface_x2.elem_2= 0;
//             mesh3D.face_name.push_back("top");              //strcpy (tmpface_x2.name,"top");
//
//             //             mesh3D.faces.push_back(tmpface_x2);
//           }
//         }
//       }
//
//
//
//       //     mesh3D.removeDoubleFaces();
//
//     }
//     else
//       cout<<"Inside: buildMesh3DNew(): The 2D faces objects of the mesh were not already created! Therefore the face connectivity for the *.msh was not created" <<endl;


     return(mesh3D);
   }

void MESHTOOLS::writeFluentMsh_Interface(const char *filename, const UNSTRUCTMESH meshL, const UNSTRUCTMESH meshR ,  deque<string> bc_names)
{
  FILE * fp;

  printf(" > starting to write fluent case file with an INTERFACE: %s\n", filename);

  if ( (fp=fopen(filename,"w"))==NULL ) {
    printf(" > could not open file %s\n", filename);
    exit(-1);
  }

  fprintf(fp,"(0 \"Project to Fluent\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(1 \"In-house Mesh Generator\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(0 \"Dimensions:\")\n");
  fprintf(fp,"(2 3)\n");
  fprintf(fp,"\n");

  //
  // header and global dimensions
  //

  int verticesR = meshR.nodes.size();
  int elementsR = meshR.elem_i.size()-1;
  int nfacesR = meshR.face_i.size();
  int nzoneR = meshR.n_zone;


  int verticesL = meshL.nodes.size();
  int elementsL = meshL.elem_i.size()-1;
  int nfacesL = meshL.face_i.size();
  int nzoneL = meshL.n_zone;

  int vertices =verticesL + verticesR;
  int elements = elementsL + elementsR;
  int nfaces = nfacesL + nfacesR;
  int nzone =nzoneL + nzoneL - 3;   // top, bottom and fluid are already in both meshes!

  fprintf(fp,"(0 \"Grid:\")\n");
//  fprintf(fp,"\n");
//  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

  fprintf(fp,"\n");

  fprintf(fp,"(0 \"Nodes:\")\n");
  fprintf(fp,"(10 (0 %x %x 1 3))\n",1, vertices);
  fprintf(fp, "(10 (1 %x %x 1 3)(\n",1, vertices);
  for (int i=0; i<meshL.nodes.size(); i++)
    fprintf(fp, "%20.11le%20.11le%20.11le\n", meshL.nodes[i].pt.x, meshL.nodes[i].pt.y, meshL.nodes[i].pt.z);
  for (int i=0; i<meshR.nodes.size(); i++)
      fprintf(fp, "%20.11le%20.11le%20.11le\n", meshR.nodes[i].pt.x, meshR.nodes[i].pt.y, meshR.nodes[i].pt.z);
  fprintf(fp, "))\n");
  fprintf(fp,"\n");


  //
  // start of faces
  //
  int fa_element_type = 4; // hexahedral = 4
  int fa_type = 0; //1; // active = 1 - inactive = 32
//  nfaces=106;
  fprintf(fp,"(0 \"Faces:\")\n");
  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

  int ind_faces[nzone+1];      //Index of the current zone

  bool zone_noname=false;
  for(int i=0; i<nfacesL; i++)
    if(meshL.face_name[i].compare("noname")==0) //if(strcmp (meshL.face_name[i], "noname") == 0)
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces in the Left side was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }
  for(int i=0; i<nfacesR; i++)
    if(meshR.face_name[i].compare("noname")==0) //if(strcmp (meshR.faces[i].name, "noname") == 0)
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces in the Right side was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }


  if(zone_noname)
  {

    int faces_wall=0;
    for(int i=0; i<nfacesL; i++)
      if(meshL.face_elem[(2*i)+1]==0)  faces_wall++;
    for(int i=0; i<nfacesR; i++)
      if(meshR.face_elem[(2*i)+1]==0)  faces_wall++;

    //WALL FACES
    ind_faces[0]=3;

    fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[0] ,1, faces_wall);
    for(int i=0; i<nfacesL; i++)
      if(meshL.face_elem[(2*i)+1]==0)  fprintf(fp,"4 %x %x %x %x %x %x \n",  meshL.face_v[meshL.face_i[i]]+1 ,           meshL.face_v[meshL.face_i[i]+1]+1 ,               meshL.face_v[meshL.face_i[i]+2]+1 ,             meshL.face_v[meshL.face_i[i]+3]+1 ,               meshL.face_elem[(2*i)],             meshL.face_elem[(2*i)+1]);
    for(int i=0; i<nfacesR; i++)
      if(meshR.face_elem[(2*i)+1]==0)  fprintf(fp,"4 %x %x %x %x %x %x \n",  meshR.face_v[meshR.face_i[i]]+1+verticesL , meshR.face_v[meshR.face_i[i]+1]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]+elementsL);


    fprintf(fp, "))\n");

    //INNER FACES
    ind_faces[1]=5;
    fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[1] ,faces_wall+1, nfaces);
    for(int i=0; i<nfacesL; i++)
      if(meshL.face_elem[(2*i)+1]>0)   fprintf(fp,"4 %x %x %x %x %x %x \n",  meshL.face_v[meshL.face_i[i]]+1 ,           meshL.face_v[meshL.face_i[i]+1]+1 ,               meshL.face_v[meshL.face_i[i]+2]+1 ,             meshL.face_v[meshL.face_i[i]+3]+1 ,               meshL.face_elem[(2*i)],             meshL.face_elem[(2*i)+1]);
    for(int i=0; i<nfacesR; i++)
      if(meshR.face_elem[(2*i)+1]>0)   fprintf(fp,"4 %x %x %x %x %x %x \n",  meshR.face_v[meshR.face_i[i]]+1+verticesL , meshR.face_v[meshR.face_i[i]+1]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]+elementsL);



    fprintf(fp, "))\n");
  }
  else
  {
    int nfa_i[nzone];        //Number of faces in the current zone
    int total_faces=0;

    for(int n=0 ; n<bc_names.size(); n++)
    {
      // Defining the index of the face zone
      ind_faces[n+1]=n+4;                                 // the index starts in 3
      //Counting the number of faces of the present face zone, needed in the declaration for the face
      nfa_i[n]=0;
      for(int i=0; i<nfacesL; i++)
        if(meshL.face_name[i].compare(bc_names[n])==0)
          nfa_i[n]++;
      for(int i=0; i<nfacesR; i++)
        if(meshR.face_name[i].compare(bc_names[n])==0)
          nfa_i[n]++;


      // Declaration of the face zone
      if(n!=bc_names.size()-1) fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[n+1] ,total_faces+1, total_faces+nfa_i[n]);
      else fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[n+1]+1 ,total_faces+1, total_faces+nfa_i[n]);
      // Giving the face connectivity
      for(int i=0; i<nfacesL; i++)
        if(meshL.face_name[i].compare(bc_names[n])==0)
          fprintf(fp,"4 %x %x %x %x %x %x \n",                                  meshL.face_v[meshL.face_i[i]]+1 ,             meshL.face_v[meshL.face_i[i]+1]+1 ,               meshL.face_v[meshL.face_i[i]+2]+1 ,               meshL.face_v[meshL.face_i[i]+3]+1 ,               meshL.face_elem[(2*i)],             meshL.face_elem[(2*i)+1]);

      for(int i=0; i<nfacesR; i++)
        if(meshR.face_name[i].compare(bc_names[n])==0)
          if (meshR.face_elem[(2*i)+1]==0) fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.face_v[meshR.face_i[i]]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+1]+1 +verticesL ,    meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]);
          else                          fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.face_v[meshR.face_i[i]]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+1]+1 +verticesL ,    meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]+elementsL);

      fprintf(fp, "))\n");
      total_faces+=nfa_i[n];
      cout<<"Number of faces in boundary "<<bc_names[n] <<" "<<nfa_i[n]<<endl;

    }

    cout<<"Total number of faces in the numerical domain "<<total_faces<<endl;


  }


  fprintf(fp,"\n");

  //
  // cells
  //
  int cv_element_type = 4; // hexahedral = 4
  int cv_type = 1; //1; // active = 1 - inactive = 32
  char this_zone[128];
  sprintf(this_zone,"fluid");
  fprintf(fp,"(0 \"Cells:\")\n");
  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 3,elementsR+1,elements,cv_type,cv_element_type);
  fprintf(fp,"\n");

  if(zone_noname)
  {
    fprintf(fp,"(0 \"Zones:\")\n");
    fprintf(fp,"(45 (2 fluid fluid))\n");
//    fprintf(fp,"(45 (3 fluid fluid 0))\n");
    fprintf(fp,"(45 (%i wall boundaries))\n", ind_faces[0]);
    fprintf(fp,"(45 (%i interior default-interior)\n", ind_faces[1]);
  }
  else
  {
    fprintf(fp,"(0 \"Zones:\")\n");
    fprintf(fp,"(45 (2 fluid fluid)())\n");
//    fprintf(fp,"(45 (3 fluid fluid 0))\n");
    for(int n=0 ; n<bc_names.size()-1; n++)
    {
      fprintf(fp,"(45 (%i wall %s)())\n", ind_faces[n+1], bc_names[n].c_str());
    }
    fprintf(fp,"(45 (%i interior default-interior)())\n", ind_faces[bc_names.size()]+1);
  }
//  fprintf(fp,"(12 (%x %x %x %x))\n", 0,1,elements,cv_type);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"(45 (%d fluid %s)())\n",  4001,this_zone);
  fprintf(fp,"\n");



  fclose(fp);

  printf(" > Fluent case file written for INTERFACE: %s\n", filename);
}

void MESHTOOLS::writeFluentMsh_Interface_obj(const char *filename, const UNSTRUCTMESH meshL, const UNSTRUCTMESH meshR ,  deque<string> bc_names)
{
  FILE * fp;

  printf(" > starting to write fluent case file with an INTERFACE: %s\n", filename);

  if ( (fp=fopen(filename,"w"))==NULL ) {
    printf(" > could not open file %s\n", filename);
    exit(-1);
  }

  fprintf(fp,"(0 \"Project to Fluent\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(1 \"In-house Mesh Generator\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(0 \"Dimensions:\")\n");
  fprintf(fp,"(2 3)\n");
  fprintf(fp,"\n");

  //
  // header and global dimensions
  //

  int verticesR = meshR.nodes.size();
  int elementsR = meshR.elem_i.size()-1;
  int nfacesR = meshR.faces.size();         //face_i.size();
  int nzoneR = meshR.n_zone;


  int verticesL = meshL.nodes.size();
  int elementsL = meshL.elem_i.size()-1;
  int nfacesL = meshL.faces.size();         //face_i.size();
  int nzoneL = meshL.n_zone;

  int vertices =verticesL + verticesR;
  int elements = elementsL + elementsR;
  int nfaces = nfacesL + nfacesR;
  int nzone =nzoneL + nzoneL - 3;   // top, bottom and fluid are already in both meshes!

  fprintf(fp,"(0 \"Grid:\")\n");
//  fprintf(fp,"\n");
//  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

  fprintf(fp,"\n");

  fprintf(fp,"(0 \"Nodes:\")\n");
  fprintf(fp,"(10 (0 %x %x 1 3))\n",1, vertices);
  fprintf(fp, "(10 (1 %x %x 1 3)(\n",1, vertices);
  for (int i=0; i<meshL.nodes.size(); i++)
    fprintf(fp, "%20.11le%20.11le%20.11le\n", meshL.nodes[i].pt.x, meshL.nodes[i].pt.y, meshL.nodes[i].pt.z);
  for (int i=0; i<meshR.nodes.size(); i++)
      fprintf(fp, "%20.11le%20.11le%20.11le\n", meshR.nodes[i].pt.x, meshR.nodes[i].pt.y, meshR.nodes[i].pt.z);
  fprintf(fp, "))\n");
  fprintf(fp,"\n");


  //
  // start of faces
  //
  int fa_element_type = 4; // hexahedral = 4
  int fa_type = 0; //1; // active = 1 - inactive = 32
//  nfaces=106;
  fprintf(fp,"(0 \"Faces:\")\n");
  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

  int ind_faces[nzone+1];      //Index of the current zone

  bool zone_noname=false;
  for(int i=0; i<nfacesL; i++)
    if(strcmp (meshL.faces[i].name, "noname") == 0)   //(meshL.face_name[i].compare("noname")==0) //if(strcmp (meshL.face_name[i], "noname") == 0)
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces in the Left side was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }
  for(int i=0; i<nfacesR; i++)
    if(strcmp (meshR.faces[i].name, "noname") == 0)   //(meshR.face_name[i].compare("noname")==0) //if(strcmp (meshR.faces[i].name, "noname") == 0)
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces in the Right side was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }
  cout<<"Segmentation before"<<endl;
  if(zone_noname)
  {

    int faces_wall=0;
    for(int i=0; i<nfacesL; i++)
      if(meshL.faces[i].elem_2==0) faces_wall++;   //(meshL.face_elem[(2*i)+1]==0)  faces_wall++;
    for(int i=0; i<nfacesR; i++)
      if(meshR.faces[i].elem_2==0) faces_wall++;   //(meshR.face_elem[(2*i)+1]==0)  faces_wall++;

    //WALL FACES
    ind_faces[0]=3;

    fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[0] ,1, faces_wall);
    for(int i=0; i<nfacesL; i++)
      if(meshL.faces[i].elem_2==0)           fprintf(fp,"4 %x %x %x %x %x %x \n",    meshL.faces[i].node[0]+1,              meshL.faces[i].node[1]+1,           meshL.faces[i].node[2]+1,           meshL.faces[i].node[3]+1,           meshL.faces[i].elem_1,              meshL.faces[i].elem_2);
                    //if(meshL.face_elem[(2*i)+1]==0)  fprintf(fp,"4 %x %x %x %x %x %x \n",  meshL.face_v[meshL.face_i[i]]+1 ,           meshL.face_v[meshL.face_i[i]+1]+1 ,               meshL.face_v[meshL.face_i[i]+2]+1 ,             meshL.face_v[meshL.face_i[i]+3]+1 ,               meshL.face_elem[(2*i)],             meshL.face_elem[(2*i)+1]);
    for(int i=0; i<nfacesR; i++)
      if(meshR.faces[i].elem_2==0)           fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.faces[i].node[0]+1+verticesL,    meshR.faces[i].node[1]+1+verticesL, meshR.faces[i].node[2]+1+verticesL, meshR.faces[i].node[3]+1+verticesL, meshR.faces[i].elem_1+elementsL,    meshR.faces[i].elem_2);
                    //if(meshR.face_elem[(2*i)+1]==0)  fprintf(fp,"4 %x %x %x %x %x %x \n",  meshR.face_v[meshR.face_i[i]]+1+verticesL , meshR.face_v[meshR.face_i[i]+1]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]+elementsL);


    fprintf(fp, "))\n");

    //INNER FACES
    ind_faces[1]=5;
    fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[1] ,faces_wall+1, nfaces);
    for(int i=0; i<nfacesL; i++)
      if(meshL.faces[i].elem_2>0)           fprintf(fp,"4 %x %x %x %x %x %x \n",    meshL.faces[i].node[0]+1,              meshL.faces[i].node[1]+1,           meshL.faces[i].node[2]+1,           meshL.faces[i].node[3]+1,           meshL.faces[i].elem_1,              meshL.faces[i].elem_2);
                          //if(meshL.face_elem[(2*i)+1]>0)   fprintf(fp,"4 %x %x %x %x %x %x \n",  meshL.face_v[meshL.face_i[i]]+1 ,           meshL.face_v[meshL.face_i[i]+1]+1 ,               meshL.face_v[meshL.face_i[i]+2]+1 ,             meshL.face_v[meshL.face_i[i]+3]+1 ,               meshL.face_elem[(2*i)],             meshL.face_elem[(2*i)+1]);
    for(int i=0; i<nfacesR; i++)
      if(meshR.faces[i].elem_2>0)           fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.faces[i].node[0]+1+verticesL,    meshR.faces[i].node[1]+1+verticesL, meshR.faces[i].node[2]+1+verticesL, meshR.faces[i].node[3]+1+verticesL, meshR.faces[i].elem_1+elementsL,    meshR.faces[i].elem_2+elementsL);
                          //if(meshR.face_elem[(2*i)+1]>0)   fprintf(fp,"4 %x %x %x %x %x %x \n",  meshR.face_v[meshR.face_i[i]]+1+verticesL , meshR.face_v[meshR.face_i[i]+1]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]+elementsL);



    fprintf(fp, "))\n");
  }
  else
  {
    int nfa_i[nzone];        //Number of faces in the current zone
    int total_faces=0;


    for(int n=0 ; n<bc_names.size(); n++)
    {
      // Defining the index of the face zone
      ind_faces[n+1]=n+4;                                 // the index starts in 3
      //Counting the number of faces of the present face zone, needed in the declaration for the face
      nfa_i[n]=0;
      for(int i=0; i<nfacesL; i++)
        if(strcmp (meshL.faces[i].name, bc_names[n].c_str()) == 0)    //((meshL.face_name[i].compare(bc_names[n])==0)
          nfa_i[n]++;
      for(int i=0; i<nfacesR; i++)
        if(strcmp (meshR.faces[i].name, bc_names[n].c_str()) == 0)    //(meshR.face_name[i].compare(bc_names[n])==0)
          nfa_i[n]++;
;
      // Declaration of the face zone
      if(n!=bc_names.size()-1) fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[n+1] ,total_faces+1, total_faces+nfa_i[n]);
      else fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[n+1]+1 ,total_faces+1, total_faces+nfa_i[n]);
      // Giving the face connectivity
      for(int i=0; i<nfacesL; i++)
        if(strcmp (meshL.faces[i].name, bc_names[n].c_str()) == 0)    //(meshL.face_name[i].compare(bc_names[n])==0)
          fprintf(fp,"4 %x %x %x %x %x %x \n",                                          meshL.faces[i].node[0]+1,              meshL.faces[i].node[1]+1,           meshL.faces[i].node[2]+1,           meshL.faces[i].node[3]+1,           meshL.faces[i].elem_1,              meshL.faces[i].elem_2);
      //meshL.face_v[meshL.face_i[i]]+1 ,             meshL.face_v[meshL.face_i[i]+1]+1 ,               meshL.face_v[meshL.face_i[i]+2]+1 ,               meshL.face_v[meshL.face_i[i]+3]+1 ,               meshL.face_elem[(2*i)],             meshL.face_elem[(2*i)+1]);

      for(int i=0; i<nfacesR; i++)
        if(strcmp (meshR.faces[i].name, bc_names[n].c_str()) == 0)    //(meshR.face_name[i].compare(bc_names[n])==0)
          if(meshR.faces[i].elem_2==0)          fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.faces[i].node[0]+1+verticesL,    meshR.faces[i].node[1]+1+verticesL, meshR.faces[i].node[2]+1+verticesL, meshR.faces[i].node[3]+1+verticesL, meshR.faces[i].elem_1+elementsL,    meshR.faces[i].elem_2);
      //fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.face_v[meshR.face_i[i]]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+1]+1 +verticesL ,    meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]);
          else                                  fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.faces[i].node[0]+1+verticesL,    meshR.faces[i].node[1]+1+verticesL, meshR.faces[i].node[2]+1+verticesL, meshR.faces[i].node[3]+1+verticesL, meshR.faces[i].elem_1+elementsL,    meshR.faces[i].elem_2+elementsL);
      //fprintf(fp,"4 %x %x %x %x %x %x \n",    meshR.face_v[meshR.face_i[i]]+1+verticesL ,   meshR.face_v[meshR.face_i[i]+1]+1 +verticesL ,    meshR.face_v[meshR.face_i[i]+2]+1+verticesL ,     meshR.face_v[meshR.face_i[i]+3]+1+verticesL ,     meshR.face_elem[(2*i)]+elementsL,   meshR.face_elem[(2*i)+1]+elementsL);

      fprintf(fp, "))\n");
      total_faces+=nfa_i[n];
      cout<<"Number of faces in boundary "<<bc_names[n] <<" "<<nfa_i[n]<<endl;

    }

    cout<<"Total number of faces in the numerical domain "<<total_faces<<endl;


  }

  fprintf(fp,"\n");
  cout<<"Segmentation after"<<endl;
  //
  // cells
  //
  int cv_element_type = 4; // hexahedral = 4
  int cv_type = 1; //1; // active = 1 - inactive = 32
  char this_zone[128];
  sprintf(this_zone,"fluid");
  fprintf(fp,"(0 \"Cells:\")\n");
  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 3,elementsR+1,elements,cv_type,cv_element_type);
  fprintf(fp,"\n");

  if(zone_noname)
  {
    fprintf(fp,"(0 \"Zones:\")\n");
    fprintf(fp,"(45 (2 fluid fluid))\n");
//    fprintf(fp,"(45 (3 fluid fluid 0))\n");
    fprintf(fp,"(45 (%i wall boundaries))\n", ind_faces[0]);
    fprintf(fp,"(45 (%i interior default-interior)\n", ind_faces[1]);
  }
  else
  {
    fprintf(fp,"(0 \"Zones:\")\n");
    fprintf(fp,"(45 (2 fluid fluid)())\n");
//    fprintf(fp,"(45 (3 fluid fluid 0))\n");
    for(int n=0 ; n<bc_names.size()-1; n++)
    {
      fprintf(fp,"(45 (%i wall %s)())\n", ind_faces[n+1], bc_names[n].c_str());
    }
    fprintf(fp,"(45 (%i interior default-interior)())\n", ind_faces[bc_names.size()]+1);
  }
//  fprintf(fp,"(12 (%x %x %x %x))\n", 0,1,elements,cv_type);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"(45 (%d fluid %s)())\n",  4001,this_zone);
  fprintf(fp,"\n");



  fclose(fp);

  printf(" > Fluent case file written for INTERFACE: %s\n", filename);
}

