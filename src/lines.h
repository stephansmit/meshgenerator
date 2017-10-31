
#ifndef LINES_H_
#define LINES_H_

#include "math.h"
#include "point.h"
#include "functions.h"

#ifdef WITH_OPENGL
  #include <GL/glu.h>
  #include <GL/glut.h>
#else
  #include "gldummy.h"
#endif

#include <deque>
using namespace std;


//------------------------------------------------------------------------------------
// basic line class, contains deque list for the points called controlPts
// method calcPoint is a simple linear interpolation
// in calcPoint it is possible to define a sub-line by giving a start and end point
// in this case t is relative to this sub-line
class LINE
{
public:
  deque<POINT> controlPt;
  double length;
  int iStart, iEnd;

public:
  virtual ~LINE() {}
  LINE() {}
  LINE(int N) {controlPt.resize(N);}
  LINE(const deque<POINT> &ctrPts) : controlPt(ctrPts)
  {
    iStart = 0;
    iEnd = controlPt.size()-1;
    init();
  }

  LINE(const POINT &Pt1, const POINT &Pt2)
  {
    controlPt.push_back(Pt1);
    controlPt.push_back(Pt2);
    iStart = 0;
    iEnd = controlPt.size()-1;
    init();
  }

  LINE(const string &name)
  {
    readPts(controlPt, name);
    iStart = 0;
    iEnd = controlPt.size()-1;
    init();
  }

  int size() const;

  virtual void init();

//  LINE operator+(const POINT point);

  POINT & operator[](int i);

  POINT const & operator[](int i) const;

  void calcS(int i1, int i2);

  virtual POINT calcPoint(double t, double offset1 = 0.0, double offset2 = 1.0) const;


  void split(deque<POINT> &line1, deque<POINT> &line2, double t);

  template <class T> deque<T> multipleSplit(deque<double> splitT)
  {
    deque<T> lines;

    deque<int> start_i;
    start_i.push_back(iStart);

    for (int i=0; i<splitT.size(); i++)
    {
      int ind=iStart;
      while ((ind<iEnd) && (splitT[i]>controlPt[ind].s)) ind++;
      controlPt.insert(controlPt.begin()+ind, calcPoint(splitT[i]));
      start_i.push_back(ind);
      iEnd++;
    }
    start_i.push_back(iEnd);

    init();

    // define the splines
    for (int i=0; i<start_i.size()-1; i++)
    {
      T tmp;
      tmp.controlPt = controlPt;
      tmp.iStart = start_i[i];
      tmp.iEnd = start_i[i+1];
      tmp.init();
      lines.push_back(tmp);
    }

    return(lines);
  };

  double intersect(const POINT &p11, const POINT &p12, const POINT &p21, const POINT &p22);

  /*! \brief Computes intersection point of two lines using their control points.
   *
   *  Intersection point is only accurate if both lines are of type LINE
   */
  POINT intersectXY(LINE& line2, bool insertPoint);

  /*! \brief Computes intersection point of two lines using their non-dimensional coordinate s.
   *
   *  To be used for SPLINE and BEZIER curves
   *  \param line2 is the line to perform intersection with
   *  \param insertPoint boolean to set if intersection point should be added
   *  \param eps is the truncation error of the iteration to compute the intersection point
   */
  POINT intersectXY(LINE& line2, bool insertPoint, double eps);

  /*! \brief Computes the orthogonal projection of a point on the line
   *
   *  \param pt is the point to be projected
   */
  POINT project2D(POINT &pt);

  void rotateDegX(double phi);

  void rotateDegZ(double phi);

  void write(const char *name, int N);

  void write(const char *name, LINE &distr, int N);

  void writeCtrPts(const char *name);

  virtual void drawLine();
};

//---------------------------------------------------------------------------
// BEZIER class inherits data from LINE and overloads calcPoint(double t, double offset1 = 0.0, double offset2 = 1.0)
class BEZIER: public LINE
{
public:
  BEZIER() : LINE() {}
  BEZIER(const deque<POINT> &ctrPts) : LINE(ctrPts) {}
  virtual ~BEZIER() {}


  virtual POINT calcPoint(double t, double offset1 = 0.0, double offset2 = 1.0) const;

  virtual void drawLine();
};


//---------------------------------------------------------------------------
// SPLINE class inherits data from LINE and overloads calcPoint(double t, double offset1 = 0.0, double offset2 = 1.0)
// it is based on Numerical recipes and uses a third order spline interpolation
class SPLINE: public LINE
{
public:
  deque<POINT> deriv;

public:
  virtual ~SPLINE() {}
  SPLINE() : LINE() {}
  SPLINE(const deque<POINT> &ctrPts) : LINE(ctrPts)  { init(); }
  SPLINE(const deque<POINT> &ctrPts, int i1, int i2) : LINE(ctrPts)
  {
    iStart = i1;
    iEnd = i2;
    init();
  }

  SPLINE(const string &name)
  {
    readPts(controlPt, name);
    iStart = 0;
    iEnd = controlPt.size()-1;
    init();
  }


//  SPLINE(LINE line)
//  {
//    controlPt = line.controlPt;
//    iStart = line.iStart;
//    iEnd = line.iEnd;
//    init();
//  }

  virtual void init();

  virtual POINT calcPoint(double t, double offset1 = 0.0, double offset2 = 1.0) const;

  void calcDerivative();

  /*! \brief spline offset
   *
   */
  SPLINE offsetRadial(const double thickness);


  POINT calcNorm2D(double t) const;

  virtual void drawLine();

};


//---------------------------------------------------------------------------
// POLYGON class (closed LINE)
class POLYGON: public LINE
{
public:
//  POLYGON() : LINE() {}
  POLYGON(const deque<POINT> &ctrPts) : LINE(ctrPts)
  {
  	if (controlPt[controlPt.size()-1]!=controlPt[0])
  		{controlPt.push_back(ctrPts[0]);	init();}
  }
  virtual ~POLYGON() {}

  /*! \brief Check if a point is inside the polygon. Works just in x,y plane.
   *
   *  \param pt point to check
   */
  bool pointInsidePolygon(POINT &pt);
};


////===================================================
////                    OPERATORS
////===================================================
template <class T>
T operator+(const T curve, const POINT point)
{
  deque<POINT> newPts;

  for (int i=0; i<curve.size(); i++)
    newPts.push_back(curve[i]+point);

  T newCurve(newPts, curve.iStart, curve.iEnd);

  return(newCurve);
};

#endif
