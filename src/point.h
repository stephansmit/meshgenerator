/*
 * point.h
 *
 *  Created on: Jan 27, 2012
 *      Author: renep
 */

#ifndef POINT_H_
#define POINT_H_

#include "stdio.h"
#include "math.h"
#include <string>
#include <list>
#include <deque>
#include <algorithm>
using namespace std;

#ifdef WITH_OPENGL
  #include <GL/glu.h>
  #include <GL/glut.h>
#else
  #include "gldummy.h"
#endif

#define PI (4.*atan(1.0))


class POINT
{
public:
  double x, y, z, s, dummy;

public:
  POINT(): x(0.0), y(0.0), z(0.0), s(0.0) {}
  POINT(double value): x(value), y(value), z(value) {}
  POINT(double x_, double y_): x(x_), y(y_), z(0) {}
  POINT(double x_, double y_, double z_): x(x_), y(y_), z(z_) {}

  double radX() const;
  double radZ() const;
  double phiRadX() const;
  double phiRadZ() const;

  double phiDegX() const;
  double phiDegZ() const;

  void rotateDegX(double phi);

  void rotateDegZ(double phi);

  double dot(const POINT &p);

  void drawPoint();

  void print();

  bool operator<(const POINT& p1) const{
    return x < p1.x;     // GUSTAVO 20.06.16 comparison
  }

};


#define VEC POINT

double mag(const POINT &pt);

POINT crossProd(const POINT &p1, const POINT &p2);


POINT operator+(const double val, const POINT &pt);

POINT operator+(const POINT &pt, const double val);

POINT operator+(const POINT &p1, const POINT &p2);

POINT operator-(const POINT &p1, const POINT &p2);

POINT operator-(const double val, const POINT &pt);

POINT operator-(const POINT &pt, const double val);

POINT operator+=(POINT &p1, const POINT &p2);

POINT operator-=(POINT &p1, const POINT &p2);


bool operator==(POINT &p1, const POINT &p2);

bool operator!=(POINT &p1, const POINT &p2);

POINT operator*(const double fact, const POINT &pt);

POINT operator*(const POINT &pt, const double fact);

POINT operator*(const POINT &pt1, const POINT &pt2);

POINT operator/(const POINT &pt1, const POINT &pt2);

POINT operator/(const POINT &pt1, const double val);

POINT operator/(const double val, const POINT &pt1);

POINT rotateDegZ(const double phi, const POINT &p1, const POINT &p2);

// read points from file
void readPts(deque<POINT> &pts, const string &name);



#endif /* POINT_H_ */



