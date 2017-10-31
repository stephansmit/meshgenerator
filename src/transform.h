
#ifndef TRANSFORM_H_
#define TRANSFORM_H_

#include "math.h"
#include "lines.h"
#include "struct.h"
#include "unstruct.h"


#include <deque>
using namespace std;

class TRANSFORM
{
public:
  LINE hubXYZ, hubMPS;
  LINE shrXYZ, shrMPS;

public:
//  TRANSFORM() {}

  void initTransform(const LINE &hub, const LINE &shr);

  void transformRZPHI_to_XYZ(deque<POINT> &pts);

  void transformHubShr(const LINE &lineXYZ, LINE &lineMPS, double span);

  double calcSpan(POINT &pts);


  void transform(POINT &pts, double span);
  void transform(deque<POINT> &pts, double span);
  void transform(STRUCTMESH &mesh, double span);
  void transform(UNSTRUCTMESH &unstr, double span);

  void transform(POINT &pts);
  void transform(deque<POINT> &pts);
  void transform(STRUCTMESH &mesh);
  void transform(UNSTRUCTMESH &unstr);


  void transformBack(POINT &pts);
  void transformBack(deque<POINT> &pts);
  void transformBack(deque<NODE> &nodes);
  void transformBack(STRUCTMESH &mesh);
  void transformBack(UNSTRUCTMESH &unstr);
};

#endif
