#ifndef UNSTRUCTMESH_H_
#define UNSTRUCTMESH_H_



#define REAL double
#define VOID int
#define ANSI_DECLARATORS


// function for the unstructured mesh generation


#ifdef WITH_OPENGL
  #include <GL/glu.h>
  #include <GL/glut.h>
#else
  #include "gldummy.h"
#endif

#include "string.h"
#include "math.h"
#include "struct.h"
#include "triangle.h"
#include "list"
using namespace std;

//#include "point.h"
//#include "lines.h"

typedef struct {
  POINT pt;
  list<int> neig;
  int marker;
} NODE;

typedef struct {
  deque<int> node;
  int elem_1;
  int elem_2;
  char name[64]="noname";     //
} FACE;          //GUSTAVO 09-05-2016

//class FACE
//{
//public:
//  deque<int> node;
//  int elem_1;
//  int elem_2;
//  char str[24]="gfhgfgh";     //
//};          //GUSTAVO 09-05-2016


typedef struct {
  deque<int> node;
} ELEMENT;

typedef struct{
  POINT insidePoints;
  deque<POINT> holesPoints;
}HOLE;


#define GAMBIT_VER "2.4.6"    // Version of Gambit used to create input file
/***************************************************************************************************************/


/* Define neutral file specific constants */
#define NGRPS 1               // Number of element groups
#define NBSETS 0              // Number of boundary condition sets
#define NDFCD 3               // Number of coordinate directions (2 for neutral file read, 3 for neutral file written)
#define NDFVL 3               // Number of velocity components (2 for neutral file read, 3 for neutral file written, but not used here)

#define NTYPE_TRI 3           // Element geometry type TRIANGLE
#define NTYPE_PRISM 5         // Element geometry type PRISM (=WEDGE)
#define NTYPE_QUAD 2          // Element geometry type QUAD
#define NTYPE_BRICK 4         // ELement geometry type BRICK (=HEX)

#define NDP_TRI 3             // Number of nodes that define element TRIANGLE (node order is 0 -> 1 -> 2)
#define NDP_PRISM 6           // Number of nodes that define element PRISM (node order is 0 -> 1 -> 2 -> 3 -> 4 -> 5)
#define NDP_QUAD 4            // Number of nodes that define element QUAD (node order is 0 -> 1 -> 2 -> 3)
#define NDP_BRICK 8           // Number of nodes that define element BRICK (node order is 0 -> 1 -> 3 -> 2 -> 4 -> 5 -> 7 -> 6)

class UNSTRUCTMESH
{
public:
  deque<NODE> nodes;         ///< list of nodes

  int nno_i;                 ///< number of internal nodes
  deque<int> elem_v;         ///< list of elements nodes
  deque<int> elem_i;         ///< list of elements starting index in elem_v
  int nel_i;                 ///< number of internal elements

  deque<int>    face_v;        ///< list of face nodes                         size: nface_i*nodes/per/face
  deque<int>    face_i;        ///< list of face starting index in face_v      size: nface_i
  deque<int>    face_elem;     ///< list of elements that share the face       size: nface_i*2
  deque<string> face_name;     ///< list of names of each face                 size: nface_i
  int nface_i;                 ///< number of faces

  deque<FACE> faces;        ///< list of faces          //GUSTAVO 09-05-2016
  int n_zone; //=3;                        ///< number of zones faces, by default there are 3: top, bottom and fluid    //GUSTAVO 09-05-2016
  int nfa_i;                        ///< number of faces          //GUSTAVO 09-05-2016

  // old structure
  deque<ELEMENT> elems;

public:

  UNSTRUCTMESH() {}

  UNSTRUCTMESH(STRUCTMESH &structMesh)
  {
    nno_i = 0;
    nel_i = 0;
    // first assign values of the boundary in anti clock wise direction
    for (int j=0; j<structMesh.jmax; j++)
      for (int i=0; i<structMesh.imax; i++)
        if ((i == 0) || (i==structMesh.imax-1) || (j==0) || (j==structMesh.jmax-1))  push_back_node(structMesh.mesh[i][j], 1);
        else                                                                         push_back_node(structMesh.mesh[i][j], 0);

    int count=0;
    for (int j=0; j<structMesh.jmax-1; j++)
      for (int i=0; i<structMesh.imax-1; i++)
      {
        elem_v.push_back(structMesh.imax*j     + i  );
        elem_v.push_back(structMesh.imax*j     + i+1);
        elem_v.push_back(structMesh.imax*(j+1) + i+1);
        elem_v.push_back(structMesh.imax*(j+1) + i  );
        //--------------------------------------------//GUstavo 10/06/2016
        //              FACE tmp_face;
        //
        //              tmp_face.node.push_back(structMesh.imax*j     + i  );
        //              tmp_face.node.push_back(structMesh.imax*j     + i+1);
        //              tmp_face.node.push_back(structMesh.imax*(j+1) + i+1);
        //              tmp_face.node.push_back(structMesh.imax*(j+1) + i  );
        //              tmp_face.elem_1 = count+1;
        //              tmp_face.elem_2 = 0;
        //              faces.push_back(tmp_face);
        //--------------------------------------------
        elem_i.push_back(4*count);
        count++;
      }
    elem_i.push_back(4*count);

    //=============================================================================
    // move boundary nodes and elements (ordered) at the end of the nodes list
    //=============================================================================
    deque<POINT> boundPts;
    for (int i=0; i<structMesh.imax; i++)       boundPts.push_back(structMesh.mesh[i][0]);
    for (int j=1; j<structMesh.jmax; j++)        boundPts.push_back(structMesh.mesh[structMesh.imax-1][j]);
    for (int i=structMesh.imax-2; i>=0; i--)    boundPts.push_back(structMesh.mesh[i][structMesh.jmax-1]);
    for (int j=structMesh.jmax-2; j>0; j--)      boundPts.push_back(structMesh.mesh[0][j]);

    moveBoundaryNodesEnd(boundPts);
    moveBoundaryElementsEnd(boundPts);
    removeDoublePoints();
  }

    void push_back_node(const POINT &pt, int marker = 0)
    {
        NODE nd;
        nd.pt = pt;
        nd.marker = marker;
        nodes.push_back(nd);
    }

  // OK
  UNSTRUCTMESH & operator+=(const UNSTRUCTMESH &mesh);

//  UNSTRUCTMESH & operator-=(const UNSTRUCTMESH &mesh);

  // OK
  UNSTRUCTMESH operator+(const UNSTRUCTMESH &unstructMesh);

//  UNSTRUCTMESH operator-(const UNSTRUCTMESH &unstructMesh);

  // OK
  void moveBoundaryNodesEnd(deque<POINT> boundPts);
  // OK
  void moveBoundaryElementsEnd(deque<POINT> boundPts);
  // OK
  void moveNodeEnd(int ind);
  // OK
  void moveElementEnd(int ind);
  // OK
  void moveFaceEnd_old(int ind); //GUstavo 10/06/2016

  void moveFaceEnd(int ind); //GUstavo 27/07/2016

  void removeDoublePoints();
  // OK
  void findNodeNeighbors();
  // OK
  void elementObject(); //GUstavo 10/06/2016

  bool CompareFaceFace(FACE &face1, FACE &face2); //GUstavo 10/06/2016

  void removeDoubleFaces(); //GUstavo 10/06/2016

  void removeDoubleFaces2D_obj(); //GUstavo 10/06/2016

  void removeDoubleFaces2D(); //GUstavo 27/07/2016

  bool CompareElementFace(ELEMENT &element_, FACE &face_); //GUstavo 10/06/2016

  void findFaces2D(); //GUstavo 27/07/2016

  void findFaces2D_obj(); //GUstavo 10/06/2016
    //
  void setMarker(const deque<POINT> &points, int marker);
  // OK
  void replacePoints(const deque<POINT> &oldPoints, const deque<POINT> &newPoints);
  // OK
  // Laplacian smoother
  void smoothMesh(int maxit, double tol);
  // OK
  void smoothMesh_(int maxit, double tol, double relax = 0.6);
    // OK
  void writeMarkedNodes(const char * name, int marker);
  // OK

  void ScaleMesh(const double scaling_factor);  //Gustavo 24/06/2016

  void writeTecplot(const char * name, int ndim, bool onlyBoundary = false);
  // OK
  void writeGambitNeu(const char *name, int ndim, bool onlyBoundary = false);

  void readFluentMsh2D(const char *name); //Gustavo 29/07/2016

  void readFluentMsh(const char *name, int nvert_gambit); //Gustavo 29/07/2016

  void readFluentMsh_named(const char *name, int nvert_gambit); //Gustavo 01/08/2016    NOT FINISHED!!!

  void readGambitNeu2D(const char *name);

  void readGambitNeu(const char *name);

  void readNastranMod(const char *name);

#define FLUENT_TYPE_INTERIOR 2
#define FLUENT_TYPE_WALL 3
  //#define OUTPUTHEX

  void writeFluentMsh_Original(char * filename, int nbl);//, block *bl)

//  void writeFluentMsh_New(const char *name, int ndim, deque<string> bc_name); // With old face structure

  void writeFluentMsh(const char *name, int ndim, deque<string> bc_name); //Gustavo 28/07/16

  void writeFluentMsh_obj(const char *name, int ndim, deque<string> bc_name); //Gustavo 22/08/16

  void drawMesh();

  void drawSolid();
};





class TRIANGULATE : public UNSTRUCTMESH
{
public:
    // these are the external boundaries
    deque<POINT> extBoundary;

    // these are the holes in the mesh
    deque<HOLE> holes;

    // string of parameters for the triangulation
    char *triangParameters;


public:

    TRIANGULATE(): UNSTRUCTMESH() {}

  void unstructuredMesh2D();
};


UNSTRUCTMESH operator*(const UNSTRUCTMESH &mesh, const double val);

UNSTRUCTMESH operator*=(UNSTRUCTMESH &mesh, const double val);


#endif
