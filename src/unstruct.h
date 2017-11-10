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
  int marker; //int marker;
} NODE;

typedef struct {
  deque<int> node;
  int elem_1;
  int elem_2;
  char name[64];//="noname";     //
} FACE;          //GUSTAVO 09-05-2016



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

//  deque<int>    face_v;      ///< list of face nodes                         size: number of face *nodes/per/face
//  deque<int>    face_i;      ///< list of face starting index in face_v      size: number of face
//  deque<int>    face_elem;   ///< list of elements that share the face       size: number of face*2
//  deque<string> face_name;   ///< list of names of each face                 size: number of face


  deque<FACE> faces;         ///< list of faces          //GUSTAVO 09-05-2016
  int n_zone;     ;          ///< number of zones faces, by default there are 3: top, bottom and fluid    //GUSTAVO 09-05-2016
  int nfa_nno_i;			 ///< number of internal faces that have only internal nodes
  int nfa_i;				 ///< number of internal faces (includes faces with external nodes)   nfa_i>nfa_nno_i


public:

  UNSTRUCTMESH() {}

  UNSTRUCTMESH(STRUCTMESH &structMesh)
  {
    nno_i 	= 0;
    nel_i 	= 0;
    nfa_i = 0;
    nfa_nno_i = 0;
    n_zone=0;

    //======================================================
    // Assigning the NODES!!!!
    //------------------------------------------------------
    //First assigning the internal nodes
    for (int j=1; j<structMesh.jmax-1; j++)
      for (int i=1; i<structMesh.imax-1; i++)
        push_back_node(structMesh.mesh[i][j], 0);
    //Number of internal nodes
    nno_i=nodes.size();
    //------------------------------------------------------
    //Now assigning the boundary nodes in anti clock wise direction
    for (int i=0; i<structMesh.imax-1; i++)
      push_back_node(structMesh.mesh[i][0], 1);
    for (int j=0; j<structMesh.jmax-1; j++)
      push_back_node(structMesh.mesh[structMesh.imax-1][j], 1);
    for (int i=0; i<structMesh.imax-1; i++)
      push_back_node(structMesh.mesh[structMesh.imax-1-i][structMesh.jmax-1], 1);
    for (int j=0; j<structMesh.jmax-1; j++)
      push_back_node(structMesh.mesh[0][structMesh.jmax-1-j], 1);
    //======================================================
    // Assigning the ELEMENTS!!!!
    int count=0;
    //------------------------------------------------------
    //Internal Elements
    int imax_inside= structMesh.imax - 2;
    for (int j=0; j<structMesh.jmax-3; j++)
      for (int i=0; i<structMesh.imax-3; i++)
      {
        elem_v.push_back(imax_inside*j     + i  );
        elem_v.push_back(imax_inside*j     + i+1);
        elem_v.push_back(imax_inside*(j+1) + i+1);
        elem_v.push_back(imax_inside*(j+1) + i  );
        elem_i.push_back(4*count);
        count++;
      }
    //Number of internal elements
    nel_i=elem_i.size();

    //Now assigning the boundary nodes in anti clock wise direction
    //------------------------------------------------------
    //Boundary Elements
    int nno_i_tmp=nno_i;		//Define as the node in the "corner", bottom left
    //*** CORNER, bottom left corner element
    elem_v.push_back(nno_i);
    elem_v.push_back(nno_i+1);
    elem_v.push_back(0);
    elem_v.push_back(nodes.size()-1);
    elem_i.push_back(4*count);
    count++;
    //*** Bottom boundary elements
    for (int i=1; i<structMesh.imax-2; i++)
    {
      elem_v.push_back(nno_i    	+ i  );
      elem_v.push_back(nno_i     	+ i+1);
      elem_v.push_back(			+ i  );
      elem_v.push_back(			+ i-1);
      elem_i.push_back(4*count);
      count++;
    }
    nno_i_tmp+= structMesh.imax-1;		//Define as the node in the "corner", bottom right
    //*** CORNER, bottom right corner element
    elem_v.push_back(nno_i_tmp-1);		//nno_i    	+ structMesh.imax-2);
    elem_v.push_back(nno_i_tmp);		//nno_i    	+ structMesh.imax-1);
    elem_v.push_back(nno_i_tmp+1);		//nno_i    	+ structMesh.imax  );
    elem_v.push_back(			+ structMesh.imax-3);
    elem_i.push_back(4*count);
    count++;
    //*** Right boundary elements
    for (int j=1; j<structMesh.jmax-2; j++)
    {
      elem_v.push_back(-1			+	imax_inside*j		);
      elem_v.push_back(nno_i_tmp	+	j					);		//nno_i    	+ structMesh.imax	+	j-1 				);
      elem_v.push_back(nno_i_tmp	+	j+1					);		//nno_i    	+ structMesh.imax	+	j   				);
      elem_v.push_back(-1			+	imax_inside*(j+1)	);
      elem_i.push_back(4*count);
      count++;
    }
    nno_i_tmp+= structMesh.jmax-1;		//Define as the node in the "corner", top right
    //*** CORNER, top right corner element
    elem_v.push_back(nno_i-1	 );
    elem_v.push_back(nno_i_tmp-1);			//nno_i    	+ structMesh.imax-2	+ 	structMesh.jmax-2 +1	);
    elem_v.push_back(nno_i_tmp	 );			//nno_i    	+ structMesh.imax-2	+ 	structMesh.jmax 	);
    elem_v.push_back(nno_i_tmp+1);			//nno_i    	+ structMesh.imax-2 + 	structMesh.jmax+1	);
    elem_i.push_back(4*count);
    count++;
    //*** Top boundary elements
    for (int i=1; i<structMesh.imax-2; i++)
    {
      elem_v.push_back(nno_i    		- i-1  	);
      elem_v.push_back(nno_i     		- i		);
      elem_v.push_back(nno_i_tmp		+ i		);			//nno_i	    	+ structMesh.imax-2	+ 	structMesh.jmax			+i	);
      elem_v.push_back(nno_i_tmp		+ i+1	);			//nno_i	    	+ structMesh.imax-2	+ 	structMesh.jmax			+i+1);
      elem_i.push_back(4*count);
      count++;
    }
    nno_i_tmp+= structMesh.imax-1;		//Define as the node in the "corner", top left
    //*** CORNER, top left corner element
    elem_v.push_back(nno_i_tmp		+	1				);		//nodes.size()	- structMesh.jmax+2	 		 					);
    elem_v.push_back(nno_i    		- structMesh.imax+2	);
    elem_v.push_back(nno_i_tmp		-	1				);		//nodes.size()	- structMesh.jmax								);
    elem_v.push_back(nno_i_tmp  						);		//nodes.size()	- structMesh.jmax+1 		 					);
    elem_i.push_back(4*count);
    count++;
    //*** Left boundary elements
    for (int j=1; j<structMesh.jmax-2; j++)
    {
      elem_v.push_back(nno_i_tmp		+	j+1						);			//nodes.size()	- structMesh.jmax+2	+   j										);
      elem_v.push_back(imax_inside*(structMesh.jmax-2-j-1)		);
      elem_v.push_back(imax_inside*(structMesh.jmax-2-j)			);
      elem_v.push_back(nno_i_tmp		+	j						);			//nodes.size()	- structMesh.jmax+2	+   j-1										);
      elem_i.push_back(4*count);
      count++;
    }
    elem_i.push_back(4*count);
    //cout<<"Number of internal elements"<< nel_i <<endl;

    //======================================================
    // Assigning the FACES!!!
    int nelem_imax=structMesh.imax-3;	//Number of internal elements in i
    int nelem_jmax=structMesh.jmax-3;	//Number of internal elements in j
    int	nelem_i=nelem_imax*nelem_jmax;
    if(nelem_i-nel_i!=0)		cout<<"Error in the number of  elements: " << nelem_i-nel_i <<endl;

    //------------------------------------------------------
    //Internal Faces
    //*** From the internatl Elements
    count=0;
    for (int j=0; j<nelem_jmax; j++)
      for (int i=0; i<nelem_imax; i++)
      {
        int count_face=0;
        int fstNodeInd = 4*(count);
        int lstNodeInd = 4*(count+1) -1;
        for(int k=fstNodeInd; k<lstNodeInd; k++)
        {
          FACE tmp_face;
          tmp_face.node.push_back(elem_v[k+1]);
          tmp_face.node.push_back(elem_v[k]);
          tmp_face.elem_1 = count+1;
          sprintf(tmp_face.name,"fluid");
          count_face++;

          if((count_face==1) && (j==0))							tmp_face.elem_2 = nel_i+2+i;
          else if((count_face==2) && (i< (nelem_imax-1)))		tmp_face.elem_2 = count+2;
          else if((count_face==2) && (i==(nelem_imax-1)))		tmp_face.elem_2 = nel_i+nelem_imax+2+j+1;
          else if((count_face==3) && (j< (nelem_jmax-1)))		tmp_face.elem_2 = nelem_imax*(j+1) + i + 1;
          else if((count_face==3) && (j==(nelem_jmax-1)))		tmp_face.elem_2 = nel_i+2*(nelem_imax+1)+(nelem_jmax+1)-i;
          else				  									continue;

          faces.push_back(tmp_face);

        }
        if(i==0)		//%nelem_imax==0)
        {
          FACE tmp_face;
          tmp_face.node.push_back(elem_v[fstNodeInd]);
          tmp_face.node.push_back(elem_v[lstNodeInd]);
          tmp_face.elem_1 = count+1;
          tmp_face.elem_2 = elem_i.size()-1-j;
          sprintf(tmp_face.name,"fluid");
          faces.push_back(tmp_face);
        }
        count++;
      }

    //Number of internal faces that have only internal nodes
    nfa_nno_i=faces.size();
    //*** From the boundary Elements
    int count_bc=count;
    //Bottom elements
    for (int i=0; i<nelem_imax+1; i++)
    {	
      int fstNodeInd = 4*(count)+1;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd+1]);
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = count+2;
      sprintf(tmp_face.name,"fluid");
      faces.push_back(tmp_face);
      count++;
    }
    //Right elements
    for (int j=0; j<nelem_jmax+1; j++)
    {	
      int fstNodeInd = 4*(count)+2;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd+1]);
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = count+2;
      sprintf(tmp_face.name,"fluid");
      faces.push_back(tmp_face);
      count++;
    }
    //Top elements
    for (int i=0; i<nelem_imax+1; i++)
    {	
      int fstNodeInd = 4*(count);
      int lstNodeInd = 4*(count+1) -1;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.node.push_back(elem_v[lstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = count+2;
      sprintf(tmp_face.name,"fluid");
      faces.push_back(tmp_face);
      count++;
    }
    //Left elements
    for (int j=0; j<nelem_jmax+1; j++)
    {	
      int fstNodeInd = 4*(count);
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd+1]);
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.elem_1 = count+1;
      if(j<nelem_jmax)	tmp_face.elem_2 = count+2;
      else				tmp_face.elem_2 = nel_i+1;
      sprintf(tmp_face.name,"fluid");
      faces.push_back(tmp_face);
      count++;
    }

    //Number of internal faces
    nfa_i=faces.size();

    //------------------------------------------------------
    //External Faces
    count= count_bc;
    //Bottom elements
    for (int i=0; i<nelem_imax+2; i++)
    {	
      int fstNodeInd = 4*(count);
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd+1]);
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = 0;
      sprintf(tmp_face.name,"noname");
      faces.push_back(tmp_face);
      if(i<nelem_imax+1)	count++;
    }
    //Right elements
    for (int j=0; j<nelem_jmax+2; j++)
    {	
      int fstNodeInd = 4*(count)+1;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd+1]);
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = 0;
      sprintf(tmp_face.name,"noname");
      faces.push_back(tmp_face);
      if(j<nelem_jmax+1)	count++;
    }
    //Top elements
    for (int i=0; i<nelem_imax+2; i++)
    {	
      int fstNodeInd = 4*(count)+2;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd+1]);
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = 0;
      sprintf(tmp_face.name,"noname");
      faces.push_back(tmp_face);
      if(i<nelem_imax+1)	count++;
    }
    //Left elements
    for (int j=0; j<nelem_jmax+1; j++)
    {	
      int fstNodeInd = 4*(count);
      int lstNodeInd = 4*(count+1) -1;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.node.push_back(elem_v[lstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = 0;
      sprintf(tmp_face.name,"noname");
      faces.push_back(tmp_face);
      count++;
    }
    count= count_bc;
    {
      int fstNodeInd = 4*(count);
      int lstNodeInd = 4*(count+1) -1;
      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.node.push_back(elem_v[lstNodeInd]);
      tmp_face.elem_1 = count+1;
      tmp_face.elem_2 = 0;
      sprintf(tmp_face.name,"noname");
      faces.push_back(tmp_face);

    }

    //======================================================
    // Finding the nodes Neighbours
    // identify neighbours of every point
    findNodeNeighbors();

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
  void removeDoublePoints();
  // OK
  void Store_NodesElemFaces(const string name);	//Gustavo 31/08/2016
  // ok
  void UpdateBoundaries();	//Gustavo 25/08/2016

  void findNodeNeighbors();

  void findFaces2D();

  bool compareFaces(const FACE &f1, const FACE &f2);

  void removeDoubleFaces(); //GUstavo 10/06/2016
  //
  void setMarker(const deque<POINT> &points, int marker, bool boundarypts = false);
  // OK
  void replacePoints(const deque<POINT> &oldPoints, const deque<POINT> &newPoints);
  // OK
  // Laplacian smoother
  void smoothMesh(int maxit, double tol);
  // OK
  void smoothMesh_(int maxit, double tol, double relax = 0.6);
  // OK
  void smoothMesh_options(int maxit, double tol , SPLINE splineA, SPLINE splineB, double relax = 0.6 );
    // OK       /Gustavo 8/9/2016
  void smoothMesh_generic(int maxit, double tol , deque<SPLINE> spline, double relax=0.6 ,  bool periodic_splines=false);
    // OK       //Gustavo 9/9/2016
  void smoothMesh_freeSplines(int maxit, double tol , SPLINE original_spline_A, SPLINE original_spline_B, double relax=0.6, bool perSp_CYL=false);
    // OK       //Gustavo 9/9/2016
  void writeMarkedNodes(const char * name, int marker);
  // OK

  void ScaleMesh(const double scaling_factor);  //Gustavo 24/06/2016

  void writeTecplot(const char * name, int ndim, bool onlyBoundary = false);
  // OK
  void writeGambitNeu(const char *name, int ndim, bool onlyBoundary = false);

  void readFluentMsh2D(const char *name, bool writtenMeshGen=false); //Gustavo 29/07/2016

  void readFluentMsh(const char *name, int nvert_gambit);

  void readGambitNeu2D(const char *name);

  void readGambitNeu(const char *name);

  void readNastranMod(const char *name);

#define FLUENT_TYPE_INTERIOR 2
#define FLUENT_TYPE_WALL 3

  void writeFluentMsh(const char *name, int ndim, deque<string> bc_name);

  void writeSU2(const char *name, int ndim, deque<string> bc_names); // GUstavo 3/11/2017

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

  void moveBoundaryElementsEnd(deque<POINT> boundPts);
};


UNSTRUCTMESH operator*(const UNSTRUCTMESH &mesh, const double val);

UNSTRUCTMESH operator*=(UNSTRUCTMESH &mesh, const double val);


#endif
