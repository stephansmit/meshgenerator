/*
 * unstruct.cpp
 *
 *  Created on: Mar 19, 2012
 *      Author: enrico
 */
#include "unstruct.h"

//=============================================
//        UNSTRUCTMESH CLASS FUNCTIONS
//=============================================
UNSTRUCTMESH & UNSTRUCTMESH::operator+=(const UNSTRUCTMESH &mesh)
{

  //Sizes of the mesh this->
  int nno_i_this					    = this->nno_i;				      //Number of internal nodes
  int nodes_size_this				    = this->nodes.size();		      //Number of nodes
  int nel_i_this					    = this->nel_i;				      //Number of internal elements
  int elem_i_nel_i_this			        = this->elem_i[nel_i];            //Location of last internal elements
  int elem_i_size_this			        = this->elem_i.size();		      //Number of elements
  int elem_i_elem_i_size_this	        = this->elem_i[elem_i.size()-1];  //Location of the last element
  int elem_v_size_this			        = this->elem_v.size();            //Total of nodes in all elements
  int nfa_nno_i_this			        = this->nfa_nno_i;			      //Number of internal faces
  int nfa_i_this				        = this->nfa_i;			          //Number of internal faces
  int faces_size_this				    = this->faces.size();		      //Number of faces


  UNSTRUCTMESH tmp_mesh=mesh;

  //************************************
  //		NODES
  //Adding the list of nodes together
  nodes.insert(nodes.begin()+nno_i_this				,mesh.nodes.begin()							    ,mesh.nodes.begin()+mesh.nno_i);
  nodes.insert(nodes.end()							,mesh.nodes.begin()+mesh.nno_i					,mesh.nodes.end());

  printf("MESH + operator: Summing of the nodes done. \n");

  //************************************
  //		ELEMENTS
  //Updating ELEM_V
  for (int i=0; i<elem_v.size(); i++)
    if (elem_v[i]>=nno_i_this)                              elem_v[i]			      +=mesh.nno_i;
  for (int i=0; i<tmp_mesh.elem_v.size(); i++)
    if (tmp_mesh.elem_v[i]>=tmp_mesh.nno_i)                 tmp_mesh.elem_v[i]	+=nodes_size_this;
    else                                                    tmp_mesh.elem_v[i]	+=nno_i_this;

  //Adding the list of elements nodes together ELEM_V
  elem_v.insert(elem_v.begin()+elem_i[nel_i_this]		    ,tmp_mesh.elem_v.begin()									,tmp_mesh.elem_v.begin()+tmp_mesh.elem_i[tmp_mesh.nel_i]);
  elem_v.insert(elem_v.end()							    ,tmp_mesh.elem_v.begin()+tmp_mesh.elem_i[tmp_mesh.nel_i]	,tmp_mesh.elem_v.end());



  //Updating ELEM_I
  for(int i=0; i<tmp_mesh.nel_i; i++)                       tmp_mesh.elem_i[i]  +=elem_i_nel_i_this;
  for(int i=nel_i; i<elem_i_size_this; i++)                 elem_i[i]           +=tmp_mesh.elem_i[tmp_mesh.nel_i];
  for(int i=tmp_mesh.nel_i; i<tmp_mesh.elem_i.size(); i++)  tmp_mesh.elem_i[i]  +=elem_i_elem_i_size_this;

  //Adding the list of elements index together ELEM_I
  elem_i.insert(elem_i.begin()+nel_i_this			        ,tmp_mesh.elem_i.begin()					              ,tmp_mesh.elem_i.begin()+tmp_mesh.nel_i);
  elem_i.insert(elem_i.end()-1						        ,tmp_mesh.elem_i.begin()+tmp_mesh.nel_i                   ,tmp_mesh.elem_i.end()-1);
  elem_i[elem_i.size()-1]=elem_v.size();

  printf("MESH + operator: Summing of the elements done. \n");


  //************************************
  //		FACES
  //Updating vertices and elements of the face
  for (int i=0; i<faces.size(); i++)
  {
    //--- Nodes
    for(int j=0; j<faces[i].node.size(); j++)
      if (faces[i].node[j]>=nno_i_this)               faces[i].node[j]          +=mesh.nno_i;

    //--- Elements
    if (faces[i].elem_1 >nel_i_this)                  faces[i].elem_1				    +=mesh.nel_i;
    if (faces[i].elem_2 >nel_i_this)                  faces[i].elem_2				    +=mesh.nel_i;
  }
  for (int i=0; i<tmp_mesh.faces.size(); i++)
  {
    //--- Nodes
    for(int j=0; j<tmp_mesh.faces[i].node.size(); j++)
      if (tmp_mesh.faces[i].node[j]>=tmp_mesh.nno_i)  tmp_mesh.faces[i].node[j] +=nodes_size_this;
      else                                            tmp_mesh.faces[i].node[j] +=nno_i_this;

    //--- Elements
    if (tmp_mesh.faces[i].elem_1 >tmp_mesh.nel_i)     tmp_mesh.faces[i].elem_1	+=elem_i_size_this-1;	//-1
    else											                        tmp_mesh.faces[i].elem_1	+=nel_i_this;

    if (tmp_mesh.faces[i].elem_2 >tmp_mesh.nel_i)     tmp_mesh.faces[i].elem_2	+=elem_i_size_this-1;	//-1
    else if (tmp_mesh.faces[i].elem_2 != 0) 		      tmp_mesh.faces[i].elem_2	+=nel_i_this;
  }

  //Adding the list of faces objects together
  if((nfa_nno_i==0) || (tmp_mesh.nfa_nno_i==0))
  {
    faces.insert(faces.begin()+nfa_i_this		                    , tmp_mesh.faces.begin()					  ,tmp_mesh.faces.begin()+tmp_mesh.nfa_i);
    faces.insert(faces.end()					                    , tmp_mesh.faces.begin()+tmp_mesh.nfa_i	      ,tmp_mesh.faces.end());
  }
  else
  {
    faces.insert(faces.begin()+nfa_nno_i_this						, tmp_mesh.faces.begin()					  ,tmp_mesh.faces.begin()+tmp_mesh.nfa_nno_i);
    faces.insert(faces.begin()+nfa_i_this+tmp_mesh.nfa_nno_i		, tmp_mesh.faces.begin()+tmp_mesh.nfa_nno_i	  ,tmp_mesh.faces.begin()+tmp_mesh.nfa_i);
    faces.insert(faces.end()										, tmp_mesh.faces.begin()+tmp_mesh.nfa_i		  ,tmp_mesh.faces.end());
  }

  printf("MESH + operator: Summing of the faces done. \n");

  //Number of internal nodes, elements, faces (with internal nodes, with internal and internal nodes

  nno_i 		+= tmp_mesh.nno_i;
  nel_i 		+= tmp_mesh.nel_i;
  nfa_nno_i     += tmp_mesh.nfa_nno_i;
  nfa_i 	    += tmp_mesh.nfa_i;


  // This values are not correct, they need to be updated as nodes, elements and faces are internal now that the two meshes were summed.


  //************************************
  // Removing doubles, updating boundaries, finding neighbors
  printf("MESH + operator: Before Double points removed. \n");
  //----------------------------------
  //	REMOVING DOUBLE POINTS
  // Based on how meshes are added, double points can be on the boundaries only
  removeDoublePoints();

  printf("MESH + operator: Before Double faces removed. \n");
  //----------------------------------
  //	REMOVING DOUBLE FACES
  // Based on how meshes are added, double faces can be on the boundaries only
  removeDoubleFaces();

  //----------------------------------
  //  UPDATING NODES, ELEMENTS, FACES IN THE BOUNDARIES
  // Updating the nodes, elements and faces in the boundaries, also updating the integers above.
  UpdateBoundaries();

  printf("MESH + operator: Boundaries updated. \n");

  //----------------------------------
  // FINDING THE NEIGHBORS OF THE NODES
  //Clearing the out-dated neighbors!
  for (int i=0; i<nodes.size(); i++)
    nodes[i].neig.clear();

  // identify neighbours of every point using the faces
  findNodeNeighbors();

  printf("MESH + operator: Sum of meshes finish. \n");

  //Returning the summed mesh!
  return *this;
}

UNSTRUCTMESH UNSTRUCTMESH::operator+(const UNSTRUCTMESH &unstructMesh)
{
  UNSTRUCTMESH tmp = *this;
  tmp += unstructMesh;
  return tmp;
}

void UNSTRUCTMESH::moveBoundaryNodesEnd(deque<POINT> boundPts)
{
  for (int i=0; i<(int)boundPts.size(); i++)
  {
    int ind=0;
    while ((ind<nodes.size()) && (nodes[ind].pt != boundPts[i])) ind++;
    if (!(ind == nodes.size())) moveNodeEnd(ind);
  }
  nno_i = nodes.size()-boundPts.size();
}

void UNSTRUCTMESH::moveNodeEnd(int ind)
{
  // move node to the end of the deque
  NODE tmpNode = nodes[ind];

  nodes.erase(nodes.begin()+ind);
  nodes.push_back(tmpNode);

  // update nodes in elements;
  for (int i=0; i<elem_v.size(); i++)
    if (elem_v[i]==ind)      elem_v[i] = nodes.size()-1;
    else if (elem_v[i]>ind)  elem_v[i]--;

  // update nodes in Faces;
  for (int fa=nfa_nno_i; fa<faces.size(); fa++)
  {
    for(int j=0; j<faces[fa].node.size(); j++)
      if (faces[fa].node[j]==ind)      faces[fa].node[j] = nodes.size()-1;
      else if (faces[fa].node[j]>ind)  faces[fa].node[j]--;
  }
}


void UNSTRUCTMESH::moveBoundaryElementsEnd(deque<POINT> boundPts)
{
  int elemMoved = 0;
  for (int i=0; i<boundPts.size()-1; i++)
  {
    int ind=0;

    while (ind<elem_i.size()-1-elemMoved)
    {
      int count = 0;
      int nNodes = elem_i[ind+1]-elem_i[ind];

      for (int n=0; n<nNodes; n++)
        if (nodes[elem_v[elem_i[ind]+n]].pt==boundPts[i])    count++;

      if (count>0)
      {moveElementEnd(ind); elemMoved++;}
      else ind++;
    }
  }
  nel_i = elem_i.size()-1-elemMoved; // -1 because the last i is not an element
}

void UNSTRUCTMESH::moveElementEnd(int ind)
{
  // move node to the end of the deque
  int nNodes = elem_i[ind+1]-elem_i[ind];

  // save the element into a temp variable
  deque<int> tmpElement;
  for (int i=0; i<nNodes; i++)
    tmpElement.push_back(elem_v[elem_i[ind]+i]);

  // remove element from the list
  for (int i=0; i<nNodes; i++)
    elem_v.erase(elem_v.begin() + elem_i[ind]);
  // remove element from elem_i list
  elem_i.erase(elem_i.begin()+ind);

  // update startNewElem from ind to the end
  for (int i=ind; i<elem_i.size(); i++)
    elem_i[i] -= nNodes;

  // insert the element at the end of the list
  for (int i=0; i<nNodes; i++)
    elem_v.push_back(tmpElement[i]);

  // insert elem_i for the last element
  elem_i.push_back(elem_i[elem_i.size()-1]+nNodes);

  // updating elements in the faces
  for (int i=0; i<faces.size(); i++)
  {
    if (faces[i].elem_1==ind+1)      faces[i].elem_1 = elem_i.size()-1;
    else if (faces[i].elem_1>ind+1)  faces[i].elem_1--;

    if (faces[i].elem_2==ind+1)      faces[i].elem_2 = elem_i.size()-1;
    else if (faces[i].elem_2>ind+1)  faces[i].elem_2--;
  }

}

void UNSTRUCTMESH::removeDoublePoints()
{

  // Based on how meshes are added, double points can be on the boundaries only
  for (int i=nno_i; i<nodes.size(); i++)
    for (int j=i+1; j<nodes.size(); j++)
      if (nodes[j].pt == nodes[i].pt)
      {
        //Erasing double nodes
        nodes.erase(nodes.begin()+j);

        //Updating elements and faces

        // change nodes index in the boundary elements
        for (int e=elem_i[nel_i]; e<elem_v.size(); e++)		//elem_i[nel_i]
          if (elem_v[e]==j)    elem_v[e] = i;
          else if (elem_v[e]>j) elem_v[e]--;

        // change nodes index in the boundary faces
        for (int f=nfa_nno_i; f<faces.size(); f++)
          for(int n=0; n<faces[f].node.size(); n++)
            if (faces[f].node[n]==j)    	faces[f].node[n] = i;
            else if (faces[f].node[n]>j) 	faces[f].node[n]--;

        break;
      }

  // nodes and elements that are now internal should be put outside the boundary list
}

void UNSTRUCTMESH::Store_NodesElemFaces(const string name)
{
  //Function to store the nodes, elements and face information
  //      Mostly used for debugging

  //***********
  //Nodes

  char filename1[32];
  sprintf(filename1, "%s_nodes.dat", name.c_str());

  FILE *fp1 = fopen(filename1, "wt");
  for (int i=0; i<nodes.size(); i++)
  {
    int num_neigh= nodes[i].neig.size();
    fprintf(fp1,"%d\t%f\t%f\t%f\t%d\n", i+1, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z, num_neigh);
  }
  fclose(fp1);

  //***********
  //Elements

  char filename2[32];
  sprintf(filename2, "%s_elements.dat", name.c_str());

  FILE *fp2 = fopen(filename2, "wt");

  for (int i=0; i<elem_i.size()-1; i++)
  {
    //    int size_elements=elem_i[i+1]-elem_i[i];
    fprintf(fp2,"%d\t%d", i+1, elem_i[i]);
    for (int n=elem_i[i]; n<elem_i[i+1]; n++)
      fprintf(fp2, "%8d", elem_v[n]+1);

    fprintf(fp2, "\n");
  }
  fclose(fp2);

  //***********
  //Faces

  char filename3[32];
  sprintf(filename3, "%s_face.dat", name.c_str());

  FILE *fp_ = fopen(filename3, "wt");
  for (int i=0; i<faces.size(); i++)
  {
    int size_faces=faces[i].node.size();
    fprintf(fp_,"%d\t%s", i+1, faces[i].name);

    for (int j=0; j<size_faces; j++)
      fprintf(fp_, "%8d", faces[i].node[j]+1);

    fprintf(fp_, "%8d\t%8d", faces[i].elem_1, faces[i].elem_2);
    fprintf(fp_, "\n");
  }
  fclose(fp_);

}

void UNSTRUCTMESH::UpdateBoundaries()
{

  //----------------------------------
  //	Moving Boundary NODES, ELEMENTS, FACES to the end of the array

  // The update of the boundary nodes, elements and faces are done with the knowledge of the faces: icv0 and icv1
  // The control volumes that have icv1 (or in this case elem_2) equal to 0 (zero) then it means it is in the boundary.

  //****************************
  //Node
  deque<int>	move_node;	    // This deque will store the nodes that were already moved, so the same node is not moved twice
  move_node.push_back(-1);      // Initializing the array

  int count=0;

  for (int i=nfa_i; i<faces.size(); i++)
    if(faces[i].elem_2==0)
    {
      int ind;
      //...... Nodes of the face
      for(int j=0; j<faces[i].node.size(); j++)
      {
        ind=faces[i].node[j];
        count=0;
        //Checking if the node was already moved
        for(int j=0; j<move_node.size(); j++)
          if (move_node[j]==ind)
            count++;
        if(count==0)
        {
          //Moving node to the end (this function also updated elements and faces)
          moveNodeEnd(ind);

          // update nodes in move_node deque
          for(int k=0; k<move_node.size(); k++)
            if (move_node[k]>ind)  move_node[k]--;
          //Adding the node to the move_node deque
          move_node.push_back(nodes.size()-1);
        }
      }
    }

  nno_i=move_node[1];   //First node moved is the first non internal node.
  move_node.clear();


  //****************************
  //Element
  deque<int> 	move_elem_i;			// This deque will store the elements that were already moved, so the same element is not moved twice
  move_elem_i.push_back(-1);            // Initializing the array

  for (int i=nfa_i; i<faces.size(); i++)
    if(faces[i].elem_2==0)
    {

      int ind=faces[i].elem_1-1;
      count=0;
      //Checking if the element was already moved
      for(int j=0; j<move_elem_i.size(); j++)
        if (move_elem_i[j]==ind)
          count++;

      if(count==0)
      {
        //Moving element to the end (this function also updates faces)
        moveElementEnd(ind);

        // update elements in move_elem_i deque
        for(int k=0; k<move_elem_i.size(); k++)
          if (move_elem_i[k]>ind)  move_elem_i[k]--;
        //Adding the node to the move_elem_i deque
        move_elem_i.push_back(elem_i.size()-2);

      }
    }
  nel_i=move_elem_i[1];   //First element moved is the first non internal element.
  move_elem_i.clear();


  //****************************
  //Faces

  // Moving the internal faces that have only internal nodes
  for (int i=nfa_nno_i; i<nfa_i; i++)
  {
    count=0;
    for(int j=0; j<faces[i].node.size(); j++)
      if(faces[i].node[j]<nno_i)
        count++;

    if(count == faces[i].node.size())
    {
      deque<FACE> tmp_face;
      tmp_face.push_back(faces[i]);
      faces.erase(faces.begin()+i);
      faces.insert(faces.begin()+nfa_nno_i,tmp_face.begin(), tmp_face.end() );
      nfa_nno_i++;
    }

  }

  // Moving the internal faces that have external nodes
  for (int i=nfa_i; i<faces.size(); i++)
    if(faces[i].elem_2!=0)
    {
      deque<FACE> tmp_face;
      tmp_face.push_back(faces[i]);
      faces.erase(faces.begin()+i);

      count=0;
      for(int j=0; j<tmp_face[0].node.size(); j++)
        if(tmp_face[0].node[j]>=nno_i)
          count++;

      if(count>0)
        faces.insert(faces.begin()+nfa_i,tmp_face.begin(), tmp_face.end() );
      else
      {
        faces.insert(faces.begin()+nfa_nno_i,tmp_face.begin(), tmp_face.end() );
        nfa_nno_i++;
      }
      nfa_i++;

    }

}


void UNSTRUCTMESH::findNodeNeighbors()
{
  for (int i=0; i<faces.size(); i++)
  {
    nodes[faces[i].node[0]].neig.push_back(faces[i].node[1]);
    nodes[faces[i].node[1]].neig.push_back(faces[i].node[0]);
  }

  // sort nodes indexes and remove duplicates
  for (int i=0; i<nodes.size(); i++)
  {
    nodes[i].neig.sort();
    nodes[i].neig.unique();
    nodes[i].neig.remove(i);
  }

  // there was an old function that didn't use the faces to find the neighbors...
  // find it in Code backups before October 2017

}

void UNSTRUCTMESH::findFaces2D()          //GUSTAVO 09-05-2016
{
  //If the mesh is 2D, then the faces consists of only two nodes
  if(faces.size()==0)
  {
    for(int i=0; i<elem_i.size()-1; i++)
    {
      int fstNodeInd = elem_i[i];
      int lstNodeInd = elem_i[i+1]-1;

      for(int j=fstNodeInd; j<lstNodeInd; j++)
      {
        FACE tmp_face;
        tmp_face.node.push_back(elem_v[j+1]);
        tmp_face.node.push_back(elem_v[j]);
        tmp_face.elem_1 = i+1;
        tmp_face.elem_2 = 0;
        sprintf(tmp_face.name,"noname");
        faces.push_back(tmp_face);
      }

      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.node.push_back(elem_v[lstNodeInd]);
      tmp_face.elem_1 = i+1;
      tmp_face.elem_2 = 0;
      sprintf(tmp_face.name,"noname");
      faces.push_back(tmp_face);
  }
    nfa_i=0; nfa_nno_i=0;

    removeDoubleFaces();

    //-----------------------------------
    // Moving faces to the back
    //      internal faces, internal faces with boundary nodes, and boundary faces
    int nfaces = faces.size();
    int faceMoved = 0;
    deque<FACE> nfa_nno_iface;
    deque<FACE> nfa_i_face;


    for(int i=0; i<nfaces-faceMoved; i++)
    {
      if(faces[i].elem_2==0)
      {
        FACE tmp= faces[i];
        faces.erase(faces.begin() + i);
        nfa_i_face.push_back(tmp);

        faceMoved++;
        i--;

      }
      else if((faces[i].node[0]>=nno_i || faces[i].node[1]>=nno_i) && (faces[i].elem_2!=0))
      {
        FACE tmp= faces[i];
        faces.erase(faces.begin() + i);
        nfa_nno_iface.push_back(tmp);

        faceMoved++;
        i--;
      }
    }


    nfa_i=nfaces-nfa_i_face.size();
    nfa_nno_i=nfa_i-nfa_nno_iface.size();

    faces.insert(faces.end()    ,nfa_nno_iface.begin()  , nfa_nno_iface.end() );
    faces.insert(faces.end()    ,nfa_i_face.begin()     , nfa_i_face.end() );


    // Checking that the number of faces was mantained
    if(faces.size() - nfaces != 0 )
    {
      printf("\n Number of faces are wrong \n");
      throw(-1);
    }

    nfa_i_face.clear();
    nfa_nno_iface.clear();

  }
  else
    cout<<"Inside: findFaces2D(): The objects faces were already created!" <<endl;
}



void UNSTRUCTMESH::setMarker(const deque<POINT> &points, int marker, bool boundarypts)
{

  //If points are in the boundary then the search domain is much smaller
  int n_start=0;
  if(boundarypts)    n_start=nno_i;

  for (int i=0; i<points.size(); i++)
    for (int j=n_start; j<nodes.size(); j++)
      if (nodes[j].pt == points[i])
      {
        nodes[j].marker = marker;
        break;
      }
}


void UNSTRUCTMESH::removeDoubleFaces()
{
  // Based on how meshes are added, double faces can be on the boundaries only

  for (int i=nfa_i; i<faces.size(); i++)
    for (int j=i+1; j<faces.size(); j++)
      if (compareFaces(faces[i],faces[j]))
      {
        faces[i].elem_2 = faces[j].elem_1;
        sprintf(faces[i].name,"fluid");
        faces.erase(faces.begin()+j);
        break;
      }


}


bool UNSTRUCTMESH::compareFaces(const FACE &f1, const FACE &f2)
{
  if(f1.node.size()<4)
    return ((f1.node[0] == f2.node[0] && f1.node[1] == f2.node[1]) || (f1.node[0] == f2.node[1] && f1.node[1] == f2.node[0]));
  else
  {
    for(int i=0; i<f2.node.size(); i++)
      if(f1.node[0] == f2.node[i])
      {
        if(i==0)
          if ((f1.node[1] == f2.node[i+1] && f1.node[2] == f2.node[i+2])  || (f1.node[1] == f2.node[3]    && f1.node[2] == f2.node[2]))
            return true;


        if(i==1)
          if ((f1.node[1] == f2.node[i+1] && f1.node[2] == f2.node[i+2])  || (f1.node[1] == f2.node[i-1]  && f1.node[2] == f2.node[3]))
            return true;

        if(i==2)
          if ((f1.node[1] == f2.node[i+1] && f1.node[2] == f2.node[0])    || (f1.node[1] == f2.node[i-1]  && f1.node[2] == f2.node[i-2]))
            return true;

        if(i==3)
          if ((f1.node[1] == f2.node[0]   && f1.node[2] == f2.node[1])    || (f1.node[1] == f2.node[i-1]  && f1.node[2] == f2.node[i-2]))
            return true;

      }

    return false;

  }

}



void UNSTRUCTMESH::replacePoints(const deque<POINT> &oldPoints, const deque<POINT> &newPoints)
{
  for (int i=0; i<newPoints.size(); i++)
    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].pt == oldPoints[i])
        nodes[j].pt = newPoints[i];
}



// Laplacian smoother, simplest
void UNSTRUCTMESH::smoothMesh(int maxit, double tol)
{
  deque<POINT> temp_points(nodes.size());
  for (int j=0; j<nodes.size(); j++)
    temp_points[j] = nodes[j].pt;

  int n = 0;
  double err = 1.0;

  while ((n<maxit) && (err>tol))
  {
    double diff = 0.0;

    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].marker != 1)
      {
        POINT pt(0.0, 0.0, 0.0);
        for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
        POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;
        //          if (mag(delta)>diff)
          //              diff = mag(delta);
        diff += mag(delta);
        temp_points[j] = nodes[j].pt + 0.6*delta;
      }

    n++;
    err = diff;

    if (n%100 == 0)        printf("laplace iter = %d\t %.4le\n", n, err);

    // update all the points x and y
    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].marker != 1)
        nodes[j].pt = temp_points[j];
  }

  if (err>tol)
    printf("smoother not converged \n");
  printf("iter = %d, err = %1.5le\n", n, err);
}


// Laplacian smoother upgrade1, movement in x and y
void UNSTRUCTMESH::smoothMesh_(int maxit, double tol, double relax)
{
  // Relax is the relaxation coefficient which by default is 0.6. However, the computation could get unstable and a smaller relaxation coefficient would be better
  // Important in this new function are the markers!
  //    Market 0: nodes that can be moved freely
  //    Marker 1: nodes that shouldn't be move. For example boundary layer nodes, or boundary nodes: corners for example
  //    Market 2: nodes in the boundary in y! They can be moved but only in x coordinate.
  //    Market 3: nodes in the boundary in x! They can be moved but only in y coordinate.
  deque<POINT> temp_points(nodes.size());
  for (int j=0; j<nodes.size(); j++)
    temp_points[j] = nodes[j].pt;

  int n = 0;
  double err = 1.0;


  while ((n<maxit) && (err>tol))
  {
    double diff = 0.0;

    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].marker != 1)
      {
        POINT pt(0.0, 0.0, 0.0);
        for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
        POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;
        //          if (mag(delta)>diff)
        //              diff = mag(delta);
        if(nodes[j].marker == 0)
        {
          diff += mag(delta);
          temp_points[j] = nodes[j].pt + relax*delta;
        }
        else if(nodes[j].marker == 2)       // movement of the boundary point only on that line: in this case Y
        {
          diff += delta.x;
          temp_points[j].x = nodes[j].pt.x + relax*delta.x;
        }
        else if(nodes[j].marker == 3)       // movement of the boundary point only on that line: in this case X
        {
          diff += delta.y;
          temp_points[j].y = nodes[j].pt.y + relax*delta.y;
        }
      }

    n++;
    err = diff;

    if (n%100 == 0)        printf("laplace iter = %d\t %.4le\n", n, err);

    // update all the points x and y
    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].marker != 1)
        nodes[j].pt = temp_points[j];
  }

  if (err>tol)
    printf("smoother not converged \n");
  printf("iter = %d, err = %1.5le\n", n, err);
}



// Laplacian smoother upgrade2, movement in radius, angle, splines and same movement in 2 splines
void UNSTRUCTMESH::smoothMesh_options(int maxit, double tol , SPLINE spline_A, SPLINE spline_B, double relax)
 {
   // Relax is the relaxation coefficient which by default is 0.6. However, the computation could get unstable and a smaller relaxation coefficient would be better
   // Important in this new function are the markers!
   //    Market 0: nodes that can be moved freely
   //    Marker 1: nodes that shouldn't be move. For example boundary layer nodes, or boundary nodes: corners for example
   //    Market 2: nodes in the boundary in constant y!            They can be moved but only in x direction.
   //    Market 3: nodes in the boundary in constant x!            They can be moved but only in y direction.
   //    Market 4: nodes in the boundary in constant radius!       They can be moved but only in tangential direction.
   //    Market 5: nodes in the boundary in constant angle!        They can be moved but only in radial direction.
   //    Market 6: nodes in the boundary that consist of splineA!  They can be moved but only inside the spline.
   //    Market 7: nodes in the boundary that consist of splineB!  They can be moved but only inside the spline.
   //    Market 8: nodes in the boundary that consist of splineB!  They can be moved but only in the splineB with the same offset as splineA.
  //                Only applicable for two periodic boundaries coupled to each other
   deque<POINT> temp_points(nodes.size());
   for (int j=0; j<nodes.size(); j++)
     temp_points[j] = nodes[j].pt;

   int n = 0;
   double err = 1.0;
   int n_spline=1000;

   SPLINE splineA= spline_A;
   SPLINE splineB= spline_B;


   while ((n<maxit) && (err>tol))
   {
     double diff = 0.0;

     //To move the same of set the other periodic boundary
     deque <POINT> deltaA_vec;
     deque <double> offsetA_vec;

     //To had the weight of the other side of the periodic boundary
     deque <POINT> deltaB_vec;
     deque <double> offsetB_vec;

     //THe moved points in the periodic boundary need to generate a need spline
     deque <POINT> splineA_sort;
     deque <POINT> splineB_sort;

     deque <POINT> splineA_vec;
     deque <POINT> splineB_vec;

     //To sort
     deque <double> offset_sort;

     for (int j=nno_i; j<nodes.size(); j++)
       if ((nodes[j].marker == 8))
       {
         POINT pt(0.0, 0.0, 0.0);
         for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
         POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;
         double offset= splineB.findS(nodes[j].pt);
         deltaB_vec.push_back(delta);
         offsetB_vec.push_back(offset);
       }

     for (int j=0; j<nodes.size(); j++)
       if (nodes[j].marker != 1)
       {
         POINT pt(0.0, 0.0, 0.0);
         for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
         POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;


         if(nodes[j].marker == 0)
         {
           diff += mag(delta);
           temp_points[j] = nodes[j].pt + relax*delta;
         }
         else if(nodes[j].marker == 2)       // movement of the boundary point only on that line: in this case Y
         {
           diff += delta.x;
           temp_points[j].x = nodes[j].pt.x + relax*delta.x;
         }
         else if(nodes[j].marker == 3)       // movement of the boundary point only on that line: in this case X
         {
           diff += delta.y;
           temp_points[j].y = nodes[j].pt.y + relax*delta.y;
         }
         else if(nodes[j].marker == 4)       // movement of the boundary point only on that line: in this case rotating a certain angle
         {
           POINT tmp_point_mv=nodes[j].pt + relax*delta;
           double angle=  atan2(tmp_point_mv.y, tmp_point_mv.x); //tmp_point_mv.phiRadZ();
           double radius = sqrt(nodes[j].pt.x*nodes[j].pt.x+nodes[j].pt.y*nodes[j].pt.y);//mesh.nodes[j].pt.radZ();
           temp_points[j].x = radius*cos(angle);
           temp_points[j].y = radius*sin(angle);

           diff += mag(temp_points[j]-nodes[j].pt);

         }
         else if(nodes[j].marker == 5)       // movement of the boundary point only on that line: in this case move it a radius
         {
           POINT tmp_point_mv=nodes[j].pt + relax*delta;
           double angle= atan2(nodes[j].pt.y, nodes[j].pt.x);
           double radius =  sqrt(tmp_point_mv.x*tmp_point_mv.x+tmp_point_mv.y*tmp_point_mv.y);
           temp_points[j].x = radius*cos(angle);
           temp_points[j].y = radius*sin(angle);

           diff += mag(temp_points[j]-nodes[j].pt);

         }
         else if(nodes[j].marker == 6)
         {
           //Periodic 1

           double offset= splineA.findS(nodes[j].pt);

//          Projection!!!!!
//           POINT p1=splineA.calcPoint(offset);
//           POINT p2=splineA.calcPoint(offset+0.00001);
//           POINT tangent= (p2-p1)*(1/0.00001);
//           double proj= ((tangent.x*delta.x)+(tangent.y*delta.y))/mag(tangent);
//           double delta_s=proj/splineA.length;
//
//           diff += mag(delta);
//           temp_points[j] = splineA.calcPoint( offset+ relax*delta_s);
//           //cout<<"Difference between the offset and actual point "<< mag(nodes[j].pt - spline4.calcPoint(offset))<< " and the offsetS is: "<< offset << " the point is x: "<< nodes[j].pt.x<< " y: "<< nodes[j].pt.y << " Number of iterations: "<< nn<< " tolerance in the newton solver "<< toln <<endl;
//           offset_vec.push_back(offset);
//           ds_vec.push_back(delta_s);
//           delta_vec.push_back(delta);

           //*******************************************************88
           int ind=0;
           while(fabs(offsetB_vec[ind]-offset)>0.000001 && ind<offsetB_vec.size()) ind++;

           //Average movement between the two periodic boundaries
           POINT delta_avg = (deltaB_vec[ind]+delta)*0.5;

           //Storing offset and movement to use it in the other spline points
           offsetA_vec.push_back(offset);
           deltaA_vec.push_back(delta_avg);


           diff += mag(delta_avg);
           temp_points[j] = nodes[j].pt + relax*delta_avg;

           splineA_vec.push_back(temp_points[j]);

         }
         else if(nodes[j].marker == 7)
         {
           double offset= splineB.findS(nodes[j].pt);
           POINT p1=splineB.calcPoint(offset);
           POINT p2=splineB.calcPoint(offset+0.00001);

           POINT tangent= (p2-p1)*(1/0.00001);
           double proj= ((tangent.x*delta.x)+(tangent.y*delta.y))/mag(tangent);
           double delta_s=proj/splineB.length;

           diff += mag(delta);
           temp_points[j] = splineB.calcPoint( offset+ relax*delta_s);
         }
       }

     offsetB_vec.clear();
     offset_sort=offsetA_vec;
     sort(offset_sort.begin(), offset_sort.end());
//     offset_sort.sort();

     for (int j=nno_i; j<nodes.size(); j++)
       if (nodes[j].marker == 8)
       {
         //Periodic 2

         double offset= splineB.findS(nodes[j].pt);
         int ind=0;
         while(fabs(offsetA_vec[ind]-offset)>0.000001 && ind<offsetA_vec.size()) ind++;

         diff += mag(deltaA_vec[ind]);
         temp_points[j] = nodes[j].pt + relax*deltaA_vec[ind];

         offsetB_vec.push_back(offset);
         splineB_vec.push_back(temp_points[j]);

       }
     //**************************************************************************
     //Generating the new splines on the periodic boundary.
     //Corner!
     splineA_sort.push_back(splineA.calcPoint(0));
     splineB_sort.push_back(splineB.calcPoint(0));
     for(int j=0; j<offset_sort.size() ;j++)
     {
       int indA=0;
       while(fabs(offset_sort[j]-offsetA_vec[indA])>0.000001 && indA<offsetA_vec.size()) indA++;
       splineA_sort.push_back(splineA_vec[indA]);

       int indB=0;
       while(fabs(offset_sort[j]-offsetB_vec[indB])>0.000001 && indB<offsetB_vec.size()) indB++;
       splineB_sort.push_back(splineB_vec[indB]);
     }

     //Corner!
     splineA_sort.push_back(splineA.calcPoint(1));
     splineB_sort.push_back(splineB.calcPoint(1));

     SPLINE tmp_splineA(splineA_sort);
     SPLINE tmp_splineB(splineB_sort);

     splineA= tmp_splineA;
     splineB= tmp_splineB;

     n++;
     err = diff;

     if (n%100 == 0)        printf("laplace iter = %d\t %.4le\n", n, err);

     // update all the points x and y
     for (int j=0; j<nodes.size(); j++)
       if (nodes[j].marker != 1)
         nodes[j].pt = temp_points[j];
   }

   if (err>tol)
     printf("smoother not converged \n");
   printf("iter = %d, err = %1.5le\n", n, err);
 }


void UNSTRUCTMESH::smoothMesh_generic(int maxit, double tol , deque<SPLINE> spline, double relax ,  bool periodic_splines)
{
  // Relax is the relaxation coefficient which by default is 0.6. However, the computation could get unstable and a smaller relaxation coefficient would be better
  // Important in this new function are the markers!
  //    Market 0: nodes that can be moved freely
  //    Marker 1: nodes that shouldn't be move. For example boundary layer nodes, or boundary nodes: corners for example
  //    Market 2: nodes in the boundary in constant radius!       They can be moved but only in tangential direction.
  //               Note: This marker was added because spline in tangential change was giving poor results
  //    Market larger: nodes in the boundary that consist of a spline!  They can be moved but only along the spline.
  //                   The marker value and the spline array must follow the same sequence. e.g. Marker 2 must correlate the the first spline in the spline array.
  //                   For the periodic splines: the two spline that need to have "periodic" nodes must be in the end of the spline array!!
  //                                             Moreover, both spline needs to be equal: in terms of length and control points.
  //                Only applicable for two periodic boundaries coupled to each other
  deque<POINT> temp_points(nodes.size());
  for (int j=0; j<nodes.size(); j++)
    temp_points[j] = nodes[j].pt;

  //For the periodic splines!!!
  int per_splines=spline.size()+3;
  if(periodic_splines) per_splines--;

  int n = 0;
  double err = 1.0;
  cout<<"Smoothing function: "<<per_splines <<endl;
  while ((n<maxit) && (err>tol))
  {
    double diff = 0.0;
    deque <double> offset_vec;
    deque <double> ds_vec;

    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].marker != 1)
      {
        POINT pt(0.0, 0.0, 0.0);
        for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
        POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;
        //          if (mag(delta)>diff)
        //              diff = mag(delta);
        if(nodes[j].marker == 0)
        {
          diff += mag(delta);
          temp_points[j] = nodes[j].pt + relax*delta;
        }
        else if(nodes[j].marker == 2)       // movement of the boundary point only on that line: in this case rotating a certain angle
        {
          POINT tmp_point_mv=nodes[j].pt + relax*delta;
          double angle=  atan2(tmp_point_mv.y, tmp_point_mv.x); //tmp_point_mv.phiRadZ();
          double radius = sqrt(nodes[j].pt.x*nodes[j].pt.x+nodes[j].pt.y*nodes[j].pt.y);//mesh.nodes[j].pt.radZ();
          temp_points[j].x = radius*cos(angle);
          temp_points[j].y = radius*sin(angle);

          diff += mag(temp_points[j]-nodes[j].pt);

        }
        else if((nodes[j].marker > 2) && (nodes[j].marker<per_splines))     //Movement of the node only along the spline!
        {
          int spl= nodes[j].marker-3;                      //spline location in the array must follow the same sequence as the marker value!
          double offset= spline[spl].findS(nodes[j].pt);

          //Projecting the moved point to the spline
          POINT p1=spline[spl].calcPoint(offset);
          POINT p2=spline[spl].calcPoint(offset+0.00001);

          POINT tangent= (p2-p1)*(1/0.00001);
          double proj= ((tangent.x*delta.x)+(tangent.y*delta.y))/mag(tangent);
          double delta_s=proj/spline[spl].length;

          diff += mag(delta);
          temp_points[j] = spline[spl].calcPoint(offset+ relax*delta_s);

          if(periodic_splines && (spl==spline.size()-2))
          {
            offset_vec.push_back(offset);
            ds_vec.push_back(delta_s);
          }

        }
      }
    // Moving a point in a periodic spline equal to its corresponding periodic spline!!!
    if(periodic_splines)
      for (int j=nno_i; j<nodes.size(); j++)
        if (nodes[j].marker == per_splines)
        {
//          cout<<"Inside the periodic spline thingy"<<endl;
          int spl= nodes[j].marker-3;                      //spline location in the array must follow the same sequence as the marker value!
          double offset= spline[spl].findS(nodes[j].pt);
          int ind=0;
          while(fabs(offset_vec[ind]-offset)>0.000001 && ind<offset_vec.size()) ind++;

          diff += fabs(ds_vec[ind])*spline[spl].length;
          temp_points[j] = spline[spl].calcPoint( offset_vec[ind]+ relax*ds_vec[ind]);
        }



    n++;
    err = diff;

    if (n%100 == 0)        printf("laplace iter = %d\t %.4le\n", n, err);

    // update all the points x and y
    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].marker != 1)
        nodes[j].pt = temp_points[j];
  }

  if (err>tol)
    printf("smoother not converged \n");
  printf("iter = %d, err = %1.5le\n", n, err);
}

void UNSTRUCTMESH::smoothMesh_freeSplines(int maxit, double tol , SPLINE original_spline_A, SPLINE original_spline_B, double relax, bool perSp_CYL)
  {


    // Relax is the relaxation coefficient which by default is 0.6. However, the computation could get unstable and a smaller relaxation coefficient would be better
    // Important in this new function are the markers!
    //    Market 0: nodes that can be moved freely
    //    Marker 1: nodes that shouldn't be move. For example boundary layer nodes, or boundary nodes: corners for example
    //    Market 2: nodes in the boundary in constant angle!        They can be moved but only in tangential or vertical direction (depending on how is the boundary). Last and first point are connect to the periodic spline
    //    Market 3: nodes in the boundary that consist of splineA!  They can be moved freely in space. Takes into account the laplacian of the corresponding point on the other periodic boundary
    //    Market 4: nodes in the boundary that consist of splineB!  They can be exactly the same is the once in Market 3 (splineA).
  //        Spline A and B are updated each time!

    deque<POINT> temp_points(nodes.size());
    for (int j=0; j<nodes.size(); j++)
      temp_points[j] = nodes[j].pt;

    int n = 0;
    double err = 1.0;

    //-------------------------------------------------------------------------------
    //Initializing the smoothing for the periodic boundaries (splines)
    SPLINE splineA= original_spline_A;
    SPLINE splineB= original_spline_B;

    double perSp_angle;
    double perSp_mov_y;

    //
    //For the movement of the periodic splines:
    if(perSp_CYL)   // For the case that the periodic boundaries are separated by an angle in Cylindrical coordinates.
    {

      POINT pA= splineA.calcPoint(1);
      POINT pB= splineB.calcPoint(1);

      perSp_angle= pA.phiDegZ() - pB.phiDegZ() ;
      cout<<"Angle between the periodic boundaries in deg: "<< perSp_angle<<endl;

      //Checking if the spline is equidistant in angle!
      pA= splineA.calcPoint(0);
      pB= splineB.calcPoint(0);
      if(fabs(perSp_angle-(pA.phiDegZ() - pB.phiDegZ()))>0.0000001)  {cout<<"Angle between the periodic splines is not ctte! "<< perSp_angle<<endl; throw(-1);}

    }
    else    // For the case that the periodic boundaries are separated by dY in Cartisian coordinates.
    {

      POINT pA= splineA.calcPoint(1);
      POINT pB= splineB.calcPoint(1);

      perSp_mov_y= pA.y - pB.y ;
      cout<<"Delta y between the periodic boundaries is: "<< perSp_mov_y<<endl;

      //Checking if the spline is equidistant in y!
      pA= splineA.calcPoint(0);
      pB= splineB.calcPoint(0);
      if(fabs(perSp_mov_y-(pA.y - pB.y))>0.0000001)  {cout<<"Delta y between the periodic splines is not ctte! "<< perSp_mov_y<<endl; throw(-1);}
    }
    //-------------------------------------------------------------------------------


    while ((n<maxit) && (err>tol))
    {
      double diff = 0.0;

      //----------------------------------------------
      // Arrays to transfer the data between the periodic boundaries
      //To move the same of set the other periodic boundary
      deque <POINT> deltaA;
      deque <double> offsetA;

      //To had the weight of the other side of the periodic boundary
      deque <POINT> deltaB_lap;
      deque <double> offsetB_lap;
      deque <double> offsetB_sp;

      //THe moved points in the periodic boundary need to generate a need spline
      deque <POINT> splineA_sp;
      deque <POINT> splineB_sp;

      deque <POINT> splineA_sort;
      deque <POINT> splineB_sort;

      //To sort
      deque <double> offsetA_sort;
      deque <double> offsetB_sort;

      POINT p_splineA;
      POINT p_splineB;
      //----------------------------------------------

      //To take into account movement of both periodic boundaries at the same time. Loop over one of the periodic boundaries
      for (int j=nno_i; j<nodes.size(); j++)
        if ((nodes[j].marker == 4))
        {
          // ****************** Periodic B

          POINT pt(0.0, 0.0, 0.0);
          for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
          POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;

          double offset= splineB.findS(nodes[j].pt);

          //Storing the information:
          deltaB_lap.push_back(delta);
          offsetB_lap.push_back(offset);
        }



      for (int j=0; j<nodes.size(); j++)
        if (nodes[j].marker != 1)
        {
          POINT pt(0.0, 0.0, 0.0);
          for (list<int>::iterator it = nodes[j].neig.begin(); it!=nodes[j].neig.end(); ++it)           pt += nodes[*it].pt;
          POINT delta =  pt/(double)nodes[j].neig.size() - nodes[j].pt;
          //          if (mag(delta)>diff)
          //              diff = mag(delta);
          if(nodes[j].marker == 0)
          {
            diff += mag(delta);
            temp_points[j] = nodes[j].pt + relax*delta;
          }
          else if(nodes[j].marker == 2)       // movement of the boundary point only on that line: in this case rotating a certain angle
          {
            if(perSp_CYL)
            {
              POINT tmp_point_mv=nodes[j].pt + relax*delta;
              double angle=  atan2(tmp_point_mv.y, tmp_point_mv.x); //tmp_point_mv.phiRadZ();
              double radius = sqrt(nodes[j].pt.x*nodes[j].pt.x+nodes[j].pt.y*nodes[j].pt.y);//mesh.nodes[j].pt.radZ();
              temp_points[j].x = radius*cos(angle);
              temp_points[j].y = radius*sin(angle);
            }
            else
            {
              temp_points[j].y = nodes[j].pt.y + relax*delta.y;
            }

            diff += mag(temp_points[j]-nodes[j].pt);

            //Storing the last point of the splines!
            if(nodes[j].pt==splineA.calcPoint(1)) p_splineA=temp_points[j];
            if(nodes[j].pt==splineB.calcPoint(1))
            {
              POINT p_b= p_splineA;
              if(perSp_CYL)             p_b.rotateDegZ(perSp_angle);
              else                      p_b.y=p_splineA.y-perSp_mov_y;
              temp_points[j]=p_b;
            }

          }
          else if(nodes[j].marker == 3)
          {
            // ****************** Periodic A

            double offset= splineA.findS(nodes[j].pt);
//            cout<<"In A"<<endl;
            int ind=0;
            while(fabs(offsetB_lap[ind]-offset)>0.000001 && ind<offsetB_lap.size()) ind++;

            //Average movement between the two periodic boundaries
            POINT delta_b= deltaB_lap[ind];
            if(perSp_CYL)            delta_b.rotateDegZ(-perSp_angle);
            POINT delta_avg = (delta_b+delta)*0.5;

            //Moving the point
            diff += mag(delta_avg);
            temp_points[j] = nodes[j].pt + relax*delta_avg;

            //Storing data
            splineA_sp.push_back( nodes[j].pt + relax*delta_avg);
            offsetA.push_back(offset);
            deltaA.push_back(delta_avg);

          }
//          else if(nodes[j].marker > 5)
//          {
////            cout<<" Inside the loop of large markers" <<endl;
//
//            if(n<40)
//            {
//              double move_limitation=(double)nodes[j].marker/100;
//              delta=delta*move_limitation;
//              diff += mag(delta);
//              temp_points[j] = nodes[j].pt + relax*delta;
//            }
//          }
        }


      for (int j=nno_i; j<nodes.size(); j++)
        if (nodes[j].marker == 4)
        {
          // ****************** Periodic B

          double offset= splineB.findS(nodes[j].pt);
          int ind=0;
          while(fabs(offsetA[ind]-offset)>0.000001 && ind<offsetA.size()) ind++;

          //    you need to rotate delta_avg!!!! it is a vector now!!!
          POINT delta_avg = deltaA[ind];
          if(perSp_CYL)            delta_avg.rotateDegZ(perSp_angle);

          //*****************************************************************
          ///           WorkingS
//          temp_points[j] = p1;  //nodes[j].pt + relax*delta_avg;
//
//          splineB_sp.push_back(p1);//nodes[j].pt + relax*delta_avg);
//          offsetB_sp.push_back(offset);
          //*****************************************************************
          diff += mag(delta_avg);
          temp_points[j] = nodes[j].pt + relax*delta_avg;

          splineB_sp.push_back(nodes[j].pt + relax*delta_avg);
          offsetB_sp.push_back(offset);


        }

      for(int j=0; j<offsetA.size();j++)
        offsetA_sort.push_back(offsetA[j]);

      sort(offsetA_sort.begin(), offsetA_sort.end());


      for(int j=0; j<offsetB_sp.size();j++)
        offsetB_sort.push_back(offsetB_sp[j]);

      sort(offsetB_sort.begin(), offsetB_sort.end());

      //**************************************************************************
      //Generating the new splines on the periodic boundary.
      //Corner!
      splineA_sort.push_back(splineA.calcPoint(0));
      splineB_sort.push_back(splineB.calcPoint(0));
      for(int j=0; j<offsetA_sort.size() ;j++)
      {
        int indA=0;
        while(fabs(offsetA_sort[j]-offsetA[indA])>0.000001 && indA<offsetA.size()) indA++;
        splineA_sort.push_back(splineA_sp[indA]);
      }
      for(int j=0; j<offsetB_sort.size() ;j++)
      {
        int indB=0;
        while(fabs(offsetB_sort[j]-offsetB_sp[indB])>0.000001 && indB<offsetB_sp.size()) indB++;
        splineB_sort.push_back(splineB_sp[indB]);
      }

      //Corner!
      splineA_sort.push_back(p_splineA);
      if(perSp_CYL)             p_splineA.rotateDegZ(perSp_angle);
      else                      p_splineA.y=p_splineA.y-perSp_mov_y;
//      p_splineA.rotateDegZ(-10);
      splineB_sort.push_back(p_splineA);
//      splineA_sort.push_back(splineA.calcPoint(1));
//      splineB_sort.push_back(splineB.calcPoint(1));

      SPLINE tmp_splineA(splineA_sort);
      SPLINE tmp_splineB(splineB_sort);

      //Updating Splines:
      splineA= tmp_splineA;
      splineB= tmp_splineB;


      //**************************************************************************
      n++;
      err = diff;

      if (n%100 == 0)
      {
        printf("laplace iter = %d\t %.4le\n", n, err);
      }

      // update all the points x and y
      for (int j=0; j<nodes.size(); j++)
        if (nodes[j].marker != 1)
          nodes[j].pt = temp_points[j];
    }

    if (err>tol)
      printf("smoother not converged \n");
    printf("iter = %d, err = %1.5le\n", n, err);
  }


void UNSTRUCTMESH::writeMarkedNodes(const char * name, int marker)
{
  FILE *fp = fopen(name, "wt");
  for (int i=0; i<nodes.size(); i++)
  {
    if (nodes[i].marker == marker)
      fprintf(fp, "%1.10le\t%1.10le\t%1.10le\n", nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  }
  fclose(fp);
}

void UNSTRUCTMESH::ScaleMesh(const double scaling_factor)
{
  for (int j=0;j<nodes.size();j++)
  {
    nodes[j].pt.x*= scaling_factor;
    nodes[j].pt.y*= scaling_factor;
    nodes[j].pt.z*= scaling_factor;
  }
}

//****************************************************************************************************************
//****************************************************************************************************************
//                                  Writing mesh files functions!!
//****************************************************************************************************************
//****************************************************************************************************************
void UNSTRUCTMESH::writeTecplot(const char *name, int ndim, bool onlyBoundary)
{
  FILE *fp;
  if ((fp = fopen(name, "wt")) == NULL)
  {
    printf("couldn't open %s for writing\n", name);
    return;
  }

  string elemtype;
  if      (ndim == 2) elemtype = "FEQUADRILATERAL";
  else if (ndim == 3) elemtype = "FEBRICK";
  else
  {
    printf("ERROR: NUMBER OF DIMENSIONS NOT CORRECT!\n");
    throw(-100);
  }

  fprintf(fp, "TITLE=\"mesh\"\n");
  fprintf(fp, "VARIABLES=\"X\" \"Y\" \"Z\"\n");
  if (onlyBoundary==true)
    fprintf(fp, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=%s\n", (int)nodes.size(), (int)(elem_i.size()-1-nel_i), elemtype.c_str());
  else    fprintf(fp, "ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=%s\n", (int)nodes.size(), (int)(elem_i.size()-1), elemtype.c_str());

  for (int i=0; i<nodes.size(); i++)
    fprintf(fp, "%1.10le\t%1.10le\t%1.10le\n", nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  fprintf(fp, "\n");


  int iStart = 0;
  if (onlyBoundary==true)
    iStart = nel_i;


  for (int i=iStart; i<elem_i.size()-1; i++)
  {
    int nNodes = elem_i[i+1]-elem_i[i];

    if (nNodes == 3)
      fprintf(fp, "%d\t%d\t%d\t%d\n", elem_v[elem_i[i]]+1, elem_v[elem_i[i]]+1, elem_v[elem_i[i]+1]+1, elem_v[elem_i[i]+2]+1);

    else if (nNodes == 4)
      fprintf(fp, "%d\t%d\t%d\t%d\n", elem_v[elem_i[i]]+1, elem_v[elem_i[i]+1]+1, elem_v[elem_i[i]+2]+1, elem_v[elem_i[i]+3]+1);

    else if (nNodes == 6)
      fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", elem_v[elem_i[i]]+1, elem_v[elem_i[i]]+1, elem_v[elem_i[i]+1]+1, elem_v[elem_i[i]+2]+1,
          elem_v[elem_i[i]+3]+1, elem_v[elem_i[i]+3]+1, elem_v[elem_i[i]+4]+1, elem_v[elem_i[i]+5]+1);

    else if (nNodes == 8)
      fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", elem_v[elem_i[i]]+1, elem_v[elem_i[i]+1]+1, elem_v[elem_i[i]+2]+1, elem_v[elem_i[i]+3]+1,
          elem_v[elem_i[i]+4]+1, elem_v[elem_i[i]+5]+1, elem_v[elem_i[i]+6]+1, elem_v[elem_i[i]+7]+1);

    else {printf("ERROR DIMENSIONS: %d\t i: %d\n", nNodes, i); throw(-333);}
  }
  fclose(fp);

  cout << " > " << name << " written." << endl;
}

void UNSTRUCTMESH::writeGambitNeu(const char *name, int ndim, bool onlyBoundary)
{
  //Function to generate *.neu mesh files (Gambit Neutral format) from an existing UNSTRUCTMESH mesh object.
  //          The *.neu mesh is element base.
  FILE *fp;
  if ((fp = fopen(name, "wt")) == NULL)
  {
    printf("couldn't open %s for writing\n", name);
    return;
  }

  //Header of the file
  fprintf(fp, "        CONTROL INFO 2.4.6\n");
  fprintf(fp, "** GAMBIT NEUTRAL FILE\n");
  fprintf(fp, "NONAME\n");
  fprintf(fp, "PROGRAM:                Gambit     VERSION:  2.4.6\n");
  fprintf(fp, " 1 Jan 2012    00:00:00\n");
  fprintf(fp, "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n");
  if (onlyBoundary==true)
    fprintf(fp, "%10d%10d         1         1%10d%10d\n", (int)nodes.size(), (int)(elem_i.size()-1-nel_i), ndim, ndim);
  else fprintf(fp, "%10d%10d         1         1%10d%10d\n", (int)nodes.size(), (int)(elem_i.size()-1), ndim, ndim);
  fprintf(fp, "ENDOFSECTION\n");

  //    int iStart = 0;
  //  if (onlyBoundary==true)
  //    iStart = nno_i;


  //Writing the node coordinates.
  fprintf(fp, "   NODAL COORDINATES 2.4.6\n");
  if (ndim == 2)
    for (int i=0; i<nodes.size(); i++)
      fprintf(fp, "%10d%20.11le%20.11le\n", i+1, nodes[i].pt.x, nodes[i].pt.y);
  else if (ndim == 3)
    for (int i=0; i<nodes.size(); i++)
      fprintf(fp, "%10d%20.11le%20.11le%20.11le\n", i+1, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  else
  {
    printf("ERROR: NUMBER OF DIMENSIONS NOT CORRECT!\n");
    throw(-100);
  }
  fprintf(fp, "ENDOFSECTION\n");

  //Writing the element connectivity.
  fprintf(fp, "      ELEMENTS/CELLS 2.4.6\n");
  int cStart = 0;
  if (onlyBoundary==true)
    cStart = nel_i;

  for (int c=cStart; c<elem_i.size()-1; c++)
  {
    int nNodes = elem_i[c+1]-elem_i[c];

    //Depending on the number of nodes per element:
    if (nNodes == 3)
    {
      fprintf(fp, "%8d%3d%3d ", c+1, NTYPE_TRI, NDP_TRI);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]+1);
      fprintf(fp, "\n");
    }
    else if (nNodes == 4)
    {
      fprintf(fp, "%8d%3d%3d ", c+1, NTYPE_QUAD, NDP_QUAD);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]+1);
      fprintf(fp, "\n");
    }

    else if (nNodes == 6)
    {
      fprintf(fp, "%8d%3d%3d ", c+1, NTYPE_PRISM, NDP_PRISM);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]+1);
      fprintf(fp, "\n");
    }

    else if(nNodes == 8)
    {
      fprintf(fp, "%8d%3d%3d ", c+1, NTYPE_BRICK, NDP_BRICK);
      fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", elem_v[elem_i[c]]+1, elem_v[elem_i[c]+1]+1, elem_v[elem_i[c]+3]+1, elem_v[elem_i[c]+2]+1,
          elem_v[elem_i[c]+4]+1, elem_v[elem_i[c]+5]+1, elem_v[elem_i[c]+7]+1, elem_v[elem_i[c]+6]+1);
    }

    else {printf("ERROR DIMENSIONS: %d\t i: %d\n", nNodes, c); throw(-333);}

  }
  fprintf(fp, "ENDOFSECTION\n");
  fclose(fp);

  printf("\n GAMBIT file %s written.\n", name);
}


void UNSTRUCTMESH::writeFluentMsh(const char *filename, int ndim, deque<string> bc_names)
{
  FILE * fp;

  printf(" > starting to write fluent case file: %s\n", filename);

  if ( (fp=fopen(filename,"w"))==NULL ) {
    printf(" > could not open file %s\n", filename);
    exit(-1);
  }

  //Header of the file
  fprintf(fp,"(0 \"Project to Fluent\")\n");
  fprintf(fp,"\n");

  //Writing the dimension of the mesh
  fprintf(fp,"(0 \"Dimension:\")\n");
  if(ndim==2)       fprintf(fp,"(2 2)\n");
  else              fprintf(fp,"(2 3)\n");

  //Extra information, Useful when reading the file by the mesh Generator.
  //                    In Gambit this will not be read
  fprintf(fp,"\n");
  fprintf(fp,"(0 \"Information:\")\n");
  fprintf(fp,"(0 (0 %x %x %x %x))\n",nno_i, nfa_nno_i, nel_i, 1000);
  fprintf(fp,"\n");
  //
  // header and global dimensions
  //

  int vertices = nodes.size(), elements = elem_i.size()-1, nfaces =  faces.size();


  fprintf(fp,"(0 \"Grid:\")\n");

  fprintf(fp,"\n");

  fprintf(fp,"(0 \"Nodes:\")\n");
  if(ndim==2)
  {
    fprintf(fp,"(10 (0 %x %x 1 2))\n",1, vertices);
    fprintf(fp, "(10 (1 %x %x 1 2)(\n",1, vertices);
    for (int i=0; i<nodes.size(); i++)
      fprintf(fp, "%20.11le%20.11le\n", nodes[i].pt.x, nodes[i].pt.y);
  }
  else
  {
    fprintf(fp,"(10 (0 %x %x 1 3))\n",1, vertices);
    fprintf(fp, "(10 (1 %x %x 1 3)(\n",1, vertices);
    for (int i=0; i<nodes.size(); i++)
      fprintf(fp, "%20.11le%20.11le%20.11le\n", nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  }
  fprintf(fp, "))\n");
  fprintf(fp,"\n");
  cout<<"Number of nodes: "<< vertices <<endl;

  //
  // start of faces
  //
  fprintf(fp,"(0 \"Faces:\")\n");
  fprintf(fp,"(13(0 %x %x 0))\n",1, nfaces);

  if(n_zone==0)
    n_zone = bc_names.size();

  int ind_faces[n_zone];      //Index of the current zone

  cout<<"Number of faces: "<< nfaces <<endl;

  //Checking if all the faces have a name!
  bool zone_noname=false;
  for(int i=0; i<nfaces; i++)
    if(strcmp (faces[i].name, "noname") == 0)
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }

  cout<<"Number of elements: "<< elements <<endl;

  if(zone_noname)
  {
    int faces_wall=0;
    for(int i=0; i<nfaces; i++)
      if(faces[i].elem_2==0) faces_wall++;

    //WALL FACES
    ind_faces[0]=3;

    fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[0] ,1, faces_wall);
    for(int i=0; i<nfaces; i++)
    {

      if(faces[i].elem_2==0)
      {
        int size_faces=faces[i].node.size();
        fprintf(fp,"%d ", size_faces);

        for (int j=0; j<size_faces; j++)
          fprintf(fp, "%x ", faces[i].node[j]+1);

        fprintf(fp, "%x %x \n", faces[i].elem_1, faces[i].elem_2);
      }
    }

    fprintf(fp, "))\n");

    //INNER FACES
    ind_faces[1]=5;
    fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[1] ,faces_wall+1, nfaces);
    for(int i=0; i<nfaces; i++)
    {
      if(faces[i].elem_2>0)
      {
        int size_faces=faces[i].node.size();
        fprintf(fp,"%d ", size_faces);

        for (int j=0; j<size_faces; j++)
          fprintf(fp, "%x ", faces[i].node[j]+1);

        fprintf(fp, "%x %x \n", faces[i].elem_1, faces[i].elem_2);
      }
    }

    fprintf(fp, "))\n");
  }
    else
    {
      cout<<"All faces have name, number of boundaries: "<< bc_names.size()<<endl;


      int nfa_i[n_zone];        //Number of faces in the current zone
      int total_faces=0;


      for(int n=0 ; n<bc_names.size(); n++)
      {
        cout<<"loop over the boundaries boundaries: "<<bc_names[n].c_str()<<endl;
        // Defining the index of the face zone
        ind_faces[n]=n+3;                                 // the index starts in 3

        //Counting the number of faces of the present face zone, needed in the declaration for the face
        nfa_i[n]=0;
        for(int i=0; i<nfaces; i++)
          if(strcmp (faces[i].name, bc_names[n].c_str()) == 0)    //(face_name[i].compare(bc_names[n])==0) //
            nfa_i[n]++;

        // Declaration of the face zone
        if(n!=bc_names.size()-1)fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[n] ,total_faces+1, total_faces+nfa_i[n]);
        else fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[n]+1 ,total_faces+1, total_faces+nfa_i[n]);

        // Giving the face connectivity
        for(int i=0; i<nfaces; i++)
          if(strcmp (faces[i].name, bc_names[n].c_str()) == 0)    //(face_name[i].compare(bc_names[n])==0) //
          {
            int size_faces=faces[i].node.size();
            fprintf(fp,"%d ", size_faces);

            for (int j=0; j<size_faces; j++)
              fprintf(fp, "%x ", faces[i].node[j]+1);

            fprintf(fp, "%x %x \n", faces[i].elem_1, faces[i].elem_2);

          }
        fprintf(fp, "))\n");

        total_faces+=nfa_i[n];
        cout<<"Number of faces in boundary "<<bc_names[n] <<" "<<nfa_i[n]<<endl;
      }

      cout<<"Total number of faces in the numerical domain "<<total_faces<<endl;

    }


  fprintf(fp,"\n");

  //
  // cells/ elements
  //

  int nNodes_old = elem_i[2]-elem_i[1];
  int cv_element_type;
  if(nNodes_old==3)
    cv_element_type   = 1;          // triangular = 1
  else if(nNodes_old==4)
    cv_element_type   = 3;          // quadrilateral = 3
  else if(nNodes_old==6)
    cv_element_type   = 6;          // wedge = 6
  else if(nNodes_old==8)
      cv_element_type   = 4;        // hexahedral = 6
  else
  {
    printf(" > elements are of unknown size %s\n", filename);
    exit(-1);
  }

  for(int i=1; i<elem_i.size()-1; i++)
  {
    int nNodes = elem_i[i+1]-elem_i[i];
    if(nNodes-nNodes_old==0)
      nNodes_old = nNodes;
    else
    {
      // not all elements are the same, therefore is mixed
      cv_element_type= 0;          // mixed = 0
      break;
    }
  }


  int cv_type = 1; //1; // active = 1 - inactive = 32
  char this_zone[128];
  sprintf(this_zone,"fluid");
  fprintf(fp,"(0 \"Cells:\")\n");
  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
  fprintf(fp,"(12 (%x %x %x %x %x)", 2,1,elements,cv_type,cv_element_type);
  if(cv_element_type==0)    //in case it is mixed elements...
  {
    fprintf(fp," (\n");
    int count=1;
    for(int i=0; i<elem_i.size()-1; i++)
    {
      int nNodes = elem_i[i+1]-elem_i[i];
      if(nNodes==3)
        fprintf(fp,"1 ");
      else if(nNodes==4)
        fprintf(fp,"3 ");
      else if(nNodes==6)
        fprintf(fp,"6 ");
      else if(nNodes==8)
        fprintf(fp,"4 ");

      if(count%8==0)
      {
        count=1;
        if(i<elem_i.size()-2)
          fprintf(fp," \n");
      }
      else
        count++;
    }
    fprintf(fp," \n");
    fprintf(fp," ))\n");
  }
  else
    fprintf(fp,")\n");

  fprintf(fp,"\n");

  //
  // boundaries
  //
  if(zone_noname)
  {
    fprintf(fp,"(0 \"Zones:\")\n");
    fprintf(fp,"(45 (2 fluid fluid 0))\n");
    fprintf(fp,"(45 (%i wall boundaries 0))\n", ind_faces[0]);
    fprintf(fp,"(45 (%i interior default-interior 0))\n", ind_faces[1]);
  }
  else
  {
    fprintf(fp,"(0 \"Zones:\")\n");
    fprintf(fp,"(45 (2 fluid fluid 0))\n");
    for(int n=0 ; n<bc_names.size()-1; n++)
    {
      fprintf(fp,"(45 (%i wall %s 0))\n", ind_faces[n], bc_names[n].c_str());
    }
    fprintf(fp,"(45 (%i interior default-interior 0))\n", ind_faces[bc_names.size()-1]+1);
  }

  fprintf(fp,"\n");

  fclose(fp);

  printf(" > Fluent case file written: %s\n", filename);
}

void UNSTRUCTMESH::writeSU2(const char *name, int ndim, deque<string> bc_names)
{
  FILE *fp;

  printf(" > starting to write SU2 mesh file: %s\n", name);


  if ((fp = fopen(name, "wt")) == NULL)
  {
    printf("couldn't open %s for writing\n", name);
    return;
  }

  //Number of
  int num_elem =elem_i.size()-1;    // elements
  int num_nodes=nodes.size();       // elements
  int num_faces=faces.size();    // boundaries
  int num_bc   =bc_names.size();    // boundaries

  // Types of elementes
  int triangle        = 5;     // 2D triangle
  int quadrilateral   = 9;     // 2D quadrilateral
  int wedge           = 13;    // triangle extruded
  int hexahedral      = 12;    // quadrilateral extruded



  // Dimension of the mesh
  fprintf(fp,"NDIME= %d\n",ndim);

  //**********************************
  // Writing elements


  fprintf(fp,"NELEM= %d\n",num_elem);
  for (int c=0; c<num_elem; c++)
  {
    int nNodes = elem_i[c+1]-elem_i[c];
    if (nNodes == 3)
    {
      fprintf(fp, "%3d ", triangle);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]);
      fprintf(fp, "%8d ", c);
      fprintf(fp, "\n");
    }
    else if (nNodes == 4)
    {
      fprintf(fp, "%3d ", quadrilateral);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]);
      fprintf(fp, "%8d ", c);
      fprintf(fp, "\n");
    }

    else if (nNodes == 6)
    {
      fprintf(fp, "%3d ", wedge);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]);
      fprintf(fp, "%8d ", c);
      fprintf(fp, "\n");
    }

    else if(nNodes == 8)
    {
      fprintf(fp, "%3d ", hexahedral);
      for (int n=elem_i[c]; n<elem_i[c+1]; n++)
        fprintf(fp, "%8d", elem_v[n]);
      fprintf(fp, "%8d ", c);
      fprintf(fp, "\n");
    }

    else {printf("ERROR DIMENSIONS: %d\t i: %d\n", nNodes, c); throw(-333);}


  }

  //**********************************
  // Writing nodes

  if (ndim == 2)
  {
    fprintf(fp,"NPOIN= %d\n",num_nodes);
    for (int i=0; i<num_nodes; i++)
      fprintf(fp, "%20.11le%20.11le%10d\n",          nodes[i].pt.x, nodes[i].pt.y, i);
  }
  else if (ndim == 3)
  {
    fprintf(fp,"NPOIN= %d\n",num_nodes);
    for (int i=0; i<num_nodes; i++)
      fprintf(fp, "%20.11le%20.11le%20.11le%10d\n",  nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z, i);
  }
  else
  {
    printf("ERROR: NUMBER OF DIMENSIONS NOT CORRECT!\n");
    throw(-100);
  }

  //**********************************
  // Writing boundaries, face boundaries....
  fprintf(fp,"NMARK= %d\n",num_bc); //, deque<string> bc_names
  for (int i=0; i<num_bc; i++)
  {
    fprintf(fp,"MARKER_TAG= %s\n",bc_names[i].c_str());

    //Counting the number of faces of the present face zone, needed in the declaration for the face
    int nfaces = 0;      //Number of faces in the current zone
    for(int j=0; j<num_faces; j++)
      if(strcmp (faces[j].name, bc_names[i].c_str()) == 0)    //(face_name[i].compare(bc_names[n])==0) //
        nfaces++;

    printf("Number of faces in %s is: %d !\n ",bc_names[i].c_str() ,nfaces);

    fprintf(fp,"MARKER_ELEMS= %8d\n",nfaces); // number of elements in the boundary, in reality is the face!
    for(int j=nfa_i; j<faces.size(); j++)
      if(strcmp (faces[j].name, bc_names[i].c_str()) == 0)    //(face_name[i].compare(bc_names[n])==0) //
      {
        int faces_size = faces[j].node.size();

        if      (faces_size==2)         fprintf(fp,"3 %8d %8d\n",        faces[j].node[0], faces[j].node[1]);                                     //line
        else if (faces_size == 3)       fprintf(fp,"5 %8d %8d %8d \n",   faces[j].node[0], faces[j].node[1], faces[j].node[2]);                   //triangle
        else if (faces_size == 4)       fprintf(fp,"9 %8d %8d %8d %8d\n",faces[j].node[0], faces[j].node[1], faces[j].node[2], faces[j].node[3]); //quadrilateral
        else
        {
          printf("ERROR: NUMBER OF DIMENSIONS NOT CORRECT!\n");
          throw(-100);
        }

      }
  }
  fclose(fp);
  printf(" > SU2 mesh file written: %s\n", name);

}

//****************************************************************************************************************
//****************************************************************************************************************
//                                  Reading mesh files functions!!
//****************************************************************************************************************
//****************************************************************************************************************

void UNSTRUCTMESH::readFluentMsh2D(const char *name, bool writtenMeshGen)
{

  FILE *fp;
  if ((fp = fopen(name, "rt")) == NULL)
  {
    printf("couldn't open %s for reading\n", name);
    exit(-1);
  }
  else printf("***** Opening %s for reading\n", name);

  int dummy1, dummy2, dummy3;
  int ndim;
  int nvert, nfaces_bc, nfaces_fluid, nelem, nfaces;
  string str1, str2, str3;
  char char1[64], char2[64], char3[64];
  bool wittenGambit=true;

  //**************************************
  //****** Read header
  char line[1000];
  fscanf(fp, "%s\n", line);


  //**************************************
  //****** Reading the dimensions of the mesh!!
  int n=0;
  while ((strcmp(line, "\"Dimension:\")") != 0) && (n<100000))    {fscanf(fp, "%s\n", line); n++;}
//  fscanf(fp, "%s %d %s \n", str1.c_str(), &ndim, str3.c_str());
  fscanf(fp, "%s %d %s \n", char1, &ndim, char3);
  //Checking if it is a 2D mesh!
  if (ndim != 2)
  {
    printf("couldn't read the mesh because it is of dimension %d, and the function only reads 2D meshes. \n", ndim);
    exit(-1);
  }

  //**************************************
  //****** Reading Information that will speed up the creation of the mesh
  //****** Extra information: number of internal nodes, number of internal faces with internal nodes and number of internal elements
  if(writtenMeshGen)
  {
    n=0;
    while ((strcmp(line, "\"Information:\")") != 0) && (n<100000))    {fscanf(fp, "%s\n", line); n++;}
    fscanf(fp, "%s %s %x %x %x %s \n", char1, char2, &nno_i, &nfa_nno_i, &nel_i, char3);
    wittenGambit=false;
  }

  //**************************************
  //****** Reading the nodes
  n=0;
  while ((strcmp(line, "2))") != 0)  && (n<100000))    {fscanf(fp, "%s\n", line); n++;}

  fscanf(fp, "%s %s %x %x %i %s \n", char1, char2, &dummy1, &nvert, &dummy3, char3);
  nodes.resize(nvert);
  // Read nodal coordinates
  for (int i = 0; i < nvert; i++)
  {
    fscanf(fp, "%lf%lf\n", &nodes[i].pt.x, &nodes[i].pt.y); nodes[i].pt.z = 0.0;
    if (i%1000 == 0) printf("%10d\t%lf\t%lf\t%lf\n", i, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  }

  //**************************************
  //****** Reading the Faces
  int fir_face;
  n=0;
  while ((strcmp(line, "\"Faces:\")") != 0) && (n<100000))    {fscanf(fp, "%s\n", line) ; n++;}


  fscanf(fp, "%s %x %x %s \n", char1, &fir_face, &nfaces , char3);


  //Internal array for reading the faces
  int tmps[4];
  deque<FACE> faces_bc;
  deque<FACE> faces_fluid;

  int nface_total=0;
  int n1;
  int n2;
  do
  {
    fscanf(fp, "%s %x %x %x %s \n", char1, &n1, &n2, &dummy2, char3);

    for (int i = n1; i <= n2; i++)
    {
      fscanf(fp, "%10d%x%x%x%x\n", &dummy1, &tmps[0], &tmps[1], &tmps[2], &tmps[3]);

      FACE tmp_face;

      if(tmps[2]!=0)
      {
        tmp_face.node.push_back(tmps[1]-1);
        tmp_face.node.push_back(tmps[0]-1);
        tmp_face.elem_1 = tmps[2];
        tmp_face.elem_2 = tmps[3];
        if(( (dummy2 != 2) && (tmps[3]!=0) ) || ( (dummy2 == 2) && (tmps[3]==0) ))
        {printf("Problem when reading faces, elem_2 doesn't match! \n");  exit(-1);}

        if(tmps[3]==0)		{sprintf(tmp_face.name,"noname");		faces_bc.push_back(tmp_face);}
        else				{sprintf(tmp_face.name,"fluid"); 		faces_fluid.push_back(tmp_face);}
      }
      else
      {
        tmp_face.node.push_back(tmps[0]-1);
        tmp_face.node.push_back(tmps[1]-1);
        tmp_face.elem_1 = tmps[3];
        tmp_face.elem_2 = tmps[2];
        if(( (dummy2 != 2) && (tmps[2]!=0) ) || ( (dummy2 == 2) && (tmps[2]==0) ))
        {printf("Problem when reading faces, elem_2 doesn't match! \n");  exit(-1);}

        if(tmps[2]==0)		{sprintf(tmp_face.name,"noname");		faces_bc.push_back(tmp_face);}
        else				{sprintf(tmp_face.name,"fluid"); 		faces_fluid.push_back(tmp_face);}
      }
      nface_total++;
    }
    fscanf(fp, "%s\n", line);
  }while( (dummy2 != 2) && (nface_total<nfaces));

  //**************************************
  //****** Reading the number of ELEMENTS
  while (strcmp(line, "\"Cells:\")") != 0)    fscanf(fp, "%s\n", line);
  fscanf(fp, "%s %s %x %x %i %s \n", char1, char2, &dummy1, &nelem, &dummy2, char3);

  printf("\n Fluent input file (%s format)! \n", name);




  //**************************************
  //**************************************
  //****** Generating the actual mesh

  printf("\n Starting to sort the nodes and generate the elements in %s \n", name);

  //------ Faces: internal and boundary
  nfa_i         =faces_fluid.size();
  int nfa_bc    =faces_bc.size();
  int nface     = faces.size();
  for(int i=0; i<nfa_i; i++)        faces.push_back(faces_fluid[i]);
  for(int i=0; i<nfa_bc; i++)       faces.push_back(faces_bc[i]);

  int count=0;

  if(wittenGambit)
  {
    printf(" Inside when it is a Gambit written file@@@ %s \n", name);
    //**************************************
    //******* Node: sorting nodes
    deque<int>	move_node;			// This deque will store the nodes that were already moved, so the same node is not moved twice
    move_node.push_back(-1);


    for (int i=nfa_i; i<nface; i++)
      if(faces[i].elem_2==0)
      {
        int ind;
        //...... Node 1 of the face
        ind=faces[i].node[0];
        count=0;
        for(int j=0; j<move_node.size(); j++)
          if (move_node[j]==ind)
            count++;
        if(count==0)
        {
          // move node to the end of the deque
          NODE tmpNode = nodes[ind];
          nodes.erase(nodes.begin()+ind);
          nodes.push_back(tmpNode);

          // update nodes in Faces;
          for (int fa=0; fa<nface; fa++)
          {
            if (faces[fa].node[0]==ind)      faces[fa].node[0] = nodes.size()-1;
            else if (faces[fa].node[0]>ind)  faces[fa].node[0]--;

            if (faces[fa].node[1]==ind)      faces[fa].node[1] = nodes.size()-1;
            else if (faces[fa].node[1]>ind)  faces[fa].node[1]--;

          }
          // update nodes in move_node deque
          for(int k=0; k<move_node.size(); k++)
            if (move_node[k]>ind)  move_node[k]--;

          move_node.push_back(nodes.size()-1);
        }

        //...... Node 2 of the face
        ind=faces[i].node[1];
        count=0;
        for(int j=0; j<move_node.size(); j++)
          if (move_node[j]==ind)
            count++;
        if(count==0)
        {
          // move node to the end of the deque
          NODE tmpNode = nodes[ind];
          nodes.erase(nodes.begin()+ind);
          nodes.push_back(tmpNode);

          // update nodes in Faces;
          for (int fa=0; fa<nface; fa++)
          {
            if (faces[fa].node[0]==ind)      faces[fa].node[0] = nodes.size()-1;
            else if (faces[fa].node[0]>ind)  faces[fa].node[0]--;

            if (faces[fa].node[1]==ind)      faces[fa].node[1] = nodes.size()-1;
            else if (faces[fa].node[1]>ind)  faces[fa].node[1]--;

          }
          // update nodes in move_node deque
          for(int k=0; k<move_node.size(); k++)
            if (move_node[k]>ind)  move_node[k]--;

          move_node.push_back(nodes.size()-1);
        }
      }

    nno_i=move_node[1];



    //**************************************
    //****** Faces: Moving the internal faces that have only internal nodes to the front of the array faces[i]
    nfa_nno_i=0;
    printf(" Moving the internal faces that have only internal nodes in %s \n", name);
    for (int i=0; i<nfa_i; i++)
      if((faces[i].node[0]<nno_i) && (faces[i].node[1]<nno_i))
      {
        FACE tmp_face;
        tmp_face=faces[i];
        faces.erase(faces.begin()+i);
        faces.push_front(tmp_face);
        nfa_nno_i++;

      }
  }


  //**************************************
  //****** Elements: generating the elements
  deque< deque<int> > ind_elem_int; //internal elements
  deque< deque<int> > ind_elem_ext; //external elemnts

  // Initializing the element deques
  for(int e=1; e<=nelem; e++)
  {
    deque <int> tmp;
    tmp.push_back(e);
    ind_elem_ext.push_back(tmp);
    ind_elem_int.push_back(tmp);
  }
//  int nelem_i;
  printf(" Finding the faces of each element in %s \n", name);
  for (int i=0; i<faces.size(); i++)
  {
    bool bc_elem=false;
    int e1 = faces[i].elem_1;
    int e2 = faces[i].elem_2;

    if(e2==0)  bc_elem=true;

    e1--;
    e2--;

    if  (bc_elem)
    {
      ind_elem_ext[e1].push_back(i);
      ind_elem_int[e1].push_front(-1);  //Flag to know that element e1 is in the boundary!
    }
    else
    {
      ind_elem_int[e1].push_back(i);
      ind_elem_int[e2].push_back(i);
      ind_elem_ext[e1].push_back(i);
      ind_elem_ext[e2].push_back(i);
    }

  }

  // All the elements that are not Boundary are therefore internal!
  for(int e=nelem-1; e>=0; e--)
    if(ind_elem_int[e][0]>0)
      ind_elem_ext[e].push_front(-1);  //Flag to know that element e is internal!

  int newElem [nelem] ;		//This array will be used later to update the elements in the face object!!
  count=0;
  elem_i.push_back(0);


  //------------------
  //	Generating the internal elements
  printf(" Generating the internal elements in %s \n", name);
  for(int i=0; i<ind_elem_int.size(); i++)
    if(ind_elem_int[i][0]>0)
    {
      //Array where the nodes of the element will be stored
      deque <int> nodes_elem;
      nodes_elem.push_back(-1);     //initializing the array

      int tmp_elem= ind_elem_int[i][0];   //the first element of the array is the element index
      newElem[tmp_elem-1]=elem_i.size();  //array used later to update the element in the faces

      // Element can be generated with faces-1
      int ind_elem_int_size=ind_elem_int[i].size()-2;
      for(int j=0; j<ind_elem_int_size; j++)
      {
        FACE tmp_face=faces[ind_elem_int[i][j+1]];
        deque <int> tmp_nodes_elem;

        //The face are generated in a certain order, right-hand rule!
        if(tmp_face.elem_1==tmp_elem)	{tmp_nodes_elem.push_back(tmp_face.node[1]);	tmp_nodes_elem.push_back(tmp_face.node[0])	;}
        else							            {tmp_nodes_elem.push_back(tmp_face.node[0]);	tmp_nodes_elem.push_back(tmp_face.node[1])	;}

        // The nodes of the elements most be sorted in a certain way: counter-clockwise, the following if conditions make sure of this!
        if (j==0)
        {
          nodes_elem.clear();
          nodes_elem.push_back(tmp_nodes_elem[0]);
          nodes_elem.push_back(tmp_nodes_elem[1]);
        }
        else if (nodes_elem[nodes_elem.size()-1]==tmp_nodes_elem[0])				nodes_elem.push_back(tmp_nodes_elem[1]);
        else if (nodes_elem[0]==tmp_nodes_elem[1])									        nodes_elem.push_front(tmp_nodes_elem[0]);
        else																		                            ind_elem_int_size=ind_elem_int[i].size()-1;  //In case all faces are needed to generate the element

      }
      for(int j=0; j<nodes_elem.size(); j++)	{elem_v.push_back(nodes_elem[j]);	count++;}
      elem_i.push_back(count);
      nodes_elem.clear();
    }
  nel_i=elem_i.size()-1;

  //------------------
  //	Generating the boundary elements
  printf(" Generating the boundary elements in %s \n", name);
  for(int i=0; i<ind_elem_ext.size(); i++)
    if(ind_elem_ext[i][0]>0)
    {
      deque <int> nodes_elem;
      nodes_elem.push_back(-1);

      int tmp_elem= ind_elem_ext[i][0];
      newElem[tmp_elem-1]=elem_i.size();
      int ind_elem_ext_size=ind_elem_ext[i].size()-2;
      for(int j=0; j<ind_elem_ext_size; j++)
      {
        FACE tmp_face=faces[ind_elem_ext[i][j+1]];
        deque <int> tmp_nodes_elem;

        if(tmp_face.elem_1==tmp_elem)	{tmp_nodes_elem.push_back(tmp_face.node[1]);	tmp_nodes_elem.push_back(tmp_face.node[0])	;}
        else							            {tmp_nodes_elem.push_back(tmp_face.node[0]);	tmp_nodes_elem.push_back(tmp_face.node[1])	;}

        // The nodes of the elements most be sorted in a certain way: counter-clockwise
        if (j==0)
        {
          nodes_elem.clear();
          nodes_elem.push_back(tmp_nodes_elem[0]);
          nodes_elem.push_back(tmp_nodes_elem[1]);
        }
        else if (nodes_elem[nodes_elem.size()-1]==tmp_nodes_elem[0])					nodes_elem.push_back(tmp_nodes_elem[1]);
        else if(nodes_elem[0]==tmp_nodes_elem[1])										          nodes_elem.push_front(tmp_nodes_elem[0]);
        else																			                            ind_elem_ext_size=ind_elem_ext[i].size()-1;

      }
      for(int j=0; j<nodes_elem.size(); j++)	{elem_v.push_back(nodes_elem[j]);	count++;}
      elem_i.push_back(count);
      nodes_elem.clear();
    }

  //Updating the elements in the face!
  for (int i=0; i<nface; i++)
  {
    faces[i].elem_1=newElem[faces[i].elem_1-1];
    if(faces[i].elem_2!=0)	faces[i].elem_2=newElem[faces[i].elem_2-1];
  }

  ind_elem_int.clear();
  ind_elem_ext.clear();

  //******************************************************************************************************************



  printf(" Finding the Nodes Neighbors \n");
  //Finding the Neighbors of each node
  findNodeNeighbors();
  printf("**** Mesh %s has been generated!\n", name);
}


void UNSTRUCTMESH::readFluentMsh(const char *name, int nvert_gambit)
{
  //THis function needs to be combined with readGambitNeu2D, because this function only reads the faces: its connectivity!!
  //The elements needs to be know before handed (related to the nodes).

  printf("\n Starting to read the Fluent input file (*.msh format)! \n");

  FILE *fp;
  if ((fp = fopen(name, "rt")) == NULL)
  {
    printf("couldn't open %s for reading\n", name);
    exit(-1);
  }

  int ndim;
  int dummy1, dummy2, dummy3;
  int nvert, nfaces_bc, nfaces_fluid, nfaces, nnodes;
  string str1, str2, str3;
  char line[1000];
  fscanf(fp, "%s\n", line);

  // Reading the dimensiones of the mesh!! header
  while (strcmp(line, "\"Dimension:\")") != 0)    fscanf(fp, "%s\n", line);
  fscanf(fp, "%s %d %s \n", str1.c_str(), &ndim, str3.c_str());



  while (strcmp(line, "\"Nodes:\")") != 0)    fscanf(fp, "%s\n", line);


  //  cout<<"Before reading the vertices"<<endl;

  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nvert, &dummy2, str3.c_str());

  printf(" %dD Grid --> Vertices: %i      \n", ndim, nvert);

  if(nvert-nvert_gambit == 0)       cout<<"Number of vertices of mesh read in gambit and fluend format the same: "<< nvert<< endl;
  else {
    cout<<"Check that the mesh being read from Gambit format and Fluent format are the same"<< nvert<< endl;
  }

  while (strcmp(line, "\"Faces:\")") != 0)    fscanf(fp, "%s\n", line);
  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nfaces, &dummy2, str3.c_str());
  cout<<"Starting to read bc faces: Number of total faces: "<< nfaces<< endl;
  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nfaces_bc, &dummy2, str3.c_str());
  //  fscanf(fp, "%s\n", line);

  int nodes[6];
  int count=0;

  //  cout<<"Starting to read bc faces: Number of boundary faces: "<< nfaces_bc<< endl;

  //  cout<<"Number of boundary faces: "<< nfaces_bc<< endl;

  //======================================================
  //Boundary faces!!!
  for (int i = 0; i < nfaces_bc; i++)
  {
    if      (ndim==2)    fscanf(fp, "%10d%x%x%x%x\n", &nnodes, &nodes[0], &nodes[1], &nodes[2], &nodes[3]);
    else if (ndim==3)    fscanf(fp, "%10d%x%x%x%x%x%x\n", &nnodes, &nodes[0], &nodes[1], &nodes[2], &nodes[3], &nodes[4], &nodes[5]);
    else
    {
      printf("\n Dimension of the mesh is wrong!! \n");
      throw(-1);
    }

    FACE tmp_face;
    for(int nf=0; nf<nnodes; nf++)    tmp_face.node.push_back(nodes[nf]-1);

    if(nodes[nnodes]==0)
      tmp_face.elem_1 = nodes[nnodes+1]; //face_elem.push_back(nodes[nnodes]);
    else
      tmp_face.elem_1 = nodes[nnodes]; //face_elem.push_back(nodes[nnodes]);

    tmp_face.elem_2 = 0; //face_elem.push_back(0);
    strcpy (tmp_face.name,"noname");
    faces.push_back(tmp_face);

//    for(int nf=0; nf<nnodes; nf++)    face_v.push_back(nodes[nf]-1);
//    face_elem.push_back(nodes[nnodes]);
//    face_elem.push_back(0);
//    face_i.push_back(count*nnodes);
//    face_name.push_back("noname");



    count++;

  }
  fscanf(fp, "%s\n", line);

  //======================================================
  //Fluid faces!!!
  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nfaces_fluid, &dummy2, str3.c_str());
  //  fscanf(fp, "%s\n", line);
  nfaces_fluid=nfaces-nfaces_bc;

  for (int i = 0; i < nfaces_fluid; i++)
  {
    if      (ndim==2)    fscanf(fp, "%10d%x%x%x%x\n",       &nnodes, &nodes[0], &nodes[1], &nodes[2], &nodes[3]);
    else if (ndim==3)    fscanf(fp, "%10d%x%x%x%x%x%x\n",   &nnodes, &nodes[0], &nodes[1], &nodes[2], &nodes[3], &nodes[4], &nodes[5]);
    else
    {
      printf("\n Dimension of the mesh is wrong!! \n");
      throw(-1);
    }


    FACE tmp_face;
    for(int nf=0; nf<nnodes; nf++)    tmp_face.node.push_back(nodes[nf]-1);
    tmp_face.elem_1 = nodes[nnodes]; //face_elem.push_back(nodes[nnodes]);
    tmp_face.elem_2 = nodes[nnodes+1]; //face_elem.push_back(0);
    strcpy (tmp_face.name,"fluid");
    faces.push_back(tmp_face);

//    for(int nf=0; nf<nnodes; nf++)    face_v.push_back(nodes[nf]-1);
//    face_elem.push_back(nodes[nnodes]);
//    face_elem.push_back(nodes[nnodes+1]);
//    face_i.push_back(count*nnodes);
//
//    face_name.push_back("fluid");
    count++;

  }
  printf("\n Fluent input file (*.msh format) was read! \n");
  //  cout<<"Number of vertices: "<< nvert<<endl;

}

//***********************************

void UNSTRUCTMESH::readGambitNeu2D(const char *name)
{
  FILE *fp;
  if ((fp = fopen(name, "rt")) == NULL)
  {
    printf("couldn't open %s for reading\n", name);
    exit(-1);
  }

  // Read header
  char line[1000];
  fscanf(fp, "%s\n", line);
  while (strcmp(line, "NDFVL") != 0)    fscanf(fp, "%s\n", line);

  int ngrps, nbsets, ndfcd, ndfvl, dummy;
  int nvertices, nelements;
  fscanf(fp, "%i %i %i %i %i %i \n", &nvertices, &nelements, &ngrps, &nbsets, &ndfcd, &ndfvl);

  printf(" %dD Grid --> Vertices: %i        Elements: %i \n", ndfcd, nvertices, nelements);

  // Allocate memory to save vertices and elements
  nodes.resize(nvertices);


  while (strcmp(line, GAMBIT_VER) != 0)    fscanf(fp, "%s\n", line);

  // Read nodal coordinates
  for (int i = 0; i < nvertices; i++)
  {
    fscanf(fp, "%10d%lf%lf\n", &dummy, &nodes[i].pt.x, &nodes[i].pt.y); nodes[i].pt.z = 0.0;
    if (i%100 == 0) printf("%10d\t%lf\t%lf\t%lf\n", dummy, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  }

  fscanf(fp, "%s\n", line);
  while (strcmp(line, GAMBIT_VER) != 0)    fscanf(fp, "%s\n", line);

  // Read elements
  elem_i.push_back(0);

  for (int c = 0; c < nelements; c++)
  {
    int nvert, type;
    fscanf(fp, "%8i%2i%2i\n", &dummy, &type, &nvert);


    int nodes[10];

    if (type == NTYPE_TRI)
    {
      fscanf(fp, "%8i%8i%8i\n", &nodes[0], &nodes[1], &nodes[2]);
    }
    else if (type == NTYPE_QUAD)
    {
      fscanf(fp, "%8i%8i%8i%8i\n", &nodes[0], &nodes[1],
          &nodes[2], &nodes[3]);
    }
    else if (type == NTYPE_BRICK)
    {
      fscanf(fp, "%8i%8i%8i%8i%8i%8i%8i\n%10d\n",
          &nodes[0], &nodes[1], &nodes[2], &nodes[3],
          &nodes[4], &nodes[5], &nodes[6], &nodes[7]);
    }


    for (int ne=0; ne<nvert; ne++) elem_v.push_back(nodes[ne]-1);
    elem_i.push_back(elem_i[elem_i.size()-1]+nvert);

  }
  fclose(fp);

  nel_i = 0;
  nno_i = 0;

  findFaces2D();
  findNodeNeighbors();

  printf("\n input file in Gambit Neutral format read! \n");
}

void UNSTRUCTMESH::readGambitNeu(const char *name)
{
  FILE *fp;
  if ((fp = fopen(name, "rt")) == NULL)
  {
    printf("couldn't open %s for reading\n", name);
    exit(-1);
  }

  // Read header
  char line[1000];
  fscanf(fp, "%s\n", line);
  while (strcmp(line, "NDFVL") != 0)    fscanf(fp, "%s\n", line);

  int ngrps, nbsets, ndfcd, ndfvl, dummy;
  int nvertices, nelements;
  fscanf(fp, "%i %i %i %i %i %i \n", &nvertices, &nelements, &ngrps, &nbsets, &ndfcd, &ndfvl);
  printf(" %dD Grid --> Vertices: %i        Elements: %i \n", ndfcd, nvertices, nelements);

  // Allocate memory to save vertices and elements
  nodes.resize(nvertices);

  while (strcmp(line, GAMBIT_VER) != 0)    fscanf(fp, "%s\n", line);

  // Read nodal coordinates
  for (int i = 0; i < nvertices; i++)
  {
    if(ndfcd==2)    fscanf(fp, "%10d%lf%lf\n",      &dummy, &nodes[i].pt.x, &nodes[i].pt.y); nodes[i].pt.z = 0.0;

    if(ndfcd==3)    fscanf(fp, "%10d%lf%lf%lf\n",   &dummy, &nodes[i].pt.x, &nodes[i].pt.y, &nodes[i].pt.z);

    if (i%100 == 0) printf("%10d\t%lf\t%lf\t%lf\n", dummy, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  }

  fscanf(fp, "%s\n", line);
  while (strcmp(line, GAMBIT_VER) != 0)    fscanf(fp, "%s\n", line);

  // Read elements
  elem_i.push_back(0);

  for (int c = 0; c < nelements; c++)
  {
    int nvert, type;
    fscanf(fp, "%8i%2i%2i\n", &dummy, &type, &nvert);

    int nodes[10];

    if (type == NTYPE_TRI)
    {
      fscanf(fp, "%8i%8i%8i\n", &nodes[0], &nodes[1], &nodes[2]);
    }
    else if (type == NTYPE_QUAD)
    {
      fscanf(fp, "%8i%8i%8i%8i\n", &nodes[0], &nodes[1],
          &nodes[2], &nodes[3]);
    }
    else if (type == NTYPE_BRICK)
    {
      fscanf(fp, "%8i%8i%8i%8i%8i%8i%8i\n%10d\n",
          &nodes[0], &nodes[1], &nodes[3], &nodes[2],
          &nodes[4], &nodes[5], &nodes[7], &nodes[6]);
    }


    for (int ne=0; ne<nvert; ne++)
    {
      elem_v.push_back(nodes[ne]-1);
    }
    elem_i.push_back(elem_i[elem_i.size()-1]+nvert);

  }
  fclose(fp);

  nel_i = 0;
  nno_i = 0;

//  findFaces2D();
//  findNodeNeighbors();

  printf("\n input file in Gambit Neutral format read! \n");
}


void UNSTRUCTMESH::readNastranMod(const char *name)
{
  FILE *fp;
  if ((fp = fopen(name, "rt")) == NULL)
  {
    printf("couldn't open %s for reading\n", name);
    exit(-1);
  }

  // Read header
  //  char line[1000];
  //  fscanf(fp, "%s\n", line);
  //  while (strcmp(line, "NDFVL") != 0)    fscanf(fp, "%s\n", line);


  int nvertices, nelements;
  fscanf(fp, "%d\n%d\n", &nvertices, &nelements);

  printf(" NASTRAN Grid --> Vertices: %i        Elements: %i \n", nvertices, nelements);

  // Allocate memory to save vertices and elements
  nodes.resize(nvertices);

  //  while (strcmp(line, GAMBIT_VER) != 0)    fscanf(fp, "%s\n", line);

  int dummy;
  char line[1000];

  // Read nodal coordinates
  for (int i = 0; i < nvertices; i++)
  {
    fscanf(fp, "%d%d%lf%lf\n", &dummy, &dummy, &nodes[i].pt.x, &nodes[i].pt.y);
    fscanf(fp, "%s\n", line);
    fscanf(fp, "%7s%lf\n", line, &nodes[i].pt.z);
  }

  // Read elements
  elem_i.push_back(0);

  for (int c = 0; c < nelements; c++)
  {
    int nodes[10];

    fscanf(fp, "%18s%d%d%d%d%d%d%d%d", line, &dummy, &dummy, &nodes[0], &nodes[1], &nodes[2], &nodes[3],&nodes[4], &nodes[5]);
    fscanf(fp, "%s\n", line);
    fscanf(fp, "%7s%d%d\n", line, &nodes[6], &nodes[7]);

    for (int ne=0; ne<8; ne++) elem_v.push_back(nodes[ne]-1);
    elem_i.push_back(elem_i[elem_i.size()-1]+8);
  }
  fclose(fp);

  findFaces2D();
  findNodeNeighbors();

  printf("\n input file in NASTRAN format read! \n");
}

void UNSTRUCTMESH::drawMesh()
{
  //***************** with Elements *************************
//  for (int e=0; e<elem_i.size()-1; e++)
//  {
//    //        glLineWidth(0.2);
//    glBegin(GL_LINE_LOOP);
//    {
//      double grey = 0.5;
//      glColor3d(0.0, 0.0, 0.0);
//
//
//      for (int ne=elem_i[e]; ne<elem_i[e+1]-1; ne++)
//      {
//        glVertex3d(nodes[elem_v[ne]].pt.x,   nodes[elem_v[ne]].pt.y,   nodes[elem_v[ne]].pt.z);
//        glVertex3d(nodes[elem_v[ne+1]].pt.x, nodes[elem_v[ne+1]].pt.y, nodes[elem_v[ne+1]].pt.z);
//      }
//    }
//    glEnd();
//    //      glLineWidth(1.0);
//  }
  //***************** with Faces ****************************
  for (int f=0; f<faces.size(); f++)
  {
    //        glLineWidth(0.2);
    glBegin(GL_LINE_LOOP);
    {
      double grey = 0.5;
      glColor3d(0.10, 0.0, 0.0);

      int nnodes = faces[f].node.size();
      for (int fe=0; fe<nnodes-1; fe++)
      {
        POINT p1=nodes[faces[f].node[fe]].pt;
        POINT p2=nodes[faces[f].node[fe+1]].pt;
        glVertex3d(p1.x,   p1.y,   p1.z);
        glVertex3d(p2.x,   p2.y,   p2.z);

      }
//      POINT p1=nodes[faces[f].node[nnodes-1]].pt;
//      POINT p2=nodes[faces[f].node[0]].pt;
//      glVertex3d(p1.x,   p1.y,   p1.z);
//      glVertex3d(p2.x,   p2.y,   p2.z);


    }
    glEnd();
    //      glLineWidth(1.0);
  }



}

void UNSTRUCTMESH::drawSolid()
{
  //***************** with Elements *************************
//  for (int e=0; e<elem_i.size()-1; e++)
//  {
//    glBegin(GL_TRIANGLES);
//    glColor3d(0.85, 0.85, 0.85);
//
//    POINT t1 = nodes[elem_v[elem_i[e]]].pt;
//    for (int ne=elem_i[e]; ne<elem_i[e+1]-1; ne++)
//    {
//      POINT t2 = nodes[elem_v[ne]].pt;
//      POINT t3 = nodes[elem_v[ne+1]].pt;
//
//      glVertex3d(t1.x, t1.y, t1.z);
//      glVertex3d(t2.x, t2.y, t2.z);
//      glVertex3d(t3.x, t3.y, t3.z);
//    }
//    glEnd();
//  }
//  //***************** with Faces ****************************
  for (int f=0; f<faces.size(); f++)
  {
    glBegin(GL_TRIANGLES);
    glColor3d(0.85, 0.85, 0.85);

    POINT t1 = nodes[faces[f].node[0]].pt;
    for (int fe=0; fe<faces[f].node.size()-1; fe++)
    {
      POINT t2 = nodes[faces[f].node[fe]].pt;
      POINT t3 = nodes[faces[f].node[fe+1]].pt;

      glVertex3d(t1.x, t1.y, t1.z);
      glVertex3d(t2.x, t2.y, t2.z);
      glVertex3d(t3.x, t3.y, t3.z);
    }
    glEnd();
  }

}




//=============================================
//        UNSTRUCTMESH CLASS FUNCTIONS
//=============================================

void TRIANGULATE::unstructuredMesh2D()
{
  struct triangulateio in, out, vorout;

  int Next = extBoundary.size();
  // ***********************
  // Define number of points
  in.numberofpoints = Next;

  if (holes.size()>0)
  {
    for (int i=0; i<holes.size(); i++)
      in.numberofpoints += holes[i].holesPoints.size();
  }
  //        printf("number of points: %d\n",in.numberofpoints);


  // ****************************
  // Assign points [x y x y ....]
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

  for (int i=0; i<2*extBoundary.size(); i=i+2)
  {
    // assign all the points
    in.pointlist[i]   = extBoundary[i/2].x;
    in.pointlist[i+1] = extBoundary[i/2].y;
  }
  //        printf("External boundary points assigned\n");

  // HOLES
  if (holes.size()>0)
  {
    int N=0;
    for (int j=0; j<holes.size(); j++)
    {
      for (int i=0; i<holes[j].holesPoints.size()*2; i=i+2)
      {
        in.pointlist[Next*2+N+i]   = holes[j].holesPoints[i/2].x;
        in.pointlist[Next*2+N+i+1] = holes[j].holesPoints[i/2].y;
      }
      N += holes[j].holesPoints.size()*2;
    }
  }
  //        printf("Hole assigned\n");


  // *****************************************************
  // Assign markers to points (1 means they are a boundary
  in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
  for (int i=0; i<in.numberofpoints; i++)
    in.pointmarkerlist[i] = 1;
  //        printf("point markers assigned\n");



  // ******************************************************************
  // Assign points attributes (can specify area for example). note used
  // point attributes
  in.numberofpointattributes = 1;
  in.pointattributelist = (REAL *) malloc(in.numberofpoints *
      in.numberofpointattributes *
      sizeof(REAL));
  for (int i=0; i<in.numberofpoints*in.numberofpointattributes; i++)
    in.pointattributelist[i] = 0.0;
  //        printf("point attributes assigned\n");



  // *****************************************
  // Define segments [start end start end ...]
  in.numberofsegments = Next;
  if (holes.size()>0)
    for (int i=0; i<holes.size(); i++)
      in.numberofsegments += holes[i].holesPoints.size();

  printf("number of segments %d\n", in.numberofsegments);
  in.segmentlist = (int *) malloc(in.numberofsegments*2 * sizeof(int));

  int vert = 1;
  for (int i=0; i<Next*2; i=i+2)
  {
    in.segmentlist[i] = vert;
    in.segmentlist[i+1] = vert+1;
    vert++;
  }
  in.segmentlist[2*Next-1] = 1;

  in.numberofholes = holes.size();

  // HOLES
  if (holes.size()>0)
  {
    int N=0;
    int sum=0;
    for (int j=0; j<holes.size(); j++)
    {
      sum += holes[j].holesPoints.size();
      vert = Next+N+1;
      for (int i=0; i<holes[j].holesPoints.size()*2; i=i+2)
      {
        in.segmentlist[2*Next+N*2+i] = vert;
        in.segmentlist[2*Next+N*2+i+1] = vert+1;
        vert++;
      }
      in.segmentlist[2*Next+sum*2-1] = Next+N+1;
      N = holes[j].holesPoints.size();
    }
  }
  //    printf("Holes segments assigned\n");


  // *********************************************
  // Define 1 point inside each hole [x y x y ...]
  in.holelist = (REAL *) malloc(in.numberofholes*2 * sizeof(REAL));
  for (int i=0; i<2*holes.size(); i=i+2)
  {
    in.holelist[i]   = holes[i/2].insidePoints.x;
    in.holelist[i+1] = holes[i/2].insidePoints.y;
  }


  // ******************************************
  // Segment markers (1 means it is a boundary)
  in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
  for (int i=0; i<in.numberofsegments; i++)
    in.segmentmarkerlist[i] = 1;
  printf("Segments markers assigned\n");

  in.numberofregions = 1;
  in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));



  // **********************************
  // Initialize out and vorout as NULL
  out.pointlist                 = (double *) NULL;
  out.pointattributelist        = (double *) NULL;
  out.pointmarkerlist           = (int *) NULL;
  out.trianglelist              = (int *) NULL;
  out.triangleattributelist     = (double *) NULL;
  out.neighborlist              = (int *) NULL;
  out.segmentlist               = (int *) NULL;
  out.segmentmarkerlist         = (int *) NULL;
  out.edgelist                  = (int *) NULL;
  out.edgemarkerlist            = (int *) NULL;

  vorout.pointlist              = (double *) NULL;
  vorout.pointattributelist     = (double *) NULL;
  vorout.edgelist               = (int *) NULL;
  vorout.normlist               = (double *) NULL;


  // *******************************************************
  // Call to external function to generate unstructured mesh
  triangulate(triangParameters,&in, &out, &vorout);

  // *****************************************
  // Save data in POINT and TRIANGLE structure

  // assign nodes
  for (int i=0; i<2*out.numberofpoints; i=i+2)
    push_back_node(POINT(out.pointlist[i], out.pointlist[i+1], extBoundary[0].z));

  for (int i=0; i<3*out.numberoftriangles; i++)
    elem_v.push_back(out.trianglelist[i]-1);

  for (int i=0; i<=out.numberoftriangles; i++)
    elem_i.push_back(3*i);

  // points are listed such that boundary points (marker = 1) are first in the array
  for (int i=0; i<out.numberofpoints; i++)
    nodes[i].marker = out.pointmarkerlist[i];


  if (holes.size()<=0)
  {
    moveBoundaryNodesEnd(extBoundary);
    moveBoundaryElementsEnd(extBoundary);
  }
  else
  {
    deque<POINT> completeBoundary;
    completeBoundary = extBoundary;
    for (int h=0; h<holes.size(); h++)
    {
      completeBoundary.insert(completeBoundary.end(), holes[h].holesPoints.begin(), holes[h].holesPoints.end());
      completeBoundary.push_back(holes[h].holesPoints[0]);
    }

    moveBoundaryNodesEnd(completeBoundary);
    moveBoundaryElementsEnd(completeBoundary);
  }

//  findNodeNeighbors();
  findFaces2D();
  findNodeNeighbors();
}

void TRIANGULATE::moveBoundaryElementsEnd(deque<POINT> boundPts)
{
  int elemMoved = 0;
  for (int i=1; i<boundPts.size()-1; i++)
  {
    int ind=0;

    while (ind<elem_i.size()-1-elemMoved)
    {
      int count = 0;
      int nNodes = elem_i[ind+1]-elem_i[ind];

      if(nNodes == 3)
      {
        if ((nodes[elem_v[elem_i[ind]]].pt==boundPts[i-1])      && (nodes[elem_v[elem_i[ind]+1]].pt==boundPts[i]))          count++;
        if ((nodes[elem_v[elem_i[ind]]].pt==boundPts[i])        && (nodes[elem_v[elem_i[ind]+1]].pt==boundPts[i+1]))        count++;
        if ((nodes[elem_v[elem_i[ind]]].pt==boundPts[i+1])      && (nodes[elem_v[elem_i[ind]+1]].pt==boundPts[i-1]))        count++;

        if ((nodes[elem_v[elem_i[ind]+1]].pt==boundPts[i-1])    && (nodes[elem_v[elem_i[ind]+2]].pt==boundPts[i]))          count++;
        if ((nodes[elem_v[elem_i[ind]+1]].pt==boundPts[i])      && (nodes[elem_v[elem_i[ind]+2]].pt==boundPts[i+1]))        count++;
        if ((nodes[elem_v[elem_i[ind]+1]].pt==boundPts[i+1])    && (nodes[elem_v[elem_i[ind]+2]].pt==boundPts[i-1]))        count++;

        if ((nodes[elem_v[elem_i[ind]+2]].pt==boundPts[i-1])    && (nodes[elem_v[elem_i[ind]]].pt==boundPts[i]))            count++;
        if ((nodes[elem_v[elem_i[ind]+2]].pt==boundPts[i])      && (nodes[elem_v[elem_i[ind]]].pt==boundPts[i+1]))          count++;
        if ((nodes[elem_v[elem_i[ind]+2]].pt==boundPts[i+1])    && (nodes[elem_v[elem_i[ind]]].pt==boundPts[i-1]))          count++;

      }



      if (count>0)
      {
        moveElementEnd(ind);
        elemMoved++;
        cout<<"element moved!! "<<elemMoved<<endl;
      }
      else ind++;
    }
  }
  nel_i = elem_i.size()-1-elemMoved; // -1 because the last i is not an element
}

UNSTRUCTMESH operator*(const UNSTRUCTMESH &mesh, const double val)
{
  UNSTRUCTMESH res = mesh;
  return (res*=val);
}

UNSTRUCTMESH operator*=(UNSTRUCTMESH &mesh, const double val)
    {
  for (int i=0; i<mesh.nodes.size(); i++)
    mesh.nodes[i].pt = mesh.nodes[i].pt*val;
  return (mesh);
    }
