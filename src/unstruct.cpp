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
    // TO DO:
    // CAN BE RE-WRITTEN USING ALREADY nodes, elem_i AND elem_v (NO NEED FOR *_new)

    // NODES
    deque<NODE> nodes_new;
    for (int i=0; i<nno_i; i++)                         nodes_new.push_back(this->nodes[i]);
    for (int i=0; i<mesh.nno_i; i++)                    nodes_new.push_back(mesh.nodes[i]);
    for (int i=nno_i; i<nodes.size(); i++)              nodes_new.push_back(this->nodes[i]);
    for (int i=mesh.nno_i; i<mesh.nodes.size(); i++)    nodes_new.push_back(mesh.nodes[i]);


//    nodes.insert(nodes.begin()+nno_i, mesh.nodes.begin(),mesh.nodes.begin()+mesh.nno_i);
//    nodes.insert(nodes.end(), mesh.nodes.begin()+mesh.nno_i,mesh.nodes.end());

    // ELEM_I
    deque<int> elem_i_new;
    for (int i=0; i<nel_i; i++)                             elem_i_new.push_back(elem_i[i]);
    for (int i=0; i<mesh.nel_i; i++)                        elem_i_new.push_back(mesh.elem_i[i]+elem_i[nel_i]);
    for (int i=nel_i; i<elem_i.size()-1; i++)               elem_i_new.push_back(elem_i[i]+mesh.elem_i[mesh.nel_i]);
    for (int i=mesh.nel_i; i<=mesh.elem_i.size()-1; i++)    elem_i_new.push_back(mesh.elem_i[i]+(int)elem_v.size());

    // ELEM_V
    deque<int> elem_v_new;
    for (int i=0; i<elem_i[nel_i]; i++)
        if (elem_v[i]>=nno_i)                              elem_v_new.push_back(elem_v[i]+mesh.nno_i);
        else                                               elem_v_new.push_back(elem_v[i]);

    for (int i=0; i<mesh.elem_i[mesh.nel_i]; i++)
        if (mesh.elem_v[i]>=mesh.nno_i)                    elem_v_new.push_back(mesh.elem_v[i]+nodes.size());
        else                                               elem_v_new.push_back(mesh.elem_v[i]+nno_i);

    for (int i=elem_i[nel_i]; i<elem_v.size(); i++)
        if (elem_v[i]>=nno_i)                              elem_v_new.push_back(elem_v[i]+mesh.nno_i);
        else                                               elem_v_new.push_back(elem_v[i]);

    for (int i=mesh.elem_i[mesh.nel_i]; i<mesh.elem_v.size(); i++)
        if (mesh.elem_v[i]>=mesh.nno_i)                    elem_v_new.push_back(mesh.elem_v[i]+nodes.size());
        else                                               elem_v_new.push_back(mesh.elem_v[i]+nno_i);

    nno_i += mesh.nno_i;
    nel_i += mesh.nel_i;
    nodes = nodes_new;
    elem_i = elem_i_new;
    elem_v = elem_v_new;

//    face.nodes = elem_v_new;

    removeDoublePoints();

    return *this;
}

UNSTRUCTMESH UNSTRUCTMESH::operator+(const UNSTRUCTMESH &unstructMesh)
{
  UNSTRUCTMESH tmp = *this;
  tmp += unstructMesh;
  return tmp;
}

//UNSTRUCTMESH & UNSTRUCTMESH::operator-=(const UNSTRUCTMESH &mesh)
//{
//    // Define the POLYGON boundary
//    POLYGON extPoly(mesh.extBoundary);
//
//    extPoly.writeCtrPts("test.dat");
//
//    deque<int> nodesInside;
//
//    // list all the nodes which are inside the external boundary of "mesh"
//    for (int i=0; i<nodes.size(); i++)
//        if (extPoly.pointInsidePolygon(nodes[i].pt)==true)
//            nodesInside.push_back(i);
//
//    // remove all the elements which contain at least one of the nodesInside
//    int e=0; bool del=false;
//    while (e<elems.size())
//    {
//        for (int i=0; i<nodesInside.size(); i++)
//            for (int j=0; j<elems[e].node.size(); j++)
//                if (elems[e].node[j]==nodesInside[i])
//                {
//                    elems.erase(elems.begin() + e);
//                    del = true;
//                }
//
//        if (del==false)    e++;
//        del=false;
//    }
//
//    return *this;
//}

//UNSTRUCTMESH UNSTRUCTMESH::operator-(const UNSTRUCTMESH &unstructMesh)
//{
//    UNSTRUCTMESH tmp = *this;
//    tmp -= unstructMesh;
//    return tmp;
//}


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
  if(faces.size()>0)
    for (int i=0; i<faces.size(); i++)
      for (int j=0; j<faces[i].node.size()-1; j++)
        if (faces[i].node[j]==ind)      faces[i].node[j] = nodes.size()-1;
        else if (faces[i].node[j]>ind)  faces[i].node[j]--;



  // update nodes in faces; //GUSTAVO 09-13-2016
//  for (int i=0; i<faces.size(); i++)
//    for (int j=0; i<faces[i].node.size(); j++)
//      if (faces[i].node[j]==ind)      faces[i].node[j] = nodes.size()-1;
//      else if (faces[i].node[j]>ind)   faces[i].node[j]--;
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
}

void UNSTRUCTMESH::moveFaceEnd(int ind)
{
    // move node to the end of the deque
    int nNodes = face_i[ind+1]-face_i[ind];

    // save the face into a temp variable
    deque<int> tmpFace;
    for (int i=0; i<nNodes; i++)
        tmpFace.push_back(face_v[face_i[ind]+i]);

    int tmpFace_elem1=face_elem[2*ind];
    int tmpFace_elem2=face_elem[(2*ind)+1];

    string tmpFace_name;
    //strcpy (tmpFace_name,face_name[ind]);
    tmpFace_name= face_name[ind];

    // remove face from the list
    for (int i=0; i<nNodes; i++)
        face_v.erase(face_v.begin() + face_i[ind]);
    // remove face from elem_i list
    face_i.erase(face_i.begin()+ind);

    face_elem.erase(face_elem.begin() + (2*ind));
    face_elem.erase(face_elem.begin() + (2*ind) + 1);

    face_name.erase(face_name.begin()+ind);

    // update startNewElem from ind to the end
    for (int i=ind; i<face_i.size(); i++)
        face_i[i] -= nNodes;

    // insert the face at the end of the list
    for (int i=0; i<nNodes; i++)
        face_v.push_back(tmpFace[i]);

    face_elem.push_back(tmpFace_elem1);
    face_elem.push_back(tmpFace_elem2);

    face_name.push_back(tmpFace_name);

    // insert elem_i for the last element
    face_i.push_back(face_i[face_i.size()-1]+nNodes);
}


void UNSTRUCTMESH::moveFaceEnd_old(int ind)
{
    // move node to the end of the deque
    int nNodes = elem_i[ind+1]-elem_i[ind];

    // save the face into a temp variable
    FACE tmpface;
    tmpface = faces[ind];

    // remove face from the list
    faces.erase (faces.begin()+ind);;

    // insert the element at the end of the list
    faces.push_back(tmpface);

}


void UNSTRUCTMESH::removeDoublePoints()
{
    // Based on how meshes are added, double points can be on the boundaries only
    for (int i=nno_i; i<nodes.size(); i++)
        for (int j=i+1; j<nodes.size(); j++)
            if (nodes[j].pt == nodes[i].pt)
            {
                nodes.erase(nodes.begin()+j);
                // change nodes index in the boundary elements
                for (int e=elem_i[nel_i]; e<elem_v.size(); e++)
                    if (elem_v[e]==j)    elem_v[e] = i;
                    else if (elem_v[e]>j) elem_v[e]--;
            }

//  for (int i=0; i<nodes.size(); i++)
//    for (int j=i+1; j<nodes.size(); j++)
//      if (nodes[j].pt == nodes[i].pt)
//      {
//        nodes.erase(nodes.begin()+j);
//        // change nodes index in the boundary elements
//        for (int e=0; e<elem_v.size(); e++)
//          if (elem_v[e]==j) elem_v[e] = i;
//          else if (elem_v[e]>j) elem_v[e]--;
//      }

    // nodes and elements that are now internal should be put outside the boundary list
}



void UNSTRUCTMESH::findNodeNeighbors()
{
    // identify neighbours of every point
    for (int i=0; i<nodes.size(); i++)
        for (int j=0; j<elem_i.size()-1; j++)
        {
            for (int n=0; n<elem_i[j+1]-elem_i[j]; n++)
                if (elem_v[elem_i[j]+n] == i)
                    for (int nn=0; nn<elem_i[j+1]-elem_i[j]; nn++)
                        nodes[i].neig.push_back(elem_v[elem_i[j]+nn]);
        }

    // sort nodes indexes and remove duplicates
    for (int i=0; i<nodes.size(); i++)
    {
        nodes[i].neig.sort();
        nodes[i].neig.unique();
        nodes[i].neig.remove(i);
    }

    //With old elemnt objects:
//    for (int i=0; i<nodes.size(); i++)
//        for (int j=0; j<elems.size(); j++)
//        {
//            for (int n=0; n<elems[j].node.size(); n++)
//                if (elems[j].node[n] == i)
//                    for (int nn=0; nn<elems[j].node.size(); nn++)
//                        nodes[i].neig.push_back(elems[j].node[nn]);
//        }
//
//    // sort nodes indexes and remove duplicates
//    for (int i=0; i<nodes.size(); i++)
//    {
//        nodes[i].neig.sort();
//        nodes[i].neig.unique();
//        nodes[i].neig.remove(i);
//    }
}

void UNSTRUCTMESH::setMarker(const deque<POINT> &points, int marker)
{
  for (int j=0; j<nodes.size(); j++)
    nodes[j].marker = 0;

  for (int i=0; i<points.size(); i++)
    for (int j=0; j<nodes.size(); j++)
      if (nodes[j].pt == points[i])
      {
        nodes[j].marker = marker;
        break;
      }
}

void UNSTRUCTMESH::elementObject()          //GUSTAVO 09-05-2016
{
  if(elems.size()==0)
  {
    for(int i=0; i<elem_i.size()-1; i++)
    {
      ELEMENT tmp_element;

      for(int j=elem_i[i]; j<elem_i[i+1] ; j++)
        tmp_element.node.push_back(elem_v[j]);

      elems.push_back(tmp_element);
    }
  }
  else
    cout<<"Inside: elementObject: The objects elements were already created!" <<endl;
}

bool UNSTRUCTMESH::CompareElementFace(ELEMENT &element_, FACE &face_)          //GUSTAVO 09-05-2016
{
    int count = 0;

    for(int j=0; j<element_.node.size() ; j++)
      for(int i=0; i<face_.node.size() ; i++)
        if(element_.node[j] == face_.node[i])   count++;

    if(count>3) return(true);           //rectangular faces and tetrahedral elements
    else return(false);

}

bool UNSTRUCTMESH::CompareFaceFace(FACE &face1, FACE &face2)          //GUSTAVO 09-05-2016
{
    int count = 0;

    for(int j=0; j<face1.node.size() ; j++)
      for(int i=0; i<face2.node.size() ; i++)
        if(face1.node[j] == face2.node[i])   count++;

    if(count>3) return(true);           //rectangular faces and tetrahedral elements
    else return(false);

}

void UNSTRUCTMESH::removeDoubleFaces()
{


  for(int i=0; i<faces.size(); i++)
    for(int j=i+1; j<faces.size(); j++)
    {
      int count = 0;
      for(int k=0; k<faces[i].node.size(); k++)
        for(int kk=0; kk<faces[j].node.size(); kk++)
          if(faces[i].node[k] == faces[j].node[kk])
          {
            count++;
            if(count>3)
            {
              //              faces[j].elem_2=faces[i].elem_1;
              faces[i].elem_2=faces[j].elem_1;
              faces.erase (faces.begin()+j);
              count = 0;
              //              if(faces[i].elem_1==faces[i].elem_2) faces.erase (faces.begin()+i);
            }
          }
    }


  for(int i=0; i<faces.size(); i++)
     if(faces[i].elem_1==faces[i].elem_2)
       faces.erase (faces.begin()+i);

//------------------------------------------------------------------------------------------------------------------------
//      Using the CompareFaces function, Check if it is ok!!
//  for(int i=0; i<faces.size(); i++)
//  {
//    for(int j=i+1; j<faces.size(); j++)
//      if(CompareFaceFace(faces[i], faces[j]))
//      {
//        faces[i].elem_2=faces[j].elem_1;
//        faces.erase (faces.begin()+j);
//      }
//
//    if(faces[i].elem_1==faces[i].elem_2)
//      faces.erase (faces.begin()+i);
//  }

  //------------------------------------------------------------------------------------------------------------------------
}


bool compareFaces(const FACE &f1, const FACE &f2)
{
  return ((f1.node[0] == f2.node[0] && f1.node[1] == f2.node[1]) || (f1.node[0] == f2.node[1] && f1.node[1] == f2.node[0]));
}


void UNSTRUCTMESH::removeDoubleFaces2D_obj()
{
  deque<int> face_erase;

  clock_t begin = clock();

    for (int i=0; i<faces.size(); i++)
        for (int j=i+1; j<faces.size(); j++)
            if (compareFaces(faces[i],faces[j]))
            {
                faces[i].elem_2 = faces[j].elem_1;
                sprintf(faces[i].name,"fluid");
                faces.erase(faces.begin()+j);
                break;
            }

    clock_t end = clock();                                                                //for profiling
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;                           //for profiling
    cout<<"Took old struct: "<< elapsed_secs<< " seconds to do the find the double faces"<<endl;     //for profiling

  //------------------------------------------------------------------------------------------------------------------------
}

void UNSTRUCTMESH::removeDoubleFaces2D()
{
//  cout<<"Total of number of faces before "<< face_i.size() <<endl;
//  cout<<"Total of number of faces_elem before  "<< face_elem.size() <<endl;
//  cout<<"Total of number of face name before  "<< face_name.size() <<endl;


  deque<int> face_erase;

  clock_t begin = clock();                                                              //for profiling

  int nNodes=face_i[1]-face_i[0];
  int nfaces=face_i.size();

  //Coopering each face with the rest of the faces (2xloop)
  for(int i=0; i<nfaces-1; i++)
    for(int j=i+1; j<nfaces; j++)
    {
      int count = 0;
      int index_i=face_i[i];
      int index_j=face_i[j];

      int facei= face_v[index_i];
      int facej1= face_v[index_j];
      int facej2= face_v[index_j+1];

      if(facei == facej1 || facei == facej2)//(face_v[index_i] == face_v[index_j] || face_v[index_i] == face_v[index_j+1])
      {
        facei= face_v[index_i+1];
        if(facei == facej1 || facei == facej2) //(face_v[index_i+1] == face_v[index_j] || face_v[index_i+1] == face_v[index_j+1])
        {
          face_elem[(2*j)+1]=face_elem[(2*i)];
          face_name[j]="fluid";
          face_erase.push_back(i);
          break;
        }
      }
    }
  clock_t end = clock();                                                                //for profiling
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;                           //for profiling
  cout<<"Took new struct: "<< elapsed_secs<< " seconds to do the find the double faces"<<endl;     //for profiling
//  cout<<"Before erasing the extra faces"<<endl;

  cout<<"Number of face before removing old: "<< face_i.size()<<endl;
  //Erasing double faces!!!

  for(int i=face_erase.size()-1; i>=0; i--)
  {
      int j= face_erase[i];
      int index_j=face_i[j];

      // remove face from the list
      for (int k=0; k<nNodes; k++)
        face_v.erase(face_v.begin() + index_j);

      face_i.erase(face_i.begin()+j);
      face_name.erase(face_name.begin()+j);

      // update startNewElem from ind to the end
      for (int ii=j; ii<face_i.size(); ii++)
        face_i[ii] -= nNodes;

      face_elem.erase(face_elem.begin() + (2*j)+1);
      face_elem.erase(face_elem.begin() + (2*j) );




  }
}

void UNSTRUCTMESH::findFaces2D()
{
  int count=0;
  for(int i=0; i<elem_i.size()-1; i++)
  {
    //Each element has 4 faces;
    // As this is done in 2D: each faces consists of 2 vertices.
    // FACE 1
    face_v.push_back(elem_v[elem_i[i]]);
    face_v.push_back(elem_v[elem_i[i]+1]);
    face_i.push_back(int(count*2));
    face_elem.push_back(i+1);
    face_elem.push_back(int(0));
    face_name.push_back("noname");
    count++;
    // FACE 2
    face_v.push_back(elem_v[elem_i[i]+1]);
    face_v.push_back(elem_v[elem_i[i]+2]);
    face_i.push_back(int(count*2));
    face_elem.push_back(i+1);
    face_elem.push_back(int(0));
    face_name.push_back("noname");
    count++;
    // FACE 3
    face_v.push_back(elem_v[elem_i[i]+2]);
    face_v.push_back(elem_v[elem_i[i]+3]);
    face_i.push_back(int(count*2));
    face_elem.push_back(i+1);
    face_elem.push_back(int(0));
    face_name.push_back("noname");
    count++;
    // FACE 4
    face_v.push_back(elem_v[elem_i[i]+3]);
    face_v.push_back(elem_v[elem_i[i]]);
    face_i.push_back(int(count*2));
    face_elem.push_back(i+1);
    face_elem.push_back(int(0));
    face_name.push_back("noname");
    count++;
  }

  removeDoubleFaces2D();

}

void UNSTRUCTMESH::findFaces2D_obj()          //GUSTAVO 09-05-2016
{
  //  findNodeNeighbors();
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
        faces.push_back(tmp_face);
      }

      FACE tmp_face;
      tmp_face.node.push_back(elem_v[fstNodeInd]);
      tmp_face.node.push_back(elem_v[lstNodeInd]);
      tmp_face.elem_1 = i+1;
      tmp_face.elem_2 = 0;
      faces.push_back(tmp_face);
  }

    cout<<"Number of face before removing old: "<< faces.size()<<endl;
/*    


      //Each element has 4 faces;
      FACE tmp_face1;    //Face1, in 2D the face has 2 nodes;
      FACE tmp_face2;    //Face2, in 2D the face has 2 nodes;
      FACE tmp_face3;    //Face3, in 2D the face has 2 nodes;
      FACE tmp_face4;    //Face4, in 2D the face has 2 nodes;

      for(int j=elem_i[i]; j<elem_i[i+1]-2 ; j++)
      {
        tmp_face1.node.push_front(elem_v[j]);
        tmp_face2.node.push_front(elem_v[j+1]);
        tmp_face3.node.push_front(elem_v[j+2]);
      }
      tmp_face4.node.push_front(elem_v[elem_i[i+1]-1]);
      tmp_face4.node.push_front(elem_v[elem_i[i]]);

      tmp_face1.elem_1= i+1;
      tmp_face1.elem_2= 0;
      tmp_face2.elem_1= i+1;
      tmp_face2.elem_2= 0;
      tmp_face3.elem_1= i+1;
      tmp_face3.elem_2= 0;
      tmp_face4.elem_1= i+1;
      tmp_face4.elem_2= 0;

      faces.push_back(tmp_face1);
      faces.push_back(tmp_face2);
      faces.push_back(tmp_face3);
      faces.push_back(tmp_face4);
      //------------------------------------------------------------------
//      int num_nodes=elem_i[i+1]-elem_i[i];
//
//      for(int j=0; j<num_nodes; j++)
//      {
//        FACE tmp_face;
//        for(int k=elem_i[i]; k<elem_i[i+1]-2; k++)
//          if(j+k<elem_i[i+1])
//            tmp_face.node.push_back(elem_v[j+k]);
//          else
//            tmp_face.node.push_back(elem_v[elem_i[i]]);
//
//        tmp_face.elem_1= i+1;
//        tmp_face.elem_2= 0;
//        faces.push_back(tmp_face);
//      }
    }

*/

    removeDoubleFaces2D_obj();

//    for(int i=0; i<faces.size(); i++)
//      if(faces[i].elem_2 != 0)
//        strcpy (faces[i].name,"fluid");

  }
  else
    cout<<"Inside: findFaces2D(): The objects faces were already created!" <<endl;
}


void UNSTRUCTMESH::replacePoints(const deque<POINT> &oldPoints, const deque<POINT> &newPoints)
{
    for (int i=0; i<newPoints.size(); i++)
        for (int j=0; j<nodes.size(); j++)
                if (nodes[j].pt == oldPoints[i])
                    nodes[j].pt = newPoints[i];
}



// Laplacian smoother
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
// Laplacian smoother
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

void UNSTRUCTMESH::ScaleMesh(const double scaling_factor)  //Gustavo 24/06/2016
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

//  int iStart = 0;
//  if (onlyBoundary==true)
//    iStart = nno_i;

    for (int i=0; i<nodes.size(); i++)
        fprintf(fp, "%1.10le\t%1.10le\t%1.10le\n", nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
    fprintf(fp, "\n");


    int iStart = 0;
    if (onlyBoundary==true)
        iStart = nel_i;


//  if (onlyBoundary==true)
//  {
//    printf("number: %d\n",nno_i);
//    for (int i=elem_i[iStart]; i<elem_v.size(); i++)
//      elem_v[i] -= nno_i;
//  }

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
  FILE *fp;
  if ((fp = fopen(name, "wt")) == NULL)
  {
    printf("couldn't open %s for writing\n", name);
    return;
  }

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

  fprintf(fp, "      ELEMENTS/CELLS 2.4.6\n");

  int cStart = 0;
  if (onlyBoundary==true)
    cStart = nel_i;

  for (int c=cStart; c<elem_i.size()-1; c++)
  {
    int nNodes = elem_i[c+1]-elem_i[c];

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

//void UNSTRUCTMESH::writeFluentMshNew(const char *filename, int ndim, deque<string> bc_names)
//{
//  FILE * fp;
//
//  printf(" > starting to write fluent case file: %s\n", filename);
//
//  if ( (fp=fopen(filename,"w"))==NULL ) {
//    printf(" > could not open file %s\n", filename);
//    exit(-1);
//  }
//
//  fprintf(fp,"(0 \"Project to Fluent\")\n");
//  fprintf(fp,"\n");
//  fprintf(fp,"(0 \"Dimensions:\")\n");
//  if(ndim==2)       fprintf(fp,"(2 2)\n");
//  else              fprintf(fp,"(2 3)\n");
//  fprintf(fp,"\n");
//
//  //
//  // header and global dimensions
//  //
//
//  int vertices = nodes.size(), elements = elem_i.size()-1, nfaces = faces.size();
//
////  for (int n=0; n<nbl; ++n)
////  {
////    vertices += 0;//(bl[n].imax+1)*(bl[n].jmax+1)*(bl[n].kmax+1);
////    elements += 0;//bl[n].imax*bl[n].jmax*bl[n].kmax;
////    nfaces += 0;//bl[n].imax*bl[n].jmax*(bl[n].kmax+1)+bl[n].imax*(bl[n].jmax+1)*bl[n].kmax+(bl[n].imax+1)*bl[n].jmax*bl[n].kmax;
////  }
//
//  fprintf(fp,"(0 \"Grid:\")\n");
////  fprintf(fp,"\n");
////  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);
//
//  fprintf(fp,"\n");
//
//  fprintf(fp,"(0 \"Nodes:\")\n");
//  if(ndim==2)
//  {
//    fprintf(fp,"(10 (0 %x %x 1 2))\n",1, vertices);
//    fprintf(fp, "(10 (1 %x %x 1 2)(\n",1, vertices);
//    for (int i=0; i<nodes.size(); i++)
//      fprintf(fp, "%20.11le%20.11le\n", nodes[i].pt.x, nodes[i].pt.y);
//  }
//  else
//  {
//    fprintf(fp,"(10 (0 %x %x 1 3))\n",1, vertices);
//    fprintf(fp, "(10 (1 %x %x 1 3)(\n",1, vertices);
//    for (int i=0; i<nodes.size(); i++)
//      fprintf(fp, "%20.11le%20.11le%20.11le\n", nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
//  }
//  fprintf(fp, "))\n");
//  fprintf(fp,"\n");
//
//
//  //
//  // start of faces
//  //
//  int fa_element_type = 4; // hexahedral = 4
//  int fa_type = 0; //1; // active = 1 - inactive = 32
////  nfaces=106;
//  fprintf(fp,"(0 \"Faces:\")\n");
//  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);
//
//  int ind_faces[n_zone];      //Index of the current zone
//
//  //==================================================================
//  //This part of the script is old, before all the names of the boundaries are known
////  deque<string> bc_names;     //String with the names of the boundary zones
////
////
////  for(int i=0; i<bc_name.size(); i++)
////    bc_names.push_back(bc_name[i]);
////
////  for(int i=0; i<bc_names.size(); i++)
////      cout<<"Name of the boundaries: "<< bc_names[i] <<endl;
//  //==================================================================
//
//  //Checking if all the faces have a name!
//  bool zone_noname=false;
//  for(int i=0; i<nfaces; i++)
//    if(strcmp (faces[i].name, "noname") == 0)
//    {
//      zone_noname=true;
//      cout<<"WARNING: At least one of the faces was not named! All the face zones will be named together as: boundaries"<<endl;
//      break;
//    }
//
//  if(zone_noname)
//  {
//    int faces_wall=0;
//    for(int i=0; i<nfaces; i++)
//      if(faces[i].elem_2==0) faces_wall++;
//
//    //WALL FACES
//    ind_faces[0]=3;
//
//    fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[0] ,1, faces_wall);
//    for(int i=0; i<nfaces; i++)
//    {
//      //      fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 , faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
//     if(ndim==2)    if(faces[i].elem_2==0)      fprintf(fp,"2 %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].elem_1, faces[i].elem_2);
//     else           if(faces[i].elem_2==0)      fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 ,faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
//
//    }
//
//    fprintf(fp, "))\n");
//
//    //INNER FACES
//    ind_faces[1]=5;
//    fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[1] ,faces_wall+1, nfaces);
//    for(int i=0; i<nfaces; i++)
//    {
////            fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 , faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
//      if(ndim==2)    if(faces[i].elem_2>0)      fprintf(fp,"2 %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].elem_1, faces[i].elem_2);
//      else           if(faces[i].elem_2>0)      fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 ,faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
//
//    }
//
//    fprintf(fp, "))\n");
//  }
//  else
//  {
//    int nfa_i[n_zone];        //Number of faces in the current zone
//
//    /*
//    //==================================================================
//    //This part of the script is old, before all the names of the boundaries are known
//    if(ndim!=2)
//    {
//      bc_names.push_back("top");
//      bc_names.push_back("bottom");
//    }
//    else
//      bc_names.push_back("inlet");
//
////    bc_names.push_back("inlet");
//    bc_names.push_back("fluid");
//    cout<<"Size of BC names: "<< bc_names.size() <<endl;
//
//
//    //Finding the names of the rest of the boundaries!!!!
//    for(int i=0; i<nfaces; i++)
//    {
//      int count=0;
//      for(int j=0; j<bc_names.size(); j++)
//        if(strcmp (faces[i].name, bc_names[j].c_str()) != 0)
//        {
//          for(int k=0; k<bc_names.size(); k++)
//            if(strcmp (faces[i].name , bc_names[k].c_str()) != 0) //if(strcmp (bc_names[j].c_str() , bc_names[k].c_str()) != 0)
//            {
//              count++;
//              break;
//            }
//          if(count==bc_names.size())
//            break;
//        }
//
//      if(count==bc_names.size())
//        bc_names.push_back(faces[i].name);
//    }
//    //Making the fluid boundary the last of the deque!!!!
//    for(int k=0; k<bc_names.size(); k++)
//      if(strcmp (bc_names[k].c_str(), "fluid") == 0)
//      {
//        bc_names.erase (bc_names.begin()+k);
//        bc_names.push_back("fluid");
//        break;
//      }
//        //==================================================================
//      */
//
//
//    int total_faces=0;
//
//
//    for(int n=0 ; n<bc_names.size(); n++)
//    {
//      // Defining the index of the face zone
//      ind_faces[n]=n+3;                                 // the index starts in 3
//
//      //Counting the number of faces of the present face zone, needed in the declaration for the face
//      nfa_i[n]=0;
//      for(int i=0; i<nfaces; i++)
//        if(strcmp (faces[i].name, bc_names[n].c_str()) == 0)
//          nfa_i[n]++;
//
//      // Declaration of the face zone
//      if(n!=bc_names.size()-1)fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[n] ,total_faces+1, total_faces+nfa_i[n]);
//      else fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[n]+1 ,total_faces+1, total_faces+nfa_i[n]);
//
//      // Giving the face connectivity
//      for(int i=0; i<nfaces; i++)
//        if(strcmp (faces[i].name, bc_names[n].c_str()) == 0)
//        {
//          if(ndim==2)   fprintf(fp,"2 %x %x %x %x  \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].elem_1, faces[i].elem_2);
//          else   fprintf(fp,"4 %x %x %x %x %x %x  \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 ,  faces[i].node[2]+1 ,faces[i].node[3]+1 ,  faces[i].elem_1, faces[i].elem_2);
//        }
//      fprintf(fp, "))\n");
//
//      total_faces+=nfa_i[n];
//      cout<<"Number of faces in boundary "<<bc_names[n] <<" "<<nfa_i[n]<<endl;
//    }
//
//    cout<<"Total number of faces in the numerical domain "<<total_faces<<endl;
//
//
//  }
//
//
//  fprintf(fp,"\n");
//
//  //
//  // cells
//  //
//  int cv_element_type = 4; // hexahedral = 4
//  int cv_type = 1; //1; // active = 1 - inactive = 32
//  char this_zone[128];
//  sprintf(this_zone,"fluid");
//  fprintf(fp,"(0 \"Cells:\")\n");
//  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"\n");
//
//  if(zone_noname)
//  {
//    fprintf(fp,"(0 \"Zones:\")\n");
//    fprintf(fp,"(45 (2 fluid fluid 0))\n");
//    fprintf(fp,"(45 (%i wall boundaries 0))\n", ind_faces[0]);
//    fprintf(fp,"(45 (%i interior default-interior 0))\n", ind_faces[1]);
//  }
//  else
//  {
//    fprintf(fp,"(0 \"Zones:\")\n");
//    fprintf(fp,"(45 (2 fluid fluid 0))\n");
//    for(int n=0 ; n<bc_names.size()-1; n++)
//    {
//      fprintf(fp,"(45 (%i wall %s 0))\n", ind_faces[n], bc_names[n].c_str());
//    }
//    fprintf(fp,"(45 (%i interior default-interior 0))\n", ind_faces[bc_names.size()-1]+1);
//  }
////  fprintf(fp,"(12 (%x %x %x %x))\n", 0,1,elements,cv_type);
////  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
////  fprintf(fp,"(45 (%d fluid %s)())\n",  4001,this_zone);
//  fprintf(fp,"\n");
//
//  fclose(fp);
//
//  printf(" > Fluent case file written: %s\n", filename);
//}


void UNSTRUCTMESH::writeFluentMsh(const char *filename, int ndim, deque<string> bc_names)
{
  FILE * fp;

  printf(" > starting to write fluent case file: %s\n", filename);

  if ( (fp=fopen(filename,"w"))==NULL ) {
    printf(" > could not open file %s\n", filename);
    exit(-1);
  }

  fprintf(fp,"(0 \"Project to Fluent\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(0 \"Dimensions:\")\n");
  if(ndim==2)       fprintf(fp,"(2 2)\n");
  else              fprintf(fp,"(2 3)\n");
  fprintf(fp,"\n");

  //
  // header and global dimensions
  //

  int vertices = nodes.size(), elements = elem_i.size()-1, nfaces =  face_i.size();

  cout<<"Number of faces: "<< nfaces <<endl;
  fprintf(fp,"(0 \"Grid:\")\n");
//  fprintf(fp,"\n");
//  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

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


  //
  // start of faces
  //
  int fa_element_type = 4; // hexahedral = 4
  int fa_type = 0; //1; // active = 1 - inactive = 32
//  nfaces=106;
  fprintf(fp,"(0 \"Faces:\")\n");
  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

  int ind_faces[n_zone];      //Index of the current zone
  //Checking if all the faces have a name!
  bool zone_noname=false;
  for(int i=0; i<nfaces; i++)
    if(face_name[i].compare("noname")==0)//(strcmp (faces[i].name, "noname") == 0)
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }

  if(zone_noname)
  {
    int faces_wall=0;
    for(int i=0; i<nfaces; i++)
      if(face_elem[(2*i)+1]==0) faces_wall++;


//    cout<<"Number of wall faces:  "<<faces_wall<< " Dimension of the mesh "<< ndim<< " Number of faces: "<< nfaces <<endl;

    //WALL FACES
    ind_faces[0]=3;

    fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[0] ,1, faces_wall);
    for(int i=0; i<nfaces; i++)
    {
      //      fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 , faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
     if(ndim==2)
       if(face_elem[(2*i)+1]==0)      fprintf(fp,"2 %x %x %x %x \n",           face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

     if(ndim==3)
       if(face_elem[(2*i)+1]==0)      fprintf(fp,"4 %x %x %x %x %x %x \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

    }

    fprintf(fp, "))\n");

    //INNER FACES
    ind_faces[1]=5;
    fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[1] ,faces_wall+1, nfaces);
    for(int i=0; i<nfaces; i++)
    {
//            fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 , faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
      if(ndim==2)
        if(face_elem[(2*i)+1]>0)      fprintf(fp,"2 %x %x %x %x \n",           face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);
      if(ndim==3)
        if(face_elem[(2*i)+1]>0)      fprintf(fp,"4 %x %x %x %x %x %x \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

    }

    fprintf(fp, "))\n");
  }
  else
  {
    int nfa_i[n_zone];        //Number of faces in the current zone
    int total_faces=0;


    for(int n=0 ; n<bc_names.size(); n++)
    {
      // Defining the index of the face zone
      ind_faces[n]=n+3;                                 // the index starts in 3

      //Counting the number of faces of the present face zone, needed in the declaration for the face
      nfa_i[n]=0;
      for(int i=0; i<nfaces; i++)
        if(face_name[i].compare(bc_names[n])==0) //(strcmp (faces[i].name, bc_names[n].c_str()) == 0)
          nfa_i[n]++;

      // Declaration of the face zone
      if(n!=bc_names.size()-1)fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[n] ,total_faces+1, total_faces+nfa_i[n]);
      else fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[n]+1 ,total_faces+1, total_faces+nfa_i[n]);

      // Giving the face connectivity
      for(int i=0; i<nfaces; i++)
        if(face_name[i].compare(bc_names[n])==0) //if(strcmp (faces[i].name, bc_names[n].c_str()) == 0)
        {
          if(ndim==2)   fprintf(fp,"2 %x %x %x %x  \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);
          else   fprintf(fp,"4 %x %x %x %x %x %x  \n",      face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);
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
  int cv_element_type= 4;    // hexahedral = 4
  if(ndim==2)   cv_element_type = 3;    // quadrilateral = 3


  int cv_type = 1; //1; // active = 1 - inactive = 32
  char this_zone[128];
  sprintf(this_zone,"fluid");
  fprintf(fp,"(0 \"Cells:\")\n");
  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
  fprintf(fp,"\n");

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
//  fprintf(fp,"(12 (%x %x %x %x))\n", 0,1,elements,cv_type);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"(45 (%d fluid %s)())\n",  4001,this_zone);
  fprintf(fp,"\n");

  fclose(fp);

  printf(" > Fluent case file written: %s\n", filename);
}


void UNSTRUCTMESH::writeFluentMsh_obj(const char *filename, int ndim, deque<string> bc_names)
{
  FILE * fp;

  printf(" > starting to write fluent case file: %s\n", filename);

  if ( (fp=fopen(filename,"w"))==NULL ) {
    printf(" > could not open file %s\n", filename);
    exit(-1);
  }

  fprintf(fp,"(0 \"Project to Fluent\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(0 \"Dimensions:\")\n");
  if(ndim==2)       fprintf(fp,"(2 2)\n");
  else              fprintf(fp,"(2 3)\n");
  fprintf(fp,"\n");

  //
  // header and global dimensions
  //

  int vertices = nodes.size(), elements = elem_i.size()-1, nfaces =  faces.size();  //face_i.size();

  cout<<"Number of faces: "<< nfaces <<endl;
  fprintf(fp,"(0 \"Grid:\")\n");
//  fprintf(fp,"\n");
//  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

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


  //
  // start of faces
  //
  int fa_element_type = 4; // hexahedral = 4
  int fa_type = 0; //1; // active = 1 - inactive = 32
//  nfaces=106;
  fprintf(fp,"(0 \"Faces:\")\n");
  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);

  int ind_faces[n_zone];      //Index of the current zone
  //Checking if all the faces have a name!
  bool zone_noname=false;
  for(int i=0; i<nfaces; i++)
    if(strcmp (faces[i].name, "noname") == 0)   //(face_name[i].compare("noname")==0)//
    {
      zone_noname=true;
      cout<<"WARNING: At least one of the faces was not named! All the face zones will be named together as: boundaries"<<endl;
      break;
    }

  if(zone_noname)
  {
    int faces_wall=0;
    for(int i=0; i<nfaces; i++)
      if(faces[i].elem_2==0) faces_wall++;   //(face_elem[(2*i)+1]==0) faces_wall++;


//    cout<<"Number of wall faces:  "<<faces_wall<< " Dimension of the mesh "<< ndim<< " Number of faces: "<< nfaces <<endl;

    //WALL FACES
    ind_faces[0]=3;

    fprintf(fp,"(13 (%x %x %x 3 0)(\n", ind_faces[0] ,1, faces_wall);
    for(int i=0; i<nfaces; i++)
    {

     if(ndim==2)
       if(faces[i].elem_2==0)           fprintf(fp,"2 %x %x %x %x \n",          faces[i].node[0]+1, faces[i].node[1]+1, faces[i].elem_1, faces[i].elem_2);
       //(face_elem[(2*i)+1]==0)      fprintf(fp,"2 %x %x %x %x \n",           face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

     if(ndim==3)
       if(faces[i].elem_2==0)           fprintf(fp,"4 %x %x %x %x %x %x \n",    faces[i].node[0]+1, faces[i].node[1]+1, faces[i].node[2]+1, faces[i].node[3]+1, faces[i].elem_1, faces[i].elem_2);
              //(face_elem[(2*i)+1]==0)      fprintf(fp,"4 %x %x %x %x %x %x \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

    }

    fprintf(fp, "))\n");

    //INNER FACES
    ind_faces[1]=5;
    fprintf(fp,"(13 (%x %x %x 2 0)(\n", ind_faces[1] ,faces_wall+1, nfaces);
    for(int i=0; i<nfaces; i++)
    {
//    with Obj          fprintf(fp,"4 %x %x %x %x %x %x \n",  faces[i].node[0]+1 ,faces[i].node[1]+1 , faces[i].node[2]+1 , faces[i].node[3]+1 , faces[i].elem_1, faces[i].elem_2);
//    with list         fprintf(fp,"4 %x %x %x %x %x %x \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

      if(ndim==2)
        if(faces[i].elem_2>0)           fprintf(fp,"2 %x %x %x %x \n",          faces[i].node[0]+1, faces[i].node[1]+1, faces[i].elem_1, faces[i].elem_2);
               //(face_elem[(2*i)+1]>0)      fprintf(fp,"2 %x %x %x %x \n",           face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);
      if(ndim==3)
        if(faces[i].elem_2>0)           fprintf(fp,"4 %x %x %x %x %x %x \n",    faces[i].node[0]+1, faces[i].node[1]+1, faces[i].node[2]+1, faces[i].node[3]+1, faces[i].elem_1, faces[i].elem_2);
                      //(face_elem[(2*i)+1]>0)      fprintf(fp,"4 %x %x %x %x %x %x \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);

    }

    fprintf(fp, "))\n");
  }
  else
  {
    int nfa_i[n_zone];        //Number of faces in the current zone
    int total_faces=0;


    for(int n=0 ; n<bc_names.size(); n++)
    {
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
          if(ndim==2)   fprintf(fp,"2 %x %x %x %x \n",          faces[i].node[0]+1, faces[i].node[1]+1, faces[i].elem_1, faces[i].elem_2);
          //fprintf(fp,"2 %x %x %x %x  \n",     face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);
          else          fprintf(fp,"4 %x %x %x %x %x %x \n",    faces[i].node[0]+1, faces[i].node[1]+1, faces[i].node[2]+1, faces[i].node[3]+1, faces[i].elem_1, faces[i].elem_2);
          //fprintf(fp,"4 %x %x %x %x %x %x  \n",      face_v[face_i[i]]+1 ,face_v[face_i[i]+1]+1 , face_v[face_i[i]+2]+1 ,face_v[face_i[i]+3]+1 , face_elem[(2*i)], face_elem[(2*i)+1]);
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
  int cv_element_type= 4;    // hexahedral = 4
  if(ndim==2)   cv_element_type = 3;    // quadrilateral = 3


  int cv_type = 1; //1; // active = 1 - inactive = 32
  char this_zone[128];
  sprintf(this_zone,"fluid");
  fprintf(fp,"(0 \"Cells:\")\n");
  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
  fprintf(fp,"\n");

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
//  fprintf(fp,"(12 (%x %x %x %x))\n", 0,1,elements,cv_type);
//  fprintf(fp,"(12 (%x %x %x %x %x))\n", 2,1,elements,cv_type,cv_element_type);
//  fprintf(fp,"(45 (%d fluid %s)())\n",  4001,this_zone);
  fprintf(fp,"\n");

  fclose(fp);

  printf(" > Fluent case file written: %s\n", filename);
}

void UNSTRUCTMESH::writeFluentMsh_Original(char * filename, int nbl)//, block *bl)
{
  FILE * fp;

  printf(" > starting to write fluent case file: %s\n", filename);

  if ( (fp=fopen(filename,"w"))==NULL ) {
    printf(" > could not open file %s\n", filename);
    exit(-1);
  }

  fprintf(fp,"(0 \"Project to Fluent\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(0 \"Dimensions:\")\n");
  fprintf(fp,"(2 3)\n");
  fprintf(fp,"\n");

  //
  // header and global dimensions
  //

  int vertices = 0, elements = 0, nfaces = 0;

  for (int n=0; n<nbl; ++n)
  {
    vertices += 0;//(bl[n].imax+1)*(bl[n].jmax+1)*(bl[n].kmax+1);
    elements += 0;//bl[n].imax*bl[n].jmax*bl[n].kmax;
    nfaces += 0;//bl[n].imax*bl[n].jmax*(bl[n].kmax+1)+bl[n].imax*(bl[n].jmax+1)*bl[n].kmax+(bl[n].imax+1)*bl[n].jmax*bl[n].kmax;
  }

  fprintf(fp,"(0 \"Grid:\")\n");
  fprintf(fp,"\n");
  fprintf(fp,"(12 (0 %x %x 0))\n",1, elements);
  fprintf(fp,"(13 (0 %x %x 0))\n",1, nfaces);
  fprintf(fp,"(10 (0 %x %x 0 3))\n",1, vertices);
  fprintf(fp,"\n");

  //
  // cells
  //
  int cv_element_type = 7; // polyhedra = 7
  int cv_type = 1; // active = 1 - inactive = 32
  char this_zone[128];
  sprintf(this_zone,"fluid");
  fprintf(fp,"(0 \"Cells:\")\n");
  fprintf(fp,"(12 (%x %x %x %x %x))\n", 4001,1,elements,cv_type,cv_element_type);
  fprintf(fp,"(45 (%d fluid %s)())\n",  4001,this_zone);
  fprintf(fp,"\n");

  //
  // start of faces
  //
  fprintf(fp,"(0 \"Faces:\")\n");
  fprintf(fp,"\n");

/*
  //
  // faces
  //
  int fa_element_type = 5; // 5 is polygonal; 4 is quadrilateral

  int n = 0, i, j, k;
  int this_zone_id = 11;
  int fa_offset = 1;
  int fa_number = bl[n].imax*bl[n].jmax;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, fa_number, (int)FLUENT_TYPE_WALL, fa_element_type);
  k = 0;
  for (i=0; i<bl[n].imax; ++i)
    for (j=0; j<bl[n].jmax; ++j)
    {
      fprintf(fp,"%x ", 4); // number of nodes of face
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+1) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 1);
      fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j + k + 1, 0);
    }
  fprintf(fp, "))\n");

  this_zone_id++;
  fa_offset = fa_number+1;
  fa_number = fa_offset + bl[n].imax*bl[n].jmax-1;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, fa_number, (int)FLUENT_TYPE_WALL, fa_element_type);
  k = bl[n].kmax;
  for (i=0; i<bl[n].imax; ++i)
    for (j=0; j<bl[n].jmax; ++j)
    {
      fprintf(fp,"%x ", 4); // number of nodes of face
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+1) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j + (k-1) + 1, 0);
    }
  fprintf(fp, "))\n");

  //################################################################################################
  this_zone_id++;
  fa_offset = fa_number+1;
  fa_number = fa_offset + bl[n].imax*bl[n].kmax-1;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, fa_number, (int)FLUENT_TYPE_WALL, fa_element_type);
  j = 0;
  for (i=0; i<bl[n].imax; ++i)
    for (k=0; k<bl[n].kmax; ++k)
    {
      fprintf(fp,"%x ", 4); // number of nodes of face
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 2);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 2);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j + k + 1, 0);
    }
  fprintf(fp, "))\n");

  this_zone_id++;
  fa_offset = fa_number+1;
  fa_number = fa_offset + bl[n].imax*bl[n].kmax-1;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, fa_number, (int)FLUENT_TYPE_WALL, fa_element_type);
  j = bl[n].jmax;
  for (i=0; i<bl[n].imax; ++i)
    for (k=0; k<bl[n].kmax; ++k)
    {
      fprintf(fp,"%x ", 4); // number of nodes of face
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 2);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 2);
      fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*(j-1) + k + 1, 0);
    }
  fprintf(fp, "))\n");

  //################################################################################################
  this_zone_id++;
  fa_offset = fa_number + 1;
  fa_number = fa_offset + bl[n].jmax*bl[n].kmax-1;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, fa_number, (int)FLUENT_TYPE_WALL, fa_element_type);
  i = 0;
  for (j=0; j<bl[n].jmax; ++j)
    for (k=0; k<bl[n].kmax; ++k)
    {
      fprintf(fp,"%x ", 4); // number of nodes of face
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 2);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 2);
      fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j + k + 1, 0);
    }
  fprintf(fp, "))\n");

  this_zone_id++;
  fa_offset = fa_number+1;
  fa_number = fa_offset + bl[n].jmax*bl[n].kmax-1;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, fa_number, (int)FLUENT_TYPE_WALL, fa_element_type);
  i = bl[n].imax;
  for (j=0; j<bl[n].jmax; ++j)
    for (k=0; k<bl[n].kmax; ++k)
    {
      fprintf(fp,"%x ", 4); // number of nodes of face
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 2);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 2);
      fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 1);
      fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*(i-1) + bl[n].kmax*j + k + 1, 0);
    }
  fprintf(fp, "))\n");

  //################################################################################################
  //################################################################################################
  this_zone_id++;
  fa_offset = fa_number+1;
  //  fa_number += bl[n].jmax*bl[n].kmax;
  fprintf(fp,"(13 (%x %x %x %x %x)(\n", this_zone_id, fa_offset, nfaces, (int)FLUENT_TYPE_INTERIOR, fa_element_type);

  for (i=0; i<bl[n].imax; ++i)
    for (j=0; j<bl[n].jmax; ++j)
      for (k=1; k<bl[n].kmax; ++k)
      {
        fprintf(fp,"%x ", 4); // number of nodes of face
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 1);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+1) + k + 1);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 1);
        fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j + k,
            bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j + k+1);
      }

  for (i=0; i<bl[n].imax; ++i)
    for (j=1; j<bl[n].jmax; ++j)
      for (k=0; k<bl[n].kmax; ++k)
      {
        fprintf(fp,"%x ", 4); // number of nodes of face
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 1);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+1) + (bl[n].kmax+1)*(j+0) + k + 2);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 2);
        fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*i + bl[n].kmax*(j-1) + k+1,
            bl[n].jmax*bl[n].kmax*i + bl[n].kmax*j     + k+1);
      }

  for (i=1; i<bl[n].imax; ++i)
    for (j=0; j<bl[n].jmax; ++j)
      for (k=0; k<bl[n].kmax; ++k)
      {
        fprintf(fp,"%x ", 4); // number of nodes of face
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 1);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+0) + k + 2);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 2);
        fprintf_(fp,"%x ", (bl[n].jmax+1)*(bl[n].kmax+1)*(i+0) + (bl[n].kmax+1)*(j+1) + k + 1);
        fprintf(fp,"%x %x\n", bl[n].jmax*bl[n].kmax*(i-1) + bl[n].kmax*j + k+1,
            bl[n].jmax*bl[n].kmax*i     + bl[n].kmax*j + k+1);
      }
  fprintf(fp, "))\n");

  //################################################################################################
  //################################################################################################

  fprintf(fp, "(39 (%d wall z0 1)())\n", 11);
  fprintf(fp, "(39 (%d wall z1 1)())\n", 12);
  fprintf(fp, "(39 (%d wall y0 1)())\n", 13);
  fprintf(fp, "(39 (%d wall y1 1)())\n", 14);
  fprintf(fp, "(39 (%d wall x0 1)())\n", 15);
  fprintf(fp, "(39 (%d wall x1 1)())\n", 16);
  fprintf(fp, "(39 (%d default-interior default-interior 1)())\n", 17);



  // start of nodes...
  fprintf(fp,"(0 \"Nodes:\")\n");
  fprintf(fp,"(10 (1 %x %x 1 3 )(\n",1,vertices);

  //
  // nodes
  //
  double scale = 1.0;//0.0254;//e-3;

  int indexBlock = 0;
  for (n=0; n<nbl; ++n)
  {
    for (i=0; i<=bl[n].imax; ++i)
      for (j=0; j<=bl[n].jmax; ++j)
        for (k=0; k<=bl[n].kmax; ++k)
          fprintf(fp, "%15.8e %15.8e %15.8e\n", bl[n].pt[i][j][k][0]*scale, bl[n].pt[i][j][k][1]*scale, bl[n].pt[i][j][k][2]*scale);

    indexBlock += (bl[n].imax+1)*(bl[n].jmax+1)*(bl[n].kmax+1);
  }

  fprintf(fp,"))\n");*/

  fclose(fp);

  printf(" > Fluent case file written: %s\n", filename);
}


//****************************************************************************************************************
//****************************************************************************************************************
//                                  Reading mesh files functions!!
//****************************************************************************************************************
//****************************************************************************************************************

void UNSTRUCTMESH::readFluentMsh2D(const char *name)
{
  //THis function needs to be combined with readGambitNeu2D, because this function only reads the faces: its connectivity!!
  //The elements needs to be know before handed (related to the nodes).

  FILE *fp;
  if ((fp = fopen(name, "rt")) == NULL)
  {
    printf("couldn't open %s for reading\n", name);
    exit(-1);
  }

  // Read header
  char line[1000];
  fscanf(fp, "%s\n", line);
  while (strcmp(line, "\"Nodes:\")") != 0)    fscanf(fp, "%s\n", line);

  int dummy1, dummy2, dummy3;
  int nvert, nfaces_bc, nfaces_fluid;
  string str1, str2, str3;

  cout<<"Before reading the vertices"<<endl;

  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nvert, &dummy2, str3.c_str());

  cout<<"Number of vertices: "<< nvert<< endl;

  while (strcmp(line, "\"Faces:\")") != 0)    fscanf(fp, "%s\n", line);
  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nfaces_bc, &dummy2, str3.c_str());
  cout<<"Starting to read bc faces: Number of total faces: "<< nfaces_bc<< endl;
  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nfaces_bc, &dummy2, str3.c_str());
//  fscanf(fp, "%s\n", line);

  int nodes[4];
  int count=0;

//  cout<<"Starting to read bc faces: Number of boundary faces: "<< nfaces_bc<< endl;

//  cout<<"Number of boundary faces: "<< nfaces_bc<< endl;


  for (int i = 0; i < nfaces_bc; i++)
  {
    fscanf(fp, "%10d%x%x%x%x\n", &dummy1, &nodes[0], &nodes[1], &nodes[2], &nodes[3]);
//    if (i%100 == 0) printf("%10d\t%lf\t%lf\t%lf\n", i, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);

//    for (int j=0; j<2; j++) face_v.push_back(nodes[j]-1);

    face_v.push_back(nodes[0]-1);
    face_v.push_back(nodes[1]-1);
    face_elem.push_back(nodes[2]);
    face_elem.push_back(0);//nodes[3]-1);
    face_name.push_back("noname");
//    if(nodes[3]==0)     face_name.push_back("noname");
//    else                face_name.push_back("fluid");
    face_i.push_back(count*2);
    count++;

  }
  fscanf(fp, "%s\n", line);
  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &nfaces_fluid, &dummy2, str3.c_str());
//  fscanf(fp, "%s\n", line);
  nfaces_fluid=nfaces_fluid-nfaces_bc;
//  count=0;

  for (int i = 0; i < nfaces_fluid; i++)
  {
    fscanf(fp, "%10d%x%x%x%x\n", &dummy1, &nodes[0], &nodes[1], &nodes[2], &nodes[3]);
    //    if (i%100 == 0) printf("%10d\t%lf\t%lf\t%lf\n", i, nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);

    //    for (int j=0; j<2; j++) face_v.push_back(nodes[j]-1);

    face_v.push_back(nodes[0]-1);
    face_v.push_back(nodes[1]-1);
    face_elem.push_back(nodes[2]);
    face_elem.push_back(nodes[3]);
    face_name.push_back("fluid");
//    if(nodes[3]==0)     face_name.push_back("noname");
//    else                face_name.push_back("fluid");
    //    face_i.push_back(face_i[face_i.size()-1]+2);
    face_i.push_back(count*2);
    count++;

  }
  printf("\n Fluent input file (*.msh format) read! \n");
//  cout<<"Number of vertices: "<< nvert<<endl;

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
  while (strcmp(line, "\"Dimensions:\")") != 0)    fscanf(fp, "%s\n", line);
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

    for(int nf=0; nf<nnodes; nf++)    face_v.push_back(nodes[nf]-1);
    face_elem.push_back(nodes[nnodes]);
    face_elem.push_back(0);
    face_i.push_back(count*nnodes);

//    face_v.push_back(nodes[0]-1);
//    face_v.push_back(nodes[1]-1);
//    if      (ndim==2)
//    {
//      face_elem.push_back(nodes[2]);
//      face_elem.push_back(0);//nodes[3]-1);
//      face_i.push_back(count*2);
//    }
//    else
//    {
//      face_v.push_back(nodes[2]-1);
//      face_v.push_back(nodes[3]-1);
//      face_elem.push_back(nodes[4]);
//      face_elem.push_back(0);//nodes[3]-1);
//      face_i.push_back(count*4);
//    }

    face_name.push_back("noname");
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

    for(int nf=0; nf<nnodes; nf++)    face_v.push_back(nodes[nf]-1);
    face_elem.push_back(nodes[nnodes]);
    face_elem.push_back(nodes[nnodes+1]);
    face_i.push_back(count*nnodes);

//    face_v.push_back(nodes[0]-1);
//    face_v.push_back(nodes[1]-1);
//    if      (ndim==2)
//    {
//      face_elem.push_back(nodes[2]);
//      face_elem.push_back(nodes[3]);//nodes[3]-1);
//      face_i.push_back(count*2);
//    }
//    else
//    {
//      face_v.push_back(nodes[2]-1);
//      face_v.push_back(nodes[3]-1);
//      face_elem.push_back(nodes[4]);
//      face_elem.push_back(nodes[5]);//nodes[3]-1);
//      face_i.push_back(count*4);
//    }

    face_name.push_back("fluid");
    count++;

  }
  printf("\n Fluent input file (*.msh format) was read! \n");
//  cout<<"Number of vertices: "<< nvert<<endl;

}

void UNSTRUCTMESH::readFluentMsh_named(const char *name, int nvert_gambit)  //NOT FINISHED!!
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
  int nvert, nfaces, nnodes;
  deque <int> nfaces_bc;
  string str1, str2, str3, str4, str5;
  char line[1000];
  fscanf(fp, "%s\n", line);

  // Reading the dimensiones of the mesh!! header
  while (strcmp(line, "\"Dimensions:\")") != 0)    fscanf(fp, "%s\n", line);
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
  cout<<"Starting to read faces: Number of total faces: "<< nfaces<< endl;

  int num_bc=0;
  int total_faces=0;
  int tmp_nfaces_bc;

//
//  fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str() ,str2.c_str(), &dummy2, &tmp_nfaces_bc, &dummy2, str3.c_str());
//  printf(" %s is what is said with the number of faces %i      \n", str1, tmp_nfaces_bc);
//
////  cout<<"Dummy: "<< dummy2<<" second dumm " <<tmp_nfaces_bc<<" third dumm " <<dummy3<<endl;
//  while (false)
//  {
//    cout<<"Does it enter even once?"<<endl;
//
////    fscanf(fp, "%s %s %x %x %i %s \n", str1.c_str(), str2.c_str(), &dummy1, &tmp_nfaces_bc, &dummy2, str3.c_str());
//    int nodes[6];
//    int count=0;
//    int bc_num_faces= tmp_nfaces_bc -total_faces;
//    cout<<"Read a section of boundary faces: Number of bc faces: "<< bc_num_faces<<" this is boundary number: " << num_bc<< endl;
//
//    nfaces_bc.push_back(bc_num_faces);
//
//    total_faces+=bc_num_faces;
//
//
//    for (int i = 0; i < bc_num_faces; i++)
//     {
//       if      (ndim==2)    fscanf(fp, "%10d%x%x%x%x\n", &nnodes, &nodes[0], &nodes[1], &nodes[2], &nodes[3]);
//       else if (ndim==3)    fscanf(fp, "%10d%x%x%x%x%x%x\n", &nnodes, &nodes[0], &nodes[1], &nodes[2], &nodes[3], &nodes[4], &nodes[5]);
//       else
//       {
//         printf("\n Dimension of the mesh is wrong!! \n");
//         throw(-1);
//       }
//
//       for(int nf=0; nf<nnodes; nf++)    face_v.push_back(nodes[nf]-1);
//       face_elem.push_back(nodes[nnodes]);
//       face_elem.push_back(nodes[nnodes+1]);
//       face_i.push_back(count*nnodes);
//
////       face_name.push_back("noname");
//       count++;
//
//     }
////     fscanf(fp, "%s\n", line);
////     int tmp_nfaces_bc;
//    fscanf(fp, "%s %i %s %x %x %i %s \n", str1.c_str(), &dummy1 ,str2.c_str(), &dummy2, &tmp_nfaces_bc, &dummy2, str3.c_str());
//    num_bc++;
//  }
//
//  while (strcmp(line, "\"Zones:\")") != 0)    fscanf(fp, "%s\n", line);
//  for(int i = 0; i < nfaces_bc.size()-1; i++)
//  {
//    fscanf(fp, "%s %s %s %s %s \n", str1.c_str(), str2.c_str(), str3.c_str(), str4.c_str(), str5.c_str());
//    for(int j=0 ; j<nfaces_bc[i]; j++)
//      face_name.push_back(str4);
//
//  }
//
//  fscanf(fp, "%s %i %s %x %x %i %s \n", str1.c_str(), &dummy1 ,str2.c_str(), &dummy2, &tmp_nfaces_bc, &dummy2, str3.c_str());
//  for(int j=0 ; j<nfaces_bc[nfaces_bc.size()-1]; j++)
//    face_name.push_back("fluid");

  printf("\n Fluent input file (*.msh format) was read! \n");
//  cout<<"Number of vertices: "<< nvert<<endl;

}

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
//  if (ndfcd != 2)  // Check if it is indeed a 2D grid
//  {
//    printf(" The input file should be 3D but ndfcd = %i \n", ndfcd);
//    printf(" Exiting! \n");
//    exit(-1);
//  }
  printf(" %dD Grid --> Vertices: %i        Elements: %i \n", ndfcd, nvertices, nelements);

  // Allocate memory to save vertices and elements
  nodes.resize(nvertices);
  elems.resize(nelements);


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
//  if (ndfcd != 2)  // Check if it is indeed a 2D grid
//  {
//    printf(" The input file should be 3D but ndfcd = %i \n", ndfcd);
//    printf(" Exiting! \n");
//    exit(-1);
//  }
  printf(" %dD Grid --> Vertices: %i        Elements: %i \n", ndfcd, nvertices, nelements);

  // Allocate memory to save vertices and elements
  nodes.resize(nvertices);
  elems.resize(nelements);


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
  elems.resize(nelements);


//  while (strcmp(line, GAMBIT_VER) != 0)    fscanf(fp, "%s\n", line);

  int dummy;
  char line[1000];

  // Read nodal coordinates
  for (int i = 0; i < nvertices; i++)
  {
    fscanf(fp, "%d%d%lf%lf\n", &dummy, &dummy, &nodes[i].pt.x, &nodes[i].pt.y);
    fscanf(fp, "%s\n", line);
    fscanf(fp, "%7s%lf\n", line, &nodes[i].pt.z);
//    printf("%lf\t%lf\t%lf\n", nodes[i].pt.x, nodes[i].pt.y, nodes[i].pt.z);
  }

//  throw(-112);


  // Read elements
  elem_i.push_back(0);

  for (int c = 0; c < nelements; c++)
  {
    int nodes[10];

    fscanf(fp, "%18s%d%d%d%d%d%d%d%d", line, &dummy, &dummy, &nodes[0], &nodes[1], &nodes[2], &nodes[3],&nodes[4], &nodes[5]);
    fscanf(fp, "%s\n", line);
    fscanf(fp, "%7s%d%d\n", line, &nodes[6], &nodes[7]);

//    for (int ne=0; ne<8; ne++) printf("%d\t", nodes[ne]);
//    printf("\n");

    for (int ne=0; ne<8; ne++) elem_v.push_back(nodes[ne]-1);
    elem_i.push_back(elem_i[elem_i.size()-1]+8);
  }
  fclose(fp);

  printf("\n input file in NASTRAN format read! \n");
}



void UNSTRUCTMESH::drawMesh()
{
  for (int e=0; e<elem_i.size()-1; e++)
    {
//        glLineWidth(0.2);
    glBegin(GL_LINE_LOOP);
    {
      double grey = 0.5;
      glColor3d(0.0, 0.0, 0.0);


      for (int ne=elem_i[e]; ne<elem_i[e+1]-1; ne++)
      {
        glVertex3d(nodes[elem_v[ne]].pt.x,   nodes[elem_v[ne]].pt.y,   nodes[elem_v[ne]].pt.z);
        glVertex3d(nodes[elem_v[ne+1]].pt.x, nodes[elem_v[ne+1]].pt.y, nodes[elem_v[ne+1]].pt.z);
      }
    }
    glEnd();
//      glLineWidth(1.0);
  }

}

void UNSTRUCTMESH::drawSolid()
{
  for (int e=0; e<elem_i.size()-1; e++)
    {
    glBegin(GL_TRIANGLES);
    glColor3d(0.85, 0.85, 0.85);

    POINT t1 = nodes[elem_v[elem_i[e]]].pt;
    for (int ne=elem_i[e]; ne<elem_i[e+1]-1; ne++)
    {
      POINT t2 = nodes[elem_v[ne]].pt;
      POINT t3 = nodes[elem_v[ne+1]].pt;

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
  out.pointlist = (double *) NULL;
  out.pointattributelist = (double *) NULL;
  out.pointmarkerlist = (int *) NULL;
  out.trianglelist = (int *) NULL;
  out.triangleattributelist = (double *) NULL;
  out.neighborlist = (int *) NULL;
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int *) NULL;
  out.edgelist = (int *) NULL;
  out.edgemarkerlist = (int *) NULL;

  vorout.pointlist = (double *) NULL;
  vorout.pointattributelist = (double *) NULL;
  vorout.edgelist = (int *) NULL;
  vorout.normlist = (double *) NULL;


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

  findNodeNeighbors();
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
