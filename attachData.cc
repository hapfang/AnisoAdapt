#include "attachData.h"
#include "phastaIO.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

// attaches array to mesh entities
// `dataID' is the MeshDataId
// `nVar' is the no. of variables at each dof
// e.g., `nVar'=5 (for flow problems) or 27 (for hessians)
// `poly' is the polynomial order 
// this routine attaches  "valueArray"
// the incoming "valueArray" which contains data
// for ALL vertices (general dofs) is split into 
// local entity-level arrays to handle the memory
// during local mesh modifications 
void attachArray( double *valueArray, 
		  pMesh mesh, 
		  pMeshDataId dataID,
		  int nVar, 
		  int poly ) {
  int count,i;
  int nem = (poly > 1) ? (poly - 1) : 0;
  int nfm = (poly > 3) ? ((poly-3)*(poly-2)/2) : 0;
  int nrm = (poly > 5) ? ((poly-4)*(poly-5)*(poly-3)/6) : 0;
  if(poly==3) nfm =1;

  /* attach the vertex coefficients */
  count = 0;
  pVertex vertex;
  VIter vIter = M_vertexIter( mesh );
  while (vertex = VIter_next( vIter ) ) {
    double *data = new double[nVar];
    int vtxID;
    vtxID = EN_id(vertex);

    for (i=0; i < nVar; i++) data[i] = valueArray[count++];
    EN_attachDataPtr( (pEntity)vertex, dataID, (void *) data );
  }
  VIter_delete( vIter );

  /* attach the edge coefficients */
  if (nem > 0){
    pEdge edge;
    EIter eIter = M_edgeIter( mesh );
    while (edge = EIter_next( eIter ) ) {
      double* data = new double[nVar*nem];
      for (i=0; i < nVar*nem; i++) data[i] = valueArray[count++];
      EN_attachDataPtr( (pEntity)edge, dataID, (void *) data );
    }
    EIter_delete( eIter );
  }
  
  /* attach face coefficients */
  if (nfm > 0){
    pFace face;
    FIter fIter = M_faceIter( mesh );
    while (face = FIter_next( fIter ) ) {
      double* data = new double[nVar*nfm];
      for ( i=0; i < nVar*nfm; i++ ) data[i] = valueArray[count++];
      EN_attachDataPtr( (pEntity)face, dataID, (void *) data );
    }
    FIter_delete( fIter );
  }

  /* attach region coefficients */
  if (nrm > 0) {
     cerr << " No code to attach Region Modes " << endl;
     exit ( -1 );
  }

}
// get data (previously attached) from mesh
// `dataID' is the MeshDataId
// `nVar' is the no. of variables at each dof
// e.g., `nVar'=5 (for flow problems) or 27 (for hessians)
// `poly' is the polynomial order (ONLY order 1 is supported as of now) 
// this routine gets attached data array from mesh   
// in restart-writable format 
// memory is allocated within the function
// user has to delete the memory
void getAttachedArray( double *&valueArray,
                       pMesh mesh,
                       pMeshDataId dataID,
                       int nVar,
                       int poly)
{


  int i;
  if(poly!=1) {
      cerr << "\nError in getAttachedData() [in attachData.cc]" << endl;
      cerr << "Polynomial order [" << poly << "] NOT supported" << endl;
      exit(-1);
  }

  int nshg = M_numVertices(mesh);
  valueArray = new double[nshg*nVar];

  /* attach the vertex coefficients */
  int vCount = 0;
  pVertex vertex;
  VIter vIter = M_vertexIter( mesh );
  while (vertex = VIter_next( vIter ) ) {
    double *data;
    if(!EN_getDataPtr( (pEntity)vertex, dataID, (void **)&data )) {
        cerr << "\nError in getAttachedData() [in attachData.cc]" << endl;
        cerr << "Data not attached to vertex [" << EN_id((pEntity)vertex) << "]" << endl;
        exit(0);
    }
    for (i=0; i < nVar; i++) valueArray[vCount+i*nshg] = data[i];
    vCount++;

  }
  VIter_delete( vIter );
}



// cleans data attached to mesh entities
// `dataID' is the MeshDataId
// `en_type' is the entity type on which data is atached
// e.g., 0 (for vertex), 3 (for regions), 4 (for all)
// can use polynomial order to delete data on entities
// i.e., for different modes, instead user should call routine
// for different entities depending on poly. order/mode
// with this attached data doesn't have to be solution
void cleanAttachedData(pMesh mesh,
		       pMeshDataId dataID,
		       int en_type,
		       int array) {
  switch(en_type) {
  case 0:
    {      
      pVertex vtx;
      VIter vt_iter=M_vertexIter(mesh);
      while(vtx = VIter_next(vt_iter)) {
//	double *data;
//	if(EN_getDataPtr((pEntity)vtx,dataID,(void**)&data)) {
	if(EN_getDataPtr((pEntity)vtx,dataID,NULL)) {

//	  if(array)
//	    delete [] data;
//	  else
//	    free(data);

	  EN_deleteData((pEntity)vtx,dataID);
	}
	else {
	  printf("\nError in cleanAttachedData() : Data not attached to vertex\n");
	  exit(1);
	}
      }
      VIter_delete(vt_iter);
    }
    break;
  case 1:
    {
      pEdge edge;
      EIter eg_iter=M_edgeIter(mesh);
      while(edge = EIter_next(eg_iter)) {
	double *data;
	if(EN_getDataPtr((pEntity)edge,dataID,(void**)&data)) {

	  if(array)
	    delete [] data;
	  else
	    delete data;

	  EN_deleteData((pEntity)edge,dataID);
	}
	else {
	  printf("\nError in cleanAttachedData() : Data not attached to edge\n");
	  exit(1);
	}
      }
      EIter_delete(eg_iter);
    }
    break;
  case 2:
    {
      pFace face;
      FIter face_iter=M_faceIter(mesh);
      while(face = FIter_next(face_iter)) {
	double *data;
	if(EN_getDataPtr((pEntity)face,dataID,(void**)&data)) {

	  if(array)
	    delete [] data;
	  else
	    delete data;

	  EN_deleteData((pEntity)face,dataID);
	}
	else {
	  printf("\nError in cleanAttachedData() : Data not attached to face\n");
	  exit(1);
	}
      }
      FIter_delete(face_iter);
    }
    break;
  case 3:
    {
      pRegion rg;
      RIter rg_iter=M_regionIter(mesh);
      while(rg = RIter_next(rg_iter)) {
	double *data;
	if(EN_getDataPtr((pEntity)rg,dataID,(void**)&data)) {

	  if(array)
	    delete [] data;
	  else
	    delete data;

	  EN_deleteData((pEntity)rg,dataID);
	}
	else {
	  printf("\nError in cleanAttachedData() : Data not attached to region\n");
	  exit(1);
	}
      }
      RIter_delete(rg_iter);
    }
    break;
  case 4://all
    {
      for(int iENType=0; iENType<4; iENType++)
	cleanAttachedData(mesh,dataID,iENType,array);
    }
    break;
  default :
    printf("Check the en_type in cleanAttachedData(...)\n");
    break;
  }
}
