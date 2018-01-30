#ifndef __attachData_H_
#define __attachData_H_
#ifdef SIM
#include "MeshSim.h"
#endif
#ifdef FMDB
#include "AOMD.h"
#include "AOMDInternals.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
  
  // attaches array to mesh entities
  // `dataID' is the MeshDataId
  // `nVar' is the no. of variables at each dof
  // e.g., `nVar'=5 (for flow problems) or 27 (for hessians)
  // `poly' is the polynomial order 
  // this routine attaches  valueArray
  // the incoming "valueArray" which contains data
  // for ALL vertices (general dofs) is split into 
  // local entity-level arrays to handle the memory
  // during local mesh modifications 
  void attachArray( double *valueArray,
		    pMesh mesh,
		    pMeshDataId dataID,
		    int nVar,
		    int poly);


  // get data (previously attached) from mesh
  // `dataID' is the MeshDataId
  // `nVar' is the no. of variables at each dof
  // e.g., `nVar'=5 (for flow problems) or 27 (for hessians)
  // `poly' is the polynomial order 
  // this routine gets attached data  array from mesh   
  // in restart-writable format 
  // memory is allocated within the function
  //  user has to delete the memory
  void getAttachedArray( double *&valueArray,
                         pMesh mesh,
                         pMeshDataId dataID,
                         int nVar,
                         int poly);

  // cleans data attached to mesh entities
  // `dataID' is the MeshDataId
  // `en_type' is the entity type on which data is atached
  // e.g., 0 (for vertex), 3 (for regions), 4 (for all)
  // can use polynomial order to delete data on entities
  // i.e., for different modes, instead user should call routine
  // for different entities depending on poly. order/mode
  // with this attached data doesn't have to be solution
  void cleanAttachedData( pMesh mesh,
			  pMeshDataId dataID,
			  int en_type,
			  int array=1);
  
#ifdef __cplusplus
}
#endif


#endif
