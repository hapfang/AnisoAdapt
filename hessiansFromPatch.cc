#include "phParAdapt.h"
#include <stdio.h>
#ifdef SIM
#include "SimMeshTools.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif


extern pMeshDataId nodalGradientID;
extern pMeshDataId localHessianID;
      
// recover hessian using a patch (from gradients)
// hessian computed : 6-component (symmetric)
// u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
// attaches a nodal Hessian  value to ea/ vertex via nodalHessianID
void
hessiansFromPatch(pMesh mesh)
{
  VIter vIter;
  int nvtx = M_numVertices(mesh);
  vIter = M_vertexIter(mesh);
  pVertex vertex;
  
  double elemHessian[6];
  
  // loop over vertices and get patch of elements
  while(vertex = VIter_next(vIter)) {
    pPList elementPatch;
    elementPatch = V_regions(vertex);


    double* patchHessian = new double[6];
    for(int i=0; i<6; i++) {
      patchHessian[i]=0.0;
    }
    
#if  ( defined  DEBUG ) 
//     printf("\n=============================\n");
//     printf("[%i]:in hessiansFromPatch: vertex:",PMU_rank());
//     double c[3];
//     V_coord(vertex,c);
//     printf("coords %f %f %f\n",c[0],c[1],c[2]);
//     printf("\n=============================\n");
#endif    


    // loop over elements and recover hessian at vertex
    for (int i=0; i < PList_size(elementPatch); i++) {

      elementHessian( (pRegion)PList_item(elementPatch,i), elemHessian);
      
#if  ( defined  DEBUG )
//       printf("[%i]:\n",PMU_rank());
//       printf("\nELEM hess :\n %f\n %f\n %f\n",elemHessian[0],elemHessian[1],elemHessian[2]);
//       printf("\n  %f\n %f\n %f\n",elemHessian[3],elemHessian[4],elemHessian[5]);
//       pPList regionVerts;
//       regionVerts = R_vertices( (pRegion)PList_item(elementPatch,i),0);
//       double co[3];
//       printf("\nand region:\n");
//       for(int k=0;k<4;k++){
//           V_coord((pVertex)PList_item(regionVerts,k),co);
//           printf(" %f %f %f\n",co[0],co[1],co[2]);
//       }
//       PList_delete(regionVerts);   
#endif      
#ifdef SIM
      double vol =  R_VolumeSim((pRegion) PList_item(elementPatch,i));
#endif      
#ifdef FMDB 
      double vol =  R_Volume2((pRegion) PList_item(elementPatch,i));
#endif      
      // use element hessian (each component) to recover hessian at vertex
      for(int j=0; j<6; j++) {
	patchHessian[j]= patchHessian[j] + 
	  elemHessian[j]*vol;
      }
    }//loop over element patch ends
    
    // attach local contribs
    // attach the recovered hessian to the vertex    
    EN_attachDataPtr( (pEntity)vertex, localHessianID, (void *)
		      patchHessian );

    PList_delete(elementPatch);
  }//while (vertex = VIter_next(vIter))      
  VIter_delete(vIter);


}

#ifdef __cplusplus
}
#endif
