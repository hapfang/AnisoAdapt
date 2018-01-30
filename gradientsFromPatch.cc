#ifdef SIM
#include "SimMeshTools.h"
#endif
#include "phParAdapt.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


extern pMeshDataId localPatchVolID;
extern pMeshDataId localGradientID;
// recover gradient at a vertex using patch of elements around it
//  attaches:
//  Sum of weighted (local) gradients to vertices via localGradientID
//  Sum of (local) patch Volumes to vertices via localPatchVolID
// this is a LOCAL operation
void
gradientsFromPatch(pMesh mesh)
{
  VIter vIter;
  vIter = M_vertexIter(mesh);
  pVertex vertex;
  pRegion region;

  double elemGradient[3];
  int count = 0;
  
  // loop over vertices and get patch of elements
  while(vertex = VIter_next(vIter)) {
      pPList elementPatch;
      elementPatch = V_regions(vertex);
      
      
      double* patchVolume = new double[1];
      patchVolume[0] = 0;
      
      double* patchGradient = new double[3];
      for(int i=0; i<3; i++){
          patchGradient[i]= 0.0  ;
      }

#if  ( defined  DEBUG ) 
//       printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
//       printf("[%i]:in gradientsFromPatch: vertex:",PMU_rank());
//       double c[3];
//       V_coord(vertex,c);
//       printf("coords %f %f %f\n",c[0],c[1],c[2]);
//       printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif      

      // loop over elements and recover gradient at vertex   
      for (int i=0; i<PList_size(elementPatch); i++) {

	region = (pRegion) PList_item(elementPatch,i);
#ifdef SIM   
	double vol = R_VolumeSim(region);
#endif   
#ifdef FMDB   
	double vol = R_Volume2(region);
#endif   
	*patchVolume+=vol;

        // retrieve element grad. for each element in the patch at vertex
        elementGradient(region, elemGradient);

#if  ( defined  DEBUG )
//           printf("\n[%i]: ELEM grad :\n %f\n %f\n %f\n",PMU_rank(),elemGradient[0],elemGradient[1],elemGradient[2]);
//           pPList regionVerts;
//           regionVerts = R_vertices( (pRegion)PList_item(elementPatch,i),0);
//           double co[3];
//           printf("\nand region for vertex %i:\n",count);
//           for(int k=0;k<4;k++){
//               V_coord((pVertex)PList_item(regionVerts,k),co);
//               printf(" %f %f %f\n",co[0],co[1],co[2]);
//           }
//           PList_delete(regionVerts);   
#endif      


          // use element grad. (each component) to recover gradient at vertex
          for(int j=0; j<3; j++) {
              patchGradient[j]= patchGradient[j] + 
                  elemGradient[j]*vol;
          }
      }//loop over element patch ends
      
      
#if  ( defined  DEBUG )
//       printf("\n[%i]:in gradientsFromPatch  local vals:\n %f\n %f\n %f\n",PMU_rank(),patchGradient[0],patchGradient[1],patchGradient[2]);
//       double ci[3];
//       printf("\nand local vol %f:\n",patchVolume[0]);
//       V_coord(vertex,ci);
//       printf(" for vertex %i   %f %f %f\n",count,ci[0],ci[1],ci[2]);
#endif       

      // attach the recovered gradient to the vertex
      EN_attachDataPtr( (pEntity)vertex, localGradientID, (void *)
                        patchGradient );
      
      EN_attachDataPtr( (pEntity)vertex, localPatchVolID, (void *)
                        patchVolume);
      
      PList_delete(elementPatch);
      count++;
  }
  VIter_delete(vIter);

}

#ifdef __cplusplus
}
#endif
