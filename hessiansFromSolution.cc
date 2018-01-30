#include "phParAdapt.h"
#include "mpi.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalHessianID;

// hessian  returned : 6-component (symmetric)
// u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
// the nodal data later can be retrieved via
// nodalHessianID
void
hessiansFromSolution(pParMesh pmesh, pMesh mesh,int stepNumber)
{  
  // compute the hessian field from the solution
  
  // recover gradients using a patch
  // attaches:
  //  Sum of weighted (local) gradients to vertices via localGradientID
  //  Sum of (local) patch Volumes to vertices via localPatchVolId
  
//  for (int k=0; k<10; k++){
//    SmoothErrorIndicators(mesh, 1);
//  }

  gradientsFromPatch(mesh);
  MPI_Barrier(MPI_COMM_WORLD);


  // communicates missing contributions over procs/=parts
  // each part p's  vertex has
  // (double localPatchVolume, double localPatchGradient) , i.e 
  // V_p = Sum_numSurrEles(VolE) 
  // G_p = Sum_numSurrEles(VolE * grad(E))
  //
  // data are gathered via
  // localGradientID
  // localPatchVolId
  // The commu attaches the combined data  via nodalGradientID
  commuGradientsFromPatch(pmesh, mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // recover hessians using a patch (from gradients)
  // attaches local hessian to vertices via localHessianId
  // (local) patch Volumes is Still there via localPatchVolId 
  // hessian  attached : 6-component (symmetric)
  // u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
  hessiansFromPatch(mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // communicates missing contributions over procs/=parts
  // analogous to commuGradientsFromPatch
  commuHessiansFromPatch(pmesh,mesh);  
  MPI_Barrier(MPI_COMM_WORLD);

  // simple average over a patch surrounding the vertex
  // needs:
  // (local) average Hessian
  // (local) number of contributing vertices
  SmoothHessians(mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // communicates missing contributions over procs/=parts
  // attaches the averaged Hessian via nodalHessianId
  commuSmoothHessians(pmesh,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // writeRestartHessians( mesh ) ; 

}

#ifdef __cplusplus
}
#endif
