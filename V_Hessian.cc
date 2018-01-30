#include "phParAdapt.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalHessianID;

// hessian returned : 6-component (symmetric)
// u_xx, u_xy, u_xz, u_yy, u_yz, u_zz
// called in setSizeFieldUsingHessians (sizefield.cc)
void
V_Hessian(pVertex v, double T[3][3])
{
  double* nodalHessian;
  
  if(EN_getDataPtr((pEntity)v, nodalHessianID,(void**)&nodalHessian) == NULL){
    cout<<"\nerror in :V_Hessian no data attached to vertex\n";
    exit(0);
  }
  
  T[0][0]=nodalHessian[0];
  T[0][1]=T[1][0]=nodalHessian[1];
  T[0][2]=T[2][0]=nodalHessian[2];
  T[1][1]=nodalHessian[3];
  T[1][2]=T[2][1]=nodalHessian[4];
  T[2][2]=nodalHessian[5];  
}
    
#ifdef __cplusplus
}
#endif
