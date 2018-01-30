#include "phParAdapt.h"

#ifdef __cplusplus
extern "C" {
#endif

// build the linear system
void
buildSystem(pRegion region, double* eltMatrix)
{

  pPList regionVerts;
  regionVerts = R_vertices(region, 1);
  double X1[3];
  double X2[3];
  double X3[3];
  double X4[3];
  
  V_coord((pVertex)PList_item(regionVerts,0), X1);
  V_coord((pVertex) PList_item(regionVerts,1), X2);
  V_coord((pVertex) PList_item(regionVerts,2), X3);
  V_coord((pVertex) PList_item(regionVerts,3), X4);
  
  PList_delete(regionVerts);
  
  ///////////////////////////////////////////////////////////////////////
  // matrix used for reconstruction of the polynomial: only geometric  //
  // facts are needed                                                 //
  // The right hand side is given by the field values                  //
  ///////////////////////////////////////////////////////////////////////
  
  // coordinates
  double x1 = X1[0];
  double x2 = X2[0];
  double x3 = X3[0];
  double x4 = X4[0];
  double y1 = X1[1]; 
  double y2 = X2[1];
  double y3 = X3[1];
  double y4 = X4[1];
  double z1 = X1[2];
  double z2 = X2[2];
  double z3 = X3[2];
  double z4 = X4[2];
    
  /////////////////////////////////////////////////////////////////////////
  // we are starting with the zeroth order polynomial                    //
  // that is, m(1,1) = 1 , the coefficient in front of the constant term //
  /////////////////////////////////////////////////////////////////////////
  // first row: field at node 1 
  eltMatrix[0] = 1;  // const

  eltMatrix[1] = x1; // x

  eltMatrix[2] = y1; // y

  eltMatrix[3] = z1; // z
  // second row: field at node 2 
  eltMatrix[4] = 1;  // const

  eltMatrix[5] = x2; // x

  eltMatrix[6] = y2; // y

  eltMatrix[7] = z2; // z

  // third row: field at node 3 
  eltMatrix[8] = 1;   // const

  eltMatrix[9] = x3;  // x

  eltMatrix[10] = y3; // y

  eltMatrix[11] = z3; // z

  // forth row: field at node 4 
  eltMatrix[12] = 1;  // const

  eltMatrix[13] = x4; // x

  eltMatrix[14] = y4; // y

  eltMatrix[15] = z4; // z

}

void
buildSystemXYZ(dArray *xyz, double* eltMatrix) {
  int counter = 0;

  for(int iVert=0; iVert<4; iVert++) {
    eltMatrix[counter++] = 1.;
    for(int iComp=0; iComp<3; iComp++) {
      eltMatrix[counter++] = xyz[iVert][iComp];
    }
  }

}
   
#ifdef __cplusplus
}
#endif
