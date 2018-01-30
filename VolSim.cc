#include <stdlib.h>
#include <stdio.h>
#ifdef SIM
#include "SimPartitionedMesh.h" 
#include "MeshSim.h"
#include "MeshSimInternal.h"
#include "Eigen.h"
typedef class MVertex * pVertex;
typedef class MEdge * pEdge;
typedef class MFace * pFace;
typedef class MRegion * pRegion;
#include <vector>


double XYZ_volumeSim (double xyz[4][3])
{
   double v01[3],v02[3],v03[3];
   double normal[3];
   diffVt(xyz[1],xyz[0],v01);
   diffVt(xyz[2],xyz[0],v02);
   crossProd(v01,v02,normal);
   diffVt(xyz[3],xyz[0],v03);
   
   return(dotProd(normal,v03)/6);
}

//to calculate the volumes of the pyramids
double R_VolumeSim(pRegion region) {
   int topoType = R_topoType(region);
   int numChildTets;
   int elemType;
   // fit two tets. in a pyramid
   // (diagonal edge of a pyramid is : V0-V2)
   // fit three tets. in a prism/wedge
   // (diagonal edges of a prism are : V0-V4, V1-V5 and V5-V0)
   //  int tetVerts[3][4] = {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}};
   int mapNodes[3][3][4] ={{{0,1,2,3}, {-1,-1,-1,-1}, {-1,-1,-1,-1}},
      {{0,1,2,4}, {0,2,3,4}, {-1,-1,-1,-1}},
      {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}}};
   
   switch(topoType) {
      case Rtet :
         numChildTets = 1;
         elemType = 0;
         break;
      case Rpyramid :
         numChildTets = 2;
         elemType = 1;
         break;
      case Rwedge :
         numChildTets = 3;
         elemType = 2;
         break;
      //  case HEX :
      //     numChildTets = 6; 
     //     elemType = 3; 
     //     break; 
      default:
        printf("\nError in R_volume2()...");
        printf("region topology NOT supported");
        exit(0);
      }

  dArray *X_Y_Z;
  pPList rverts = R_vertices(region,1);
  int numRVerts = PList_size(rverts);
  X_Y_Z = new dArray[numRVerts];
  pVertex vtx;
  for(int iRVert=0; iRVert<numRVerts; iRVert++) {
       vtx = (pVertex)PList_item(rverts,iRVert);
       V_coord(vtx,X_Y_Z[iRVert]);
  }
  PList_delete(rverts);

  double volume = 0.;
  dArray sub_xyz[4];
  for(int iTet=0; iTet<numChildTets; iTet++) {
       for(int iTVert=0; iTVert<4; iTVert++)
              for(int iComp=0; iComp<3; iComp++)
                  sub_xyz[iTVert][iComp] = X_Y_Z[mapNodes[elemType][iTet][iTVert]][iComp];

         volume += XYZ_volumeSim(sub_xyz);
         }

  delete [] X_Y_Z;

  return volume;
}

#endif //endif SIM

