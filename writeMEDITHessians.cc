#include "phParAdapt.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <strings.h>
#include "Eigen.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId errorIndicatorID;

void
writeMEDITSolution(pMesh mesh)
{
  ofstream fout;
  fout.open("sizefield-sol.bb");
  
  fout<<"3 1 "<<M_numVertices(mesh)<<" 2\n";

  double *nodalData;
  pVertex vertex;
  VIter vIter = M_vertexIter( mesh );
  while (vertex = VIter_next( vIter ) ) {
    if(!EN_getDataPtr((pEntity)vertex, errorIndicatorID,
		     (void**)&nodalData)) {
      cout<<"\nerror in writeMEDITSolution : no data attached to vertex\n";
      exit(0);
    }
    fout<<nodalData[4]<<endl;
  }
  VIter_delete(vIter);
  fout.close();
}

////////////////////////////////////////////
// write in MEDIT format                  //
// for visualization of mesh-metric field //
// hess is the meshSim formatted sizefield
////////////////////////////////////////////
void
writeMEDITSizeField(Hessian* hess, pMesh mesh,int currentTimeStep, int gprocID)
{
    pVertex vertex;
    VIter vIter = M_vertexIter( mesh );

    int counter =0;
    
    ofstream fout;
    char sizefield_file[50];
    sprintf(sizefield_file,"sizefield.%i.%i.bb",currentTimeStep,gprocID+1);
    fout.open(sizefield_file);

    // file header 
    fout<<"3 6 "<<M_numVertices(mesh)<<" 2\n";

    while ( vertex = VIter_next( vIter ) ) {

      int factor=1;
      double n[3];
      crossProd(hess[counter].dir[0],hess[counter].dir[1],n);

      if( dotProd(hess[counter].dir[2],n)<0. )
	factor=-1;

      double h[3];
      h[0]=hess[counter].h[0];
      h[1]=hess[counter].h[1];
      h[2]=hess[counter].h[2];
      
      double xx  = hess[counter].dir[0][0]*hess[counter].dir[1][0];
      xx=xx/(h[0]*h[1]);
      double xy  = hess[counter].dir[0][0]*hess[counter].dir[1][1];
      xy=xy/(h[0]*h[1]);
      double xz  = hess[counter].dir[0][0]*hess[counter].dir[1][2];
      xz=xz/(h[0]*h[1]);
      double yx  = hess[counter].dir[0][1]*hess[counter].dir[1][0];
      yx=yx/(h[0]*h[1]);
      double zx  = hess[counter].dir[0][2]*hess[counter].dir[1][0];
      zx=zx/(h[0]*h[1]);
      double oyx = hess[counter].dir[1][1]*factor*hess[counter].dir[2][0];
      oyx=oyx/(h[1]*h[2]);

      double xxy = xx*factor*hess[counter].dir[2][1];
      xxy=xxy/(h[2]);
      double xxz = xx*factor*hess[counter].dir[2][2];
      xxz=xxz/(h[2]);
      double xyx = xy*factor*hess[counter].dir[2][0];
      xyx=xyx/(h[2]);
      double xyy = xy*factor*hess[counter].dir[2][1];
      xyy=xyy/(h[2]);
      double xyz = xy*factor*hess[counter].dir[2][2];
      xyz=xyz/(h[2]);
      double xzx = xz*factor*hess[counter].dir[2][0];
      xzx=xzx/(h[2]);
      double xzy = xz*factor*hess[counter].dir[2][1];
      xzy=xzy/(h[2]);

      double yxx = yx*factor*hess[counter].dir[2][0];
      yxx=yxx/(h[2]);
      double yxy = yx*factor*hess[counter].dir[2][1];  
      yxy=yxy/(h[2]);
      double yxz = yx*factor*hess[counter].dir[2][2];
      yxz=yxz/(h[2]);
      double yzx = hess[counter].dir[0][1]*hess[counter].dir[1][2]*factor*hess[counter].dir[2][0];
      yzx=yzx/(h[0]*h[1]*h[2]);
      double yyx = hess[counter].dir[0][1]*oyx;
      yyx=yyx/(h[0]);

      double zxx = zx*factor*hess[counter].dir[2][0];
      zxx=zxx/(h[2]);
      double zxy = zx*factor*hess[counter].dir[2][1];
      zxy=zxy/(h[2]);
      double zyx = hess[counter].dir[0][2]*oyx;
      zyx=zyx/(h[0]);

      double vp[3];
      for(int i=0;i<3;i++) {
	vp[i]=1./(hess[counter].h[i]*hess[counter].h[i]);
      }
      
      double metric[6];

      metric[0]=  vp[0]*(xyz-xzy)+vp[1]*(zxy-yxz)+vp[2]*(yzx-zyx);
      
      metric[1]= -vp[0]*(xxz-xzx)-vp[1]*(zxx-xxz)-vp[2]*(xzx-zxx);

      metric[2]=  vp[0]*(xxy-xyx)+vp[1]*(yxx-xxy)+vp[2]*(xyx-yxx);

      metric[3]= -vp[0]*(yxz-yzx)-vp[1]*(zyx-xyz)-vp[2]*(xzy-zxy);

      metric[4]=  vp[0]*(yxy-yyx)+vp[1]*(yyx-xyy)+vp[2]*(xyy-yxy);

      metric[5]=  vp[0]*(zxy-zyx)+vp[1]*(yzx-xzy)+vp[2]*(xyz-yxz);

      fout<<metric[0]<<" ";
      fout<<metric[1]<<" "<<metric[3]<<" ";
      fout<<metric[2]<<" "<<metric[4]<<" "<<metric[5]<<"\n";
      
      counter++;
    }
    VIter_delete(vIter); 

    fout.close();
}

#ifdef __cplusplus
}
#endif
