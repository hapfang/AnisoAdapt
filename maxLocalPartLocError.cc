#include "phParAdapt.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif


extern pMeshDataId nodalHessianID;
extern pMeshDataId locMaxInterpolErrorID;

// local max interpolaton error local to each partition
// determined via the attached nodalHessian
// and all the edges of the vertex that are LOCAL to 
// this partition. 
// The max local interpolation error (LOCAL) on this partition
// is the max over ALL local edges 
//
// for each vertex, the overall max has to be communicated by 
// an additional 'max-operation' where off-proc/part vertices are also included
// (done in  commuMaxLocalPartLocError.cc)
void
maxLocalPartLocError(pMesh mesh)
{

    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;
    double H[3][3];

    while(vertex = VIter_next(vIter)) {
        double* nodalHessian;
        
        if(!EN_getDataPtr((pEntity)vertex, 
			  nodalHessianID,(void**)&nodalHessian)){
            cout<<"\nerror in :V_Hessian no data attached to vertex\n";
            exit(0);
        }
        int k = 0;
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                if(j>=i){
                    H[i][j] = nodalHessian[k];
                    H[j][i] = nodalHessian[k];
                    k++;
                }
            }
        }
                
        double* maxLocErr = new double[1];
        
        // takes edges of vertex that are LOCAL to this partition
        maxLocErr[0] =  maxLocalError( vertex,H);

        // attach that
        EN_attachDataPtr( (pEntity)vertex, locMaxInterpolErrorID, (void *)
                           maxLocErr);


#if  ( defined  DEBUG )
//             printf("\nattaching  LOCAL  Max loc error: %f \n",maxLocErr[0]);
//             printf("for vertex:");
//             double c[3];
//             V_coord(vertex,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif 

    }
    VIter_delete(vIter);
}    
#ifdef __cplusplus
}
#endif
