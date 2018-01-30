#include "phParAdapt.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalHessianID;
extern pMeshDataId numSurroundNodesID;
  
// simple average over a patch surrounding the vertex    
void 
SmoothHessians(pMesh mesh) {

  pVertex v;
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;
  int* numSurroundingVerts = new int[M_numVertices(mesh)];
  // keep the vals in memory before finally setting them
  double** averageHessian;
  averageHessian = new double*[M_numVertices(mesh)];  
  
#if  ( defined  DEBUG )
//   printf("[%i] in SmoothHessians:\n",PMU_rank());
#endif 

  while(v = VIter_next(vIter)) {
    averageHessian[vCount] = new double[6];
    
    for(int iHess=0;iHess<6;iHess++) {
      averageHessian[vCount][iHess]=0.0;
    }
    // can be different from V_numEdges
    // will include the actual vertex, too
    numSurroundingVerts[vCount]=0;
    
    // exclude non-interor vertex from contributing to
    // its own average
    // also make sure this vertex is owned by this proc
    if(V_whatInType(v)== 3 && EN_isOwnerProc((pEntity)v) ) {
      // first the "central" vertex
      double *nodalHessian;
      if(!EN_getDataPtr((pEntity)v, nodalHessianID,(void**)&nodalHessian)){
	
	cout<<"\nerror in SmoothHessians: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      for(int iHess=0; iHess<6; iHess++) {
	averageHessian[vCount][iHess]= averageHessian[vCount][iHess] + nodalHessian[iHess];
      }
      numSurroundingVerts[vCount]++;

#if  ( defined  DEBUG )
//       double c[3];
//       V_coord(v,c);
//       printf("\n[%i]contributing to vert: coords %f %f %f\n",PMU_rank(),c[0],c[1],c[2]);
#endif     


    }//if(V_whatInType(v)== 3 && ..

    int numEdges = V_numEdges(v);
    pEdge edge;

    // now for the surrounding vertices
    // now this vertex v also can be non-owned
    // but EDGES to surrounding verts MUST be owned and verts also be off model-bdry
    for(int i=0; i<numEdges ; i++){
        edge = V_edge(v,i);

        if(EN_isOwnerProc((pEntity)edge)){

            pVertex otherVertex;
            otherVertex = E_otherVertex(edge,v);
            
        
#if  ( defined  DEBUG )
//             double cfq[3], cvq[3];
//             V_coord(otherVertex,cfq);
//             V_coord(v,cvq);
//             printf("\n[%i]considering other vertex (ORIG) (ownflag %i) %f %f %f:\n",PMU_rank(),EN_isOwnerProc((pEntity)v),cvq[0],cvq[1],cvq[2]);
//             printf("coords %f %f %f (orig) %f %f %f\n",cfq[0],cfq[1],cfq[2],cvq[0],cvq[1],cvq[2]);
#endif   

            // if neighbor vert is NOT on bdry NOR is it NOT owned (neighbor  HAS to be owned)
            if(V_whatInType(otherVertex)== 3){//is in interior
            
                double *nodalHessianOther;
                if(!EN_getDataPtr((pEntity)otherVertex, nodalHessianID,(void**)&nodalHessianOther)){
                    
                    cout<<"\nerror in SmoothHessians: no data attached to OTHER vertex\n";
                    V_info(otherVertex);
                    exit(0);
                }
                // add values up
                for(int iHess=0;iHess<6;iHess++) {
                    averageHessian[vCount][iHess] = averageHessian[vCount][iHess] + nodalHessianOther[iHess];
                }
                numSurroundingVerts[vCount]++;
                
#if  ( defined  DEBUG )
//                 printf("\n[%i]is also vert (orig):\n",PMU_rank());
//                 double cf[3], cv[3];
//                 V_coord(otherVertex,cf);
//                 V_coord(v,cv);
//                 printf("coords %f %f %f (orig) %f %f %f\n",cf[0],cf[1],cf[2],cv[0],cv[1],cv[2]);
                
#endif  
            }//if(V_whatInType(otherVertex)== 3
        }// if(EN_isOwnerProc(edge){
            
    }//for( numEdges

    // At this moment we cannnot decide whether a vertex' neighbors are
    // exclusively on the bdry --> also look at adjacent partitions
    // for now, we rely that the initial mesh hase been prepared and got 
    // all its bdry-only elements removed

    // if no surronding vertex found 
    // this is not a good point to take care of this for partitioned mesh
    if(!numSurroundingVerts[vCount]) {
      double *nodalHessian;
      if(!EN_getDataPtr((pEntity)v, nodalHessianID,(void**)&nodalHessian)){
	cout<<"\nerror in SmoothHessians: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      
      // just add its own values up
      // mostly the case when a region has no intr. node/vert.
      for(int iHess=0;iHess<6;iHess++) {
       averageHessian[vCount][iHess] = averageHessian[vCount][iHess] + nodalHessian[iHess];
      }
      numSurroundingVerts[vCount]++;
    }

//      if(!numSurroundingVerts[vCount]){
//        cout<<"\nerror in SmoothHessians: there is a boundary vertex whose\n"
//  	  <<"neighbors are exclusively classsified on modelfaces/edges/vertices\n"
//  	  <<"and NOT in the interior\n";
      
//        cout<<"For the following vertex : "<<endl;
//        V_info(v);
//        exit(0);
//      }
    vCount++;
  }// while v = vIter
  VIter_reset(vIter); 
  
  vCount=0;
  while(v = VIter_next(vIter)) {
    // delete the old Hessian data
    double* oldHessian;
    if(!EN_getDataPtr((pEntity)v, nodalHessianID,(void**)&oldHessian)){
      
      cout<<"\nerror in SmoothHessians: no data attached to OTHER vertex\n";
      V_info(v);
      exit(0);
    }
//      delete [] oldHessian;
    
//      EN_deleteData((pEntity)v, nodalHessianID);
   
    
    double* smoothHess = new double[6];

    for(int k=0;k<6;k++){
        smoothHess[k] = averageHessian[vCount][k];
    }
    delete [] oldHessian;
    EN_deleteData( (pEntity)v, nodalHessianID);
    EN_attachDataPtr( (pEntity)v, nodalHessianID, (void *)
		       smoothHess);

    // 
    double* nSurrNodes = new double[1];
    nSurrNodes[0] = 1.0*numSurroundingVerts[vCount] ; 
    EN_attachDataPtr( (pEntity)v, numSurroundNodesID, (void *)
                       nSurrNodes);
    
#if  ( defined  DEBUG )
//             printf("\n[%i]attaching in  SmoothHessians:\n %f %f %f\n",PMU_rank(),smoothHess[0],smoothHess[1],smoothHess[2]);
//             printf(" %f %f %f\n",smoothHess[3],smoothHess[4],smoothHess[5]);
//             printf("\nand local numSurrVertex Num %f\n",nSurrNodes[0]);
//             printf("for vertex:");
//             double c[3];
//             V_coord(v,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif     

    vCount++;
  }
  VIter_delete(vIter);

  for(int i=0;i<M_numVertices(mesh);i++){
      delete [] averageHessian[i];
  }
  delete [] averageHessian;
  delete [] numSurroundingVerts;

}
  
#ifdef __cplusplus
}
#endif
