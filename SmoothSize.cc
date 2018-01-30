#include "phParAdapt.h"
//#include "MA.h"
#include "Eigen.h"
#include <stdlib.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <stdio.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalSizeID;
extern pMeshDataId numSurroundNodesID;
extern pMeshDataId nodalDirectionID;
  
// simple average over a patch surrounding the vertex    
void 
SmoothSize(pMesh mesh, int num) {

  pVertex v;
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;
  int* numSurroundingVerts = new int[M_numVertices(mesh)];
  // keep the vals in memory before finally setting them
  double** averageSize;
  averageSize = new double*[M_numVertices(mesh)];  
  
#if  ( defined  DEBUG )
//   printf("[%i] in SmoothHessians:\n",PMU_rank());
#endif 

  while(v = VIter_next(vIter)) {
    averageSize[vCount] = new double[num];
    
    for(int i=0;i<num;i++) {
      averageSize[vCount][i]=0.0;
    }
    // can be different from V_numEdges
    // will include the actual vertex, too
    numSurroundingVerts[vCount]=0;
    
    // exclude non-interor vertex from contributing to
    // its own average
    // also make sure this vertex is owned by this proc
//    if(V_whatInType(v)== 3 && EN_isOwnerProc((pEntity)v)) {
    if(EN_isOwnerProc((pEntity)v)) {
      // first the "central" vertex
      double *nodalSize;
      if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&nodalSize)){
	
	cout<<"\nerror in SmoothSize: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      for(int i=0; i<num; i++) {
	averageSize[vCount][i]= nodalSize[i];
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
//            if(V_whatInType(otherVertex)== 3){//is in interior
            
                double *nodalSizeOther;
                if(!EN_getDataPtr((pEntity)otherVertex, nodalSizeID,(void**)&nodalSizeOther)){
                    
                    cout<<"\nerror in SmoothSize: no data attached to OTHER vertex\n";
                    V_info(otherVertex);
                    exit(0);
                }
                // add values up
                for(int i=0;i<num;i++) {
                    averageSize[vCount][i] = averageSize[vCount][i] + nodalSizeOther[i];
                }
                numSurroundingVerts[vCount]++;
                
#if  ( defined  DEBUG )
//                 printf("\n[%i]is also vert (orig):\n",PMU_rank());
//                 double cf[3], cv[3];
//                 V_coord(otherVertex,cf);
//                 V_coord(v,cv);
//                 printf("coords %f %f %f (orig) %f %f %f\n",cf[0],cf[1],cf[2],cv[0],cv[1],cv[2]);
                
#endif  
//            }//if(V_whatInType(otherVertex)== 3
        }// if(EN_isOwnerProc(edge){
            
    }//for( numEdges

    // At this moment we cannnot decide whether a vertex' neighbors are
    // exclusively on the bdry --> also look at adjacent partitions
    // for now, we rely that the initial mesh hase been prepared and got 
    // all its bdry-only elements removed

    // if no surronding vertex found 
    // this is not a good point to take care of this for partitioned mesh
    
    if(!numSurroundingVerts[vCount]) {
      double *nodalSize;
      if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&nodalSize)){
	cout<<"\nerror in SmoothSize: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      
      // just add its own values up
      // mostly the case when a region has no intr. node/vert.
      for(int i=0;i<num;i++) {
       averageSize[vCount][i] = nodalSize[i];
      }
      numSurroundingVerts[vCount]++;
    }
/*
    //divide by number of Surr vertices
  for(int i=0;i<3;i++) {
     averageSize[vCount][i] = averageSize[vCount][i]/(numSurroundingVerts[vCount]);
  }
*/

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
    double* oldSize;
    if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&oldSize)){
      
      cout<<"\nerror in SmoothSize: no data attached to OTHER vertex\n";
      V_info(v);
      exit(0);
    }
//      delete [] oldHessian;
    
//      EN_deleteData((pEntity)v, nodalHessianID);
   
    
    double* smoothSize = new double[num];

    for(int k=0;k<num;k++){
        smoothSize[k] = averageSize[vCount][k];
    }

//    delete [] oldSize;
//    EN_modifyData( (pEntity)v, nodalSizeID);
    EN_modifyDataPtr( (pEntity)v, nodalSizeID, (void *)
		       smoothSize);

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
      delete [] averageSize[i];
  }
  delete [] averageSize;
  delete [] numSurroundingVerts;

//  int numExtraSmoothPassesBdryFaces = 1;
//  int numExtraSmoothPassesBdryEdges = 1;
//  int numExtraSmoothPassesBdryVerts = 1;

//  SmoothSizeOnBdry(mesh,numExtraSmoothPassesBdryFaces,2);
//  SmoothSizeOnBdry(mesh,numExtraSmoothPassesBdryEdges,1);
//  SmoothSizeOnBdry(mesh,numExtraSmoothPassesBdryVerts,0);

}

void SmoothSizeOnBdry(pMesh mesh, int numPasses, int genDim) {

  pMeshDataId avgSizeID = MD_newMeshDataId("avg. at bdry.");
//  for(int iPass=0; iPass<numPasses; iPass++) {
    pVertex v;
    VIter vit = M_vertexIter(mesh);
    while(v = VIter_next(vit)) {
      // vertex on domain bdry. of dimension "genDim"
      if(V_whatInType(v)==genDim  && EN_isOwnerProc((pEntity)v))  {

   double *avgSize = new double[3];
	int numSurroundingVertsBdry = 1;

      if( EN_isOwnerProc((pEntity)v)) { //add own value only if owned
	for(int i=0;i<3;i++)
	  avgSize[i] = 0.;

	// start from 1 due to its own contri.
	double *nodalSize;
	if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&nodalSize)) {
	  cout<<"\nerror in SmoothSizeOnBdry: no data attached to bdry. vertex\n";
	  V_info(v);
	  exit(0);
	}

	for(int i=0;i<3;i++)
	  avgSize[i] += nodalSize[i];
      } //OwnweProc

   int numEdges = V_numEdges(v);
	pEdge edge;
	for(int i=0; i<numEdges; i++) {
	  edge = V_edge(v,i);
     if(EN_isOwnerProc((pEntity)edge)){
   	  pVertex otherVertex;
	     otherVertex = E_otherVertex(edge,v);
      
	     // do not account for vertices
	     // classified on model ents. with 
	     // dimension lower than "genDim"
	    if(V_whatInType(otherVertex)<genDim)
       continue;

	    double *nodalSizeOther;
       if(!EN_getDataPtr((pEntity)otherVertex, nodalSizeID,(void**)&nodalSizeOther)){
	     cout<<"\nerror in SmoothSizeOnBdry: no data attached to OTHER vertex of bdry. vertex\n";
	     V_info(otherVertex);
	     exit(0);
	    } 

	    // add values up
	    for(int i=0;i<3;i++)
	      avgSize[i] += nodalSizeOther[i];
	    numSurroundingVertsBdry++;
     } //if owned edge 
	}//for( ..i<numEdges...)

/*
	for(int i=0;i<3;i++)
	  avgSize[i] = avgSize[i]/numSurroundingVertsBdry;
*/
	EN_attachDataPtr((pEntity)v,avgSizeID,(void *)avgSize);
      }
    }
    
    VIter_reset(vit);
    while(v = VIter_next(vit)) {
      if(V_whatInType(v)==genDim) {
	double *nodalSize, *avgSize;

	if(!EN_getDataPtr((pEntity)v, nodalSizeID,
            (void**)&nodalSize)) {
	  cout<<"\nerror in SmoothSizeOnBdry: no data attached to bdry. vertex\n";
	  V_info(v);
	  exit(0);
	}
	
	if(!EN_getDataPtr((pEntity)v, avgSizeID,(void**)&avgSize)) {
	  cout<<"\nerror in SmoothSizeOnBdry: no data attached to bdry. vertex with avghessianID\n";
	  V_info(v);
	  exit(0);
	}

//	EN_modifyData((pEntity)v, nodalSizeID);
	// attach the new data ptr. (i.e., averaged hessian)
	EN_modifyDataPtr((pEntity)v, nodalSizeID, (void *)avgSize);
	// delete the old data ptr. (i.e., old nodal hessian)
//	delete [] nodalSize;
	
	// clean the ptr. attached with avghessianID for other passes 
	EN_deleteData((pEntity)v, avgSizeID);
      }
    }
    VIter_delete(vit);
//  }

  MD_deleteMeshDataId(avgSizeID);
}

  
void 
SmoothDir(pMesh mesh) {

  pVertex v;
  int vCount=0, Indx=0;
  double eOwn[3][3];
  double eOther[3][3];
  double eNorm[3][3];
  int* numSurroundingVerts = new int[M_numVertices(mesh)];
  // keep the vals in memory before finally setting them
  double tol = 1e-3;
  double dprod;
  double** averageDir;
  averageDir = new double*[M_numVertices(mesh)];  
  int Count=0; 
  VIter vIter=M_vertexIter(mesh);
  while(v = VIter_next(vIter)) {
    averageDir[vCount] = new double[9];
    
    for(int iDir=0;iDir<9;iDir++) {
      averageDir[vCount][iDir]=0.0;
    }
    // can be different from V_numEdges
    // will include the actual vertex, too
    numSurroundingVerts[vCount]=0;
    
    //  make sure this vertex is owned by this proc
    if(EN_isOwnerProc((pEntity)v) ) {
    // first the "central" vertex
      double *nodalDir;
      if(!EN_getDataPtr((pEntity)v, nodalDirectionID,(void**)&nodalDir)){
	
	   cout<<"\nerror in SmoothDirs: no data attached to  vertex\n";
	   V_info(v);
	   exit(0);
      }
      
      //build the unit vectors (only for owned) (called e) 
      //stored as row vectors
      for(int iRow=0; iRow<3; iRow++) {
         for(int iCol=0; iCol<3; iCol++) {
            eOwn[iRow][iCol]=nodalDir[iRow*3+iCol];
         }
      }
      numSurroundingVerts[vCount]++;

    }//if(V_whatInType(v)== 3 && ..
    else { //not owned so need to zero out eOwn
      for(int iRow=0; iRow<3; iRow++) {
         for(int iCol=0; iCol<3; iCol++) {
            eOwn[iRow][iCol]=0.0;
         }
      } 
    }

    int numEdges = V_numEdges(v);
    pEdge edge;

    // now for the surrounding vertices
    // now this vertex v also can be non-owned
    // but EDGES to surrounding verts MUST be owned
    for(int i=0; i<numEdges ; i++){
        edge = V_edge(v,i);

        if(EN_isOwnerProc((pEntity)edge)){
            pVertex otherVertex;
            otherVertex = E_otherVertex(edge,v);
            
            double *nodalDirOther;
            if(!EN_getDataPtr((pEntity)otherVertex, nodalDirectionID,(void**)&nodalDirOther)){
                    
            cout<<"\nerror in SmoothDirs: no data attached to OTHER vertex\n";
            V_info(otherVertex);
            exit(0);
            }
      
            for(int iRow=0; iRow<3; iRow++) {
               for(int iCol=0; iCol<3; iCol++) {
                  eOther[iRow][iCol]=nodalDirOther[iRow*3+iCol];
               }
            }
            numSurroundingVerts[vCount]++;

            //compute dot product and and the vector averages
            for(int iRow=0; iRow<3; iRow++) {
               dprod = dotProd(eOwn[iRow], eOther[iRow]);
               if(dprod < 0.0) {
                  for(int iCol=0; iCol<3; iCol++) {
                     eOwn[iRow][iCol]=eOwn[iRow][iCol]-eOther[iRow][iCol];
                  }
               }
               else {
                  for(int iCol=0; iCol<3; iCol++) {
                     eOwn[iRow][iCol]=eOwn[iRow][iCol]+eOther[iRow][iCol];
                  }
               }    
            }   
                
        }// if(EN_isOwnerProc(edge){
            
    }//for( numEdges

    // if no surronding vertex found 
    // this is not a good point to take care of this for partitioned mesh
    if(!numSurroundingVerts[vCount]) {
      double *nodalDir;
      if(!EN_getDataPtr((pEntity)v, nodalDirectionID,(void**)&nodalDir)){
	cout<<"\nerror in SmoothDirs: no data attached to  vertex\n";
	V_info(v);
	exit(0);
      }
      // just add its own values up
      // mostly the case when a region has no intr. node/vert.
      Indx=0;
      for(int iRow=0; iRow<3; iRow++) {
         for(int iCol=0; iCol<3; iCol++) {
            eOwn[iRow][iCol]=nodalDir[Indx];
            Indx++;
         }
      }
      numSurroundingVerts[vCount]++;
    } //no surr vertex


   for(int iRow=0; iRow<3; iRow++) {
       normVt(eOwn[iRow], eNorm[iRow]);
       for(int iCol=0; iCol<3; iCol++) { 
         averageDir[vCount][iRow*3+iCol]=eOwn[iRow][iCol];
 //        printf("%lf  ",averageDir[vCount][iRow*3+iCol]);
       }
    }
/*    printf("\n");
      //check orthogonality of the resultant vectors            
      int factor;
      if(!checkUnitaryOthoganal(eNorm, factor)) {
//        printf("error in check orthogonality\n");
      Count++;
      }
*/
    vCount++;
  }// while v = vIter
  VIter_reset(vIter); 
//  printf("number of non orthogonal vectors: %i\n",Count); 
  vCount=0;
  while(v = VIter_next(vIter)) {
    // delete the old Dir data
    double* oldDir;
    if(!EN_getDataPtr((pEntity)v, nodalDirectionID,(void**)&oldDir)){
      
      cout<<"\nerror in SmoothDirs: no data attached to OTHER vertex\n";
      V_info(v);
      exit(0);
    }
//      delete [] oldDir;
    
//      EN_deleteData((pEntity)v, nodalDirID);
    
    double* smoothDir = new double[9];

    for(int k=0;k<9;k++){
        smoothDir[k] = averageDir[vCount][k];
    }

    delete [] oldDir;
    EN_deleteData( (pEntity)v, nodalDirectionID);
    EN_attachDataPtr( (pEntity)v, nodalDirectionID, (void *)
		       smoothDir);

    // 
    double* nSurrNodes = new double[1];
    nSurrNodes[0] = 1.0*numSurroundingVerts[vCount] ; 
    EN_attachDataPtr( (pEntity)v, numSurroundNodesID, (void *)
                       nSurrNodes);

    vCount++;
  }
  VIter_delete(vIter);

  delete [] averageDir;
  delete [] numSurroundingVerts;

}

#ifdef __cplusplus
}
#endif
