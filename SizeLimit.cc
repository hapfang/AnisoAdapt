#include "phParAdapt.h"
#ifdef FMDB
#include "BLUtil.h"
#endif
#ifdef SIM
#include "SimAdvMeshing.h"
#endif
#include <stdlib.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <stdio.h>

using namespace std;

extern int MaxLimitFact;
extern int MinLimitFact;
extern pMeshDataId nodalSizeID;
extern pMeshDataId nodalDirectionID;

#ifdef __cplusplus
extern "C" {
#endif
  
void ReadOrgSize(double* size) {
   printf("Lets migrate the sizes!\n");
   char sizefile[256];
   sprintf(sizefile,"OrgSize%d.dat",PMU_rank());
   FILE *sf;
   sf=fopen(sizefile, "rb");
   int nshg2;
   int count=0;
   fscanf(sf,"%i",&nshg2);
   printf("%i\n",nshg2);
   while(1) {
      fscanf(sf,"%lf",&size[count]);
      count++;
      if(feof(sf)) break;
   }
   //        attachArray(size, pmesh, migrate_size,3, poly);
   //        delete[] value;
   fclose(sf);
}

//The idea with the following routine is to copy sizes from upwind direction
//on to the sizes of a model edge since it has been seen that the model edges
//particularly on the outflow get really bad size requests
void SizePropogate(pMesh mesh) {
  
  pVertex vert;
  VIter vIter = M_vertexIter(mesh);

  while(vert = VIter_next(vIter)) {
    
    pPList BaseVtxList;
    BaseVtxList = PList_new();

    if(EN_isBLEntity(vert)) { 

      pPList VGrowth;
      pVertex BaseVert;

      VGrowth = V_growthCurveVertices(vert);
      int Size = PList_size(VGrowth);
     
      BaseVert = (pVertex)PList_item(VGrowth, 0);
      int FoundVtx = PList_inList(BaseVtxList, BaseVert);
      if(!FoundVtx) {
        PList_append(BaseVtxList, BaseVert);

        pVertex LastVert; 
        LastVert = (pVertex)PList_item(VGrowth, Size-1);

        double* h;
        double* hLast;
        double* hCopy = new double[3];

        if(!EN_getDataPtr((pEntity)BaseVert, nodalSizeID, (void **)&h)) {
           printf("Error in size propogate\n");
           exit(0);
        }

        if(!EN_getDataPtr((pEntity)LastVert, nodalSizeID, (void **)&hLast)) {
           printf("Error in size propogate\n");
           exit(0);
        }

        EN_deleteData((pEntity)LastVert, nodalSizeID);

        for(int i=0; i<3; i++) {
            hCopy[i] = h[i];
        }

        EN_attachDataPtr((pEntity)LastVert, nodalSizeID, (void *)hCopy);
      }
     
      PList_delete(VGrowth);
    } //BLEntity
    PList_delete(BaseVtxList);
  }//vIter 
   VIter_delete(vIter);
/*   
   pGModel model = M_model(mesh);
   GEIter geIter = GM_edgeIter(model);
   pGEdge gedge;
   int NodeCount=0; int NumEdge=0;
   while(gedge = GEIter_next(geIter))  {

//      printf("Number of vertices on edge %d\n", M_numClassifiedVertices(mesh, (pGEntity)gedge));

       VIter vIter  = M_classifiedVertexIter(mesh,(pGEntity)gedge, 0);
       pVertex vert;       
       while (vert = VIter_next(vIter))  {

          double *h = new double[3];
          double *hCopy = new double[3];
          hCopy[0] = 0.0; hCopy[1] = 0.0; hCopy[2] = 0.0;
          if(!EN_getDataPtr((pEntity)vert, nodalSizeID,(void**)&h)){
             cout<<"\nerror in SizePropogate: no data attached to  vertex\n";
             V_info(vert);
             exit(0);
          }

          int numEdges = V_numEdges(vert);
          pEdge edge;
          int vCount=0;
          for(int i=0; i<numEdges ; i++){
             edge = V_edge(vert,i);
             pVertex otherVertex;
             otherVertex = E_otherVertex(edge,vert);

             if(V_whatInType(otherVertex)!=2) 
                continue;
             
             if(E_typeInBL(edge)==eLAYER) {
//                cout<<"\n I get here in eLayer\n";
                double *hOther;
                if(!EN_getDataPtr((pEntity)otherVertex, nodalSizeID,(void**)&hOther)) {
                      cout<<"\nerror in SizePropogate: no data attached to OTHER vertex\n";
                      V_info(otherVertex);
                      exit(0);
                }
                for(int i=0; i<3; i++) {
                   hCopy[i]=hCopy[i]+hOther[i];
                }
                vCount++;
             }
           }

           for(int i=0; i<3; i++) {
               hCopy[i]=hCopy[i]/vCount;
           }
           
           vector<pVertex> GC;
           pVertex v;
           BL_getGrowthCurveNodes(vert, GC);
           int GCSize = GC.size();
           for(int i=0; i<GC.size(); i++) {
              v=GC[i];
              EN_deleteData( (pEntity)v, nodalSizeID);
              EN_attachDataPtr((pEntity)v, nodalSizeID, (void*)
                    hCopy);
           }

           delete[] h;
           EN_deleteData( (pEntity)vert, nodalSizeID);
           EN_attachDataPtr( (pEntity)vert, nodalSizeID, (void *)
                                     hCopy);
           NodeCount++;
       }
       VIter_delete(vIter);
       NumEdge++;
    }
    GEIter_delete(geIter);
//    printf("total edges: %d, total nodes: %d\n", NumEdge, NodeCount);
*/
}

void SizeLimit(pMesh mesh){

  pMeshDataId limitSizeID = MD_newMeshDataId("avg. at bdry.");
  pVertex v;

  char sizefile[256];
  int nshg;
  
  sprintf(sizefile,"OrgSize%i.dat",PMU_rank());
  FILE *sf;
  sf=fopen(sizefile, "r");
  fscanf(sf,"%d",&nshg);
  printf("%d\n",nshg);
  
  printf("Limiting the size field request by factor of original size\n");
  VIter vIter=M_vertexIter(mesh);
  int vCount=0;
  int count=0;
//  double value[nshg*3]; 
  double *value;
  value = new double[nshg*3]; 

  while(1){
     fscanf(sf,"%lf",&value[count]);
//     printf("%i %lf\n",count, value[count]);
     count++;
     if(feof(sf)) break;
  }
  
//  printf("%lf %lf %lf\n", value[0], value[1],value[2]);
  fclose(sf);
  
  while(v = VIter_next(vIter)) {
//       if(EN_isOwnerProc((pEntity)v)){
     
     double *h;
     double anisoSize[3];
     if(!EN_getDataPtr((pEntity)v, nodalSizeID,(void**)&h)){

      cout<<"\nerror in SmoothSize: no data attached to  vertex\n";
      V_info(v);
      exit(0);
     }
   
//     MA_Vtx_Size (v, anisoSize);

//put solution from value array in anisoSize 3x3 matrix for each vertex and then
//attach it to an     
     int vtxID;
     vtxID = EN_id(v);

    double* h2 = new double[3];
    double* hCopy = new double[3];

    for(int i=0; i<3; i++) {
       h2[i]=0.0;
       hCopy[i]=h[i];
    }
     
    for(int i=0; i<3; i++){
      h2[i]=value[vCount*3+i];
    }
     
    for (int i=0;i<3;i++){
       if(hCopy[i]<h2[i]/MinLimitFact)
          hCopy[i]=h2[i]/MinLimitFact;
       if(hCopy[i]>h2[i]*MaxLimitFact)
          hCopy[i]=h2[i]*MaxLimitFact;
    }
       
//           printf("%f\n",h[i]);

//   EN_deleteData( (pEntity)v, nodalSizeID);
//   EN_attachDataPtr( (pEntity)v, nodalSizeID, (void *)
//                  hCopy);
   EN_modifyDataPtr( (pEntity)v, nodalSizeID, (void *)
                  hCopy);
   EN_attachDataPtr( (pEntity)v, limitSizeID, (void *)
                         h2);
//    } // EN_isOwner
    vCount++;   
  }//vertex loop
  VIter_delete(vIter);
#ifdef DEBUG  
//  M_writeVTKFile(mesh, "OrgSize", limitSizeID, 3);
#endif
  MD_deleteMeshDataId(limitSizeID);

  delete[] value;
}

void WriteSizeField(pMesh mesh) {

  char wsizefile[256];
  sprintf(wsizefile,"SizeField%d.dat",PMU_rank());
  FILE *wsf;
  wsf=fopen(wsizefile, "wb");
  
  pVertex vertex;
  VIter vIter = M_vertexIter(mesh);
  while(vertex = VIter_next(vIter)) {
    
    double* h;
    double* dir;
    EN_getDataPtr((pEntity)vertex, nodalSizeID, (void**)&h);
    EN_getDataPtr((pEntity)vertex, nodalDirectionID, (void**)&dir);

    double dir2[3][3]; 

    for(int i=0; i<3; i++) {
       for(int j=0; j<3; j++) {
          dir2[i][j]=dir[i*3+j]*h[i];
          fprintf(wsf, "%lf ",dir2[i][j]);
       }
    }
    fprintf(wsf, "\n"); 
  }
  VIter_delete(vIter);
  fclose(wsf);

}

void ReadSizeField(double* SizeField) {
   printf("Lets migrate the sizes!\n");
   char rsizefile[256];
   sprintf(rsizefile,"SizeField%d.dat",PMU_rank());
   FILE *rsf;
   rsf=fopen(rsizefile, "rb");
   int count=0;
   while(1) {
      fscanf(rsf,"%lf",&SizeField[count]);
      count++;
      if(feof(rsf)) break;
   }
   //        attachArray(size, pmesh, migrate_size,3, poly);
   //        delete[] value;
   fclose(rsf);
}

#ifdef __cplusplus
}
#endif

