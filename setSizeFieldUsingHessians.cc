#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef SIM
#include "SimMeshTools.h"
#endif
#include "MeshSimAdapt.h"
#include "SimAdvMeshing.h"

#include "phParAdapt.h"
#include "Eigen.h"
#include "attachData.h"
#include "phReadWrite.h"
#include "mpi.h"


#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalHessianID;
extern pMeshDataId nodalSizeID;
extern pMeshDataId nodalDirectionID;

extern pMeshDataId localGradientID; 
extern pMeshDataId localPatchVolID;
extern pMeshDataId localHessianID;
extern pMeshDataId numSurroundNodesID;
extern pMeshDataId locMaxInterpolErrorID;
extern pMeshDataId globMaxInterpolErrorID;

extern int nErrorVars;
extern int timeStepNumber;
extern int isThickAdapt;
extern int isSizeLimit;

extern pMeshDataId meshSizeID;

extern int preLBforAdaptivity;

void setSizeFieldUsingHessians(pParMesh pmesh,
                               pMesh mesh,
			       pMSAdapt simAdapter,
			       double factor,
			       double hmax,
			       double hmin,
			       int option)
{

  meshSizeID = MD_newMeshDataId("mesh size ID");


  int nshg=M_numVertices(mesh);  // only true for linear elements
  double tol=1.e-12;
  
  double T[3][3];
  double maxEigenval=0;
  int i,j;
  
  pVertex vertex;
  double eloc;  	  // local error at a vertex
  double etot=0.;	  // total error for all vertices
  double etotLoc=0.;	  // total error for all vertices local to partition
  double emean; 	  // emean = etot / nv
  double elocmax=0.;	  // max local error
  double elocmin=1.e20;   // min local error
  pProgress prog;

  char adaptSimLogFile[128];
  sprintf(adaptSimLogFile,"phAdapt.%i.log",PMU_gid(PMU_rank(),0)+1);
  // ofstream adaptSimLog(adaptSimLogFile);
  FILE *adaptSimLog=fopen(adaptSimLogFile,"w");

  // to get an idea refer Appendix A in Li's thesis
  // attach the max local interpolation error LOCAL
  // to each partition via locMaxInterpolErrorID;
  maxLocalPartLocError(mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // attach the max local interpolation error GLOBAL
  // to the parallel mesh via globMaxInterpolErrorID;
  commuMaxLocalPartLocError(pmesh,mesh);
  MPI_Barrier(MPI_COMM_WORLD);

  // struct Hessian contains decomposed values
  // mesh sizes and directional information
  Hessian *hess = new Hessian[nshg];

  VIter vit=M_vertexIter(mesh);
  i=0;
  while(vertex=VIter_next(vit)) { 
    double xyz[3];
    V_coord(vertex,xyz);

    // get hessian of appropriate variable at a `vertex'
    if(option==1 || option==2 || option==9) 
      V_Hessian(vertex,T);
    else{

//      V_AnalyticHessian(vertex,T,option);
        printf("in setSizeFieldUsingHessians: only option=1 works at the moment \n");
        exit(0);
    } 

    // `T' is the actual Hessian 
    // of any solution variable (e.g., T, u, v, w etc.)
    // `hess' has the information related to Mesh-Metric (eventually)
    // used for mapping the ellipsoid into a unit sphere 
    
    // compute eigen values and eigen vectors
    int nbEigen = eigen(T,hess[i].dir,hess[i].h);

    if(nbEigen == 1) {
      printf(" [%d] ",PMU_rank());
      printf("WARNING : nbEigen is 1\n");
      printf(" Vertex Number : %d\n",i);
      printf(" xyz : %f %f %f \n",xyz[0],xyz[1],xyz[2]);
      printf("       %f %f %f %f %f %f\n\n",T[0][0], 
             T[0][1],T[0][2],T[1][1],T[1][2],T[2][2]);
      hess[i].dir[0][0]=hess[i].dir[1][1]=hess[i].dir[2][2]=1.0;
      hess[i].dir[0][1]=hess[i].dir[0][2]=hess[i].dir[1][0]= 0.0;
      hess[i].dir[1][2]=hess[i].dir[2][0]=hess[i].dir[2][1]= 0.0;
    }

    hess[i].h[0] = ABS(hess[i].h[0]);
    hess[i].h[1] = ABS(hess[i].h[1]);
    hess[i].h[2] = ABS(hess[i].h[2]);
    
    if( hess[i].h[0] < tol ) {
      printf(" [%d] ",PMU_rank());
      printf("Error: zero maximum eigenvalue !!!\n");
      printf(" xyz : %f %f %f \n",xyz[0],xyz[1],xyz[2]);
      printf("       %f %f %f\n", hess[i].h[0],
             hess[i].h[1],hess[i].h[2]);
      printf("       %f %f %f %f %f %f\n\n",T[0][0], 
             T[0][1],T[0][2],T[1][1],T[1][2],T[2][2]);

//       adaptSimLog<<" ["<<PMU_rank()<<"]"<<endl;
//       adaptSimLog<<"Error: zero maximum eigenvalue !!!"<<endl;
//       adaptSimLog<<" xyz : "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n"<<endl;
//       adaptSimLog<<hess[i].h[0]<<" "<<hess[i].h[1]<<" "<<hess[i].h[2]<<endl;
//       adaptSimLog<<T[0][0]<<" "<<T[0][1]<<" "<<T[0][2]<<" ";
//       adaptSimLog<<T[1][1]<<" "<<T[1][2]<<" "<<T[2][2]<<endl;

// fix if nbeigen=1, for now assuming happens only at far outfield 
//  we want to set hmax for this 
    printf("For the above point, setting unit vectors as directions and hmax as sizes. If the point is not farfield, this is a serious error!\n\n");
      hess[i].dir[0][0]=hess[i].dir[1][1]=hess[i].dir[2][2]=1.0;
      hess[i].dir[0][1]=hess[i].dir[0][2]=hess[i].dir[1][0]= 0.0;
      hess[i].dir[1][2]=hess[i].dir[2][0]=hess[i].dir[2][1]= 0.0;
//      hess[i].h[0] = 1.0; hess[i].h[1]=1.0; hess[i].h[2]=1.0;
//      continue;
    }

    // estimate relative interpolation error
    // needed for scaling metric field (mesh size field)



    // retrieve the local interpolation error valid globally via 
    // globMaxInterpolErrorID
    double* globMaxInterpolError;
    if(!EN_getDataPtr((pEntity)vertex,globMaxInterpolErrorID ,
		      (void**)&globMaxInterpolError)){
      printf("\nerror in setSizeFieldUsingHEssians: no data attached to vertex\n");
                    exit(0);
    }
    eloc=globMaxInterpolError[0];


    // summing up contribs over verts
    // ONLY if vertex is OWNED
    // otherwise we would have double contributions
    if(EN_isOwnerProc((pEntity)vertex)){
        etotLoc += eloc;
    }
    if( eloc>elocmax )  elocmax=eloc;
    if( eloc<elocmin )  elocmin=eloc;

    ++i;
  }//while(vertex
  printf("Info [%d]: Reading hessian... done...\n\n",PMU_rank());

  MPI_Barrier(MPI_COMM_WORLD);
  // communicate the total
  MPI_Allreduce(&etotLoc, &etot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  // get the total number of verts, overal procs/parts
  int nshgTOT = getNSHGTOT(mesh);
  // cout<<"\n["<<PMU_rank()<<"] nshgTOT : "<<nshgTOT<<"\n"<<endl;
  
  emean =  etot / nshgTOT;

  if(PMU_rank()==0) {
    printf("\n Info on relative interpolation error: ");
    printf("   total: %f\n",etot);
    printf("   mean:  %f\n",emean);
    printf("   max local: %f\n",elocmax);
    printf("   min local: %f\n",elocmin);
  }

  eloc = emean*factor;

  if(PMU_rank()==0) {
    printf("\n towards uniform local error distribution of %f\n", eloc);
    printf("   with max edge length=%f; min edge length=%f\n\n",hmax,hmin);
  }

//   adaptSimLog<<"MeshSim   Library build ID : "<<SimMeshing_buildID()<<endl;
//   adaptSimLog<<"PMesh     Library build ID : "<<SimPMesh_buildID()<<endl;
//   adaptSimLog<<"MeshTools Library build ID : "<<SimMeshTools_buildID()<<endl;
//   adaptSimLog<<endl;
//  fprintf(adaptSimLog,"MeshSim   Library build ID : %s\n",SimMeshing_buildID());
//  fprintf(adaptSimLog,"PMesh     Library build ID : %s\n",SimPartitionedMesh_buildID());
//  fprintf(adaptSimLog,"MeshTools Library build ID : %s\n",SimMeshTools_buildID());
  fprintf(adaptSimLog,"Simmetrix Library build ID : %s\n",Sim_buildID());
  fprintf(adaptSimLog,"\n");

  if(option==1 || option==2) {
    // adaptSimLog<<"Strategy chosen is anisotropic adaptation, i.e., size-field driven"<<endl;
    fprintf(adaptSimLog,"Strategy chosen is anisotropic adaptation, i.e., size-field driven\n");
  }
  if(option==9) {
    // adaptSimLog<<"Strategy chosen is isotropic adaptation (based on local hmin of mesh metric tensor from hessians), i.e., size-field driven"<<endl;
    fprintf(adaptSimLog,"Strategy chosen is isotropic adaptation (based on local hmin of mesh metric tensor from hessians), i.e., size-field driven\n");
  }
  if(option==1 || option==2 || option==9) {
    // adaptSimLog<<"Using numerical/computed hessian... done..."<<endl;
    fprintf(adaptSimLog,"Using numerical/computed hessian... done...\n");

  }
  else {
    // adaptSimLog<<"Using analytic hessian... done..."<<endl;
    fprintf(adaptSimLog,"Using analytic hessian... done...\n");
  }

//   adaptSimLog<<"Info on relative interpolation error :"<<endl;
//   adaptSimLog<<"total : "<<etot<<endl;
//   adaptSimLog<<"mean : "<<emean<<endl;
//   adaptSimLog<<"factor : "<<factor<<endl;
//   adaptSimLog<<"min. local : "<<elocmin<<endl;
//   adaptSimLog<<"max. local : "<<elocmax<<endl;
//   adaptSimLog<<"towards uniform local error distribution of "<<eloc<<endl;
//   adaptSimLog<<"with min. edge length : "<<hmin<<endl;
//   adaptSimLog<<"with max. edge length : "<<hmax<<endl;

  fprintf(adaptSimLog,"Info on relative interpolation error :\n");
  fprintf(adaptSimLog,"total : %f\n",etot);
  fprintf(adaptSimLog,"mean : %f\n",emean);
  fprintf(adaptSimLog,"factor : %f\n",factor);
  fprintf(adaptSimLog,"min. local : %f\n",elocmin);
  fprintf(adaptSimLog,"max. local : %f\n",elocmax);
  fprintf(adaptSimLog,"towards uniform local error distribution of %f\n",eloc);
  fprintf(adaptSimLog,"with min. edge length : %f\n",hmin);
  fprintf(adaptSimLog,"with max. edge length : %f\n",hmax);

  // scale the hessian
  // double errorSqRoot = hmin*sqrt(maxEigenval);
  VIter_reset(vit);
  i=0;

  int foundHmin = 0;
  int foundHmax = 0;
  int hminCount = 0;
  int hmaxCount = 0;
  int bothHminHmaxCount = 0;
  double tol2 = 0.01*hmax, tol3 = 0.01*hmin;

#ifdef FMDB  
  if(isThickAdapt)
     DistOffWall(mesh, simAdapter);
#endif  

  while ( vertex=VIter_next(vit)) {

    foundHmin = 0;
    foundHmax = 0;

    for( j=0; j<3; j++ ) {
      if( hess[i].h[j] < tol )
        hess[i].h[j] = hmax;
      else {
        hess[i].h[j] = sqrt(eloc/hess[i].h[j]);
        if( hess[i].h[j] > hmax )
          hess[i].h[j] = hmax;
        if( hess[i].h[j] < hmin )
          hess[i].h[j] = hmin;
      }
    }

    double localHmin = hmax;
    double* h2 = new double[3]; 
    double* dir2 = new double[9];

    for(j=0; j<3; j++) {
      // for option 9
      if(localHmin > hess[i].h[j])
	localHmin = hess[i].h[j];

      if(ABS(hess[i].h[j]-hmax) <= tol2)
	foundHmax = 1;
      if(ABS(hess[i].h[j]-hmin) <= tol3)
	foundHmin = 1;
    
         h2[j]=hess[i].h[j];
         for(int k=0; k<3; k++) {
             dir2[j*3+k] = hess[i].dir[j][k];
         }
    }

    EN_attachDataPtr( (pEntity)vertex, nodalSizeID, (void *)h2);
    EN_attachDataPtr( (pEntity)vertex, nodalDirectionID, (void *)
                                             dir2);


    if(foundHmin)
      hminCount++;
    if(foundHmax)
      hmaxCount++;
    if(foundHmin && foundHmax)
      bothHminHmaxCount++;
  
    ++i;  
  } //vertex loop
  
  VIter_reset(vit);
  i=0;

  int numSmooth=3;
  for (int k=0; k<20; k++){
    SmoothSize(mesh,numSmooth); //Size field smoothing similar to hessians    
    commuSmoothSize(pmesh, mesh,numSmooth);
  }
/*  
  for (int k=0; k<0; k++) {
     SmoothDir(mesh);
     commuSmoothDir(pmesh, mesh); 
  }
*/
  if(isSizeLimit)
     SizeLimit(mesh);
//  SizePropogate(mesh);
/*
  int ReadSF=1;
  double* SizeField = new double[nshg*9];
  if(ReadSF=1) {
    WriteSizeField(mesh);
//    ReadSizeField(SizeField);
  }
*/
  while ( vertex=VIter_next(vit)) {

    // like in parallel plates (i.e, w=f(y)) 
    // if a user wants different hmax in x-direction
    // other examples can be situations where
    // user don't want to have one mesh edge 
    // connecting two geometric model edges (like no dofs)
    // ModifyMetric(vertex,hess[i].dir,hess[i].h);

    // set the data in MeshSim format
    // M'= R*Lambda   (R:eigenvecs, Lambda:eigenVals (diag)
    double* h2 = new double[3];  

    EN_getDataPtr((pEntity)vertex,nodalSizeID ,
                      (void**)&h2);
/*    
    double* SizeVert = new double[9];
    for(int j=0; j<9; j++) {
       SizeVert[j]=SizeField[i*9+j];
    }
*/
    for (int iRow=0; iRow<3; iRow++) {

       if(h2[iRow] < hmin) h2[iRow] = hmin;
       if(h2[iRow] > hmax) h2[iRow] = hmax;

       hess[i].h[iRow] = h2[iRow];
//      if(option==9)
//	     hess[i].h[iRow] = localHmin;
      for(int iDir=0; iDir<3; iDir++) {
        hess[i].dir[iRow][iDir]=hess[i].dir[iRow][iDir]*hess[i].h[iRow];
//         hess[i].dir[iRow][iDir]=SizeVert[iRow*3+iDir];
      }
    }
//    EN_deleteData((pEntity)vertex, nodalSizeID);
    EN_modifyDataPtr( (pEntity)vertex, nodalSizeID, (void *)h2);
   
#if  ( defined  DEBUG )
    // double c[3];
    // V_coord(vertex,c);
    // printf("%i %f %f %f ",EN_id((pEntity)vertex),c[0],c[1],c[2]);
    // printf("%f %f %f ",hess[i].dir[0][0],hess[i].dir[0][1],hess[i].dir[0][2]);
    // printf("%f %f %f ",hess[i].dir[1][0],hess[i].dir[1][1],hess[i].dir[1][2]);
    // printf("%f %f %f\n",hess[i].dir[2][0],hess[i].dir[2][1],hess[i].dir[2][2]);
#endif   

    if(!preLBforAdaptivity) {
      if(option==1 || option==2) {      
   	MSA_setAnisoVertexSize(simAdapter,
			       vertex,
			       hess[i].dir);
      }
      else if(option==9) {
//	MSA_setVertexSize(simAdapter,
//			  vertex,
//			  localHmin);
      }
    }
    else {
      if(option==1  || option==2) {
	double *metric = new double[9];

	for(int iRow=0; iRow<3; iRow++)
	  for(int iDir=0; iDir<3; iDir++)
	    metric[iRow*3+iDir] = hess[i].dir[iRow][iDir];

	EN_attachDataPtr((pEntity)vertex,meshSizeID,metric);
      }
      else if(option==9) {
	double *size = new double;
//	*size = localHmin;

	EN_attachDataPtr((pEntity)vertex,meshSizeID,size);
      }
    }

    ++i;
/*    if(EN_isBLEntity(vertex)) {
       printf("BL entity found!\n");
    }
*/    
  }
  VIter_delete(vit);
//  M_writeVTKFile(mesh, "nodalSize", nodalSizeID, 3);
#ifdef DEBUG  
//  M_writeVTKFile(mesh, "nodalSize", nodalSizeID, 3);
#endif  
  if(PMU_size()==1) {
     pMesh meshMerge;
     meshMerge = M_createFromParMesh(pmesh, 3, prog);
     M_write(meshMerge, "mesh_size.sms", 0, prog);
     M_release(meshMerge);
  } else {
    PM_write(pmesh, "mesh_size.sms", prog);
  }

  printf(" [%d] Nodes with hmin into effect : %d (%4.2f%%)\n",PMU_rank(),hminCount,((double)hminCount/nshg)*100);
  printf(" [%d] Nodes with hmax into effect : %d (%4.2f%%)\n",PMU_rank(),hmaxCount,((double)hmaxCount/nshg)*100);
  printf(" [%d] Nodes with both hmin/hmax into effect : %d (%4.2f%%)\n",PMU_rank(),bothHminHmaxCount,((double)bothHminHmaxCount/nshg)*100);
  printf("\n");

//   adaptSimLog<<"Nodes with hmin into effect : "<<hminCount<<endl;
//   adaptSimLog<<"Nodes with hmax into effect : "<<hmaxCount<<endl;
//   adaptSimLog<<"Nodes with both hmin/hmax into effect : "<<bothHminHmaxCount<<endl;
//   adaptSimLog.close();

  fprintf(adaptSimLog,"Nodes with hmin into effect : %d (%4.2f%%)\n",hminCount,((double)hminCount/nshg)*100);
  fprintf(adaptSimLog,"Nodes with hmax into effect : %d (%4.2f%%)\n",hmaxCount,((double)hmaxCount/nshg)*100);
  fprintf(adaptSimLog,"Nodes with both hmin/hmax into effect : %d (%4.2f%%)\n",bothHminHmaxCount,((double)bothHminHmaxCount/nshg)*100);
  fclose(adaptSimLog);

  // assuming one partition per proc
  int gid = PMU_gid(PMU_rank(),0);

  // hess is a Hessian struct
  // 3rd arg is timestep
  writeMEDITSizeField(hess,mesh,timeStepNumber,gid);

  delete [] hess;

  if(option==1 || option==2 || option==9) {
    cleanAttachedData(mesh,nodalGradientID,0);
    cleanAttachedData(mesh,nodalHessianID,0);
    cleanAttachedData(mesh,localPatchVolID,0);
    cleanAttachedData(mesh,localGradientID,0);
    cleanAttachedData(mesh,localHessianID,0);
    cleanAttachedData(mesh,numSurroundNodesID,0);
    cleanAttachedData(mesh,locMaxInterpolErrorID,0);
    cleanAttachedData(mesh,globMaxInterpolErrorID,0);  
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if(preLBforAdaptivity) {
#ifndef ibm
//    printf("\n[%2d] memory usage before partioning: %d (KB)\n\n",PMU_rank(),phParAdaptProcSize());
#endif
    int partLB4AptOption;
    if(option==1 || option==2)
       partLB4AptOption = 1;
    else if(option==9)
      partLB4AptOption = 9;
    partitionMeshToLoadBalanceForAdaptivity(pmesh,mesh,partLB4AptOption,nErrorVars);

#ifndef ibm
//    printf("\n[%2d] memory usage after partioning: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
#endif
    pMesh LBmesh = PM_mesh(pmesh,0);
    mesh = LBmesh;

    VIter vIter = M_vertexIter(LBmesh);
    while(vertex = VIter_next(vIter)) {

      double *metric;
      if(!EN_getDataPtr((pEntity)vertex,meshSizeID,(void**)&metric)) {
	printf("\nerror in setSizeFieldUsingHEssians: no data attached with meshSizeID to vertex\n");
	exit(0);
      }

      if(option==1 || option==2) {
	double simMetric[3][3];

	for(int iRow=0; iRow<3; iRow++)
	  for(int iDir=0; iDir<3; iDir++)
	    simMetric[iRow][iDir] = metric[iRow*3+iDir];

	MSA_setAnisoVertexSize(simAdapter,vertex,simMetric);
      }
      else if(option==9) {

	MSA_setVertexSize(simAdapter,vertex,*metric);
      }
    }
    VIter_delete(vIter);

    if(option==1 || option==2) {
      cleanAttachedData(LBmesh,meshSizeID,0);
    }
    else if(option==9) {
      cleanAttachedData(LBmesh,meshSizeID,0,0);
    }
  }
  MD_deleteMeshDataId(meshSizeID);

  MPI_Barrier(MPI_COMM_WORLD);
}

// max relative interpolation error at a vertex
// H is the hessian as it comes via nodalHessianID
// in parallel:
// only takes edges that are local to the partition/proc
// Commu required for overall max at a particular vertex
double maxLocalError(pVertex vertex, double H[3][3])
{
  pEdge edge;
  double locE;
  double maxLocE=0;

  for(int i=0; i<V_numEdges(vertex); i++ ) {
    edge=V_edge(vertex,i);
    locE=E_error(edge,H);
    if( locE > maxLocE )
      maxLocE=locE;
  }

  return maxLocE;
}

// relative interpolation error along an edge
double E_error(pEdge edge, double H[3][3])
{
  double xyz[2][3], vec[3];
  double locE=0;
  int i,j;

  V_coord(E_vertex(edge,0),xyz[0]);
  V_coord(E_vertex(edge,1),xyz[1]);
  diffVt(xyz[0],xyz[1],vec);

  for( i=0; i<3; i++ )
    for( j=0; j<3; j++ ) 
      locE += H[i][j]*vec[i]*vec[j];

  return ABS(locE);
}



#ifdef __cplusplus
}
#endif
