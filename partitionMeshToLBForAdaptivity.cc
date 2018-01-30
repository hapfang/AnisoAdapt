#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef SIM
#include "SimPartitionedMesh.h"
#include "SimMeshTools.h"
#endif
#include "phParAdapt.h"
#include "mpi.h"



using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

double optVolInv = sqrt(72.);

extern double masterProcWgt;

extern pMeshDataId errorIndicatorID;

extern pMeshDataId meshSizeID;

extern pProgress prog;

void partitionMeshToLoadBalanceForAdaptivity(pParMesh pmesh,
					     pMesh mesh,
					     int option,
					     int nCurrentErrorVars)
{

  printf("-- Before partitioning for LB...\n");
  printf(" [%2d]: (BP) local # of vertices: %d\n", PMU_rank(),M_numVertices(mesh));
  printf(" [%2d]: (BP) local # of edges   : %d\n", PMU_rank(),M_numEdges(mesh));
  printf(" [%2d]: (BP) local # of faces   : %d\n", PMU_rank(),M_numFaces(mesh));
  printf(" [%2d]: (BP) local # of regions : %d\n", PMU_rank(),M_numRegions(mesh));

  pMeshDataId newRgnsID = MD_newMeshDataId("new rgns");

  double nNewLowerLimit = 0.01, nNewUpperLimit = 100;
  double minNumNewRgnsLoc = 1e12, maxNumNewRgnsLoc = 0.;
  double newRgnsCountLoc = 0.;
  pRegion rgn;
  RIter rit = M_regionIter(mesh);
  while(rgn = RIter_next(rit)) {

    double numNewRgns = estimateNumNewRegions(rgn,option);

    newRgnsCountLoc += numNewRgns;

    if(minNumNewRgnsLoc>numNewRgns)
      minNumNewRgnsLoc = numNewRgns;
    if(maxNumNewRgnsLoc<numNewRgns)
      maxNumNewRgnsLoc = numNewRgns;

    EN_attachDataDbl((pEntity)rgn,newRgnsID,numNewRgns);
  }
  RIter_reset(rit);

  MPI_Barrier(MPI_COMM_WORLD);

  double minNumNewRgnsGlobal;
  MPI_Allreduce(&minNumNewRgnsLoc, &minNumNewRgnsGlobal, 1, MPI_DOUBLE, 
		MPI_MIN, MPI_COMM_WORLD);

  double maxNumNewRgnsGlobal;
  MPI_Allreduce(&maxNumNewRgnsLoc, &maxNumNewRgnsGlobal, 1, MPI_DOUBLE, 
		MPI_MAX, MPI_COMM_WORLD);

  double newRgnsCountGlobal;
  MPI_Allreduce(&newRgnsCountLoc, &newRgnsCountGlobal, 1, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);

  printf("\n[%d] local units of estimated load (before partitioning): %d\n",PMU_rank(),(int)newRgnsCountLoc);

  if(PMU_rank()==0) {
    printf("\nTotal units of ESTIMATED load (before partitioning): %d (roughly)\n",(int)newRgnsCountGlobal);

    printf("\nMin. (global) units of load for an old region: %f",minNumNewRgnsGlobal);
    printf("\nMax. (global) units of load for an old region: %f\n",maxNumNewRgnsGlobal);
    printf("\n(lower limit on load units to apply weight for partitioning: %f)",nNewLowerLimit);
    printf("\n(upper limit on load units to apply weight for partitioning: %f)\n\n",nNewUpperLimit);
  }

  double scaleWgtNumNewRgns = minNumNewRgnsGlobal;
  if(minNumNewRgnsGlobal<nNewLowerLimit)
    scaleWgtNumNewRgns = nNewLowerLimit;
  double wgtUpperLimit = nNewUpperLimit/scaleWgtNumNewRgns;

  pMeshDataId rWtID = EntityWeightId_new("rgn weight id");
  RIter_reset(rit);
  while(rgn = RIter_next(rit)) {
    double numNewRgns;
    if(!EN_getDataDbl((pEntity)rgn,newRgnsID,&numNewRgns)) {
      printf("\nError in partitionMeshForLB : no data attached to region with newRgnsID");
      exit(1);
    }

    int weight;
    double weightDbl = numNewRgns/scaleWgtNumNewRgns;
    weight = (int)weightDbl;

    if(weightDbl<1)
      weight = 1;
    else if(weightDbl>wgtUpperLimit)
      weight = (int)wgtUpperLimit;

    EN_setWeight(rgn,rWtID,weight);

    EN_deleteDataDbl((pEntity)rgn,newRgnsID);
  }
  RIter_delete(rit);
  MD_deleteMeshDataId(newRgnsID);

  MPI_Barrier(MPI_COMM_WORLD);

  int nVarsForSize = 9;
  if(option==9)
    nVarsForSize = 1;

  pAttachDataCommu adcSize = AttachDataCommu_new(sizeof(double)/sizeof(int),0,nVarsForSize);
  MD_setMeshCallback(meshSizeID, CBmigrateOut, pm_sendDblArray, adcSize);
  MD_setMeshCallback(meshSizeID, CBmigrateIn,  pm_recvDblArray, adcSize);
  if(nVarsForSize>1)
    MD_setMeshCallback(meshSizeID, CBdelete,  delDblArray, NULL);
  else
    MD_setMeshCallback(meshSizeID, CBdelete,  delDbl, NULL);
  PM_setMigrId(pmesh, meshSizeID);

  pAttachDataCommu adcEI = AttachDataCommu_new(sizeof(double)/sizeof(int),0,nCurrentErrorVars);
  MD_setMeshCallback(errorIndicatorID, CBmigrateOut, pm_sendDblArray, adcEI);
  MD_setMeshCallback(errorIndicatorID, CBmigrateIn,  pm_recvDblArray, adcEI);
  if(nCurrentErrorVars>1)
    MD_setMeshCallback(errorIndicatorID, CBdelete,  delDblArray, NULL);
  else
    MD_setMeshCallback(errorIndicatorID, CBdelete,  delDbl, NULL);
  PM_setMigrId(pmesh, errorIndicatorID);

  MPI_Barrier(MPI_COMM_WORLD);

  pPartitionOpts popts = PartitionOpts_new();
  PartitionOpts_setTotalNumParts(popts,PM_totalNumParts(pmesh));
  if(masterProcWgt>0.) {
    int tnp = PM_totalNumParts(pmesh);
    double *partwts = new double[tnp];
    
    partwts[0] = masterProcWgt;
    for(int iPart=1; iPart<tnp; iPart++)
      partwts[iPart] = 1.;

    if(tnp>1)
      PartitionOpts_setProcWt(popts, partwts); 
    else
      PartitionOpts_setPartWtEqual(popts);

    delete [] partwts;
  }
  else
    PartitionOpts_setPartWtEqual(popts);
  PartitionOpts_setNumEntityWeightIds(popts,1); 
  PartitionOpts_setEntityWeightId(popts,0,rWtID);

  PM_partition(pmesh,popts,prog);

  PM_removeMigrId(pmesh,errorIndicatorID);
  PM_removeMigrId(pmesh,meshSizeID);

  MD_removeMeshCallback(errorIndicatorID,CBmigrateOut);
  MD_removeMeshCallback(errorIndicatorID,CBmigrateIn);
  MD_removeMeshCallback(errorIndicatorID,CBdelete);

  MD_removeMeshCallback(meshSizeID,CBmigrateOut);
  MD_removeMeshCallback(meshSizeID,CBmigrateIn);
  MD_removeMeshCallback(meshSizeID,CBdelete);

  PartitionOpts_delete(popts);
  AttachDataCommu_delete(adcSize);
  AttachDataCommu_delete(adcEI);
  EntityWeightId_delete(rWtID);

  pMesh LBmesh = PM_mesh(pmesh,0);
  mesh = LBmesh;

  newRgnsCountLoc = 0.;
  RIter rIter = M_regionIter(LBmesh);
  while(rgn = RIter_next(rIter)) {
    double numNewRgns = estimateNumNewRegions(rgn,option);

    newRgnsCountLoc+=numNewRgns;
  }
  RIter_delete(rIter);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&newRgnsCountLoc, &newRgnsCountGlobal, 1, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);

  printf("\n[%d] local units of estimated load (after partitioning): %d\n",PMU_rank(),(int)newRgnsCountLoc);

  if(PMU_rank()==0)
    printf("\nTotal units of ESTIMATED load (after partitioning): %d (roughly)\n\n",(int)newRgnsCountGlobal);

  printf("-- After partitioning for LB...\n");
  printf(" [%2d]: (AP) local # of vertices: %d\n", PMU_rank(),M_numVertices(LBmesh));
  printf(" [%2d]: (AP) local # of edges   : %d\n", PMU_rank(),M_numEdges(LBmesh));
  printf(" [%2d]: (AP) local # of faces   : %d\n", PMU_rank(),M_numFaces(LBmesh));
  printf(" [%2d]: (AP) local # of regions : %d\n", PMU_rank(),M_numRegions(LBmesh));
}

double estimateNumNewRegions(pRegion rgn, int option) {

  double avgOfDetM = 0.;

  pVertex vtx;
  pPList rverts = R_vertices(rgn,1);
  int numVerts = PList_size(rverts);

  if(option==1) {
    for(int iRVert=0; iRVert<numVerts; iRVert++) {
      vtx = (pVertex)PList_item(rverts,iRVert);

      double detMInv = 1.;
      double *metric;
      if(!EN_getDataPtr((pEntity)vtx,meshSizeID,(void**)&metric)) {
	printf("\nerror in estimateNumNewRegions: no data attached with meshSizeID to vertex (option==1)\n");
	exit(0);
      }

      for(int iRow=0; iRow<3; iRow++) {
	double magSq = 0.;
	for(int iDir=0; iDir<3; iDir++) {
	  magSq += metric[iRow*3+iDir]*metric[iRow*3+iDir];
	}
	detMInv = detMInv*magSq;
      }
      avgOfDetM += 1./(detMInv);
    }
    avgOfDetM /= numVerts;
  }
  else if(option==9) {
    double avgSize = 0.;
    for(int iRVert=0; iRVert<numVerts; iRVert++) {
      vtx = (pVertex)PList_item(rverts,iRVert);

      double *metric;
      if(!EN_getDataPtr((pEntity)vtx,meshSizeID,(void**)&metric)) {
	printf("\nerror in estimateNumNewRegions: no data attached with meshSizeID to vertex (option==9)\n");
	exit(0);
      }
    
      double inv = 1.;
      double temp = (*metric)*(*metric);
      for(int iDir=0; iDir<3; iDir++)
	inv = inv*temp;

      avgOfDetM += 1/inv;
    }
    avgOfDetM /= numVerts;
  }

  PList_delete(rverts);

  return R_volume(rgn)*sqrt(avgOfDetM)*optVolInv;
}


#ifdef __cplusplus
}
#endif
