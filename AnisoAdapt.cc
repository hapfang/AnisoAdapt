#include "SimPartitionedMesh.h"
#include "SimParasolidKrnl.h"
#include "SimAdvMeshing.h"
#include "MeshSimAdapt.h"
#include "SimModel.h"
#include "SimUtil.h"
#include <stdio.h>
#include <set>

#include "phReadWrite.h"
#include "attachData.h"

pMeshDataId phasta_solution;
pMeshDataId errorIndicatorID;
pMeshDataId ybarID;
pMeshDataId nodalGradientID;
pMeshDataId nodalHessianID;
pMeshDataId nodalSizeID;
pMeshDataId wallStressID;
pMeshDataId wallDistID;
pMeshDataId SeparatedID;
pMeshDataId ShockID;
pMeshDataId isOrgNodeID;
pMeshDataId nodalVorticityID;
pMeshDataId OrgSizeID;
pMeshDataId nodalDirectionID;
pMeshDataId QualitySizeID;
pMeshDataId modes;
pMeshDataId incorp;
pMeshDataId localGradientID;
pMeshDataId localPatchVolID;
pMeshDataId localHessianID;
pMeshDataId numSurroundNodesID;
pMeshDataId locMaxInterpolErrorID;
pMeshDataId globMaxInterpolErrorID;
pMeshDataId meshSizeID;

//Used in elementGradient.cc
double refthreshold = 0.5; //this value is found at line 25 of Input.cc
double rhoinp = 1.225; //this value is found at line 52 of assignGlobalVars.cc
//Used in setSizeFieldUsingHessians.cc
int isSizeLimit = 0;
int preLBforAdaptivity = 0; //line 28 of Input.cc
int timeStepNumber = 0;
int nErrorVars;
//Used in SizeLimit.cc
int MinLimitFact = 2;
int MaxLimitFact = 2; 
//Used in localInfo.cc
int numParts=1;
int globalP=1;
int isReorder = 0;
//Used in setPeriodic.cc
int prCd=0;
//Used in callback.cc
int CBcounter = 0;
int delDblArraycounter = 0;
int delDblcounter = 0;
int DisplacementMigration = 0;
int dwalMigration = 0; //May enable it later
int numVars; //Caution
//
int ensa_dof=6;
//Used in partitionMeshToLBForAdaptivity.cc
double masterProcWgt = 0.0;
pProgress prog;

extern "C" void hessiansFromSolution(pParMesh pmesh, pMesh mesh,int stepNumber);
extern "C" void setSizeFieldUsingHessians(pParMesh pmesh,
                               pMesh mesh,
                               pMSAdapt simAdapter,
                               double factor,
                               double hmax,
                               double hmin,
                               int option);

int main(int argc, char **argv)
{
  SimPartitionedMesh_start(&argc, &argv);
  SimParasolid_start(1);
  SimAdvMeshing_start();
  Sim_readLicenseFile(NULL);     // Fix this

  prog = Progress_new();
  Progress_setDefaultCallback(prog);
  pNativeModel nmodel = ParasolidNM_createFromFile("geom_nat.x_t", 0);
  pGModel model = GM_load("geom.smd", nmodel, prog);
  pParMesh pmesh = PM_load("geom.sms", model, prog);
  pMesh mesh = PM_mesh(pmesh, 0);

// Task to do, load the adapt.inp parameters


  pMSAdapt msa = MSA_new(pmesh, 1); //0 - tag-based; 1 - size-based
// Here are parameters we set within adapt.cc
  int isBLAdapt=1;
  int localAdapt=0;
  int coarsenMode=0;
  double factor=2e5; //need to double check the numbers below
  double hmin=1e-5;
  double hmax=1e6;
  int option=11;
//  MSA_setMaxIterations seems to be not included in the libraries in this script
//  Let us not worry about this parameter now
//  int numSplit=80; 
  char solution_file[256];

  MSA_setAdaptBL(msa, isBLAdapt);
  MSA_setLocal(msa, localAdapt);     
  MSA_setCoarsenMode(msa, coarsenMode);
  
  if(isBLAdapt==1) {
     MSA_setExposedBLBehavior(msa, BL_DisallowExposed);
//     MSA_setExposedBLBehavior(simAdapter, bl_DisallowExposed);
     MSA_setBLMinLayerAspectRatio(msa, 1.0);
  }

  if(PMU_size()>1 && isBLAdapt==1) {
//     MSA_setCoarsenMode(simAdapter, 0);
//     MSA_setBLSnapping(simAdapter, 0);
     MSA_setBLMinLayerAspectRatio(msa, 0.0);
     MSA_setExposedBLBehavior(msa, BL_DisallowExposed);
   }

//   MSA_setMaxIterations(msa, numSplit);
   MSA_setBoundaryMeshModification(msa, isBLAdapt);
// End of paste of parameter settings copied from adapt.cc

// reading solution extracted from phParAdapt.....still needs PHASTAIO link
   int lstep=500;
   int ndof=6;
   int poly=1;
   sprintf(solution_file,"restart.%i.%i",lstep,PMU_rank()+1);
   printf(" ...%s (for \"solution\")\n",solution_file);
// attaching the solution to the original mesh
   double *sol;
   readArrayFromFile(solution_file,"solution",sol);
   attachArray(sol,mesh,phasta_solution,ndof,poly);
   delete [] sol;

   hessiansFromSolution(pmesh,mesh,lstep);
   setSizeFieldUsingHessians(pmesh,mesh,msa,factor,hmax,hmin,option);

// size-based adapt
//  pMSAdapt msa = MSA_new(pmesh, 1);
  MSA_adapt(msa, prog);
  MSA_delete(msa);

  PM_write(pmesh, "adapted.sms", prog);

  M_release(pmesh);
  GM_release(model);
  NM_release(nmodel);
  Progress_delete(prog);

  Sim_unregisterAllKeys();
  SimAdvMeshing_stop();
  SimParasolid_stop(1);
  SimPartitionedMesh_stop();
  return 0;
}
