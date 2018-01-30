#ifndef _H_CCFUNC_
#define _H_CCFUNC_

#include "parallel.h"

/* #ifndef SIM */
/* using namespace SCOREC_mesh; */
/* using namespace SCOREC_model; */
/* #endif */
#include "mesh_interface.h"

#ifdef __cplusplus
void printPeriodicBC(pMesh mesh);
void mdb2phasta(char* fname,char* mname,pGModel model,std::vector<pParMesh> pmesh,int nshgTot );

int attachEssentialBC(pMesh mesh, globalInfo* info, double *qtot, int *nBC, int *iBC, double **BC);
void attachNaturalBC(pMesh mesh, int ***iBCB, double ***BCB, int numNBC);
void attachInitialCondition(pMesh mesh, globalInfo*, double *qTot,
                            double **q);
void setPeriodic(pGModel model);



#ifdef PARALLEL
void setup(pParMesh pmesh, globalInfo *info, int ipart);

// not to be used in phParAdapt
extern "C" void allpart(pParMesh pmesh, pPartitioner up, int chkwgt,void*);
extern "C" void partMeshBypass(pParMesh pmesh, pPartitioner up, int chkwgt,void*);


void writeEnsaFiles(pParMesh pmesh, globalInfo *info, pPList bdry, int ipart);
int partitionMesh(pParMesh pmesh);
#else
void setup(pMesh mesh, globalInfo *info);
void writeEnsaFiles(pMesh mesh, globalInfo *info, pPList bdry);
#endif

//pPList return_vertices(pRegion region);
void lin2quad(pMesh mesh, pMesh lmesh, double* qtmpl, double* qtmp, int nshg);
void genblock(pMesh mesh, globalInfo *info, pPList bdry);
void face_coordinate_extraction(pMesh mesh, pGModel model);

void restart(double* q, int nshgTot, int nshgReqst, int nvReqst, int* lstep, pMesh mesh, int ipart);

////////////////////////////////////////////////////////////////////////////////
// handles the fork to adapt/preprocess
////////////////////////////////////////////////////////////////////////////////
extern "C" int
switchAdapt_Preproc(int argc, char *argv[]);
#endif

#endif
