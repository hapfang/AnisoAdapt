#ifndef H_MeshSimInternal1
#define H_MeshSimInternal1
#include "MeshSimAdapt.h"

//from PMops.h
#ifdef B64
  typedef long P_int;
#else
  typedef int P_int;
#endif


typedef double dArray[3];

#ifdef SIM
#include "MeshSim.h"
#include "SimPartitionedMesh.h"
/* #ifdef FMDB */
/* #include "AOMD.h" */
/* #include "AOMD_cint.h" */
/* #endif */


/* These functions are no longer supported by simmetix 6.2 */
/* void diffVt(const double*,const double*,double*); */
/* int normVt(const double*,double*); */
/* double dotProd(const double*,const double*); */
/* void crossProd(const double*,const double*,double*); */
void F_coord(pFace,double(*)[3],int);
void M_setTolerance(pMesh);
void MSA_writePartData(pMSAdapt adapter, int rank, int iter, const char *optsfile, const char *meshfile );

//int GEN_dataI(pGEntity , const char *tag);
//void * GEN_dataP(pGEntity , const char *tag);
pEdge E_exists(pVertex v1, pVertex v2);
pEdge R_gtOppEdg(pRegion region, pEdge edge);
pFace   R_vtOpFc(pRegion ,pVertex);
pVertex R_fcOpVt(pRegion rgn, pFace face);

void R_info(pRegion);
void F_info(pFace);
void E_info(pEdge);

void F_setWhatIn(pFace  , pGEntity what);
void E_setWhatIn(pEdge  , pGEntity what);

pVertex F_oppositeVertex(pFace face, pEdge edge);

pVertex M_nextVertex(pMesh, void **restart);

void F_normalVector(pFace face, int dir, double* normal);
void F_chDir(pFace);

/* void GEN_attachDataP(pGEntity , const char *tag, void *data); */
/* void GEN_attachDataI(pGEntity , const char *tag, int data); */
/* int GEN_modifyDataP(pGEntity, const char *tag, void * data); */
/* void MSA_setPartWts(pMSAdapt, const double *partwts); */
pFace F_exists( eType type, pEntity e1, pEntity e2, pEntity e3, pEntity e4 );

// left overs
// JM 8/9/04
pPList R_verticesLeft(pRegion);
void V_info(pVertex);
pFace M_nextFaceCancel(pMesh, void **restart);
pRegion M_nextRegionCancel(pMesh, void **restart);


typedef struct BdryProcIter *PBEntProcIter;
typedef struct BdryProcIter *PBFProcIter;
typedef struct BdryProcIter *PBEProcIter;
typedef struct BdryProcIter *PBVProcIter;

typedef struct BdryPartIter *PBEntPartIter;
typedef struct BdryPartIter *PBFPartIter;
typedef struct BdryPartIter *PBEPartIter;
typedef struct BdryPartIter *PBVPartIter;

typedef struct BdryNPartsIter *PBEntNPartsIter;
typedef struct BdryNPartsIter *PBFNPartsIter;
typedef struct BdryNPartsIter *PBENPartsIter;
typedef struct BdryNPartsIter *PBVNPartsIter;

#endif /*ifdef SIM*/
#endif /* not H_MeshSimInternal */
