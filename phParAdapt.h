#ifndef _PHADAPT_H_
#define _PHADAPT_H_

#include <string>
#include <iostream>
#include <map>
#include "mpi.h"
#include "MeshSimInternal.h"
#include "MeshSimAdapt.h"
#include "mesh_interface.h"
#include "parallel.h" //include definition of globalInfo
using namespace std;

#ifndef ibm
#define ludcmp ludcmp_
#define lubksb lubksb_
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x)<(y) ? (y) : (x))
  
  struct Hessian {
    double h[3];
  double dir[3][3];
  };
  typedef struct Hessian Hessian;  

  int phParAdaptProcSize();

  // sim call back  
  // typedef int (*MSA_CallbackFunc)(pPList , pPList, int);

 
  // main routine to control mesh adaptation, solution transfer and communication with PHASTA
  int adapt(  // parallel mesh
            pParMesh pmesh,
            //serial mesh (BLAdapt)
            pMesh mesh,
            // model
            pGModel model,
            //time step
	    int timeStepNumber,
	    // strategy is to specify 
	    // how to do adaptation (i.e., size-field or tag driven)
	    // 1-2 : size-field driven (for anisotropic) 
	    // 3-4 : tag driven (for isotropic)
	    // 5-6 : size-field driven (for isotropic)
	    // < 0 sets a manual mesh size-field
	    int strategy,
	    // factor is the constant appearing in the error expression
	    // for tag driven it is used to define threshold for refinement
	    // for size-field driven it is used to define the error tolerance
	    double factor,
	    // number of solution variables (5 for incompressible)
	    int ndof,
	    // number of variables for error indicators (EI)
	    // (e.g., 5 for ybar & 10 for residual-based)
	    int nvar,
	    // the maximal mesh edge length allowed in mesh size 
	    double hmax,
	    // the minimal mesh edge length allowed in mesh size
	    double hmin,
	    // option is used to decide how to compute the error value
	    // provides different choices like analytic hessian, manual size-field etc.
	    // for isotropic (ta or size-field driven :
	    // use 3 EI for flow problem or use 1 EI for scalar problem
	    int option);


  // piecewise linear mesh size field definition in terms of Hessian
//  void setSizeFieldUsingHessians(pMesh,pMSAdapt,double,double,double,int option=1);

  void partitionMeshToLoadBalanceForAdaptivity(pParMesh pmesh, pMesh mesh, int option, int nCurrentErrorVars);

  double estimateNumNewRegions(pRegion rgn, int option);

  void getGlobalErrorInfo(pMesh mesh, double& totalError, double& sumOfError);
  void computeOldMeshSize(pMesh, int option);
  void commuOldMeshSize(pParMesh, pMesh);

  void setIsotropicSizeField(pGModel,pParMesh,pMesh,pMSAdapt,double,double,double,int option );

  void setManualSizeField(pParMesh pmesh, pMesh mesh, pMSAdapt simAdapter, int comp, int option);
    
  void setSizeFieldUsingHessians(pParMesh pmesh,
                                 pMesh mesh,
				 pMSAdapt simAdapter,
				 double factor,
				 double hmax,
				 double hmin,
                                 int option);  

  void setSizeFieldUsingCombined(pParMesh pmesh,
                                 pMesh mesh,
				 pMSAdapt simAdapter,
				 double factor,
				 double hmax,
				 double hmin,
                                 int option);   

  void SizeQuality(pMesh mesh);  

  void setSizeFieldFromFile(pMesh mesh,
			    pMSAdapt simAdapter,
			    char* sfFilename);


  // for tag driven refinement (i.e., isotropic)
  void tagEntitiesForRefinement(pMesh,pMSAdapt,double, double,double,int option);

    int applyMarkingStrategy(pMesh mesh,
                             pMSAdapt simAdapter,
                             double factor,
                             double hmax,
                             double hmin,
                             double totalError, 
                             double maxError, 
                             double minError, 
                             double threshold, 
                             int option);

  double getErrorThreshold(pMesh mesh,
                           double factor,
                           double totalError,
                           double maxError,
                           double minError,
                           int option);
    void V_reordering(pMesh mesh, globalInfo* ginfo, int ipart);
    void R_reordering(pMesh mesh, int ipart);

  double getErrorValue(double *nodalValues,int option);

  // get interpolation error values from hessians  
  double maxLocalError(pVertex vertex, double H[3][3]);

  // attach the max local interpolation error LOCAL
  // to each partition via locMaxInterpolErrorID;
  void maxLocalPartLocError(pMesh mesh);

  // attach the max local interpolation error GLOBAL
  // to the parallel mesh via globMaxInterpolErrorID;
  void commuMaxLocalPartLocError(pParMesh pmesh,pMesh mesh);  

  double E_error(pEdge edge, double H[3][3]);

  // get hessians computed from phasta
  void V_getHessians(double*,pMesh,int,int,double*);

  // for doing hessians computation outside phasta
  void V_Hessian(pVertex v, double T[3][3]);
  void V_AnalyticHessian(pVertex v, double T[3][3], int option);
  void buildSystem(pRegion region, double* eltMatrix);
  void buildSystemXYZ(dArray *xyz, double* eltMatrix);
  void elementGradient(pRegion region, double* elemGradient);
  void gradientsFromPatch(pMesh mesh);
  void elementHessian(pRegion region, double* elemHessian);
  void hessiansFromSolution(pParMesh pmesh, pMesh mesh,int stepNumber);
  void hessiansFromPatch(pMesh mesh);
  void ModifyHessiansAtBdry(pMesh mesh);
  void writeRestartHessians(pMesh mesh );
  void SmoothHessians(pMesh mesh);
  void SmoothSize(pMesh mesh, int num);
  void SmoothDir(pMesh mesh);
  void commuSmoothDir(pParMesh pmesh, pMesh mesh);
  void SmoothWallField(pMesh mesh, pMeshDataId field, int ndof);
  void SizeLimit(pMesh mesh);
  void SizePropogate(pMesh mesh);
  void SmoothSizeOnBdry(pMesh mesh, int numPasses, int genDim);
  void SmoothErrorIndicators(pMesh mesh,int option);
  void writeSmoothEIs(pMesh mesh);

  /// ad hoc things
  void WriteSizeField(pMesh mesh);
  void ReadSizeField(double* SizeField);

  void ModifyMetric(pVertex vertex, double dir[3][3], double* h);

    // nvar mis number of error indicators attached via errorIndicatorID 
    void transformToScalarErrorVal(pMesh mesh, int nvar, int option);

//    double  processErrorAG(double* nodalErrorSet,int nvar); 
double  processErrorAG(double* nodalErrorSet,double* nodalSolutionSet,int nvar, int option, double* coord);


  // for solving linear system (small) 
  void ludcmp( double*, int*, int*, int*, double* );
  // the last array is the right hand side
  // it is being passed as a reference  and overridden
  // to contain the linear system's solution !
  void lubksb( double*, int*, int*, int*, double* );

  // for viewing results in medit 
  void writeMEDITSizeField(Hessian *hess,pMesh mesh, int currentTimestep, int gprocID );
  void writeMEDITSolution(pMesh mesh);

  // void MSA_setCallback(pMSAdapt, MSA_CallbackFunc);

#ifdef SIM
  void phastaTransferTopSIM(MeshModType mtype, pMeshChanges mco, void *userData);
  typedef MeshModType modType;
  int delDblArray(void *ent, pAttachDataId id, int cb, void** data, void* c);
  int delDbl(void *ent, pAttachDataId id, int cb, void** data, void* c);
#endif

#ifdef FMDB
  void phastaTransferTopSCOREC(pPList oldEnts, pPList newRgn, void *userData, modType mtype, pEntity ent);
  int delDblArray(pAttachableData ad, pAttachDataId id, int cb, void** data, void* c);
  int delDbl(pAttachableData ad, pAttachDataId id, int cb, void** data, void* c);
  void phastaTransferBottomE(pPList old, pPList fresh, pPList oldVertices, modType mtype);
#endif

  void phastaTransferBottom(pPList old, pPList fresh, pPList oldVertices, modType mtype);

#ifdef SIM
  void phastaTransferSIMBottom(pPList Vertices, pPList Edges, pPList Faces, pPList Regions, int ndof, modType mtype);
#endif

  void fix4SolutionTransfer(pMesh mesh);
  void commuFix4SolutionTransfer(pParMesh pmesh, pMesh mesh);

  void BCInflowFaceInfo(pGModel model, pParMesh pmesh, pMesh mesh);
  int BCInflowFaceNodesInfo(pGModel model, pMesh mesh);
  int BCInflowFaceConnectInfo(pGModel model, pMesh mesh);
  void commuBCInflowFaceInfo(pParMesh pmesh, pMesh mesh);
  void BCInflowFaceGlobalInfoInVtk(int nodesTot, int facesTot);

  void printTimeStatsToFile(int);

  int
  inverseMap( pRegion region,
            double gpt[3],
            double xisol[3] ) ;
  int
  inverseMapE( pEdge edge,
            double gpt[3],
            double xisol[3] );

  double*
  InterpolateSolution( pRegion region,
                     double xi[3],
                     int ndof,
                     pMeshDataId modes);
  double*
  InterpolateSolutionE( pEdge edge,
                     double xi[3], 
                     int dof,
                     pMeshDataId modes);

  void 
  R_entitiesAdapt( pRegion region, 
            pVertex *vrts, 
            pEdge *edgs,
            pFace *fcs ) ;


  void
  display_region( pRegion region );    

  // for debug
  void check(pMesh);

    ////////////////////////////////////////////////////////////////////////////////////////////
// function that reads in erro files in restart format
// using phastaIO library
////////////////////////////////////////////////////////////////////////////////////////////
void
readErrorFiles(double* nodalErrors, int stepNumber);

////////////////////////////////////////////////////////////////////////////////////////////
// function that reads in erro files in restart format
// using phastaIO library --> directly from restart format
////////////////////////////////////////////////////////////////////////////////////////////
void
readErrorFromRestart(double* nodalErrors, int stepNumber);



////////////////////////////////////////////////////////////////////////////////
// for solution transfer onto new mesh:
// read in solution and attach it to current mesh
// also provide for inter-proc transfer
////////////////////////////////////////////////////////////////////////////////
void
readAttachSolution(pMesh mesh,int stepNumber);

////////////////////////////////////////////////////////////////////////////////
// transfer the previously attached solution data onto
// new nodes (linear case)
// another migtation ID is required
//////////////////////////////////////////////////////////////////////////////// 
void
transferSolution(pMesh mesh,
                 pMeshDataId nodalAverageID,
                 int stepNumber);





////////////////////////////////////////////////////////////////////////////////
// in localSolutionContrib.cc
//
// retrieve the solution values on vertices that surround vertex.
// If surrounding vertex DOES carry an old solution AND is owned
// by proc/part the solution vals are added up and stored in
//  nodalValues
////////////////////////////////////////////////////////////////////////////////
/*  void */
/*  thisPartsSolutionContribVal(pVertex vertex, */
/*                              pMeshDataId nodalSolutionID, */
/*                              double* nodalValues); */


////////////////////////////////////////////////////////////////////////////////
// get number of vertices that contribute to this vertex' average
// counted are all the vertices which actually DO carry a solution
// and if they are on a partition bdry this proc must be owner
//////////////////////////////////////////////////////////////////////////////// 
/*  int  */
/*  thisPartsSolutionContribNN(pVertex vertex, */
/*                             pMeshDataId nodalSolutionID); */


////////////////////////////////////////////////////////////////////////////////
// gather solution values for vertices on partition bdries
// this function operates on 
// the data ptr attached to vertices via 
// nodalAverageID
// and is migrated tru partitions
// for each part, the data ptr carrying num of contribution nodes and their
// values is modifies (added to)
////////////////////////////////////////////////////////////////////////////////
/*  void */
/*  gatherNodalSolution(pParMesh pmesh,pMeshDataId nodalSolutionID,pMeshDataId nodalAverageID); */

////////////////////////////////////////////////////////////////////////////////
// write out the restart files
////////////////////////////////////////////////////////////////////////////////
void
writeRestartFiles(double* q, int nshg_fine, int stepNumber);


////////////////////////////////////////////////////////////////////////////////////////////
// get the global node number
// via incorp
// the mapping in the files that were created during partitioning
// returns the global number of DOFs
////////////////////////////////////////////////////////////////////////////////////////////
/*  int */
/*  getGlobalNodeId(pMesh mesh, pMeshDataId incorp, char* meshDirName); */


////////////////////////////////////////////////////////////////////////////////////////////
// just get the global (total) number of nodes
////////////////////////////////////////////////////////////////////////////////////////////
int
getNSHGTOT(pMesh mesh);

////////////////////////////////////////////////////////////////////////////////////////////
// assign new global IDs to the new mesh
// old node nums are kept
// the IDs for the new nodes are set according to their
// frequency in each partition
// also needs the total number of DOFs
// which has been determined in getGlobalNodeId
// (which has to be broadcasted)
// argument list : nshgTot : the old nshgTot
// returns new nshgTot 
////////////////////////////////////////////////////////////////////////////////////////////
int
assignNewGlobalNodeId(pParMesh pmesh, pMesh mesh,pMeshDataId incorp, int nshgTot);


void 
assignGlobalVars();

// replacement for the Simmetrix model attach data functions
// why replacement???? Min Zhou
#ifdef SIM
int
GEN_dataP(pGEntity g, char* str, void** data);
void
GEN_attachDataP(pGEntity g, char* str, void* data);
int
GEN_dataI(pGEntity g, char* str, int* data);
void
GEN_attachDataI(pGEntity g, char* str, int data);
#endif

////////////////////////////////////////////////////////////////////////////////////////////
// do the solution transfer
////////////////////////////////////////////////////////////////////////////////////////////
/*  int */
/*  solutionTransfer(pMesh mesh,int stepNumber); */

// parallel communication tasks for the adaptor
void
commuSmoothHessians(pParMesh pmesh, pMesh mesh);
void
commuSmoothSize(pParMesh pmesh, pMesh mesh, int num);

void
commuGradientsFromPatch(pParMesh pmesh, pMesh mesh);

void
commuHessiansFromPatch(pParMesh pmesh, pMesh mesh);
#ifdef __cplusplus
}
#endif

///////////////////////////////////////////////////////////////////////////////////
//This functions need to be overloaded, do NOT move them inside above
///////////////////////////////////////////////////////////////////////////////////
  void PMU_commuArr(void **s, int *ns, void **r, int *nr, MPI_Datatype type, int mult);
  void PMU_commuArr(void ***s, int **ns, void ***r, int **nr, map<int,int> *nsmap, map<int,int> *nrmap, MPI_Datatype type, int mult);
  void PMU_commuInt(int *ns, int *nr);
  void PMU_commuInt(int **ns, int **nr, map<int,int> *nsmap, map<int,int> *nrmap);

//////////////////////////////////////////////////////////////////////////////////
//Adding function from M_writeVTKFile
/////////////////////////////////////////////////////////////////////////////////
  void exportPVTK(const char* fname, int numPtn, int numdof);
  void exportVTK(const char* fname, pMesh mesh, pMeshDataId field, int numdof);
  void M_writeVTKFile(pMesh mesh,const char * fname, pMeshDataId field, int numdof);

/////Thickness Adaptation
  pPList V_growthCurveVertices(pVertex vertex);
  int V_isOriginatingNode(pVertex vertex);
  void WallStress(pMesh mesh); 
  void SpaldingLaw(double ydist, double u, double& uTau, double& yplus);
  double DispThickness(pVertex vertex);
  void TotalBLHeight(pMesh mesh);
  void CalcYPlus(pMesh mesh); 
  double BLHeightVort(pVertex vertex, pPList VGrowth);
  double BLHeightVel(pPList VGrowth);
  double BLConstraint(pVertex vert, double distGC);
  void DetectShock(pParMesh pmesh, pMesh mesh);
  void DistOffWall(pMesh mesh, pMSAdapt simAdapter);
  double SeparatedBL(pVertex vertex);
  double GrowthCurveIntersect(pVertex vertex);
  void Vorticity(pParMesh pmesh, pMesh mesh);
  void VortelementGradient(pRegion region, double* elemGradient, int fieldIndexForGrad);
  void VortgradientsFromPatch(pMesh mesh, int fieldIndexForGrad);
 
////Simmetrix development functions
  double E_normalizedLength(pEdge);
  double E_normalizedLengthSq(pEdge);
  void MSA_setPrebalance(pMSAdapt, int);
  void MSA_setCoarsenMode(pMSAdapt, int);
  void MSA_setBLSnapping(pMSAdapt, int);
  void MSA_setBLMinLayerAspectRatio(pMSAdapt, double);
  void MSA_setExposedBLBehavior(pMSAdapt, int);
  void MSA_setMaxIterations(pMSAdapt, int);
  void MSA_setSizeBounds(pMSAdapt, double,double);
  void MSA_setSnapping(pMSAdapt, int);
  pAManager SModel_attManager(pModel model);

/////// calculating volume for pyramids in simmetrix 
double XYZ_volumeSim (double xyz[4][3]);
double R_VolumeSim(pRegion region);

//debug functions
void VolumeDebug(pMesh mesh);
void EdgeLength(pMesh pmesh);
void NormalizedEdgeLength(pMesh pmesh);

#endif
