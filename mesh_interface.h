#ifndef _MESH_INTERFACE_H_
#define _MESH_INTERFACE_H_

#include <stdlib.h>
#ifdef SIM
#include "SimPartitionedMesh.h" 
#include "MeshSim.h"
#include "MeshSimInternal.h"
#include "Eigen.h"
typedef class MVertex * pVertex;
typedef class MEdge * pEdge;
typedef class MFace * pFace;
typedef class MRegion * pRegion;
#include <vector>
#endif

#ifdef FMDB
// using SCOREC tools
//#include "mTSTT.h"
#include "ParUtil.h"
#include "AOMD.h"
#include "AOMD_cint.h"
#include "MeshSimAdapt.h"
#include "MeshSimInternal.h" 
#include <math.h>
// enum sthreadType {sthreadNone = 0, sthreadDefault = 1};
#ifdef __cplusplus
// for remote copies
//#include "AOMD_cc.h"
#include <vector>
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern int numParts;
extern int numTotParts;
    struct EntCopy {
        pEntity rCopy;
        int rank;
    };
    typedef struct EntCopy EntCopy;
    typedef struct EntCopy *pEntCopy;

    typedef pMesh pParMesh;
    typedef pMeshDataId pRegionWtId;
    typedef struct pmEnt *pMatch;
    typedef pParMesh pUnstructuredMesh;
    typedef int pProgress;
//typedef struct entOrig *pEntOrig;

#ifdef __cplusplus
  typedef pMeshIterator PBEntProcIter; 
  typedef pMeshIterator pBdryProcIter; 
//  typedef class AOMD::mFullIterator *pBdryProcIter;
//  typedef class AOMD::mFullIterator *PBEntProcIter;
  typedef class PList *pEntCopies;
  typedef class PList *pEntOrig;  
  typedef class SMeshChanges * pMeshChanges;
  typedef class BasePartitioner *pPartitioner;
  typedef class PartitionOpts *pPartitionOpts;
  typedef class AttachDataCommu *pAttachDataCommu;
#else
//  typedef struct BdryProcIter *pBEntProcIter; 
//  typedef struct BdryProcIter *pBdryProcIter;
  typedef pMeshIterator pBdryProcIter;
  typedef pMeshIterator PBEntProcIter;
  typedef struct PList *pEntCopies;
  typedef struct PList *pEntOrig;
  typedef struct SMeshChanges * pMeshChanges;
  typedef struct BasePartitioner *pPartitioner;
  typedef struct PartitionOpts *pPartitionOpts;
  typedef struct AttachDataCommu *pAttachDataCommu;
#endif

    inline int SimMeshing_start(void){   
      if(P_pid()==0)
        printf("\nSimMeshing_start(): Not yet implemented\n");
        return 0;
        // exit(-1);
    }
    inline int SimMeshing_stop(void){
      if(P_pid()==0)
        printf("\nSimMeshing_stop(): Not yet implemented\n");
        return 0;
        // exit(-1);
    }
    inline int SimModel_start(void){
      if(P_pid()==0)
        printf("SimModel start(): Not yet implemented\n");
        return 0;
    }
    inline int SimModel_stop(void){return 0;}
    inline int Sim_readLicenseFile(const char *filename){return 1;}
    inline void Sim_logOff(void){};
    inline void SimPartitionedMesh_start(int *pargc, char ***pargv) {
        AOMD::ParUtil::Instance()->init(*pargc,*pargv);
    }
    inline void SimPartitionedMesh_stop() {
        AOMD::ParUtil::Instance()->Finalize(); 
    }
    inline const char* SimMeshing_buildID(void){
        printf("\nSimMeshing_buildID(): Not yet implemented\n");
        return("undergoing");
    }
    inline const char *SimPartitionedMesh_buildID(){
      if(P_pid()==0)
        printf("\nSimPartitionedMesh_buildID(): Not yet implemented\n");
        return("undergoing");
    }
    inline const char* SimMeshTools_buildID(){
      if(P_pid()==0)
        printf("\nSimMeshTools_buildID(): Not yet implemented\n");
        return("undergoing");
    }
    inline const char* SimAdvMeshing_buildID(){
      if(P_pid()==0)
        printf("\nSimAdvMeshing_buildID(): Not yet implemented\n");
        return("undergoing");
    }
    inline int MS_registerKey(const char *key){
      if(P_pid()==0)
        printf("\nMS_registerKey(): Not yet implemented\n");
        return 0;
        // exit(-1);
    }
    inline void EN_deleteDataDbl(pEntity ent, pMeshDataId id) {EN_deleteData(ent,id);}

    inline int Match_gid(pMatch pmatch) {
      if(P_pid()==0)
        printf("\nMatch_gid() : NOT yet implemented\n");
        exit(-1);
    }
    inline pEntity Match_ent(pMatch pmatch){
      if(P_pid()==0)
        printf("\nMatch_ent(): Not yet implemented\n");
        return (pEntity)0;
    }
    inline pMatch PList_matchItem(pPList pl, int i){
      if(P_pid()==0)
        printf("\nPList_matchItem(): Not yet implemented\n");
        return (pMatch)0;
    }

/* lets start from here */
    inline int M_write(pMesh mesh, const char *filename, int version)
    {
        M_writeSMS(mesh,filename,2);
        return 1;
    }
    inline void M_release(pUnstructuredMesh pm){M_delete(pm);}
    inline int PM_verify(pParMesh pm, int choice, int thropt, pProgress progess){
        // if(M_verify(pm))
      if(P_pid()==0)
        printf("\nM_verify is not implemented\n");
            return 1;
            // else
            //return 0;
    }
#ifdef FMDB
    inline pMesh PM_mesh(pParMesh pm, int p){return pm;}
#endif
    inline pGModel PM_model(pParMesh pm){return M_model(pm);}
    inline int PM_numParts(pParMesh pm) {return numParts;}
    inline int PM_totalNumParts(pParMesh pm) {return numTotParts;}
    inline int PMU_rank() {return P_pid();}
    inline int PMU_size() {return P_size();}
    inline int PMU_gid(int proc, int p) {return P_pid()*numParts+p;}
    inline int PMU_proc(int gid) {return gid/numParts;}
    inline int PMU_totalPartFromGid(int tnparts, int gid) {return 0;}

    inline pBdryProcIter PM_bdryProcIter(pParMesh pm, int etype) {
        switch(etype) {
        case 0 :
            return(pBdryProcIter)M_vertexIter(pm); 
            break;
        case 1 :
            return (pBdryProcIter)M_edgeIter(pm);
            break;
        case 2 :
            return (pBdryProcIter)M_faceIter(pm);
            break;
        default:      
	  if(P_pid()==0)
            printf("\nPM_bdryProcIter() : etype=%d NOT supported\n",etype);
            exit(-1);
        }
    }

    inline pEntity BdryProcIter_next(pBdryProcIter iter) {
        pEntity ent = 0;
        // use RIter_next to get the next ent
        // even if the iter is for vertex etc. we use RIter_next()
        // see AOMD.cc for more details
        while(1) {
            ent = RIter_next(iter);
            if(ent) {
                if(EN_duplicate(ent))
                    break;
            }
            else
                return 0;
        }
        return ent;
    }

    // even if the iter is for vertex etc. we use RIter_reset()
    // see AOMD.cc for more details
    inline void BdryProcIter_reset(pBdryProcIter iter) {RIter_reset(iter);}
    // even if the iter is for vertex etc. we use RIter_delete()
    // see AOMD.cc for more details
    inline void BdryProcIter_delete(pBdryProcIter iter) {RIter_delete(iter);}
    inline int EN_ownerProc(pEntity ent){return EN_owner(ent);}
    inline int EN_isOwnerProc(pEntity ent){
        int proc1, proc2;
        proc1 = P_pid();
        proc2 = EN_owner(ent);
        if(proc1==proc2) 
            return 1;
        else 
            return 0;
    }
    inline int EN_isOwnerPart(pEntity ent, int p){
        if(!EN_duplicate(ent))return 1;
        if(EN_owner(ent)==P_pid()*numParts+p) return 1;
        return 0;
    }
    inline int EN_isOnPartBdry(pEntity ent) { return EN_duplicate(ent);}

    inline int PM_numBdryEnt(pParMesh pm, int p, int etype, int proc) {
        int count = 0;
        pBdryProcIter pbIter = PM_bdryProcIter(pm,etype);
        pEntity ent;
        while(ent = BdryProcIter_next(pbIter)) {
            if(!EN_isOwnerPart(ent, p))
                continue;
            if(EN_getCopy(ent,proc))
                count++;
        }
        BdryProcIter_delete(pbIter);

        return count;
    }

#ifdef __cplusplus
    inline pEntCopies EN_copies(pEntity ent){
        pEntCopies ec = PList_new();
        std::vector<std::pair<pEntity,int> > remoteCopies;
        std::vector<std::pair<pEntity,int> >::iterator vecIter;
        EN_getCopies(ent,remoteCopies);
        for(vecIter=remoteCopies.begin(); vecIter!=remoteCopies.end(); ++vecIter) {
            pEntCopy entCopy = new EntCopy;
            entCopy->rCopy = vecIter->first;
            entCopy->rank = vecIter->second;
            PList_append(ec,entCopy);
        }
        return ec;
    }
        
    inline int EntCopies_size(pEntCopies ec){return PList_size(ec);}
    inline pEntity EntCopies_ent(pEntCopies ec, int n){return ((pEntCopy)PList_item(ec,n))->rCopy;}
    inline int EntCopies_gid(pEntCopies ec, int n){return ((pEntCopy)PList_item(ec,n))->rank;}
    inline void EntCopies_delete(pEntCopies ec) {
        void *iter = 0;
        pEntCopy entCopy;
        while(entCopy = (pEntCopy)PList_next(ec,&iter)) {
            delete entCopy;
        }
        PList_delete(ec);
    }
    
    inline pEntOrig EN_original(pEntity ent){
        void *iter =0;
        pEntOrig eo = PList_new();
        pEntCopy entOrig = new EntCopy;
        entOrig->rank = EN_owner(ent);
        pEntCopies ec = EN_copies(ent);
        pEntCopy entCopy;
        while(entCopy = (pEntCopy)PList_next(ec,&iter)){
            if(entCopy->rank == entOrig->rank){
                entOrig->rCopy = entCopy->rCopy;
                break;
            }
            else 
                continue;
        }
        EntCopies_delete(ec);
        PList_append(eo,entOrig);
        return eo;
    }
    inline void EntOrig_delete(pEntOrig eo){
        void *iter = 0;
        pEntCopy entOrig;
        while(entOrig = (pEntCopy)PList_next(eo,&iter)) {
            delete entOrig;
        }
        PList_delete(eo);
    }
    inline pEntity EntOrig_ent(pEntOrig eo){
        return ((pEntCopy)PList_item(eo,0))->rCopy;
    }
#endif

/* lets first implement above ones */ 

typedef enum MeshModType {
  MeshModType_TRef=0, MeshModType_VMov, MeshModType_ESplit, MeshModType_EColapse,
  MeshModType_ESwap, MeshModType_ESplitClps, MeshModType_FSwap, MeshModType_RColapse,
  MeshModType_CavRetri, MeshModType_General
} MeshModType;

    inline void MSA_setPartWts(pMSAdapt simAdapter, const double *partwts){
      if(P_pid()==0)
        printf("\nMSA_setPartWts(): Not yet implemented\n");
    }
/*     inline pPList MCO_getCreatedEntities(pMeshChanges mco){ */
/*         printf("\nMCO_getCreatedEntities(): Not yet implemented\n"); */
/*         return (pPList)0; */
/*     } */
/*     inline pPList MCO_getDeletedEntities(pMeshChanges mco){ */
/*      printf("\nMCO_getDeletedEntities(): Not yet implemented\n"); */
/*      return (pPList)0; */
/*  } */

    inline pPartitionOpts PartitionOpts_new(){
      if(P_pid()==0)
        printf("\nPartitionOpts_new(): Not yet implemented\n");
        return (pPartitionOpts)0;
    }
    inline void PartitionOpts_delete(pPartitionOpts popts){      
      if(P_pid()==0)
     printf("\nPartitionOpts(): Not yet implemented\n");
 }
    inline void PartitionOpts_setTotalNumParts(pPartitionOpts popts, int dtnump){
      if(P_pid()==0)
        printf("\nPartionOpts_setTotNparts(): Not yet implemented\n");
    }
    inline void PartitionOpts_setProcWt(pPartitionOpts popts, const double *wt){
      if(P_pid()==0)
        printf("\nPartionOpts_setProcWt(): Not yet implemented\n");
    }
    inline void PartitionOpts_setPartWtEqual(pPartitionOpts popts){
      if(P_pid()==0)
        printf("\nPartionOpts_setPartWtEqual(): Not yet implemented\n");
    }
    inline void PartitionOpts_setEntityWeightId(pPartitionOpts popts, int i, pRegionWtId id){
      if(P_pid()==0)
        printf("\nPartionOpts_setsetEntityWeightId(): Not yet implemented\n");
    }
    inline void PartitionOpts_setNumEntityWeightIds(pPartitionOpts popts, int nrwtid){
      if(P_pid()==0)
        printf("\nPartionOpts_setNumEntityWeightIds(): Not yet implemented\n");
    }
    inline void PM_partition(pParMesh pm, pPartitionOpts popts, int thropt, pProgress prog){
      if(P_pid()==0)
        printf("\nPM_partition(): Not yet implemented\n");
    }
    
    inline pRegionWtId EntityWeightId_new(const char *tag){
      if(P_pid()==0)
        printf("\nEntityWeightId_new(): Not yet implemented\n");
        return (pRegionWtId)0;
    }
    inline void EN_setWeight(pRegion r, pRegionWtId id, int wgt){
      if(P_pid()==0)
        printf("\nEN_setWeight(): Not yet implemented\n");
    }
    inline void EntityWeightId_delete(pRegionWtId id){
      if(P_pid()==0)
        printf("\nEntityWeightId_delete(): Not yet implemented\n");
    }
    inline void PM_setMigrId(pParMesh pm, pMeshDataId id){
      if(P_pid()==0)
        printf("\nPM_setMigrId(): Not yet implemented\n");
    }
    inline void PM_removeMigrId(pParMesh pm, pMeshDataId id){
      if(P_pid()==0)
        printf("\nPM_removeMigrId(): Not yet implemented\n");
    }
    
    inline pAttachDataCommu AttachDataCommu_new(int mult, int var,
                                                  int maxper){
      if(P_pid()==0)
        printf("\nAttachDataCommu_new(): Not yet implemented\n");
      return (pAttachDataCommu)0;
    }
    inline void AttachDataCommu_delete(pAttachDataCommu adc){
      if(P_pid()==0)
        printf("\nAttachDataCommu_delete(): Not yet implemented\n");
    }
    
    inline int pm_sendAnInt(pAttachableData ad, pAttachDataId id, int cb, void **a, void *b){
      if(P_pid()==0)
        printf("\npm_sendAnInt(): Not yet implemented\n");
        return 0;
    }
    inline int pm_recvAnInt(pAttachableData ad, pAttachDataId id, int cb, void **a, void *b){
      if(P_pid()==0)
        printf("\npm_recvAnInt(): Not yet implemented\n");
        return 0;
    }
    inline int pm_sendDblArray(pAttachableData ad, pAttachDataId id, int cb, void **a, void *b){
      if(P_pid()==0)
        printf("\npm_sendDblArray(): Not yet implemented\n");
        return 0;
    }
    inline int pm_recvDblArray(pAttachableData ad, pAttachDataId id, int cb, void **a, void *b){
      if(P_pid()==0)
        printf("\npm_recvDblArray(): Not yet implemented\n");
        return 0;
    }
    
    inline pVertex F_oppositeVertex(pFace face, pEdge edge){
        return F_edOpVt(face,edge);
    }
    inline pPList EN_matches(pEntity ent, pGEntity filter){
      if(P_pid()==0)
        printf("\nEN_matches() : Not yet implemented\n");
        return (pPList)0;
    }

    inline double E_length(pEdge edge){
        return sqrt(E_lengthSq(edge));
    }
    inline void GM_release(pGModel model){
        GM_delete(model);
    }
    inline void MD_setMeshCallback(pMeshDataId id, int event, CBfunc f, void *cbdata){
      if(P_pid()==0)
        printf("MD_setMeshCallback() not implemented yet\n");
    }
    inline void MD_removeMeshCallback(pMeshDataId id, int event){
      if(P_pid()==0)
        printf("MD_removeMeshCallback() not implemented yet\n");
    }
#ifdef __cplusplus
}
#endif
#endif
//end of SCOREC tools

#ifdef SIM
#ifdef __cplusplus
extern "C" {
#endif
inline int F_conToFace(pFace face1, pFace face2) {
  pPList fverts, fedges;
  pEdge edge;
  pVertex vtx;
  void *temp;

  temp = 0;
  fverts = F_vertices(face1,1);
  while (vtx = (pVertex)PList_next(fverts,&temp))
    if (F_inClosure(face2,vtx)) {
      PList_delete(fverts);
      return 1;
    }
  PList_delete(fverts);

  temp = 0;
  fedges = F_edges(face1,0,0);
  while (edge = (pEdge)PList_next(fedges,&temp))
    if (F_inClosure(face2,edge)) {
      PList_delete(fedges);
      return 1;
    }
  PList_delete(fedges);

  return 0;
}

#ifdef __cplusplus
   }
#endif

#endif //endif SIM

#endif
