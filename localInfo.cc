/* attach information local to each entity */

#include <stdio.h>
#include <map>
#include "func.h"
#include "parallel.h"
#include "bits.h"
#ifdef SIM
#include "SimModel.h"
#endif
#include "MeshSimInternal.h"
#include "phParAdapt.h"
using namespace std;

/*********************************************************************
* allocate and initialize a local information array to each entity   *
* in the mesh                                                        *
**********************************************************************/
//#ifdef __cplusplus
//extern "C" {
//#endif

// local information consists of:
// globalP for each entity
// the mesh is local to the calling proc

extern int globalP;

pMeshDataId MYCTID;
pMeshDataId NDOFID; 
pMeshDataId POLYID;
pMeshDataId RNENID;
pMeshDataId ibcid;
pMeshDataId switid;
extern int isReorder;
extern int numParts;
extern int numTotParts;
// only in parallel
// replaces previous MeshSim functions
// declared in SimPMesh.h at least up to
// release 5.3-040830
void PMU_commuArr(void **s, int *ns, void **r, int *nr, MPI_Datatype type,
                  int mult)
{

    int i, m, tag = 0;

      MPI_Request* req = new MPI_Request[2*(PMU_size()-1)]; 
      MPI_Status*  stat= new MPI_Status[2*(PMU_size()-1)]; 
    
    for (m=0, i=0; i < PMU_size(); i++)  {
        if (i != PMU_rank())  {
          if (nr[i] > 0)
            MPI_Irecv(r[i], mult*nr[i], type, i, tag, MPI_COMM_WORLD, &req[m++]);
          if (ns[i] > 0)
            MPI_Isend(s[i], mult*ns[i], type, i, tag, MPI_COMM_WORLD, &req[m++]);
        }
    }
    MPI_Waitall(m, req, stat);

    delete [] req; 
    delete [] stat; 

}

void PMU_commuArr(void ***s, int **ns, void ***r, int **nr, map<int,int> *nsmap, 
                  map<int,int> *nrmap, MPI_Datatype type, int mult)
{
    int i, j, ipart, m, tag = 0, RequstSize =0;
    for(m=0,ipart=0;ipart<numParts;ipart++) 
        RequstSize += nsmap[ipart].size() + nrmap[ipart].size();

    MPI_Request* req = new MPI_Request[RequstSize]; 
    MPI_Status*  stat= new MPI_Status[RequstSize]; 
    
    std::map<int,int>::iterator nsiter, nriter;
    
    for(m=0,ipart=0;ipart<numParts;ipart++) {   
        for (nsiter=nsmap[ipart].begin();nsiter!=nsmap[ipart].end();nsiter++)  {
            i = nsiter -> first;
            j = nsiter -> second;
            if(ns[ipart][j]>0) {
                int destID = i/numParts;
                tag = i - numParts*destID;
                MPI_Isend(s[ipart][j], mult*ns[ipart][j], type, destID, tag, MPI_COMM_WORLD, &req[m++]);
            }
        }
        
        for (nriter=nrmap[ipart].begin();nriter!=nrmap[ipart].end();nriter++)  {        
            i = nriter -> first;
            j = nriter -> second;
            if (nr[ipart][j] > 0)
                MPI_Irecv(r[ipart][j], mult*nr[ipart][j], type, i/numParts, ipart, MPI_COMM_WORLD, &req[m++]);
        }
    }

    MPI_Waitall(m, req, stat);
    
    delete [] req; 
    delete [] stat;      
}

//keep this function for now, just in case it is called by other functions. Take
//care of this later. Min Zhou
void PMU_commuInt(int *ns, int *nr)
{

      MPI_Request* req = new MPI_Request[2*(PMU_size()-1)];
      MPI_Status*  stat= new MPI_Status[2*(PMU_size()-1)]; 
    int i, m, tag = 0;
    for (m=0, i=0; i < PMU_size(); i++)  {
        if (i != PMU_rank())  {
            MPI_Irecv(nr+i, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &req[m++]);
            MPI_Isend(ns+i, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &req[m++]);
        }
    }
    MPI_Waitall(m, req, stat);
    
      delete [] req; 
      delete [] stat; 
}

void PMU_commuInt(int **ns, int **nr, map<int,int> *nsmap, map<int,int> *nrmap)
{
    int i, m, j, ipart, tag = 0;
    int RequstSize = 0;
    std::map<int,int>::iterator nsiter, nriter;
    map<int,int>::iterator iter1, iter2;

    for(ipart=0; ipart<numParts; ipart++) 
       RequstSize += nsmap[ipart].size() + nrmap[ipart].size();

    MPI_Request* req = new MPI_Request[RequstSize];
    MPI_Status*  stat= new MPI_Status[RequstSize]; 

    for(m=0,ipart=0;ipart<numParts;ipart++) {
      for (nsiter=nsmap[ipart].begin();nsiter!=nsmap[ipart].end();nsiter++)  {
        i = nsiter->first;
        j = nsiter->second;  
        int destID = i/numParts;           
        tag = i - destID*numParts; 
        MPI_Isend(ns[ipart]+j, 1, MPI_INT, destID, tag, MPI_COMM_WORLD, &req[m++]);
      } 
      for (nriter=nrmap[ipart].begin();nriter!=nrmap[ipart].end();nriter++)  {
       i = nriter->first;
       j = nriter->second;  
       tag = ipart;
       int fromID = i/numParts;
       MPI_Irecv(nr[ipart]+j, 1, MPI_INT, fromID, tag, MPI_COMM_WORLD, &req[m++]);
      }
    }
    
    MPI_Waitall(m, req, stat);    
    delete [] req; 
    delete [] stat; 
}

int periocheck(pEntity ent, globalInfo *info)
{
  pMatch ment;

  // isPeriodic declared in setPeriodic.cc
  if (isPeriodic(ent, &ment))  {
    if (isOnThisProc(ment)) {
// ent is slave with master is on same proc and owner by proc
      if (EN_isOwnerProc(ent)) {
        return 1;
      } else {
// ent is slave with master is on same proc but owner by another proc
        return 3;
      }
    }
// ent is slave with master on another proc
    PList_append(info->perSlv, ent);
    PList_append(info->perMst, ment);
    return 0;
  }
// ent is not a periodic slave
  return 2;
}



//////////////////////////////////////////////////////////////////////////////
// called in mdb2phasta
//////////////////////////////////////////////////////////////////////////////
#ifdef PARALLEL
void initLocalInfo(std::vector<pParMesh> pmesh, globalInfo **info)
#else
void initLocalInfo(pMesh mesh, globalInfo* ginfo)
#endif
{

  pVertex vertex;
  pEdge   edge;
  pFace   face;
  pRegion region;
  pEntity ent;
  VIter vIter;
  EIter eIter;
  FIter fIter;
  RIter rIter;

  int Biggest,i,j,tbit;
  int localP, nen;
  int p; /* typing saver for global/localP */
  int ndof;
  int mybig;
  void *tmp;
  int ipart;

  int gid, pid, count, pc;
      
  MYCTID =  MD_newMeshDataId("MYCT"); 
  NDOFID =  MD_newMeshDataId("NDOF"); 
  POLYID =  MD_newMeshDataId("POLY"); 
  RNENID =  MD_newMeshDataId("RNEN"); 
  ibcid  =  MD_newMeshDataId("ibc ");
  switid =  MD_newMeshDataId("swit");

  for(ipart=0; ipart<numParts; ipart++) {

   globalInfo *ginfo= info[ipart];
   pMesh mesh;
#ifdef FMDB
   mesh = pmesh[ipart];
#else
   mesh = PM_mesh(pmesh[0], ipart);
#endif
   count = 0;

  if(isReorder){
      ginfo->nshg = M_numVertices(mesh);
//      ginfo->nshgOwn = ginfo->nshg-numVerticesNotOwn;
      V_reordering(mesh, ginfo, ipart);
  }

  /* initializing the edge /face/ region mode info */
  ginfo->edgeson =0;
  ginfo->faceson =0;
  ginfo->regnson =0;
  ginfo->nshg = 0;
  ginfo->nshgOwn = 0;

  
  /*******************************************************************/
  /* allocate and attach a local information structure to each       */
  /* entity in the mesh                                              */
//
// what exactly is done here?
  /*******************************************************************/
  /* vertices */
   vIter = M_vertexIter(mesh);
   while (vertex = VIter_next(vIter))  {
    /* setting the number of dofs based on localP */
    localP  = globalP;
    EN_attachDataInt(vertex,POLYID,localP);
    EN_attachDataInt(vertex,NDOFID,1);

    if(!isReorder)
        EN_attachDataInt(vertex,MYCTID,count++);
    ginfo->nshg++;
    
  }//while (vertex

  VIter_delete(vIter);

  /* edges */
  eIter = M_edgeIter(mesh);
  while (edge = EIter_next(eIter))  {

    localP  = globalP;  /* P-change */
    ndof    = (localP) - 1;

    if (localP > 1) ginfo->edgeson=1;// edgemodes only for globalP=2 and higher

    EN_attachDataInt(edge,POLYID,localP);

    if (ndof) {

      EN_attachDataInt(edge,NDOFID,ndof);
      EN_attachDataInt(edge,MYCTID,count);
      count += ndof;
      ginfo->nshg += ndof;

    }//ndof!=0
  }
  EIter_delete(eIter);

  /* faces */
  fIter = M_faceIter(mesh);
  while (face = FIter_next(fIter))  {
    localP  = globalP; /* P-change */
    switch(F_numEdges(face)){
    case 3:              /* Tri Face */
      if (localP > 2) {
          ginfo->faceson = 1; // facemodes only for globalP=3 and higher
        ndof = ((localP - 1)*(localP - 2))/2;
      } else {
        ndof = 0;
      }
      break;
    case 4:            /* Quad Face */
      if (localP > 3) {
        ginfo->faceson = 1;
        ndof = ((localP - 3)*(localP - 2))/2;
      } else {
        ndof = 0;
      }
      break;
    default:
      fprintf(stderr,"Face with neither 3 nor 4 edges \n");
      exit(1);
    }
    EN_attachDataInt(face,POLYID,localP);
    if (ndof) {
      EN_attachDataInt(face,NDOFID,ndof);
      EN_attachDataInt(face,MYCTID,count); 
      count += ndof;
    }
  }
  FIter_delete(fIter);

  /* regions */
  Biggest = 0;
  rIter = M_regionIter(mesh);
  while (region =RIter_next(rIter))  {
    localP  = globalP;  /* P-change, has to be replaced with a
                         function which returns the maximum polynomial
                          order of the all the enclosed entities */
    switch(topology(region)){
    case 1:    /* tets */
      nen  = 4;
      Biggest = setbit(Biggest, 0);
      if ((p = localP) > 3)// region modes  only for globalP=4 and higher
        ndof = (p-1)*(p-2)*(p-3)/6;
      else ndof = 0;
      break;
    case 5:    /* Pyramids */
      nen  = 5;
      Biggest = setbit(Biggest, 1);
      if ((p = localP) > 5)
        ndof = (p-3)*(p-4)*(p-5)/6;
      else ndof = 0;
      break;
    case 3:   /* Wedges */
      nen = 6;
      Biggest = setbit(Biggest, 2);
      if ((p = localP) > 4)
        ndof = (p-2)*(p-3)*(p-4)/6;
      else ndof = 0;
      break;
    case 2:   /* Hexes */
      nen = 8;
      Biggest = setbit(Biggest, 3);
      if ((p = localP) > 5)
        ndof =(p-3)*(p-4)*(p-5)/6;
      else ndof = 0;
      break;
    default:
      fprintf(stderr,"Unknown Element type in LocalInfo\n");
      exit(-1);
    }
    EN_attachDataInt(region,RNENID,nen);
    if (ndof) ginfo->regnson = 1;
    EN_attachDataInt(region,POLYID,localP);
    if (ndof) {
      EN_attachDataInt(region,NDOFID,ndof);
      EN_attachDataInt(region,MYCTID,count); 
      count += ndof;
      ginfo->nshg += ndof;
    }
  }
  RIter_delete(rIter);
  
//suppose there is only one topology, This is not going to be right if there are
//mixed topologies. Need to come back to this later. Min Zhou
#ifdef PARALLEL
  mybig = Biggest;
  MPI_Allreduce(&mybig, &Biggest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  tbit = 0;
  i = 3;
  while (!(tbit = getbit(Biggest,i--)) && (i >= 0));
  if (!tbit) {
    fprintf(stderr," Bit test for Topology has problems \n");
    fprintf(stderr," Please stop in func initlocalinfo for debugging \n");
    exit(1);
  } else { tbit = ++i ;}

  switch(tbit){
  case 0:    /* tets */
    ginfo->nen  = 4;
    ginfo->nenb = 3;
    ginfo->nedges = 6;
    ginfo->nfaces = 4;
    break;
  case 1:    /* Pyramids */
    ginfo->nen  = 5;
    ginfo->nenb = 4; /* should not matter, we never let Pyramids
                        get on the boundaries */
    ginfo->nedges = 8;
    ginfo->nfaces = 5;
    break;
  case 2:   /* Wedges */
    ginfo->nen = 6;
    ginfo->nenb = 4;  /* wedge can also have Tri Boundary face,
                         keep this in mind */
    ginfo->nedges = 9;
    ginfo->nfaces = 5;
    break;
  case 3:   /* Hexes */
    ginfo->nen = 8;
    ginfo->nenb = 4;
    ginfo->nedges = 12;
    ginfo->nfaces = 6;
    break;
  default:
    fprintf(stderr,"Unknown Element type in LocalInfo\n");
    exit(-1);
  }
  }
}

//#ifdef __cplusplus
//}
//#endif
