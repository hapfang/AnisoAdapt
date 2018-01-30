// Comments irrelevant.. SST
// at least: where it is called:(JM)
// called in mdb2phasta.cc
// unclear (JM 7-30-04)

#include <fstream>
#include <math.h>
#include "parallel.h"
#include "func.h"
#include "ccfunc.h"
#include "PeriodicBC.h"
#include "NaturalBC.h"
#include "SPEBC.h"
#include "MeshSimInternal.h"
#include "phParAdapt.h"
#include "mesh_interface.h"
#include <iostream>
using namespace std;
#define ATOL 0.01

// print boundary condition codes for geometric model
extern int prCd;
extern map<pair<string,pGEntity>,int> GEntityDataContainerInt;
extern map<pair<string,pGEntity>,void *> GEntityDataContainerPtr;
string intr("intr");
string teta("teta");
string MAST("MAST");
string AXCL("AXCL");
string sper("sper");
string iper("iper");
string NoBI("NoBI");
string bc("bc  "); 

// no attribues involved ...
// no mesh involved ...
// but: case (geom) has
// been associated with model
// previously (in switch_adaptPreproc)
void setPeriodic (pGModel model)
{
  pGEntity gment;
  pGVertex gvert;
  pGEdge gedge;
  pGFace gface;
  int njumps;
  ofstream fout;

  // number of jumps probably no longer necessary .. SST

  if (PMU_rank() == 0)  {
    if (prCd)  {
      fout.open("BC.out");
      fout << "Slave entity" << "\tMaster entity" << "\tNum jumps\n";
      fout << "------------" << "\t-------------" << "\t---------\n";
    }
  }

  // model-edges and model-nodes inherit interior=intr properties from model-face
  GFIter gftr1 = GM_faceIter(model);
  while (gface = GFIter_next(gftr1))  {    // loop over model faces
    pPList gf_regions = GF_regions(gface); // regions adjacent to face
    if (PList_size(gf_regions) == 2) {     // interior model-face

      GEN_attachDataI((pGEntity)gface,"intr",1);

      // no boundary intregral, i.e. no Neumann BC
      GEN_attachDataI((pGEntity)gface,"NoBI",1);

      // inherit onto mode-edges
      pPList gf_edges = GF_edges(gface);
      void* etmp = 0;
      while(pGEdge gedge = (pGEdge)PList_next(gf_edges,&etmp))
        GEN_attachDataI((pGEntity)gedge,"intr",1);
      PList_delete(gf_edges);

      // inherit onto model-nodes
      pPList gf_vertices = GF_vertices(gface);
      void* vtmp = 0;
      while(pGVertex gvertex = (pGVertex)PList_next(gf_vertices,&vtmp))
        GEN_attachDataI((pGEntity)gvertex,"intr",1);
      PList_delete(gf_vertices);
    }
    PList_delete(gf_regions);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////
  // loop over all model entities (faces, edges, and vertices) and
  // if they are a "real" periodic slave(as in not an SPEBC), find
  // the master of all the  mesh entities classified on them.
  ////////////////////////////////////////////////////////////////////////////////////////////
  GFIter_reset(gftr1);
  while (gface = GFIter_next(gftr1)) {// over model faces
    pGEntity gent = (pGEntity) gface;
    int tmp;
    if (!GEN_dataI(gent,"intr",&tmp))  { // not interior face

      // what is that ??? first a SPEBC is created and later
      // it is tested whether it is set at all ???
      // --> test takes place in cstr.  
      SPEBC spbc(gface);

      if (spbc.isSPEBC()) { 
        GEN_attachDataI(gent, "sper", 1);
      }

      PeriodicBC bc(model, gface);
      njumps = bc.getPerMaster ((pGFace*) &gment); // master is passed onto gment

      if (njumps != 0)  {
        GEN_attachDataP(gent, "MAST", gment);

        // check for axisymmetric boundary condition
        // If the slave face makes a nonzero angle with the master
        // then we assume that we have an Axisymmetric BC case

        // Here theta represents the  rotation for the slave face to
        // the master face and hence -theta from the masterface to the
        // slave face.

        // why attach theta if not axisym ?  .. SST

        double theta = bc.getAngle((pGFace) gment, 0);  // WRONG -- SST
        if (fabs(theta) < ATOL) theta = 0.0;
        double *ttmp = new double;
        double *ttmpn = new double;
        *ttmp = theta;
        *ttmpn = -theta;
        GEN_attachDataP(gent,"teta",ttmp);
        GEN_attachDataP(gment,"teta",ttmpn);

        if (prCd && PMU_rank()==0)
          fout << "F-" << GEN_tag(gent) << "\t\tF-" << GEN_tag(gment)
               << "\t\t" << njumps << "\n";

        NaturalBC nbc(gface);
        if(!nbc.isAttSet(SID))                  // for faces with Surf ID
          GEN_attachDataI (gent, "NoBI", 1);    // flag as noBI

        GEN_attachDataI (gment, "NoBI", 1);     // wont count in bdry integral
        GEN_attachDataI(gent, "iper", 1);
      }
    }

    // why attach theta if not even periodic!!? .. SST
    double *tmpdouble;
    if (!GEN_dataP(gent,"teta",(void**)&tmpdouble)){
      double* thetat = new double;
      *thetat = 0.0;
      GEN_attachDataP(gent,"teta",thetat);
    }
  }
  GFIter_delete(gftr1);

  GEIter getr1 = GM_edgeIter(model);
  while (gedge = GEIter_next(getr1)){          // all model edges

    SPEBC spbc(gedge);
    pGEntity gent = (pGEntity) gedge;
    if (spbc.isSPEBC()) {
      GEN_attachDataI(gent, "sper", 1);
    }
    PeriodicBC bc(model, gedge);
    njumps = bc.getPerMaster((pGEdge *) &gment);

    if (njumps == 2475)    // Recognizing the axisymmetric centerline
      GEN_attachDataI(gent,"AXCL", 1);

    else if (njumps != 0)  {
      GEN_attachDataP(gent, "MAST", gment);

      double theta = bc.getAngle((pGEdge) gment);
      if(fabs(theta) < ATOL) theta = 0.0;
      double *ttmp = new double;
      *ttmp = theta;
      GEN_attachDataP(gent,"teta",ttmp);

      if (prCd && PMU_rank()==0)
        fout << "E-" << GEN_tag((GEntity *) gedge) << "\t\tE-"
                     << GEN_tag(gment) << "\t\t" << njumps << "\n";

      GEN_attachDataI(gent, "iper", 1);
    }
    double *tmpdouble;
    if (!GEN_dataP(gent,"teta",(void**)&tmpdouble)) {
      double* thetat = new double;
      *thetat = 0.0;
      GEN_attachDataP(gent,"teta",thetat);
    }
  }
  GEIter_delete(getr1);

  GVIter gvtr1 = GM_vertexIter(model);
  while (gvert = GVIter_next(gvtr1)){       // all model vertices
    pGEntity gent = (pGEntity) gvert;
    SPEBC spbc(gvert);

    if (spbc.isSPEBC())
      GEN_attachDataI(gent, "sper", 1);

    PeriodicBC bc(model, gvert);
    njumps = bc.getPerMaster((pGVertex *) &gment);

    if (njumps == 2475)          // Recognizing the axisymmetric centerline
      GEN_attachDataI(gent,"AXCL", 1);

    else if (njumps != 0)  {              // means the entity is a slave
      GEN_attachDataP(gent, "MAST", gment);
      double theta = bc.getAngle((pGVertex) gment);
      if (fabs(theta) < ATOL) theta = 0.0;
      double *ttmp = new double;
      *ttmp = theta;
      GEN_attachDataP(gent,"teta",ttmp);

      if (prCd && PMU_rank()==0)
        fout << "V-" << GEN_tag((GEntity *) gvert) << "\t\tV-"
                     << GEN_tag(gment) << "\t\t" << njumps << "\n";

      GEN_attachDataI(gent, "iper", 1);
    }
    
    double *tmpdouble;
    if (!GEN_dataP(gent,"teta",(void**)&tmpdouble)) {
      double* thetat = new double;
      *thetat = 0.0;
      GEN_attachDataP(gent,"teta",thetat);
    }
  }
  GVIter_delete(gvtr1);

  if (PMU_rank() == 0)  {
    if(prCd)
      fout.close();
  }
}
// called in localInfo,printBC 
// called in 
// int periocheck(pEntity ent, globalInfo *info)
extern "C" int isPeriodic(pEntity ent, pMatch *ment)
{

  //Returns the type of MODEL entity the MESH entity ent is  classified on
  pGEntity ge = EN_whatIn(ent);


  // MAST assigned in setPeriodic (pGModel model)
  // if a master model ent has been assigned, GET it
  double *tmpdouble;
  if(!GEN_dataP(ge, "MAST",(void**)&tmpdouble))// internal fctn
      return 0;
  pGEntity gem = (pGEntity)tmpdouble;

  // get the MESH ents ml matching for this mesh ent ent of the MASTER model ent gem
  // this is a serial function--returns MESH entities (=pEntity) that are matched to ent
  // filtered on gem
  pPList ml =  EN_matches(ent, gem);
  if (!ml)
    return 0;

  // should be one and only one item in list..
  if (PList_size(ml) != 1)  {
    cerr << "Something wrong in isPeriodic in setPeriodic.cc\n";
    PList_delete(ml);
    return 0;
  }

  // type cast to pMatch -- in DOC 5.4 it is said that PList_matchItem operates
  // on a list obtained from EN_matches(ent, gem) that ARE of type pMatch
  // no type cast would be needed ?!? 

  *ment = PList_matchItem(ml, 0);
  
//   if(EN_type(ent)==0){
//       cout<<"["<<PMU_rank()<<"]"<<"found match:\n";
//       cout<<"["<<PMU_rank()<<"]"<<"orig vertex: "<<EN_id(ent)<<"\n";
//   }

  PList_delete(ml);
  return 1;
}

extern "C" int isOnThisProc(pMatch ment)
{ 

    return PMU_rank() == PMU_proc(Match_gid(ment));
 }
