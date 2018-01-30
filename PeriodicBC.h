#ifndef H_PeriodicBC
#define H_PeriodicBC

#include "BoundaryCondition.h"

/*compatibility*/
#include "ModelTypes.h"
#ifdef SIM
#include "SimAttribute.h"
#endif
#ifdef FMDB
#include "myAttribute.h"
#endif

class PeriodicBC : public BoundaryCondition {
public:
  PeriodicBC (pGModel, pGFace);        // constructor for face
  PeriodicBC (pGModel, pGEdge);        // constructor for edge
  PeriodicBC (pGModel, pGVertex);      // constructor for vert
  int getPerMasterTag ();             // return tag of current's master entity
  int getPerMaster (pGFace *);        // returns num jumps, also entity
  int getPerMaster (pGEdge *);        //   in parentheses is the master
  int getPerMaster (pGVertex *);      // if present entity is vertex then must
                                      //   call getPerMaster(Vertex *)
  double getAngle (pGFace, pMesh);
  double getAngle (pGEdge);
  double getAngle (pGVertex);

  void GF_normal_flat(pGFace, pMesh, double *xyz);

  static int global_sanity; // what for ???

private:
  gType gtype;
  pGFace gf;
  pGEdge ge;
  pGVertex gv;
  pGModel model;
  double myangle; /* so that I don't have to pass stuff around */
  void update_inherit ();
  double getDistance (pGFace, pGFace);
  double getDistance (pGEdge, pGEdge);
  double getDistance (pGVertex, pGVertex);
  double getDistance (double *, double *, double theta=0.0);
};

#endif
