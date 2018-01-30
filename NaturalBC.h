#ifndef H_NaturalBC
#define H_NaturalBC

// Natural Boundary condition class

#include "BoundaryCondition.h"
/* #ifndef SIM */
/* #include "modeler.h"      // for gType in essentialBC */
/* #else */
#ifdef SIM
#include "MeshSim.h"
#endif

enum { MF, NP, TV, HF, F1, F2, F3, F4, SID, TW };  // for readability

class NaturalBC : public BoundaryCondition {
public:
  NaturalBC (pGFace);
  int eval (pFace , double *, int *,int nflx);
};

#endif
