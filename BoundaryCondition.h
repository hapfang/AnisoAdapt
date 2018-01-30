#ifndef H_BoundaryCondition
#define H_BoundaryCondition

// Boundary condition class

/* #ifndef SIM */
/* #include "modeler.h"      // for gType in essentialBC */
/* #include "AttributeException.h" */
/* using namespace SCOREC_mesh; */
/* using namespace SCOREC_model; */
/* using namespace SCOREC_Util; */
/* using namespace SCOREC_exp; */
/* using namespace SCOREC_att; */
/* #else */
#include "parallel.h" 
#ifdef FMDB
#include "myAttribute.h"
#endif
/* #endif */

#define MLEN 25

#define SQ(x)       ((x)*(x))
#define SMALL       1.0e-5
#define ABS(x)      ((x) < 0 ? -(x) : (x))
#define CLOSE(x,y)  (ABS((x)-(y)) < SMALL ? 1 : 0)
#define false       0
#define true        1

/* enum { false, true }; */

class BoundaryCondition {
public:
  BoundaryCondition ();
  ~BoundaryCondition ();
  int isSet ();
  int isAttSet (int i);

protected:
//#ifdef SIM
  pAttribute *AttList;
/* #else */
/*   SCOREC_att::Attribute **AttList; */
//#endif 
  int set;
  int dontinherit;            // tells you what (not) to inherit - not for face
};

#endif
