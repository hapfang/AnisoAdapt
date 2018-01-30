#ifndef H_SPEBC
#define H_SPEBC

// scalar plane extraction boundary condition class

#include "BoundaryCondition.h"

class SPEBC : public BoundaryCondition {
public:
   SPEBC (pGFace );     // Constructor for face
   SPEBC (pGEdge );     // Constructor for edge
   SPEBC (pGVertex );     // Constructor for vertex

  int  isSPEBC(void);
};

#endif
