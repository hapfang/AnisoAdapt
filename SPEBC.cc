#include <iostream>
#include <math.h>
#include "parallel.h"
#include "SPEBC.h"
#include "bits.h"

extern int ensa_dof;  // bring in the # of variables to help find # of scalars

#define numAttsS 1        // number of SPEBC attributes
char strAttS[numAttsS][MLEN] ={"spebc"};
enum { SP };

/////////////////////////////////////////////////////////////////////////////////////////
// cstr.
// construct a SPEBC by extracting whether it is set on the given face or not
/////////////////////////////////////////////////////////////////////////////////////////
SPEBC::SPEBC(pGFace gface) : BoundaryCondition ()
{
//#ifdef SIM
  if (GEN_attrib((pGEntity)gface, strAttS[SP]))  {    // if slave then everything set
//  #else
//      std::vector<Attribute*> atts;
//      try {
//         atts = SCOREC_att::retrieveAttributePList((pGEntity)gface, strAttS[SP]);
//      } catch (AttributeNotExistent) {}
//      if ( atts.size() > 0 ) {
//#endif
    this->set = true;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////
// cstr.
// construct a SPEBC for an edge by setting SPEBCs
// to the edge's faces
/////////////////////////////////////////////////////////////////////////////////////////
SPEBC::SPEBC(pGEdge gedge) : BoundaryCondition ()
{

//#ifdef SIM
  if (GEN_attrib((pGEntity)gedge, strAttS[SP]))  {     // if user sets, no inheritance
//  #else
//      std::vector<Attribute*> atts;
//      try {
//        atts = SCOREC_att::retrieveAttributePList((pGEntity)gedge, strAttS[SP]);
//      } catch (AttributeNotExistent) {}
//      if ( atts.size() > 0 ) {
//#endif
    this->set = true;
  }  
  else  {

    pPList gefaces = GE_faces(gedge);// the faces of the edge
    void* restart = 0;
    pGFace gfi;

    while (gfi = (pGFace)PList_next(gefaces,&restart))  {       // loop over faces
      SPEBC sbc(gfi);            // set bc on the face
      if (sbc.isSet())  {        // should be slave
        this->set = true;
        break;                   // break from loop => edge is slave if
      }                          //   at least one connecting face is slave
    }
    PList_delete(gefaces);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////
// cstr.
// construct a SPEBC for a vertex  by setting SPEBCs
// to the vertex's faces
/////////////////////////////////////////////////////////////////////////////////////////
SPEBC::SPEBC(pGVertex gvert) : BoundaryCondition ()
{
  /* Exactly similar to edge above so look there for comments */

//#ifdef SIM
  if (GEN_attrib((pGEntity)gvert, strAttS[SP]))  {
//  #else
//      std::vector<Attribute*> atts;
//      try {
//        atts = SCOREC_att::retrieveAttributePList((pGEntity)gvert, strAttS[SP]);
//      } catch (AttributeNotExistent) {}
//      if ( atts.size() > 0 ) {
//#endif

    this->set = true;
  }  else  {

    pPList gvfaces = GV_faces(gvert);
    void* restart= 0;
    pGFace gfi;

    while (gfi = (pGFace)PList_next(gvfaces,&restart))  {
      SPEBC sbc(gfi);
      if (sbc.isSet())  {
        this->set = true;
        break;
      }
    }
    PList_delete(gvfaces);
  }
}

int SPEBC::isSPEBC(void)
{
  if(isSet()) return 1;
  else return 0;
}
