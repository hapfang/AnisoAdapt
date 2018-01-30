#include "BoundaryCondition.h"
#include "string.h"

//  extern int ensa_dof;
// bring in the # of variables to help find # of scalars

// the following codes are now defined and used in specific boundary
// condition classes ( C.W. )

/* D = density, P = pressure, T = temperature,
 * C1 = comp1: group of magnitude and direction, when one velocity is set,
 * C3 = comp3: group of magnitude and direction, when three velocities are set,
 * MF = mass flux, NP = natural pressure, TV = traction vector, HF = heat flux
 * IM = integrated mass
 * PS = periodic slave, NOJ = number of jumps,
 */
/* strAttE contains strings for the essential BC attributes
 * strAttN contains strings for the natural BC attributes
 * strAttP contains strings for the periodic BC attributes
 * strAttS contains strings for the SPEBC attributes
 *
 * Note: The earlier form of AttList etc. is altered. The first enum is
 *       for essential bc's gui stuff, second for natural bc stuff.
 *
 * Density, pressure and temperature are the thermodynamic quantities.
 * comp1, comp3 are not exactly attributes but velocity attribute groups.
 * Both contain a vector direction and magnitude. Both cant be set at the
 *   same time. comp1 implies velocity is set only in 1 direction and is
 *   free in the other directions. comp3 means velocities are set in all
 *   the three directions and its magnitude and direction are those of the
 *   resultant.
 * Periodic slave and number of jumps is the usual stuff. So also is the
 *   natural BC stuff, which can be only set on a face.
 * Essential BC's need to be user-set over faces and optionally over edges. If
 *   any of the thermodynamic quantities are set on an edge, we dont inherit
 *   any thermodynamic quantity from the associated faces. If any velocity
 *   component is set, we dont inherit any of the velocity components from
 *   the faces. Similarly for the periodic stuff. Everything as above for
 *   vertex wrt associated edges.
 */

BoundaryCondition::BoundaryCondition()     // top level constructor
{
  AttList = NULL;
  set = false;                             // "bad" entity assumed
}

BoundaryCondition::~BoundaryCondition ()
{
  if (AttList)  delete[] AttList;
}

int BoundaryCondition::isSet()             // public : BC's set or not set ?
{
  return set;
}

int BoundaryCondition::isAttSet (int i)    // protected : ith attribute set ?
{
  // if AttList exists at all
  // AttList  
  if (AttList)  {
    // does AttList[i] exist, if yes return true else false
    return (AttList[i]) ? true : false;    // doesnt check for validity of i
  }
  // return false if AttList DOES NOT exist
  return false;
}
