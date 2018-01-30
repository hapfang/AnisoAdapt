#include <iostream>
#include <math.h>
#include "parallel.h"
#include "NaturalBC.h"
#ifndef SIM
#include <vector>
#endif
#include "bits.h"
#include "phParAdapt.h"

extern int ensa_dof;  // bring in the # of variables to help find # of scalars
extern map<pair<string,pGEntity>,int> GEntityDataContainerInt;
extern string NoBI;

#define numAttsN 10        // number of natural bc attributes

char strAttN[numAttsN][MLEN] = { "mass flux",
                                 "natural pressure",
                                 "traction vector",
                                 "heat flux",
                                 "scalar_1 flux",
                                 "scalar_2 flux",
                                 "scalar_3 flux",
                                 "scalar_4 flux",
                                 "surf ID",
                                 "turbulence wall" };

NaturalBC::NaturalBC(pGFace gface) : BoundaryCondition()     // Natural BC's
{
  int i;
  int tmp;
  if (!GEN_dataI((pGEntity)gface,"NoBI",&tmp))  {
   // model face not flagged as not in BI
    this->set = true;                            // set is true

    AttList = new pAttribute[numAttsN];          // update attlist
    for (i=0; i < numAttsN; i++)  {
//#ifdef SIM
      AttList[i] = GEN_attrib((pGEntity)gface, strAttN[i]);
// #else
//       try {
//         std::vector<Attribute*> atts = SCOREC_att::retrieveAttributePList((pGEntity)gface,strAttN[i]);
//         AttList[i] = atts[0];
//         atts.clear();
//       } catch  (AttributeNotExistent) {AttList[i] = 0 ; }
// #endif
    }
  }
}

int NaturalBC::eval(pFace face, double *BCB, int *iBCB,int nflx)
{
  pVertex vertex[4];
  int nenb = F_numEdges(face), i, j, ibcb,k;

  /* right now we are using nenb for interpolation of Natural Boundary
     conditions . This will have to be upgraded to bnshlb when we decide
     to make higher order interpolation of NBCs */

  /* as of 3/13/2001 we have decided to go with piecewise constant natural
     boundary conditions since we never seem to use anything else. we save
     trouble moving to higher order and also avoid a potential bug on linear
     ones when we have varying bcb for quadface wedges and triface pyramids

     so in short.. don't pay much attention to the last comment made
     above the current one but just keep in mind that we have natural
     boundary conditions which are specified as constant on a face.

     by the way this function evaluates them at the centroid. */

  double xyz[3];
  pPList vlist = F_vertices(face, 1);
  double centroid[3]={0.0, 0.0, 0.0};

  for (i=0; i < PList_size(vlist); i++)  {
    vertex[i] = (pVertex) PList_item(vlist, i);
  }
  PList_delete(vlist);

  for(i=0; i < nenb; i++) {
    V_coord(vertex[i], xyz);
    for(j=0; j< 3; j++) centroid[j]+= xyz[j];
  }
  for(i=0; i< 3; i++) centroid[i]/=nenb;

  if (AttList[MF])                    // mass flux
//old
//    BCB[0] = Att_evalTensorOr0(AttList[MF], centroid);
    BCB[0] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[MF], centroid);

  if (AttList[NP])                   // natural pressure
//old
//    BCB[1] = Att_evalTensorOr0(AttList[NP], centroid);
     BCB[1] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[NP], centroid); 
  if (AttList[TV])
//old                     // traction
//    Att_evalTensorOr1(AttList[TV], centroid, BCB+2);
for( int index =0; index < 3; index++)
          BCB[2+index] = AttributeTensor1_evalDS((pAttributeTensor1)AttList[TV], index,
                                                 centroid);
  if (AttList[HF])                    // heat flux
//old
//    BCB[5] = Att_evalTensorOr0(AttList[HF], centroid);
      BCB[5] = AttributeTensor0_evalDS((pAttributeTensor0)AttList[HF],
                                       centroid);

  for( k=0; k < nflx-5; k++)
    if (AttList[F1+k])                  // scalar flux
//old
//      BCB[6+k] = Att_evalTensorOr0(AttList[F1+k], centroid);
        BCB[6+k] =
            AttributeTensor0_evalDS((pAttributeTensor0)AttList[F1+k],centroid);

  ibcb = 0;                            // Updating ibcb for this one
  ibcb += (AttList[MF]) ? 1 : 0;
  ibcb += (AttList[NP]) ? 2 : 0;
  ibcb += (AttList[TV]) ? 4 : 0;
  ibcb += (AttList[HF]) ? 8 : 0;
  ibcb += (AttList[TW]) ? 16 : 0;
  for( k=0; k < nflx-5; k++)
  ibcb += (AttList[F1+k]) ? (int)(32*pow(2.0,k)) : 0;

  iBCB[0]=ibcb;
  iBCB[1]=0;
  if(AttList[SID]) {
//old
//    iBCB[1] = (int) Att_evalTensorOr0(AttList[SID], centroid);
    iBCB[1] = (int) AttributeTensor0_evalDS((pAttributeTensor0)AttList[SID], centroid);
    //
    //  above we seemingly misuse a real attribute to obtain an
    //  integer instead of the "right-use" shown below.  We do this
    //  because we would like to be able to vary the srfID on a face
    //  without actually cutting the model face (i.e. currently used
    //  for DtN BC's).
    // old old  iBCB[1] = Att_int(AttList[SID]);
    //  iBCB[1] = AttributeInt_value(AttList[SID]);
  }
  return ibcb;
}
