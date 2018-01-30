#include "phParAdapt.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#ifndef ibm
#define ludcmp ludcmp_
#define lubksb lubksb_
#endif 
#include "Eigen.h"
#ifdef SIM
#include "MeshSim.h"
#include "MeshTypes.h"
#include <math.h>
extern "C" double XYZ_volume(dArray*);
#endif
#ifdef FMDB
  #include "mEntity.h"
#endif

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId ybarID;
extern pMeshDataId nodalGradientID;
extern int adaptOption;
extern double refthreshold;
extern double rhoinp;

// reconstruct the element gradient
void
elementGradient(pRegion region, double* elementGradient)
{
  pVertex v;
  double* nodalData;

  for(int i=0; i<3; i++) {
    elementGradient[i] = 0.;
  }

  double matrix[16];
  double fieldVector[4];
  
  int fieldIndexForGrad=4;
  int four = 4;
  double rho = rhoinp;
 double scale = refthreshold;

  pPList regionVerts;
  // regionVerts are deleted at the end of the routine
  regionVerts = R_vertices(region, 1);

#ifdef FMDB
  int rTopoType = region->getType();
  switch(rTopoType) {
  case mEntity::TET : {
#endif
#ifdef SIM
  int rTopoType = R_topoType(region);
  switch(rTopoType) {
  case Rtet : {
#endif
    // build the linear system
    buildSystem(region,matrix);

    // get the field vals
    for (int i=0; i<PList_size(regionVerts); i++) {
      v = (pVertex)PList_item(regionVerts,i);
      if(EN_getDataPtr((pEntity)v, ybarID,
                       (void**)&nodalData)==NULL){
        cout<<"\nerror in elementGradient: no data attached to vertex (tet element)\n";
        exit(0);
      }
//      fieldVector[i]=sqrt(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[3]*nodalData[3]); 
     fieldVector[i]=nodalData[3]+scale*rho*(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[0]*nodalData[0]); 
//     fieldVector[i]=nodalData[fieldIndexForGrad];
//      fieldVector[i]=nodalData[3]+scale*rho*nodalData[4]*nodalData[4];
//        fieldVector[i]=nodalData[0]+nodalData[1];
//        fieldVector[i]=nodalData[4]*nodalData[4];
//     fieldVector[i]=(nodalData[3]+50*nodalData[4])*nodalData[5]/0.2;
//      fieldVector[i]=100*(nodalData[3]/9000.0+nodalData[4]/170.0+nodalData[5]/0.15);

    }

    int indx[4];
    double fnumber;
    ludcmp(matrix , &four, &four, indx, &fnumber);

    // fieldVector is going to be overridden and will
    // contain the solution
    lubksb(matrix , &four, &four, indx, fieldVector );

    for(int i=0; i<3; i++) {
      // the poly's structure is
      // a0 + a1*x + a2*y + a3*z
      elementGradient[i] = fieldVector[i+1];

#if  ( defined  DEBUG )
      //        cout<<"\ngradient["<<i<<"] "<<elementGradient[i]<<"\n";
#endif
    }
  }
  break;

#ifdef FMDB
  case mEntity::PRISM : {
    double volume;
    volume = R_Volume2(region);
#endif
#ifdef SIM
  case Rwedge : {
    double volume;
#endif
    double volTet, x_y_z[4][3];
    // fit three tets in a prism
    // (diagonal edges of a prism are : V0-V4, V1-V5 and V5-V0)
    int tetVerts[3][4] = {0,1,2,5, 0,1,5,4, 0,4,5,3};

    double v01[3],v02[3],v03[3];
    double normal[3];

    // three tets in a prism
    // compute gradient on prism as vol. averaged gradients on three child tets
    for(int iTet=0; iTet<3; iTet++) {
      for(int iVert=0; iVert<4; iVert++) {
        v = (pVertex)PList_item(regionVerts,tetVerts[iTet][iVert]);

        if(EN_getDataPtr((pEntity)v, ybarID,
                         (void**)&nodalData)==NULL) {
          cout<<"\nerror in elementGradient: no data attached to vertex (wedge element)\n";
          exit(0);
        }
//        fieldVector[iVert]=sqrt(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[3]*nodalData[3]);
        fieldVector[iVert]=nodalData[3]+scale*rho*(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[0]*nodalData[0]);        
//        fieldVector[iVert]=nodalData[fieldIndexForGrad];
//          fieldVector[iVert]=nodalData[3]+scale*rho*nodalData[4]*nodalData[4];
//         fieldVector[iVert]=nodalData[0]+nodalData[1];
//        fieldVector[iVert]=nodalData[4]*nodalData[4];
//     fieldVector[iVert]=(nodalData[3]+50*nodalData[4])*nodalData[5]/0.2;
//      fieldVector[iVert]=100*(nodalData[3]/9000.0+nodalData[4]/170.0+nodalData[5]/0.15);

        V_coord(v,x_y_z[iVert]);
      }
      // volTet = XYZ_volume(x_y_z);

      diffVt(x_y_z[1],x_y_z[0],v01);
      diffVt(x_y_z[2],x_y_z[0],v02);
      crossProd(v01,v02,normal);
      diffVt(x_y_z[3],x_y_z[0],v03);

      volTet = dotProd(normal,v03)/6;

#ifdef FMDB
//       #ifdef DEBUG
      // volume of wedge should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped (and if base tri. is thin)
      // --- probably can just check the volume of the wedge/prism
      if(volume<0.&&volTet<0.) {
        cout<<"\nError in elementGradient()..."<<endl;
        cout<<"negative volume of tet ["<<iTet<<"] in prism/wedge of negative volume"<<endl;

        R_info(region);
        exit(0);
      }
//       #endif
#endif

      // volume of wedge should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped (and if base tri. is thin)
      volTet = fabs(volTet);

      // build the linear system
      buildSystemXYZ(x_y_z,matrix);

      int indx[4];
      double fnumber;
      ludcmp(matrix , &four, &four, indx, &fnumber);

      // fieldVector is going to be overridden and will
      // contain the solution
      lubksb(matrix , &four, &four, indx, fieldVector );

      for(int i=0; i<3; i++) {
        // the poly's structure is
        // a0 + a1*x + a2*y + a3*z
        // elementGradient[i] += volTet*fieldVector[i+1];
        elementGradient[i] += fieldVector[i+1];
      }

    } // loop over three tets fitted into prism

    // gradient for prism/wedge element is
    // vol. weighted average of gradients on three fitted child tets
    for(int i=0; i<3; i++) {
      // elementGradient[i] = elementGradient[i]/volume;
      elementGradient[i] = elementGradient[i]/3;
#if  ( defined  DEBUG )
      //        cout<<"\ngradient["<<i<<"] "<<elementGradient[i]<<"\n";
#endif
    }
  }
  break;

#ifdef FMDB
  case mEntity::PYRAMID : {
    double volume;
    volume = R_Volume2(region);
#endif
#ifdef SIM
  case Rpyramid : {
    double volume;
#endif
    double volTet, x_y_z[4][3];
    // fit two tets in a pyramid
    // (diagonal edge of a pyramid is : V0-V2)
    int tetVerts[2][4] = {{0,1,2,4}, {0,2,3,4}};

    double v01[3],v02[3],v03[3];
    double normal[3];

    // two tets in a pyramid
    // compute gradient on pyramid as vol. averaged gradients on two child tets
    for(int iTet=0; iTet<2; iTet++) {
      for(int iVert=0; iVert<4; iVert++) {
        v = (pVertex)PList_item(regionVerts,tetVerts[iTet][iVert]);

        if(EN_getDataPtr((pEntity)v, ybarID,
                         (void**)&nodalData)==NULL) {
          cout<<"\nerror in elementGradient: no data attached to vertex (wedge element)\n";
          exit(0);
        }
//        fieldVector[iVert]=sqrt(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[3]*nodalData[3]);
      fieldVector[iVert]=nodalData[3]+scale*rho*(nodalData[1]*nodalData[1]+nodalData[2]*nodalData[2]+nodalData[0]*nodalData[0]);
//       fieldVector[iVert]=nodalData[fieldIndexForGrad];
//      fieldVector[iVert]=nodalData[3]+scale*rho*nodalData[4]*nodalData[4];
//        fieldVector[iVert]=nodalData[0]+nodalData[1];
//        fieldVector[iVert]=nodalData[4]*nodalData[4];
//     fieldVector[iVert]=(nodalData[3]+50*nodalData[4])*nodalData[5]/0.2;
//      fieldVector[iVert]=100*(nodalData[3]/9000.0+nodalData[4]/170.0+nodalData[5]/0.15);
        V_coord(v,x_y_z[iVert]);
      }
      // volTet = XYZ_volume(x_y_z);

      diffVt(x_y_z[1],x_y_z[0],v01);
      diffVt(x_y_z[2],x_y_z[0],v02);
      crossProd(v01,v02,normal);
      diffVt(x_y_z[3],x_y_z[0],v03);

      volTet = dotProd(normal,v03)/6;

#ifdef FMDB
//       #ifdef DEBUG
      // volume of pyramid should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped
      // --- probably can just check the volume of the pyramid
      if(volume<0.&&volTet<0.) {
        cout<<"\nError in elementGradient()..."<<endl;
        cout<<"negative volume of tet ["<<iTet<<"] in pyramid of negative volume"<<endl;

        R_info(region);
        exit(0);
      }
//       #endif
#endif

      // volume of pyramid should not be zero but
      // volume of child tet.(s) can be negative
      // if quad. face(s) are warped
      volTet = fabs(volTet);

      // build the linear system
      buildSystemXYZ(x_y_z,matrix);

      int indx[4];
      double fnumber;
      ludcmp(matrix , &four, &four, indx, &fnumber);

      // fieldVector is going to be overridden and will
      // contain the solution
      lubksb(matrix , &four, &four, indx, fieldVector );

      for(int i=0; i<3; i++) {
        // the poly's structure is
        // a0 + a1*x + a2*y + a3*z
        // elementGradient[i] += volTet*fieldVector[i+1];
        elementGradient[i] += fieldVector[i+1];
      }

    } // loop over two tets fitted into prism

    // gradient for pyramid element is
    // vol. weighted average of gradients on two fitted child tets
    for(int i=0; i<3; i++) {
      // elementGradient[i] = elementGradient[i]/volume;
      elementGradient[i] = elementGradient[i]/2;
#if  ( defined  DEBUG )
      //        cout<<"\ngradient["<<i<<"] "<<elementGradient[i]<<"\n";
#endif
    }
  }
  break;
  default: {
    cout<<"\nError in elementGradient()..."<<endl;
    cout<<"element topology ["<<rTopoType<<"] CANNOT be handelled"<<endl;
    exit(0);
  }
  }

  PList_delete(regionVerts);
}

#ifdef __cplusplus
}
#endif
