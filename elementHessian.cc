#include "phParAdapt.h"
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef ibm
#define ludcmp ludcmp_
#define lubksb lubksb_
#endif
#include "phParAdapt.h"
#include "Eigen.h"
#ifdef SIM
#include <math.h>
extern "C" double XYZ_volume(dArray *);
#include "MeshTypes.h"
#endif
#ifdef FMDB
  #include "mEntity.h" 
  using namespace AOMD;
#endif

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalGradientID;

// reconstruct the element hessian : 6-component (symmetric)
void
elementHessian(pRegion region, double* elemHessian)
{
  pVertex v;
  double* nodalGradient;

  for(int i=0; i<6; i++) {
    elemHessian[i] = 0.;
  }

  double matrix[16];

  double fieldVectorX[4];
  double fieldVectorY[4];
  double fieldVectorZ[4];

  int four = 4;

  pPList regionVerts;
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
      if(EN_getDataPtr((pEntity)v, nodalGradientID,
                       (void**)&nodalGradient)==NULL){
        cout<<"\nerror in elementHessian: no data attached to vertex (tet element)\n";
        exit(0);
      }
      fieldVectorX[i]= nodalGradient[0];
      fieldVectorY[i]= nodalGradient[1];
      fieldVectorZ[i]= nodalGradient[2];
    }

    int indx[4];
    double fnumber;
    ludcmp(matrix , &four, &four, indx, &fnumber);

    // fieldVector is going to be overridden and will
    // contain the solution
    lubksb(matrix , &four, &four, indx, fieldVectorX );
    lubksb(matrix , &four, &four, indx, fieldVectorY );
    lubksb(matrix , &four, &four, indx, fieldVectorZ );

    // each poly is
    // a0 + a1*x + a2*y + a3*z
//     elemHessian[0] = fieldVectorX[1];// xx
//     elemHessian[1] = fieldVectorX[2];// xy
//     elemHessian[2] = fieldVectorX[3];// xz
//     elemHessian[3] = fieldVectorY[2];// yy
//     elemHessian[4] = fieldVectorY[3];// yz
//     elemHessian[5] = fieldVectorZ[3];// zz

    elemHessian[0] = (fieldVectorX[1]);// xx
    elemHessian[1] = (fieldVectorX[2]+fieldVectorY[1])/2;// xy
    elemHessian[2] = (fieldVectorX[3]+fieldVectorZ[1])/2;// xz
    elemHessian[3] = fieldVectorY[2];// yy
    elemHessian[4] = (fieldVectorY[3]+fieldVectorZ[2])/2;// yz
    elemHessian[5] = fieldVectorZ[3];// zz

#if  ( defined  DEBUG )
    //     for(int i=0; i<6; i++) {
    //         cout<<"\nhessian["<<i<<"] "<<elemHessian[i]<<"\n";
    //     }
#endif
  }
  break;

#ifdef FMDB
  case mEntity::PRISM : {
    double volume;
    volume = R_Volume2(region);
#endif
#ifdef SIM
  case Rwedge : {
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

        if(EN_getDataPtr((pEntity)v, nodalGradientID,
                         (void**)&nodalGradient)==NULL) {
          cout<<"\nerror in elementHessian: no data attached to vertex (wedge element)\n";
          exit(0);
        }

        fieldVectorX[iVert]= nodalGradient[0];
        fieldVectorY[iVert]= nodalGradient[1];
        fieldVectorZ[iVert]= nodalGradient[2];

        V_coord(v,x_y_z[iVert]);
      }
      // volTet = XYZ_volume(x_y_z);

      diffVt(x_y_z[1],x_y_z[0],v01);
      diffVt(x_y_z[2],x_y_z[0],v02);
      crossProd(v01,v02,normal);
      diffVt(x_y_z[3],x_y_z[0],v03);

      volTet = dotProd(normal,v03)/6;
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
      lubksb(matrix , &four, &four, indx, fieldVectorX );
      lubksb(matrix , &four, &four, indx, fieldVectorY );
      lubksb(matrix , &four, &four, indx, fieldVectorZ );

      double tetHessian[6];
      tetHessian[0] = (fieldVectorX[1]);// xx
      tetHessian[1] = (fieldVectorX[2]+fieldVectorY[1])/2;// xy
      tetHessian[2] = (fieldVectorX[3]+fieldVectorZ[1])/2;// xz
      tetHessian[3] = fieldVectorY[2];// yy
      tetHessian[4] = (fieldVectorY[3]+fieldVectorZ[2])/2;// yz
      tetHessian[5] = fieldVectorZ[3];// zz

      for(int i=0; i<6; i++) {
        // elemHessian[i] += volTet*tetHessian[i];
        elemHessian[i] += tetHessian[i];
      }

    } // loop over three tets fitted into prism

    for(int i=0; i<6; i++) {
      // elemHessian[i] = elemHessian[i]/volume;
      elemHessian[i] = elemHessian[i]/3;
#if  ( defined  DEBUG )
      //         cout<<"\nhessian["<<i<<"] "<<elemHessian[i]<<"\n";
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

        if(EN_getDataPtr((pEntity)v, nodalGradientID,
                         (void**)&nodalGradient)==NULL) {
          cout<<"\nerror in elementHessian: no data attached to vertex (wedge element)\n";
          exit(0);
        }

        fieldVectorX[iVert]= nodalGradient[0];
        fieldVectorY[iVert]= nodalGradient[1];
        fieldVectorZ[iVert]= nodalGradient[2];

        V_coord(v,x_y_z[iVert]);
      }
      // volTet = XYZ_volume(x_y_z);

      diffVt(x_y_z[1],x_y_z[0],v01);
      diffVt(x_y_z[2],x_y_z[0],v02);
      crossProd(v01,v02,normal);
      diffVt(x_y_z[3],x_y_z[0],v03);

      volTet = dotProd(normal,v03)/6;
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
      lubksb(matrix , &four, &four, indx, fieldVectorX );
      lubksb(matrix , &four, &four, indx, fieldVectorY );
      lubksb(matrix , &four, &four, indx, fieldVectorZ );

      double tetHessian[6];
      tetHessian[0] = (fieldVectorX[1]);// xx
      tetHessian[1] = (fieldVectorX[2]+fieldVectorY[1])/2;// xy
      tetHessian[2] = (fieldVectorX[3]+fieldVectorZ[1])/2;// xz
      tetHessian[3] = fieldVectorY[2];// yy
      tetHessian[4] = (fieldVectorY[3]+fieldVectorZ[2])/2;// yz
      tetHessian[5] = fieldVectorZ[3];// zz

      for(int i=0; i<6; i++) {
        // elemHessian[i] += volTet*tetHessian[i];
        elemHessian[i] += tetHessian[i];
      }

    } // loop over two tets fitted into prism

    for(int i=0; i<6; i++) {
      // elemHessian[i] = elemHessian[i]/volume;
      elemHessian[i] = elemHessian[i]/2;
#if  ( defined  DEBUG )
      //         cout<<"\nhessian["<<i<<"] "<<elemHessian[i]<<"\n";
#endif
    }
  }
  break;
  default: {
    cout<<"\nError in elementHessian()..."<<endl;
    cout<<"element topology ["<<rTopoType<<"] CANNOT be handelled"<<endl;
    exit(0);
  }
  }

  PList_delete(regionVerts); 
}

#ifdef __cplusplus
}
#endif
