#include <stdio.h>
//#include <fstream>
#include "MeshSimInternal.h"
#ifdef SIM
#include "SimPartitionedMesh.h"
#include "MeshSim.h"
#include <math.h>
#endif
#include "func.h"
#include "MeshSimAdapt.h"
#include <unistd.h>
#include <strings.h>
#include <string>
#include <iostream>

#ifndef ibm
#define ludcmp ludcmp_
#define lubksb lubksb_
#endif

#include "phParAdapt.h"
using namespace std;

void
display_region( pRegion region ){

  double xyz[3];
  pPList vertices = R_vertices( region, 1 );
  cout << "-------------------"<< endl ;
  for( int i=0; i<PList_size( vertices ); i++) {
      V_coord( (pVertex)PList_item( vertices, i ), xyz );
      cout << xyz[0] <<" "<< xyz[1]<<" "<<xyz[2]<< endl;
  }
  cout << "-------------------"<< endl ;
}

int 
inverseMapE (pEdge edge,
           double* qpt,
           double* pt ) {
    double xyz[3];
    double node1[3], node2[3];
    double dist1, dist2, distance;

    pVertex nd1, nd2;
    xyz[0] = qpt[0];
    xyz[1] = qpt[1];
    xyz[2] = qpt[2];
    
    nd1 = E_vertex(edge, 0);
    nd2 = E_vertex(edge, 1);
    V_coord(nd1, node1);
    V_coord(nd2, node2);

    dist1 = sqrt(dist(node1, xyz));
    dist2 = sqrt(dist(node2, xyz));

    distance = dist1+dist2;
   pt[0] = dist2/distance;
   pt[1] = dist1/distance;
   return 1;  
}

int 
inverseMap( pRegion region, 
            double* qpt,
            double* pt ) {

    // This is the version of Inverse map, which uses the algorithm 
    // by ken from Mesh2Mesh  (MTMHO3)
    // This thing basically, does a linear solution and then tries to
    // get a Iterative Newton correction to it.

    // First to setup the constants of the forward transformation
    //                   x = Ax* xi
    //                   y = Ay* xi
    //                   z = Az* xi
    // Ax,Ay,Az have 8 terms each and can be obtained using the 
    // solution of an 8x8 system,( which is what I am going to do)!!

    double** A;
    double** AA;
    static double M[4][4] ={ {1, 1, 0, 0 }, 
                             {1, 0, 1, 0 },
                             {1, 0, 0, 1 },
                             {1, 0, 0, 0 } };
			   
			   
			   
			   
    int eight = 8;
    int four  = 4;
    static double Mtemp [16];
    double x = qpt[0];
    double y = qpt[1];
    double z = qpt[2];
    double xel[8],yel[8],zel[8];
    double xisol[3];
    double xyz[3];
    int indx[4];
    double fnumber;
    //FMDB does return inverse order of vertices....
    pPList verts1 = R_vertices( region ,1);
    pPList verts = PList_new();
    int mapVerts[4] = {0,2,1,3};
    for (int iVert=0; iVert<4; iVert++)
      PList_append(verts,PList_item(verts1,mapVerts[iVert]));
    PList_delete(verts1);
    A = new double* [3];
    AA = new double*[3];
    for(int i =0; i< 3;i++) A[i] = new double [4];
    for(int i =0; i< 3;i++) AA[i] = new double [4];  


    //creating the LHS
    int k=0;
    for(int i =0; i<4; i++)
        for(int j=0; j<4; j++)
            Mtemp[k++]= M[i][j]; 
  
    
  // LU decompsing the coeff matrix
    ludcmp( Mtemp, &four, &four, indx, &fnumber);
  
    // Creating the RHS
    for(int i=0; i< 4; i++){
        V_coord( ( pVertex ) PList_item( verts, i ), xyz );
        xel[i] = xyz[0];
        yel[i] = xyz[1];
        zel[i] = xyz[2];
    }

    PList_delete(verts);

    for(int i=0; i<4;i++){
        A[0][i]=xel[i];
        A[1][i]=yel[i];
        A[2][i]=zel[i];
    }

    // Now back substituting to get back the correct set of constants
    // for this element.
  
    lubksb( Mtemp, &four, &four, indx, A[0] );
    lubksb( Mtemp, &four, &four, indx, A[1] );
    lubksb( Mtemp, &four, &four, indx, A[2] );


  // Now we have Ax, Ay and Az (where A is the inverse of matrix Mtemp). 
  // Next, we try to get xi, zeta, eta for a given x, y, z.
  // Ax contains the alpha_x in the form of A[0][0] = alpha_x0, 
  // A[0][1] = alpha_x1, and so on. A[1][0] = alpha_y0, 
  // A[1][1] = alpha_y1, and so on. A[2][0] = alpha_z0, 
  // A[2][1] = alpha_z1, and so on. 

  // But first, overwrite the alphas with the solution solved for by
  // paper and pencil.

    AA[0][0] = xel[3];
    AA[0][1] = ( xel[0] - xel[3] );
    AA[0][2] = ( xel[1] - xel[3] );
    AA[0][3] = ( xel[2] - xel[3] );

    AA[1][0] = yel[3];
    AA[1][1] = ( yel[0] - yel[3] );
    AA[1][2] = ( yel[1] - yel[3] );
    AA[1][3] = ( yel[2] - yel[3] );

    AA[2][0] = zel[3];
    AA[2][1] = ( zel[0] - zel[3] );
    AA[2][2] = ( zel[1] - zel[3] );
    AA[2][3] = ( zel[2] - zel[3] );


    int indx2[3];
    int three=3;

    double MS[9];
    k =0;
    for(int i =0;i<3;i++){
        for(int j=1; j<4;j++){
            MS[k++]= A[i][j];
        }
    }

    double xl[3],dxl[3];
    xisol[0] =  x - A[0][0];
    xisol[1] =  y - A[1][0];
    xisol[2] =  z - A[2][0];
    xl[0] = x;
    xl[1] = y;
    xl[2] = z;
    // LU decompsing the coeff matrix and solving for xisol

    ludcmp( MS, &three, &three, indx2, &fnumber);
    lubksb( MS, &three, &three, indx2, xisol);

    double tol = 0.000001;
    //double tol = 0.001;
    //double tol = 0.1;
    int truth =1;                        
  
    for( int i=0; i<3 ; i++) {         
        if ( xisol[i] > 1.0+tol || xisol[i] < 0.0-tol ) {
            truth = 0; 
            //printf("xisol[ %d ] = %f\n",i,xisol[i]);
        } 
    }
    double l4 = 1 - xisol[0] - xisol[1] - xisol[2];
    if ( l4 > 1.0+tol || l4 < 0.0-tol ) {
       truth = 0;
       //printf("l4 = %f\n",l4);
    }
  
    // If a point is outside its parent region, we still continue with a less 
    // accurate interpolation (good enough for a starting solution for can be
    // an issue for dwall)
    //if (truth){                          
        pt[0] = xisol[0];                 
        pt[1] = xisol[1];                 
        pt[2] = xisol[2];  
    //}                                    

    for(int i =0; i< 3;i++) delete [] A[i] ;
    for(int i =0; i< 3;i++) delete [] AA[i] ;
    delete [] A;
    delete [] AA;


    return truth;
}
