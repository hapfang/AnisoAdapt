#include "stdio.h"
#include "func.h"
#include "ccfunc.h"
// #ifndef SIM
// #include "MSops.h"
// #include "MSopsInternal.h"
// #else
#ifdef SIM
#include "MeshSim.h"
#include "MeshSimInternal.h"
#endif
#ifdef FMDB
#include "AOMD.h"
#include "AOMDInternals.h"
#endif

/* global variables */
#define tol 1.0e-8
extern int ensa_dof;

double dist(double xyz1[3],double xyz2[3]){
  double x1,x2,x3;
  x1=xyz1[0]-xyz2[0];
  x2=xyz1[1]-xyz2[1];
  x3=xyz1[2]-xyz2[2];
  return x1*x1+x2*x2+x3*x3;
}

void lin2quad(pMesh mesh, pMesh lmesh, double* qtmpl, double* qtmp, int nshg)
{
  int i;
  int nv=ensa_dof;
  pVertex  opvertex, lastvertT;
  pEdge edge;
  double xyzT[3],xyzL[3],xyzO[3], diste;

  int iprogress;
  double bigdist,xi;

  pVertex lastvert;
  pVertex vertex;
  VIter vIter;

  //  create the list of points that we will be getting qtmp filled at
  //  (i.e. the vertex points and the mid-side points of each edge)
  double* xq = new double [ 3*nshg ];

 int i1,inv,i2;
 int numnpq=M_numVertices(mesh);

 vIter = M_vertexIter(mesh);
 while( vertex = VIter_next(vIter)) {   /* For all the vertices in the mesh */

   /* I have no idea what the next 3 lines of comments mean but I am keeping
      them anyways.. they might mean something to someone... 4/20/01 */

   // i++
   //   int vglob1=i*3
   // ??? are the previous 2 lines same as next
   int vglob1=(EN_id((pEntity)vertex))*3;
   V_coord(vertex,xyzT);
   xq[vglob1+0]=xyzT[0];
   xq[vglob1+1]=xyzT[1];
   xq[vglob1+2]=xyzT[2];
 }
 VIter_delete(vIter);

 EIter eIter =  M_edgeIter(mesh);
 while( edge = EIter_next(eIter)) {  /* for all edges in the mesh */
    lastvert = E_vertex(edge,0);     // vertex 1 for this edge
    opvertex = E_vertex(edge,1);     // vertex 2 for this edge
    int eglob=(numnpq+EN_id((pEntity)edge))*3;
    int vglob1=(EN_id((pEntity)lastvert))*3;
    int vglob2=(EN_id((pEntity)opvertex))*3;
    for(inv=0; inv< 3; inv++){
      i=eglob+inv;
      i1=vglob1+inv;
      i2=vglob2+inv;
      xq[i] =0.5*( xq[i1] + xq[i2]);
    }
  }

 EIter_delete(eIter);
 void* tmp=0;

 vIter = M_vertexIter(lmesh);
 lastvert = VIter_next(vIter);
 int iwalk=0;
 // loop over the xq list
 double mindist;
 for (int ilv=0; ilv<nshg; ilv++){
   xyzT[0]=xq[ilv*3+0];
   xyzT[1]=xq[ilv*3+1];
   xyzT[2]=xq[ilv*3+2];
   while (!0){
     V_coord(lastvert,xyzL);
     mindist=dist(xyzL,xyzT);
     if(mindist > tol){
       iwalk++;

       iprogress=0;
       int ec1 = 0;
       int vne = V_numEdges(lastvert);

       while(ec1 < vne){
         edge = V_edge(lastvert,ec1++);
         opvertex = E_otherVertex(edge,lastvert);
         V_coord(opvertex,xyzO);
         diste=dist(xyzO,xyzT);
         if(diste< tol) {
           lastvert=opvertex;
           goto endWalkEdge;
         }
         if(diste<mindist) {// found a closer vertex
           mindist=diste;
           lastvertT=opvertex;
           iprogress=1;
         }
       } // end loop over edges
       if(iprogress)
         lastvert=lastvertT; // walk across edge that was closest to
       // target
       else   // I am stuck and did not find a better vertex so
         // interpolate between these the best and me
         goto getunstuck;
     } // the vertex was a match so we took else
     else {   // I found a matching node.   Copy solution to this node
       goto endWalkEdge;
     } // close of vertex match
   } // keep walking as long as there are vertices
 endWalkEdge:
   for(inv=0; inv<nv; inv++){
     qtmp[ilv*nv+inv]=qtmpl[(EN_id((pEntity)lastvert))*nv+inv];
   }
   continue;

 getunstuck: //  we check edges again to find point that is second closest

   int vnumEdges = V_numEdges(lastvert);
   int ecount = 0;
   bigdist=1000000.0;  // big value so that we can find second closest
   while(ecount < vnumEdges){
     edge = V_edge(lastvert,ecount++);
     opvertex = E_otherVertex(edge,lastvert);
     V_coord(opvertex,xyzO);
     diste=dist(xyzO,xyzT);
     if(diste<bigdist) {// found a closer vertex
       bigdist=diste;
       lastvertT=opvertex;
     }
   } // end loop over edges
   xi=mindist/(mindist+bigdist);

   for(inv=0; inv<nv; inv++){
     qtmp[ilv*nv+inv]=(1.0-xi)*qtmpl[EN_id((pEntity)lastvert)*nv+inv] +
       xi *qtmpl[EN_id((pEntity)lastvertT)*nv+inv];
   }
 } // find the next entity value
 VIter_delete(vIter);
 double avewalk=iwalk*1.0/nshg;
 printf("Average number of edges walked, %f, \n",avewalk);

 double factor=1.000000000000;
 //  THIS WILL FAIL FOR MULTIPLE TOPOLOGY
 //  if(info->nen == 8) factor = 0.816497;
 //    else factor = 1.00000;

 EIter eIter2 = M_edgeIter(mesh);
 //  eIter = mesh->firstEdge();
 while( edge = EIter_next(eIter2)) {  /* for all edges in the mesh */
   lastvert = E_vertex(edge,0);
   opvertex = E_vertex(edge,1);
   int eglob=((EN_id((pEntity)edge))+numnpq)*nv;
   int vglob1=(EN_id((pEntity)lastvert))*nv;
   int vglob2=(EN_id((pEntity)opvertex))*nv;
   for(inv=0; inv< nv; inv++){
     i=eglob+inv;
     i1=vglob1+inv;
     i2=vglob2+inv;
     qtmp[i] =factor*( qtmp[i1] + qtmp[i2] - 2.0 * qtmp[i]);
   }
 }
 EIter_delete(eIter2);

 delete [] xq;
}
