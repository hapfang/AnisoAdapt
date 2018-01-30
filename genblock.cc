/* This routine generates the block structure for the mesh
   Fall 2000 AKK

   Rememmber to increase the allocation of NBblock from 12 to 18 when
   we get either cubic wedges or pyramids working.

*/
#include <iostream>
#include <stdio.h>
#include <map>
#include <vector>
#include <stdlib.h>
#include "parallel.h"
#include "phParAdapt.h"
// #ifndef SIM
// #include "modeler.h"
// #include "MSopsInternal.h"
// using namespace SCOREC_mesh;
// using namespace SCOREC_model;
// #else
#ifdef SIM
#include "MeshSim.h"
#include "MeshSimInternal.h"
#include "SimPartitionedMesh.h"
#endif
#ifdef FMDB
#include "AOMD.h"
#endif

int* NIblock;
int* NBblock;
int* Nshape;
int* NshapeB;

forwardblock Iblock;
forwardblock Bblock;
blockKey *RIblock;
blockKey *RBblock;

extern "C" int topology(pRegion region);

extern pMeshDataId POLYID;
extern pMeshDataId RNENID;
extern map<pair<string,pGEntity>,int> GEntityDataContainerInt;
extern string NoBI;


static int getnumIshape(blockKey bKey)
{
  int lnen = bKey.nen;
  int poly = bKey.maxpoly;
  int nedg = lnen +(lnen+1)/2;
  int numshp=lnen;

  if( poly < 2 ){ return numshp;}
  else { numshp += nedg*(poly-1);}

  if (poly < 3 ) return numshp ;
  switch(lnen){
  case 4:
  case 5:
    numshp += 4;
    break;
  case 6:
    numshp += 2;
    break;
  default:
    break;
  }
  if ( poly > 3 ){
    fprintf(stderr, " Code not ready for p > 3 \n");
    exit(-1);
  }
  return numshp;
}

int getnumBshape(blockKey bkey)
{
  int numEdges = bkey.nenbl;
  int poly = bkey.maxpoly;
  int numshp = numEdges;
  if( poly > 1) numshp += numEdges *( poly -1);
  if(( poly > 2) && (numEdges == 3))
    numshp ++;
  if ( poly > 3 ){
    fprintf(stderr, " Code not ready for p > 3 \n");
    exit(-1);
  }
  return numshp;
}

void genblock(pMesh mesh, globalInfo *info, pPList bdry)
{
  pGModel model = M_model(mesh);
  int counter;

  /*  we use counter to give unique blockids to different blocks
      on a given pid. but since we use 0 to identify nonexistent blocks we
      should start counting from 1 */
  counter = 1;
  NIblock = new int [12];
  NBblock = new int [12];
  RIblock = new blockKey [ 12 ];
  RBblock = new blockKey [ 12 ];
  Nshape = new int [12];
  NshapeB = new int [12];
  for(int j=0; j< 12; j++) {
    NIblock[j] =0;
    NBblock[j] =0;
    Nshape[j] =0;
    NshapeB[j] =0;
    RIblock[j].nen = 0;
    RIblock[j].maxpoly = 0;
    RIblock[j].nenbl = 0;
    RIblock[j].lcsyst = 0;
    RBblock[j].nen = 0;
    RBblock[j].maxpoly = 0;
    RBblock[j].nenbl = 0;
    RBblock[j].lcsyst = 0;
  }

  RIter mrIter = M_regionIter(mesh);
  pRegion region;
  blockKey bKey;
  /* First we gerenerate the interior element blocking structure */
  /* for this we loop over all the elements in the mesh and insert and
     update them into Iblock  */

  while (region = RIter_next(mrIter)){
    EN_getDataInt((pEntity)region, POLYID, &bKey.maxpoly);
    EN_getDataInt((pEntity)region, RNENID, &bKey.nen);
    bKey.nenbl = (bKey.nen == 8 ? 4 : 3);
    bKey.lcsyst = topology(region);

    if (Iblock[bKey]) {        /* Already existing block */
      NIblock[Iblock[bKey]-1]++;
      /* increment the nbr of elements in */
      /* this block */
    } else {                        /* creation of a new block */
      Iblock[bKey]=counter;  /* assign a unique ascending id */
      counter++;
      NIblock[counter-2]=1;    /* set num elements in block to 1 */
      Nshape[counter-2]  = getnumIshape(bKey);
      RIblock[counter-2] = bKey;
    }
  }

  RIter_delete(mrIter);

  /* After the interior blocks, we create the boundary blocks. */

  counter = 1;
  pGFace gface;
  pFace face;
  GFIter gfIter = GM_faceIter(model);

  while( gface = GFIter_next(gfIter)){
    int tmp;
    if(!GEN_dataI((pGEntity)gface,"NoBI",&tmp)){

      /* If the face has boundary elements */

      FIter fIter = M_classifiedFaceIter(mesh,(GEntity*) gface,0);
      while (face = FIter_next(fIter)) {
        void* rtmp = 0;
        pPList f_regions = F_regions(face);
        region = (pRegion)PList_next(f_regions,&rtmp);
        PList_delete(f_regions);
        EN_getDataInt((pEntity)region, POLYID, &bKey.maxpoly);
        EN_getDataInt((pEntity)region, RNENID, &bKey.nen);
        bKey.nenbl = F_numEdges(face);
        bKey.lcsyst = topology(region);

        /* here we try to seperate quad bface wedges from tri bface wedges*/

        if ( 3 == bKey.lcsyst) bKey.lcsyst = bKey.nenbl;

        /* we also seperate quad face pyramids from triface pyramids */
        /* quad faced pyramids retain lcsyst = 5 and tri faced pyramids
           have their lcsyst changed to 6 */

        if ( 5 == bKey.lcsyst)
          if ( 3 == bKey.nenbl ) bKey.lcsyst = 6;

        if (Bblock[bKey]) {        /* Already existing block */
          NBblock[Bblock[bKey]-1]++ ;
          /* increment the nbr of elements in */
          /* this block */
        } else {                       /* creation of a new block */
          Bblock[bKey] = counter; /* assign a unique ascending id */
          counter++;
          NBblock[counter-2] = 1; /* set num elements in block to 1 */
          NshapeB[counter-2] = getnumIshape(bKey);
          RBblock[counter-2] = bKey;
        }
        bdry = PList_append(bdry,face);
        info->numelb++;
      }
      FIter_delete(fIter);
    }
  }
  GFIter_delete(gfIter);
}
