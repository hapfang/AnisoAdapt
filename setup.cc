// This routine is being added in place of setup.c
// This will read the P order distribution from a file or some other
// source and then set all the entity level information.
//
//
// Anilkumar Karanam  ( Summer of 2000)
//

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "parallel.h"
#include "func.h"
#include "ccfunc.h"

using namespace std;

/* global variables */
extern time_t tstart[MAXNT];
extern double eltime[MAXNT];
extern int BYPASS;
extern int CUBES;
extern char meshoutname[];
extern int parted, justpart;

extern pMeshDataId NDOFID;

extern "C" int topology(pRegion rgn)
{
  // This functions determines the topology if any given region
  // 1 : Tetrahedron
  // 2 : Hexahedron
  // 3 : Wedge
  // 5 : Pyramid

  switch(R_numFaces(rgn)){
  case 4:          /* tetrahedron */
    return 1;
    break;         /* the breaks here are unnecessary, just a coding
                      habbit */
  case 6:          /* hehaxhedron */
    return 2;
    break;
  case 5:          /* could be either a Pyramid or a Wedge */
    {
      /* To resolve this issue we check the number of vertices on
         the first face of the region ... We have reached an
         agreement with the people of Mesh Database that a wedge
         always has its tri face first and a Pyramid always has its
         quad-face first ... Or was it the other way around ????

         In anycase we use the first option !

         since numEdges == numVertices for any closed face we check
         number of edges on the first face . */


      switch(F_numEdges(R_face(rgn,0))){
      case 3:
        return 3;
        break;
      case 4:
        return 5;
        break;
      default:
        cerr << "Congratulations : You have got a Mutant Face in your Mesh"
             << endl
             << "It has neither 3 nor 4 vertices "
             << endl
             << "This preprocessor decided to give up Here "
             << endl
             << "Best of Luck ;)" << endl;
        exit(1);
      }
    }
    break;

  default:
    cerr <<" Unidentified Topology in the Input mesh "<<endl;
    cerr <<" Exiting programme " << endl;
    exit(1);
  }
  return 0;
}

// some parts of setup taken out 9/8/04 JM
// pre-release issues with partitioner,
// partitioning NOT needed anyway in phParAdapt

#ifdef PARALLEL
void setup(pParMesh pmesh, globalInfo *info, int ipart)
#else
void setup(pMesh mesh, globalInfo *info)
#endif
{
#ifdef PARALLEL

//taken out 9/8/04 JM
//    // mesh is not yet partitioned -- -P option would have to be specified instead 
//    if (!parted)  {
//      tstart[11] = time(NULL);

//      // CUBES: Use multi-directional Inertial Partitioning, file partition.in must be
//      // provided
//      // BYPASS:do sliced partitioning in specified direction, needs -x also\n
//      // -C, -Y options resp.
//      if (CUBES || BYPASS)  {
//        pPartitioner uptnr;

//        // 3 arguments suficient for PM_newUPartitioner according to examples in DOCU
//        // allpart defined in allpart.cc,partMeshBypass in  partMeshBypass.cc
//        // both of these alternatives carry 4 arguments
//        uptnr = PM_newPartitioner(pmesh, CUBES ? allpart : partMeshBypass, info);

//        // where is Desired total no. of partitions, this only came as
//        // PM_uPartition(pmesh, uptnr) ==> error !!!
//        // how about using int dtnump = PM_totNParts(pmesh) as in one of the
//        // examples???
//        int dtnump = PM_totNParts(pmesh);
//        /* Partition using uptnr into dtnump partitions, 0 => unweighted */
//        PM_uPartition(pmesh, uptnr,dtnump,0);

//        UPartitioner_delete(uptnr);
//      }//if (CUBES || BYPASS) 

//      // use metis partition if no -C, -Y options were specified
//      else  {
        
//        if (PMU_rank() == 0)  {
//            cout<<"\nattempting metisPartition\n";
//        }
      
//        int chkwt = partitionMesh(pmesh);

//        // this came as  PM_metisPartition(pmesh, chkwt); ==> error
//        // how about using int dtnump = PM_totNParts(pmesh) as in one of the
//        // examples???
//        // Do a partitioning unweighted wrt weights on mesh regions and weights on
//        // partitions, using the default partition->process map. See function
//        // doc. This function is simular to  PM_partition(...)
//  //      int dtnump = PM_totNParts(pmesh);
//        int dtnump=PMU_size;
//        PM_metisPartition(pmesh, dtnump,NULL, 0, NULL, 0, NULL);
//        if (PMU_rank() == 0)  {
//            cout<<"\n PM_metisPartition  success\n";
//        }

//      }//else : NOT (CUBES || BYPASS)

//      eltime[11] = difftime(time(NULL),tstart[11]);

//      // additional intermediate partitioned files to be created
//      // -J option has been specified
//      if (justpart)  {
//  //from saurabh/nspre
//  //        PM_write(pmesh, meshoutname);
//  //        if (PMU_rank() == 0)  {
//  //          cout << "Wrote out " << meshoutname << "\n";
//  //          eltime[0]  = difftime(time(NULL), tstart[0]);
//  //          cout << "\nElapsed time: " << eltime[0]/60.0 << " min\n";
//  //        }
//  //        return;
//          PM_write(pmesh, meshoutname);
//          if (PMU_rank() == 0)  {
//              cout << "Wrote out " << meshoutname << "\n";
//              eltime[0]  = difftime(time(NULL), tstart[0]);
//              cout << "\nElapsed time: " << eltime[0]/60.0 << " min\n";
//          }
        
//          pMesh mesh = PM_mesh(pmesh, 0);
        
//          char nctmpstr[20];
//          FILE *fp;
//          void* tmp;
//          int iwillwrite,i,nd;
//          pVertex v; pEdge e; pFace f; pRegion r;
//          sprintf(nctmpstr, "%s/ncrptmp.%d", meshoutname,PMU_rank()+1);
//          fp = fopen(nctmpstr, "w"); 
//          for (tmp=0; v = M_nextVertex(mesh, &tmp);){
//              EN_getDataInt((pEntity)v,info->incorp,&iwillwrite);
//              fprintf(fp,"%d ",iwillwrite);
//          }
//          if (info->edgeson)  {
//              for (tmp=0; e = M_nextEdge(mesh, &tmp);)  {
//                    if (EN_getDataInt((pEntity)e, NDOFID, &nd))  {
//                      EN_getDataInt((pEntity)e,info->incorp,&iwillwrite);
//                      for (i=0; i < nd; i++)
//                          fprintf(fp,"%d ",iwillwrite+i);
//                  }
//              }
//          }
//          if (info->faceson)  {
//              for (tmp=0; f = M_nextFace(mesh, &tmp);)  {
//                  if (EN_getDataInt((pEntity)f, NDOFID, &nd)) {
//                      EN_getDataInt((pEntity)f,info->incorp,&iwillwrite);
//                      for (i=0; i < nd; i++)
//                          fprintf(fp,"%d ",iwillwrite+i);
//                  }
//              }
//          }
//          if (info->regnson)  {
//              for (tmp=0; r = M_nextRegion(mesh, &tmp);)  {
//                  if (EN_getDataInt((pEntity)r, NDOFID, &nd)) { 
//                      EN_getDataInt((pEntity)r,info->incorp,&iwillwrite);
//                      for (i=0; i < nd; i++)
//                          fprintf(fp,"%d ",iwillwrite+i);
//                  }
//              }
//          }
//          fclose(fp);
//          return;        
//      }// if (justpart)  
//    } //if(!parted)

  // if the mesh is already partitioned (in the previous step,
  // having used -J option there
  // now having read in the parted mesh by the -P option
  // assign to this' procs mesh again the one that is on 0th partition
  if( PM_numParts(pmesh) == 0){
      cerr<<"Error in setup(): ["<<PMU_rank()<<"]'s partition does not exist\n";
      exit(1);
  }
  
  pMesh mesh = PM_mesh(pmesh, ipart);

#endif //PARALLEL
  // used to assign following two vars
  info->numnp = M_numVertices(mesh);
  info->numel = M_numRegions(mesh);
}
