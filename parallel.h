/**********************************************************************/
/* parallel data structure definitions and function prototypes        */
/**********************************************************************/

#ifndef _Parallel
#define _Parallel

#ifdef __cplusplus
#include <map>
#endif
#ifdef FMDB
#include "AOMD.h"
#endif
#ifdef SIM
#include "MeshSim.h"
//#include "attributes.h"
// add complimentary include
#include "SimAttribute.h"
#include "mpi.h"
#endif
#include "MeshSimInternal.h"



/* This structure which uniquely identifies each block of elements
   this includes the number of vertices to identify the topology and
   the polynomial order for the block. */



struct blockKey {
   int nen;
   int maxpoly;
   int nenbl;
   int lcsyst;
};
typedef struct blockKey blockKey;

#ifdef __cplusplus

struct lessKey{
   bool operator()(const blockKey b1, const blockKey b2) const
   {
       return ( 1000*b1.nen+b1.maxpoly+10*b1.nenbl+100*b1.lcsyst <
                1000*b2.nen+b2.maxpoly+10*b2.nenbl+100*b2.lcsyst );
   }
};

typedef struct lessKey lessKey;
typedef struct pidKey pidKey;
typedef std::map<blockKey, int, lessKey> forwardblock;

#endif//__cplusplus

/**********************************************************************/
/* this structure describes a communication "task". Each processor    */
/* has as many tasks it has other processors with which it must       */
/* communicate.                                                       */
/*                                                                    */
/* A task can be described by the following data:                     */
/*     tag:       unique integer associated with this communication   */
/*     type:      = 0  this task is a send                            */
/*                = 1  this task is a recieve                         */
/*     other_pid: pid number of the other processor involved in this  */
/*                communication                                       */
/*     numSeg:    number of "segments" (entities) involved in this    */
/*                task                                                */
/*     ents:      list of entities involved in the task               */
/**********************************************************************/
struct Task{
    int tag;  // distiguishes send/recv
    int type; // denotes whether partition is master or slave
    int other_pid;//pid of partner partiiton
    int numSeg;
  pPList ents; // what for ???
  pPList ments;  // what for ???
};
typedef struct Task Task;

/**********************************************************************/
/* this structure contains global information for the processors,     */
/* such as: number of nodes, elements, etc. on each processor         */
// also carries send AND receive tasks
// 
// more or less an equivalent for partition
/**********************************************************************/

struct globalInfo{

  /*
   These  varaibles are used to see if higher order checks need to
   be done at all ??
   they will be turned on even if one entity in the mesh has those
   modes on.

   For example we will never check for face modes on any face unless
   faceson is set to 1. Even when it is set, all the faces need not
   have face modes.
   */
  int edgeson;
  int faceson;
  int regnson;

  /* These represent the maximum values of these in the entire
     mesh  is useful in writing files  */
  // number of element nodes,  number of element nodes on boundary ... 
  int nen;
  int nenb;
  int nedges;
  int nfaces;



  int numnp;                        /* # of vertices */
  int numel;                        /* # of regions */
  int numpbc;                        /* # of prescribed essential bc's */
  int numelb;                        /* # of boundary elements */
  int nshg;       /* # of shp fns on this proc, incl. not owned */
  int nshgTot;    /* # of shp fns across all processors (none counted twice) */
  int nshgOwn;    /* # of owned shp fns on this proc */

  pMeshDataId incorp; /* key to getting sms numbered data migrated */

  int nlwork;                  /* dimension of local work array */


  Task** stask; // send task
  Task** rtask; // receive task

  // what for ??? 
  pPList perSlv, perMst;

  int *ncvec;    /* This is saved now */
};
typedef struct globalInfo globalInfo;

#endif
