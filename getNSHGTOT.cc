#ifdef SIM
#include "MeshSim.h"
#endif
#ifdef FMDB
#include "AOMD.h"
#endif
#include "phParAdapt.h"
#include "mpi.h"

#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////
// just get the global (total) number of nodes
////////////////////////////////////////////////////////////////////////////////////////////
int
getNSHGTOT(pMesh mesh)
{
    int nshgTotloc = 0; //global number of Nodes (DOFs if linear)
    int nshgTot;
    
    pVertex v;
    VIter vIter = M_vertexIter( mesh );
    while ( v = VIter_next( vIter ) ) {
        
        if( EN_isOwnerPart(( pEntity )v, 0 ) ){

            nshgTotloc++;
        }
    }

    VIter_delete(vIter);
    MPI_Barrier(MPI_COMM_WORLD);// need ALL verts to be assigned to !   

    MPI_Allreduce(&nshgTotloc, &nshgTot, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);


    return nshgTot;


}
