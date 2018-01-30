#include "phParAdapt.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mpi.h"

using namespace std;

//#ifdef __cplusplus
//extern "C" {
//#endif

extern pMeshDataId locMaxInterpolErrorID;
extern pMeshDataId globMaxInterpolErrorID;

extern 
void PMU_commuInt(int* ds, int* dr);
extern
void PMU_commuArr(void **s, int *ns, void **r, int *nr, MPI_Datatype type,
                  int mult);

void
commuMaxLocalPartLocError(pParMesh pmesh, pMesh mesh)
{
    int mult = sizeof(P_int)/sizeof(int);
    int *ns = (int*) malloc(sizeof(int)*PMU_size());
    int *nr = (int*) malloc(sizeof(int)*PMU_size());
    void ***se = (void***) malloc(sizeof(void**)*PMU_size());
    void ***re = (void***) malloc(sizeof(void**)*PMU_size());
    pEntity ent;

    double** dsend = new double*[PMU_size()];
    double** drecv = new double*[PMU_size()];

    for(int i=0; i< PMU_size(); i++)  {
        ns[i] = 0;
    }

    // loop over part bdry and exchange info about local 
    // vertex id (0 FOR VERTICES)
    pBdryProcIter myBdryIter = PM_bdryProcIter(pmesh,0);
    pEntity en;
    // determine the size of *ns[i]
    while(en = BdryProcIter_next(myBdryIter)) {
         pEntCopies vCopies = EN_copies(en); 

        // how many of them on OTHER apts/=procs ?
        int vCopiesSize = EntCopies_size(vCopies );
        
        // fill the send array with ents/vertices
        for(int i=0 ; i< vCopiesSize; i++){

            // global id of partition/proc on which en's n'th copy lies.
            int vCopyGid = EntCopies_gid(vCopies,i);

            if (EntCopies_gid(vCopies,i) != PMU_rank())  {

                // increment number of copies shared between this proc and proc
                // with rank vCopyGid
                ns[vCopyGid]++;
            }
        } 
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif
    }
    BdryProcIter_reset(myBdryIter);
    
    for (int i=0; i < PMU_size(); i++)  {

        // ns[i] is   number of ent-copies shared by this proc with
        // proc i (in turn, proc i shares the same number of bdry verts with
        // this proc)
        // ns[i] has to be determined by looping over all bdry verts (of this
        // proc) and by taking all copies that are offproc
        //
        // each se[i] is of length ns[i] and
        // therefore each se[i] will be filled with (off-proc) copies of all ents
        // that are shared with proc i
        // se[i][j] will contain the j-th vertex that is shared between this
        // proc and proc i
        // along with that we communicate a double array of length ns[i]: dsend, drecv
        // 1) local interpolation error via locMaxInterpolErrorID
        se[i] = (void**) malloc(sizeof(void*)*ns[i]);
        dsend[i] = new double[ns[i]];
        ns[i] = nr[i] = 0;
    }

    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {

        // double c[3];
        // V_coord((pVertex)en,c);
        // printf("[%d]: (local) vertex id: %d \n and coords %f %f %f\n",PMU_rank(),EN_id(en),c[0],c[1],c[2]);

        // copies of the vertex on other apts/=procs
        // Does not include the entity en that vCopies  was created with
        pEntCopies vCopies = EN_copies(en); 

        // how many of them on OTHER apts/=procs ?
        int vCopiesSize = EntCopies_size(vCopies );
        

        // fill the send array with ents/vertices
        for(int i=0 ; i< vCopiesSize; i++){

            // global id of partition/proc on which en's n'th copy lies.
            int vCopyGid = EntCopies_gid(vCopies,i);

            if (EntCopies_gid(vCopies,i) != PMU_rank())  {

                // fill send array to proc of rank vCopyGid with the
                // (ns[vCopyGid])-th copy (among all copies of verts shared
                // between this proc and proc vCopyGid)
                se[vCopyGid][ ns[vCopyGid] ] = EntCopies_ent(vCopies,i);
                
                double* localErrorMax;
                // get the local data
                if(!EN_getDataPtr(en, locMaxInterpolErrorID, (void**)&localErrorMax)){
                                
                    cout<<"\nerror in commuMaxLocalPartLocError: no data attached to vertex\n";
                    exit(0);
                }

                // local error max
                dsend[vCopyGid][ ns[vCopyGid]  ] =  localErrorMax[0];

                // increment number of copies shared between this proc and proc
                // with rank vCopyGid
                ns[vCopyGid]++;
            }
        }
#ifdef FMDB
        EntCopies_delete(vCopies);
#endif   
    }
    BdryProcIter_reset(myBdryIter);

    // communicate the NUMBER of copies that are shared with other procs
    PMU_commuInt(ns, nr);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for (int i=0; i < PMU_size(); i++){
        
        re[i] = (void**) malloc(sizeof(void*)*nr[i]);
        drecv[i] = new double[nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);    
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_DOUBLE, 1); 


    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];

            // determine the global max
            // 
            double* currentMax;
            if(!EN_getDataPtr(ent, locMaxInterpolErrorID,
			      (void**)&currentMax)){
                cout<<"\nerror in commuMaxLocalPartLocError: no data attached to vertex\n";
                exit(0);
            }
            double* newMax = new double[1];  
            
            if(drecv[i][j] > currentMax[0]){
                newMax[0] =   drecv[i][j]; 
            }
            else{
                newMax[0] =  currentMax[0]; 
            }

            
            delete [] currentMax;
            EN_deleteData(ent, locMaxInterpolErrorID);
            EN_attachDataPtr( ent, locMaxInterpolErrorID, (void *)
		      newMax );
        }
    }
    // shovel the data attached via locMaxInterpolErrorID 
    // over to data attched via globMaxInterpolErrorID
    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
        
        double* localMax;
        if(!EN_getDataPtr(en, locMaxInterpolErrorID,
			  (void**)&localMax)){
            cout<<"\nerror in commuMaxLocalPartLocError: no data attached to vertex\n";
            exit(0);
        }

        double* globMax = new double[1];
        globMax[0] = localMax[0];
        EN_attachDataPtr( en, globMaxInterpolErrorID, (void *)
                          globMax) ;


#if  ( defined  DEBUG )
//             printf("\nattaching (globBDR) glob  Max loc error: %f \n",globMax[0]);
//             printf("for vertex:");
//             double c[3];
//             V_coord((pVertex)en,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif 

    }//while(en = BdryProcIter_n
    BdryProcIter_delete(myBdryIter);

    for (int i=0; i < PMU_size(); i++)  {
        free(se[i]);
        free(re[i]);
        delete [] dsend[i];
        delete [] drecv[i];
    }

    delete [] dsend;
    delete [] drecv;

    free(se);
    free(re);
    free(ns);
    free(nr);

    // loop over all interior vertices
    // and shovel local data into gobal data
    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;

    while(vertex = VIter_next(vIter)) {

        if(! EN_isOnPartBdry((pEntity)vertex)){

            double* localMax;
            if(!EN_getDataPtr((pEntity)vertex, locMaxInterpolErrorID,
			      (void**)&localMax)){
                    cout<<"\nerror in commuMaxLocalPartLocError: no data attached to vertex\n";
                    exit(0);
            }
            
            double* globMax = new double[1];
            globMax[0] = localMax[0];
            // attach these values
            double* dummy;
            if(EN_getDataPtr((pEntity)vertex, globMaxInterpolErrorID,(void**)&dummy)){
                
                cout<<"\nerror in commuMaxLocalPartLocError: data already attached via nodalHessianID\n";
                V_info(vertex);
                exit(0);
            }
            EN_attachDataPtr( (pEntity)vertex, globMaxInterpolErrorID, (void *)
                              globMax);

#if  ( defined  DEBUG )
//             printf("\nattaching (globINTR) max glob loc error: %f\n",globMax[0]);
//             printf("for vertex:");
//             double c[3];
//             V_coord(vertex,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif

        }
    }// while loop over non bdry vertices
    VIter_delete(vIter);


}

//#ifdef __cplusplus
//}
//#endif


