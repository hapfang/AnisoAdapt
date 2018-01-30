#include "phParAdapt.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "mpi.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalHessianID;
extern pMeshDataId localHessianID;
extern pMeshDataId localPatchVolID;


// communicates missing contributions over procs/=parts
// each part p's  vertex has
// (double localPatchVolume, double[6] localPatchHessian) , i.e 
// V_p = Sum_numSurrEles(VolE) 
// H_p = Sum_numSurrEles(VolE * H(E))
//
// data are gathered via
// localHessianID
// localPatchVolID
// The commu attaches the combined data  via nodalHessianID
void
commuHessiansFromPatch(pParMesh pmesh,pMesh mesh)
{

    int mult = sizeof(P_int)/sizeof(int);
    int *ns = new int [ PMU_size() ];
    int *nr = new int [ PMU_size() ];
    void ***se = new void** [ PMU_size() ];
    void ***re = new void** [ PMU_size() ];

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
        // along with that we communicate a double array of length 7: dsend, drecv
        // 2) patchVolume
        // 3) patchhessian (6x)
        se[i] = new void* [ ns[i] ];

        if(!se[i]){
            printf("[%i] couldn't allocate in commuHessiansFromPatch\n",PMU_rank());
            exit(0);
        }

        dsend[i] = new double[7*ns[i]];
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
                
                double* localHessian;
                // get the local data
                if(!EN_getDataPtr(en, localHessianID, (void**)&localHessian)){
                                
                    cout<<"\nerror in commuHessianFromPatch: no data attached to vertex\n";
                    exit(0);
                }
                double* localPatchVol;
                if(!EN_getDataPtr(en, localPatchVolID,
				  (void**)&localPatchVol)){
                    cout<<"\nerror in commuHessianFromPatch: no data attached to vertex\n";
                    exit(0);
                }
                // patchVolume
                dsend[vCopyGid][ 7*ns[vCopyGid]  ] =  localPatchVol[0];
                
                // patch Hessian
                for(int k=0;k<6;k++){
                    
                    dsend[vCopyGid][ 7*ns[vCopyGid] +1+k ] = localHessian[k];
                }
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
        
        re[i] = new void* [ nr[i] ];
        drecv[i] = new double[7*nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);    
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_DOUBLE, 7); 


    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];

            // 
            // patchVolume: this is already the global: donme in commuGradients
            double* currentVol;
            if(!EN_getDataPtr(ent, localPatchVolID,
			      (void**)&currentVol)){
                cout<<"\nerror in commuHessiansFromPatch: no data attached to vertex\n";
                exit(0);
            }
            // patch Hessian
            double* currentHess;
            if(!EN_getDataPtr(ent, localHessianID,
			      (void**)&currentHess)){
                cout<<"\nerror in commuHessianFromPatch: no data attached to vertex\n";
                exit(0);
            }
            double* newHess = new double[6];  
            for(int k=0;k<6;k++){
                newHess[k]= drecv[i][7*j+1+k];

                newHess[k] += currentHess[k];
            }
            delete [] currentHess;
            EN_deleteData(ent, localHessianID);
            EN_attachDataPtr( ent, localHessianID, (void *)
                              newHess );            

        }
    }
    // shovel the data attached via localHessianID and localPatchVolID
    // over to data attched via nodalHessianID
    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
        
        double* localHessian;
        if(!EN_getDataPtr(en, localHessianID,
			  (void**)&localHessian)){
            cout<<"\nerror in commuHessiansFromPatch: no data attached to vertex\n";
            exit(0);
        }
        double* localPatchVol;
        if(!EN_getDataPtr(en, localPatchVolID,
			  (void**)&localPatchVol)){
            cout<<"\nerror in commuHessiansFromPatch: no data attached to vertex\n";
            exit(0);
        }
        double* nodalHessian = new double[6];
        for(int k=0;k<6;k++){
            
            nodalHessian[k] = localHessian[k]/localPatchVol[0];
        }
        EN_attachDataPtr( en, nodalHessianID, (void *)
                          nodalHessian) ;


#if  ( defined  DEBUG )
//             printf("\nattaching (globBDR) patch hessian:\n %f %f %f\n",nodalHessian[0],nodalHessian[1],nodalHessian[2]);
//             printf(" %f %f %f\n",nodalHessian[3],nodalHessian[4],nodalHessian[5]);
//             printf("for vertex:");
//             double c[3];
//             V_coord((pVertex)en,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif 

    }//while(en = BdryProcIter_n
    BdryProcIter_delete(myBdryIter);

    for (int i=0; i < PMU_size(); i++)  {
        delete [] se[i];
        delete [] re[i];
        delete [] dsend[i];
        delete [] drecv[i];
    }

    delete [] dsend;
    delete [] drecv;
    delete [] se;
    delete [] re;
    delete [] ns;
    delete [] nr;


    // loop over all interior vertices
    // and shovel local data into gobal data
    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;

    while(vertex = VIter_next(vIter)) {

        
        if(! EN_isOnPartBdry((pEntity)vertex)){

            double* localHessian;
            if(!EN_getDataPtr((pEntity)vertex, localHessianID,
			      (void**)&localHessian)){
                    cout<<"\nerror in commuHessiansFromPatch: no data attached to vertex\n";
                    exit(0);
            }
            double* patchHessian = new double[6];
            for(int i=0; i<6; i++){
                patchHessian[i]= 0.0  ;
            }
            double* localPatchVol;
            if(!EN_getDataPtr((pEntity)vertex, localPatchVolID,
			      (void**)&localPatchVol)){
                    cout<<"\nerror in commuHessiansFromPatch: no data attached to vertex\n";
                    exit(0);
            }
            for(int i=0 ; i<6;i++) {

                // localPatchVol is size 1
                patchHessian[i]=localHessian[i]/localPatchVol[0];
            }
            // attach these values
            double* dummy;
            if(EN_getDataPtr((pEntity)vertex, nodalHessianID,(void**)&dummy)){
                
                cout<<"\nerror in commuHessiansFromPatch: data already attached via nodalHessianID\n";
                V_info(vertex);
                exit(0);
            }
            EN_attachDataPtr( (pEntity)vertex, nodalHessianID, (void *)
                              patchHessian);

#if  ( defined  DEBUG )
//             printf("\nattaching (globINTR) patch hessian:\n %f %f %f\n",patchHessian[0],patchHessian[1],patchHessian[2]);
//             printf("\n  %f %f %f\n",patchHessian[3],patchHessian[4],patchHessian[5]);
//             printf("for vertex:");
//             double c[3];
//             V_coord(vertex,c);
//             printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif

        }
    }// while loop over non bdry vertices
    VIter_delete(vIter);
}

#ifdef __cplusplus
}
#endif
