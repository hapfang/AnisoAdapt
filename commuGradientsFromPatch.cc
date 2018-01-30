#include "phParAdapt.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>

#include "mpi.h"

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalGradientID;
extern pMeshDataId localGradientID;
extern pMeshDataId localPatchVolID;


// communicates missing contributions over procs/=parts
// each part p's  vertex has
// (double localPatchVolume, double[3] localPatchGradient) , i.e 
// V_p = Sum_numSurrEles(VolE) 
// G_p = Sum_numSurrEles(VolE * grad(E))
//
// data are gathered via
// localGradientID
// localPatchVolID
// The commu attaches the combined data  via nodalGradientID
void
commuGradientsFromPatch(pParMesh pmesh, pMesh mesh)
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
        // along with that we communicate a double array of length 4: dsend, drecv
        // 2) patchVolume
        // 3) patchgradient (3x)
        se[i] = new void* [ ns[i] ];
        dsend[i] = new double[4*ns[i]];
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
                
                double* localGradient;
                // get the local data
                if(!EN_getDataPtr(en, localGradientID, (void**)&localGradient)){
                                
                    cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
                    exit(0);
                }
                double* localPatchVol;
                if(!EN_getDataPtr(en, localPatchVolID,
				  (void**)&localPatchVol)){
                    cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
                    exit(0);
                }
                // patchVolume
                dsend[vCopyGid][ 4*ns[vCopyGid]  ] =  localPatchVol[0];
                
                // patch Gradient
                for(int k=0;k<3;k++){
                    
                    dsend[vCopyGid][ 4*ns[vCopyGid] +1+k ] = localGradient[k];
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
        drecv[i] = new double[4*nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);    
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_DOUBLE, 4); 


    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];

            // patchVolume: add the incoming patchVolume up to what
            // is already there in ent
            double* currentVol;
            if(!EN_getDataPtr(ent, localPatchVolID,
			      (void**)&currentVol)){
                cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
                exit(0);
            }
            double* newVol = new double[1];  
            newVol[0] =   drecv[i][4*j]; 
            newVol[0] +=  currentVol[0];
            delete [] currentVol;
            EN_deleteData(ent, localPatchVolID);
            EN_attachDataPtr( ent, localPatchVolID, (void *)
		      newVol );

            // patchGradient
            double* currentGrad;
            if(!EN_getDataPtr(ent, localGradientID,
			      (void**)&currentGrad)){
                cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
                exit(0);
            }
            double* newGrad = new double[3];  
            for(int k=0;k<3;k++){
                newGrad[k]= drecv[i][4*j+1+k];

                newGrad[k] += currentGrad[k];
            }
            delete [] currentGrad;
            EN_deleteData(ent, localGradientID);
            EN_attachDataPtr( ent, localGradientID, (void *)
                              newGrad );            

        }
    }
    // shovel the data attached via localGradientID and localPatchVolID
    // over to data attched via nodalGradientID
    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
        
        double* localGradient;
        if(!EN_getDataPtr(en, localGradientID,
			  (void**)&localGradient)){
            cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
            exit(0);
        }
        double* localPatchVol;
        if(!EN_getDataPtr(en, localPatchVolID,
			  (void**)&localPatchVol)){
            cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
            exit(0);
        }
        double* nodalGradient = new double[3];
        for(int k=0;k<3;k++){
            
            nodalGradient[k] = localGradient[k]/localPatchVol[0];
        }
        EN_attachDataPtr( en, nodalGradientID, (void *)
                          nodalGradient) ;


#if  ( defined  DEBUG )
//             printf("\nattaching (globBDR) patch gradient: %f %f %f\n",nodalGradient[0],nodalGradient[1],nodalGradient[2]);
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
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loop over all interior vertices
    // and shovel local data into gobal data
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;
    
    while(vertex = VIter_next(vIter)) {

//          double cd[3];
//          V_coord(vertex,cd);
//          printf("[%d]:  INTR vertex coords %f %f %f\n",PMU_rank(),cd[0],cd[1],cd[2]);
       
        if(! EN_isOnPartBdry((pEntity)vertex)){

            double* localGradient;
            if(!EN_getDataPtr((pEntity)vertex, localGradientID,
			      (void**)&localGradient)){
                    cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
                    exit(0);
            }
            
            double patchVolume = 0;
            double* patchGradient = new double[3];
            for(int i=0; i<3; i++){
                patchGradient[i]= 0.0  ;
            }
            double* localPatchVol; 
            if(!EN_getDataPtr((pEntity)vertex, localPatchVolID,
			      (void**)&localPatchVol)){
                    cout<<"\nerror in commuGradientsFromPatch: no data attached to vertex\n";
                    exit(0);
            }
            for(int i=0 ; i<3;i++) {

                // localPatchVol is size 1
                patchGradient[i]=localGradient[i]/localPatchVol[0];
            }
            // attach these values
            double* dummy;
            if(EN_getDataPtr((pEntity)vertex, nodalGradientID,(void**)&dummy) !=
                0){
      
                cout<<"\nerror in commuGradientsFromPatch: data already attached via nodalGradientID\n";
                V_info(vertex);
                exit(0);
            }
            EN_attachDataPtr( (pEntity)vertex, nodalGradientID, (void *)
                              patchGradient);
#if  ( defined  DEBUG )
//             printf("\nattaching (globINTR) patch gradient: %f %f %f\n",patchGradient[0],patchGradient[1],patchGradient[2]);
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
