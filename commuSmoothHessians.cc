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
extern pMeshDataId numSurroundNodesID;


// communicates missing contributions over procs/=parts
// each part p's  vertex has
// average nodal Hessian via nodalHessianID
//
// data are gathered via
// nodalHessianID
// numSurroundNodesID
// The commu attaches the combined data  via nodalHessianID
// ( delete old ones on part bdry)
void
commuSmoothHessians(pParMesh pmesh, pMesh mesh)
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
      // printf("[%i]  size  %i\n",PMU_rank(),ns[i]);
        // ns[i] is   number of ent-copies shared by this proc with
        // proc i (in turn, proc i shares the same number of bdry verts with
        // this proc)
        // ns[i] has to be determined by looping over all bdry verts (of this
        // proc) and by taking all copies that are offproc
        //
        // each se[i] is of length ns[i] and
        // therefore each se[i] will be filled with (off-proc) copies of all ents
        // that are shared with proc i
        // se[i][2*j] will contain the j-th vertex that is shared between this
        // proc and proc i
        // se[i][2*j+1] will contain the list of neighs on that proc which are
        // on part-bdry
        // along with that we communicate a double array of length 7: dsend, drecv
        // 2) numSurroundNodes
        // 3) nodalhessian (6x)
        se[i] = new void* [ ns[i] ];
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

                
                double* nodalHessian;
                // get the local data
                if(!EN_getDataPtr(en, nodalHessianID, (void**)&nodalHessian)){
                                
                    cout<<"\nerror in commuSmoothHessian: no nodalHessian data attached to vertex\n";
                    exit(0);
                }
                double* numSurroundNodes;
                if(!EN_getDataPtr(en, numSurroundNodesID,
                                 (void**)&numSurroundNodes)){
                    cout<<"\nerror in commuSmoothHessian: no numSurroundNodes data attached to vertex\n";
                    exit(0);
                }
                // numSurroundNodes
                dsend[vCopyGid][ 7*ns[vCopyGid]  ] =   numSurroundNodes[0];
                
                // nodal Hessian
                for(int k=0;k<6;k++){
                    
                    dsend[vCopyGid][ 7*ns[vCopyGid] +1+k ] = nodalHessian[k];
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


    // retrieve communicated data and globalize
    for (int i=0; i < PMU_size(); i++)  {

        // nr[i] is number of ent-copies shared by this proc with
        // proc i
        for (int j=0; j < nr[i]; j++)  {
            
            // the vertex itself
            ent = (pEntity) re[i][j];

            // numSurroundNodes: add the incoming numSurroundNodes up to what
            // is already there in ent
            double* currentNumSurrNodes;
            if(!EN_getDataPtr(ent, numSurroundNodesID,
			      (void**)&currentNumSurrNodes)){
                cout<<"\nerror in commuSmoothHessians globalize: no currentNumSurrNodes data attached to vertex\n";
                exit(0);
            }
            double* newNumSurrNodes = new double[1];  
            newNumSurrNodes[0] =   drecv[i][7*j]; 
            newNumSurrNodes[0] +=  currentNumSurrNodes[0];
            delete [] currentNumSurrNodes;
            EN_deleteData(ent, numSurroundNodesID);
            EN_attachDataPtr( ent, numSurroundNodesID, (void *)
		       newNumSurrNodes);

            // get nodal Hessian (has been commuicated previous to this call)
            // and smooth = average
            double* currentHess;
            if(!EN_getDataPtr(ent, nodalHessianID,
			      (void**)&currentHess)){
                cout<<"\nerror in commuSmoothHessian globalize: no currentHess data attached to vertex\n";
                exit(0);
            }
            double* newHess = new double[6];  
            for(int k=0;k<6;k++){
                newHess[k]= drecv[i][7*j+1+k];

                newHess[k] += currentHess[k];
            }
            delete [] currentHess;
            EN_deleteData(ent, nodalHessianID);
            EN_attachDataPtr( ent, nodalHessianID, (void *)
                              newHess );      

#if  ( defined  DEBUG )
//             printf("\n[%i]attaching globalize in  commuSmoothHessians:\n %f %f %f\n",PMU_rank(),newHess[0],newHess[1],newHess[2]);
//             printf(" %f %f %f\n",newHess[3],newHess[4],newHess[5]);
//             printf("\nand local newSurrVertex Num %f\n",newNumSurrNodes[0]);
//             printf("for vertex:");
//             double cr[3];
//             V_coord((pVertex)ent,cr);
//             printf("coords %f %f %f\n",cr[0],cr[1],cr[2]);
#endif           

        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
      
    // average the data attached via nodalHessianID and numSurroundNodesID
    // over to data attched via nodalHessianID again
    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
        

        double* currentNumSurrNodes;
        if(!EN_getDataPtr(en, numSurroundNodesID,
			  (void**)&currentNumSurrNodes)){
            printf("\n[%i]error in commuSmoothHessians average: no currentNumSurrNodes data attached to vertex\n",PMU_rank());
            exit(0);
        }
        double* nodalHessian;
        if(!EN_getDataPtr(en, nodalHessianID,
			  (void**)&nodalHessian)){
            printf("\n[%i]error in commuSmoothHessians average: no nodalHessian data attached to vertex\n",PMU_rank());
            double c[3];
            V_coord((pVertex)en,c);
            printf("coords %f %f %f\n",c[0],c[1],c[2]);
            exit(0);
        }

        double* newNodalHessian = new double[6];
        for(int k=0;k<6;k++){
	  if(currentNumSurrNodes[0])
            newNodalHessian[k] = nodalHessian[k]/currentNumSurrNodes[0];
	  else
	    newNodalHessian[k] = nodalHessian[k];
        }
        delete [] nodalHessian;
        EN_deleteData(en, nodalHessianID);
        EN_attachDataPtr( en, nodalHessianID, (void *)
                          newNodalHessian) ;


#if  ( defined  DEBUG )
//             printf("\nattaching (globBDR) smooth hessian:\n %f %f %f\n",newNodalHessian[0],newNodalHessian[1],newNodalHessian[2]);
//             printf(" %f %f %f\n",newNodalHessian[3],newNodalHessian[4],newNodalHessian[5]);
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


    // loop over all verts that are NOT on part. bdry.
    // and shovel local data into gobal data
    VIter vIter;
    vIter = M_vertexIter(mesh);
    pVertex vertex;

    while(vertex = VIter_next(vIter)) {

        double* nodalHessian; 
        if(! EN_isOnPartBdry((pEntity)vertex)){

            if(!EN_getDataPtr((pEntity)vertex, nodalHessianID,
			      (void**)&nodalHessian)){
                    cout<<"\nerror in INTR commuSmoothHessians: no nodalHessian data attached to vertex\n";
                    exit(0);
            }
            
            double*  numSurrVerts;
            double* smoothHessian = new double[6];
            for(int i=0; i<6; i++){
                smoothHessian[i]= 0.0  ;
            }
            if(!EN_getDataPtr((pEntity)vertex,numSurroundNodesID ,
			      (void**)&numSurrVerts)){
                    cout<<"\nerror in INTR commuSmoothHessians: no numSurrVerts data attached to vertex\n";
                    exit(0);
            }
            for(int i=0 ; i<6;i++) {

	      if(numSurrVerts)
                // build average
		smoothHessian[i]=nodalHessian[i]/numSurrVerts[0];
	      else
		smoothHessian[i]=nodalHessian[i];
            }
            // attach these values
            delete [] nodalHessian;
            EN_deleteData((pEntity)vertex, nodalHessianID);

            EN_attachDataPtr( (pEntity)vertex, nodalHessianID, (void *)
                              smoothHessian);
        }
    }// while loop over non bdry vertices
    VIter_delete(vIter);

}

#ifdef __cplusplus
}
#endif
