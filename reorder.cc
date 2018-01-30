/* This is the function of DOF number reordering using reverse Cuthill McKee
   algorithm 
   It also reorders the mesh regions

   Min Zhou summer 2008
*/

#include <list>
#include "phParAdapt.h"
#include "func.h"
#include <iostream>
extern pMeshDataId shpfnNum;
extern pMeshDataId MYCTID;
pMeshDataId doftag;
pMeshDataId ReorderR;
extern pMeshDataId POLYID;
extern forwardblock Iblock;
extern int* NIblock;
extern pMeshDataId RNENID;
using namespace std;
void 
V_reordering(pMesh mesh, globalInfo* info, int ipart){
    pVertex vertex;
    list<pVertex> vList;
    int tag = 1,marked;
    if(!ipart) //only creat new for the first time
        doftag = MD_newMeshDataId("dof tag");
    int nnodes = info->nshg-1;
//    int nnodesOwn =info->nshgOwn-1;
 
    VIter vIter = M_vertexIter(mesh) ;
    while(nnodes!=-1){ //has disjoint part
        vertex = VIter_next(vIter) ; // the first vertex not tagged
        while(EN_getDataInt(vertex,doftag,&marked))
            vertex = VIter_next(vIter);

        vList.push_back(vertex);
        EN_attachDataInt((pEntity)vertex,doftag,tag);
        list<pVertex>::iterator list_iter = vList.begin(); //pop out the first vertex in the list

        while(list_iter!=vList.end()){
            vertex = *list_iter;
            EN_attachDataInt(vertex,MYCTID,nnodes--);
            
            int numEdge = V_numEdges(vertex);
            for(int iedge=0; iedge< numEdge; iedge++){ // looking for the adjacent nodes
                pEdge edge = V_edge(vertex, iedge);
                pVertex V0 = E_vertex(edge,0);  // the first nodes of the edge
                pVertex V1 = E_vertex(edge,1);  // the seconde nodes of the edge
                pVertex adjV = (V0==vertex?V1:V0); // take the adjacent node
                if(!EN_getDataInt((pEntity)adjV, doftag,&marked)){
                    vList.push_back(adjV);
                    EN_attachDataInt((pEntity)adjV,doftag,tag);
                }
            }
            
            vList.pop_front(); // erase the first one in the list
            list_iter = vList.begin();
        }
    }

    VIter_reset(vIter);
    while( vertex = VIter_next(vIter))
        EN_deleteData((pEntity)vertex,doftag);
    
    VIter_delete(vIter);    
}

void R_reordering(pMesh mesh, int ipart){
    if(!ipart)
        ReorderR = MD_newMeshDataId("Reorder Region");
    pVertex vertex;
    pRegion region;
    int tag = 1,marked;
    blockKey BLOCK;
    int blockid, porder;
    int Rlables[12],sum=0; //regions are stored block by block
    for(int iblock=0;iblock<12;iblock++){
        Rlables[iblock]=NIblock[iblock]-1;
        sum += Rlables[iblock];
    }
    
    VIter vIter = M_vertexIter(mesh) ;

    while(sum!=-12){
        vertex = VIter_next(vIter) ; // the first vertex
        while(EN_getDataInt(vertex,doftag,&marked))
            vertex = VIter_next(vIter);
        
        list<pVertex> vList;
        vList.push_back(vertex);
        EN_attachDataInt((pEntity)vertex,doftag,tag);

        list<pVertex>::iterator list_iter = vList.begin(); //pop out the first vertex in the list
    
        while(list_iter!=vList.end()){
            vertex = *list_iter;
            //label the adjcent regions if they are not already been labeled
            pPList regions = V_regions(vertex); //adjacent regions;
            int numRgn = PList_size(regions);
            for(int iRgn=0;iRgn<numRgn;iRgn++){
                region = (pRegion)PList_item(regions,iRgn);
                if(!EN_getDataInt(region,doftag,&marked)){//not have been labeled
                //yet
                //find the block the current mesh region belongs to
                    EN_getDataInt(region,POLYID,&porder);
                    int nen;
                    EN_getDataInt(region,RNENID,&nen);
                    BLOCK.nen = nen;
                    BLOCK.maxpoly = porder;
                    BLOCK.nenbl = nen ==8?4:3;
                    BLOCK.lcsyst=topology(region);
                    blockid = Iblock[BLOCK]-1;
                    
                    //attach the label to region, according to block
                    EN_attachDataInt(region,ReorderR,Rlables[blockid]--);
                    EN_attachDataInt(region,doftag,tag);
                }
            }
            PList_delete(regions);

            //looking for the adjacent vertices, add to the list if not tagged.        
            int numEdge = V_numEdges(vertex);
            for(int iedge=0; iedge< numEdge; iedge++){ // looking for the adjacent nodes
                pEdge edge = V_edge(vertex, iedge);
                pVertex V0 = E_vertex(edge,0);  // the first nodes of the edge
                pVertex V1 = E_vertex(edge,1);  // the seconde nodes of the edge
                pVertex adjV = (V0==vertex?V1:V0); // take the adjacent node
                if(!EN_getDataInt((pEntity)adjV, doftag,&marked)){ //not been tagged yet
                    vList.push_back(adjV);
                    EN_attachDataInt((pEntity)adjV,doftag,tag);
                }
            }
            
            vList.pop_front(); // erase the first one in the list
            list_iter = vList.begin();
        }

        sum = 0;
        for(int iblock =0; iblock<12;iblock++)
            sum+=Rlables[iblock];
        
    }

    VIter_reset(vIter);

    while( vertex = VIter_next(vIter))
        EN_deleteData((pEntity)vertex,doftag);
    
    RIter rIter = M_regionIter(mesh);
    while( region = RIter_next(rIter))
        EN_deleteData((pEntity)region,doftag);

    RIter_delete(rIter);
    VIter_delete(vIter);    
    MD_deleteMeshDataId(doftag);
}
