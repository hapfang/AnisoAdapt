#ifdef SIM
#include <iostream>
#include <map>
#include <stdlib.h>
#include "MeshSim.h"
#include "phParAdapt.h"
#include <string.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

// global visibility to all other functions that use this
map<pair<pGEntity,string>,int> GEntityDataContainerInt;
map<pair<pGEntity,string>,double> GEntityDataContainerDbl;
map<pair<pGEntity,string>,void *> GEntityDataContainerPtr;


// replacement for the Simmetrix model attach data functions
///////////////////////////////////////////////////
int
GEN_dataP(pGEntity g, char *str, void** data)
{
    pair<pGEntity,string> p(g,string(str));


    // use global map
    if(GEntityDataContainerPtr.find(p)!=GEntityDataContainerPtr.end()){
        memcpy(data, &GEntityDataContainerPtr[p],sizeof(int*));
        return 1;
    }
    return 0;
}

void
GEN_attachDataP(pGEntity g, char *str, void* data)
{
    pair<pGEntity,string> p(g,string(str));
    // use global map
    
    GEntityDataContainerPtr[p] = data;
}


int
GEN_dataI(pGEntity g, char *str, int* data)
{
    pair<pGEntity,string> p(g,string(str));

    // use global map
    if(GEntityDataContainerInt.find(p)!=GEntityDataContainerInt.end()){
        memcpy(data, &GEntityDataContainerInt[p], sizeof(int*));
        return 1;
    }
    return 0;
}

void
GEN_attachDataI(pGEntity g, char *str, int data)
{
    pair<pGEntity,string> p(g,string(str));

    // use global map
    GEntityDataContainerInt[p] = data;
}



   




#ifdef __cplusplus
}
#endif
#endif
