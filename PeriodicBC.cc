#include <iostream>
#include <math.h>
#ifndef SIM
#include <vector>
#endif
#include "parallel.h"
#include "PeriodicBC.h"
#include "bits.h"
#include "phParAdapt.h"

using namespace std;

extern int ensa_dof;  // bring in the # of variables to help find # of scalars
//old extern "C" int AttributeInt_value(pAttribute);
extern map<pair<string,pGEntity>,void *> GEntityDataContainerPtr;
extern string teta;


#define numAttsP 2        // number of periodic bc attributes

// PS: periodic slave, NOJ: number of jumps
enum { PS, NOJ, perio };     // for readability of periodic stuff
char strAttP[numAttsP][MLEN] = { "periodic slave",
                                 "number of jumps" };

// measure for reasonability
int PeriodicBC::global_sanity = 0;


////////////////////////////////////////////////////////////////////////////////////////////
// csrtr.
//
// construct a periodic BC from a model and model face and
// setup data 
// also the inherited  AttList
////////////////////////////////////////////////////////////////////////////////////////////
PeriodicBC::PeriodicBC(pGModel gm, pGFace gface) : BoundaryCondition ()
{
  this->myangle = 0.0;
//#ifdef SIM
  // retrieve attribute on gface of type "periodic slave"
  if (GEN_attrib((pGEntity)gface,strAttP[PS]))  {    // if slave then only everything set
//  #else
//    std::vector<Attribute*> atts;
//    try {
//      atts = SCOREC_att::retrieveAttributePList((pGEntity)gface,strAttP[PS]);
//    } catch  (AttributeNotExistent) {}
//    if ( atts.size() > 0) {
//#endif

    this->set = true;
    this->dontinherit = 0;
    this->gf = gface;
    this->gtype = Gface;
    this->model = gm;

    AttList = new pAttribute[numAttsP];    // allocate space for the stuff

    /* dont look at AttList[NOJ] for face coz that has to equal 1 */

    for (int i=0; i<numAttsP; i++) {

//#ifdef SIM
      AttList[i] = GEN_attrib((pGEntity)gface, strAttP[i]);
//  #else
//      try{
//        atts = SCOREC_att::retrieveAttributePList((pGEntity)gface,strAttP[i]);
//        AttList[i] = atts[0];
//      } catch (AttributeNotExistent) { AttList[i] = 0;}
//#endif
    }
  }
}

PeriodicBC::PeriodicBC(pGModel gm, pGEdge gedge) : BoundaryCondition ()
{

  this->myangle = 0.0;
//#ifdef SIM
  if (GEN_attrib((pGEntity)gedge, strAttP[PS]))  {    // if slave then only everything set
//  #else
//    std::vector<Attribute*> atts;
//    try {
//      atts = SCOREC_att::retrieveAttributePList((pGEntity)gedge,strAttP[PS]);
//    } catch  (AttributeNotExistent) {}
//    if ( atts.size() > 0) {
//#endif
    this->set = true;
    this->gtype = Gedge;
    this->ge = gedge;
    this->model = gm;
    dontinherit = setbit (dontinherit, perio);

    AttList = new pAttribute[numAttsP];
    for (int i=0; i<numAttsP; i++) {

//#ifdef SIM
      AttList[i] = GEN_attrib((pGEntity)gedge, strAttP[i]);
//  #else
//      try{
//        atts = SCOREC_att::retrieveAttributePList((pGEntity)gedge,strAttP[i]);
//        AttList[i] = atts[0];
//      } catch (AttributeNotExistent) { AttList[i] = 0;}
//#endif
    }
  } else  {

    pPList gefaces = GE_faces(gedge);
    pGFace gfi;
    void *ltmp=0;
    while (gfi = (pGFace)PList_next(gefaces, &ltmp))  { // loop over faces
      PeriodicBC pbc(gm, gfi);    // set bc on the face
      if (pbc.isSet())  {               // should be slave
        this->set = true;
        this->gtype = Gedge;
        this->ge = gedge;
        this->dontinherit = 0;
        this->model = gm;
        break;                   // break from loop => edge is slave if
      }                          //   at least one connecting face is slave
    }
    PList_delete(gefaces);
  }
}

PeriodicBC::PeriodicBC(pGModel gm, pGVertex gvert) : BoundaryCondition ()
{
  this->myangle = 0.0;
//#ifdef SIM
  if (GEN_attrib((pGEntity)gvert, strAttP[PS]))  {    // if slave then only everything set
//  #else
//    std::vector<Attribute*> atts;
//    try {
//      atts = SCOREC_att::retrieveAttributePList((pGEntity)gvert,strAttP[PS]);
//    } catch  (AttributeNotExistent) {}
//    if ( atts.size() > 0) {
//#endif
    this->set = true;
    this->gtype = Gvertex;
    this->gv = gvert;
    this->model = gm;
    dontinherit = setbit (dontinherit, perio);

    AttList = new pAttribute[numAttsP];
    for (int i=0; i<numAttsP; i++) {

//#ifdef SIM
      AttList[i] = GEN_attrib((pGEntity)gvert, strAttP[i]);
//  #else
//      try{
//        atts = SCOREC_att::retrieveAttributePList((pGEntity)gvert,strAttP[i]);
//        AttList[i] = atts[0];
//      } catch (AttributeNotExistent) { AttList[i] = 0;}
//#endif
    }

  }
  else  {

    pPList gvfaces = GV_faces(gvert);
    void* ltmp =0;
    pGFace gfi;

    while (gfi =(pGFace) PList_next(gvfaces,&ltmp))  {
      PeriodicBC pbc(gm, gfi);
      if (pbc.isSet())  {
        this->set = true;
        this->gtype = Gvertex;
        this->gv = gvert;
        this->dontinherit = 0;
        this->model = gm;
        break;
      }
    }
    PList_delete(gvfaces);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////
// gets master face of current face
// returns 1 = numjumps for success
////////////////////////////////////////////////////////////////////////////////////////////
int PeriodicBC::getPerMaster(pGFace *pgfm) 
{                                          
  global_sanity++; // reasonable size of masters ???
  if ( global_sanity > 10000) {
    cout <<"Please check periodicity BC's for possible error"<<endl;
    cout <<"Face"<<endl;
  }
  // gtype is the private model entity type the BC is specified on
  if (isAttSet(PS) && gtype == Gface)  {               // call only if face

      int tag = AttributeInt_value((pAttributeInt)AttList[PS]);
      *pgfm = (pGFace ) GM_entityByTag(model,gtype,tag); // gtype is entity type
      return 1;                                          // number of jumps
  }
  *pgfm = gf;
  return 0;                                            // failure (sort of)
}

int PeriodicBC::getPerMaster(pGEdge *pgem)
{
  global_sanity++;
  if ( global_sanity > 10000) {
    cout <<"Please check periodicity BC's for possible error"<<endl;
    cout <<"Edge"<<endl;
  }

  if (getbit(dontinherit, perio))  {        // if user has set stuff indep.ly
    if (isAttSet(PS) && gtype == Gedge)  {  // just a routine check


       int tag = AttributeInt_value((pAttributeInt)AttList[PS]);    // entity tag
  
      *pgem = (pGEdge ) GM_entityByTag(model,gtype,tag);  // this->model

      /* default num jumps = 1 */

       return isAttSet(NOJ) ? AttributeInt_value((pAttributeInt)AttList[NOJ]) : 1;
    }
  } else if (this->isSet())  {

    pPList gefaces = GE_faces(ge);
    pGFace gfconn;
    void *eftmp=0;
    int noj = 0;

    while (gfconn = (pGFace)PList_next(gefaces,&eftmp)) {  // loop over faces
      PeriodicBC pbcf(model, gfconn); // set per. bc on curr face

      if (pbcf.isSet())  {                       // the face is a slave
        pGFace gfm = 0;
        if (pbcf.getPerMaster(&gfm))  {          // get slave's master
          noj++;                                 // one jump

          pPList gfedges = GF_edges(gfm);
          pGEdge geconn;
          void* fetmp =0;

          double mindist = 1.0e16;
          pGEdge cedge;                                // closest edge
          while (geconn = (pGEdge)PList_next(gfedges,&fetmp))  {      // loop over edges
            /* first check the case for axisymmetric edge */
            if ( ge == geconn ) return 2475;
            double distee = getDistance (ge, geconn);  // dist bet edges
            if ( distee < mindist ) {
                 mindist = distee;
                 cedge = geconn;
            }
          }
          PList_delete(gfedges);

          PeriodicBC pbce(model, cedge);// set per. bc on
                                                   // closest edge
          if (pbce.isSet())  {                     // if that is slave too
            noj += pbce.getPerMaster(pgem);        // a recursive call
           }
          else  *pgem = cedge;  // else this is the GrandMaster
          break;                 // break from the loop over edges (geconn)

        }
        break;                       // break from the loop over faces
      }
    }
    PList_delete(gefaces);
    return noj;                      // return number of jumps
  }
  return 0;                          // failure
}

int PeriodicBC::getPerMaster(pGVertex *pgvm)
{
  global_sanity++;
  if ( global_sanity > 10000) {
    cout <<"Please check periodicity BC's for possible error"<<endl;
    cout <<"Vertex"<<endl;
  }

  /* Exactly similar as for edge above. So refer to that for comments */

  if (getbit(dontinherit, perio))  {
    if (isAttSet(PS) && gtype == Gvertex)  {

        int tag = AttributeInt_value((pAttributeInt)AttList[PS]);    // entity  tag

        *pgvm = (pGVertex ) GM_entityByTag(model,gtype,tag);

       return isAttSet(NOJ) ? AttributeInt_value((pAttributeInt)AttList[NOJ]) : 1;
    }
  }
  else if (this->isSet())  {
    pPList gvfaces = GV_faces(gv);
    pGFace gfconn;
    void* vftmp =0;
    int noj = 0;

    while (gfconn = (pGFace)PList_next(gvfaces, &vftmp))  {
      PeriodicBC pbcf(model, gfconn);

      if (pbcf.isSet())  {
        pGFace gfm = 0;
        if (pbcf.getPerMaster(&gfm))  {
          noj++;
          pPList gfvertices = GF_vertices(gfm);
          pGVertex gvconn;
          double mindist = 1.0e16;
          pGVertex cvert;
          void* fvtmp=0;
          while (gvconn = (pGVertex)PList_next(gfvertices,&fvtmp)){
            if ( gv == gvconn ) return 2475;
            double distvv = getDistance (gv, gvconn);
            if( distvv < mindist ){
               mindist = distvv;
               cvert = gvconn;
            }
          }
          PList_delete(gfvertices);
          PeriodicBC pbcv(model, cvert);
          if (pbcv.isSet())  {
              noj += pbcv.getPerMaster(pgvm);
           }
           else  *pgvm = cvert;
           break;                  // breaks from gvIter(gvconn)
        }
        break;                     // breaks from gfIter(gfconn)
      }
    }
    PList_delete(gvfaces);
    return noj;                    // returns number of jumps
  }
  return 0;
}

int PeriodicBC::getPerMasterTag()     // return tag of master entity
{
  GEntity *gent = 0;

  if (gtype == Gface)  getPerMaster ((pGFace *) &gent);
  else if (gtype == Gedge)  getPerMaster ((pGEdge *) &gent);
  else if (gtype == Gvertex)  getPerMaster ((pGVertex *) &gent);

  return GEN_tag(gent);  // may not be safe coz not checking if master exists
}

double PeriodicBC::getDistance(pGFace gf1, pGFace gf2)   // dist bet centroids
{
  pGVertex gvtx;
  int n1 = 0, n2 = 0, i;
  double centroid1[] = { 0.0, 0.0, 0.0 };
  double centroid2[] = { 0.0, 0.0, 0.0 };
  double x[3];

  pPList gfvertices = GF_vertices(gf1);
  void* fv=0;
  while (gvtx = (pGVertex)PList_next(gfvertices,&fv)){
    n1++;
    GV_point(gvtx, x);
    for (i=0; i<3; i++)  centroid1[i] += x[i];
  }
  PList_delete(gfvertices);
  for (i=0; i<3; i++)  centroid1[i] /= n1;

  gfvertices = GF_vertices(gf2);
  fv=0;
  while (gvtx = (pGVertex)PList_next(gfvertices,&fv)){
    n2++;
    GV_point(gvtx, x);
    for (i=0; i<3; i++)  centroid2[i] += x[i];
  }
  PList_delete(gfvertices);
  for (i=0; i<3; i++)  centroid2[i] /= n2;

  if (n1 == n2)  {     // gf1, gf2 are a periodic couple. so this must be
                       //   true. anyway, routine check.
    return getDistance (centroid1, centroid2,myangle);
  }
  return -1;                 // failure. shouldnt normally happen.
}

double PeriodicBC::getDistance(pGEdge ge1, pGEdge ge2)   // dist bet centers
{
  pGVertex gv11 = GE_vertex(ge1,0);
  pGVertex gv12 = GE_vertex(ge1,1);
  pGVertex gv21 = GE_vertex(ge2,0);
  pGVertex gv22 = GE_vertex(ge2,1);
  double center1[3], center2[3], x1[3], x2[3];
  int i;

  GV_point(gv11, x1);
  GV_point(gv12, x2);
  for (i=0; i<3; i++)  center1[i] = (x1[i] + x2[i])/2.0;

  GV_point(gv21, x1);
  GV_point(gv22, x2);
  for (i=0; i<3; i++)  center2[i] = (x1[i] + x2[i])/2.0;

  return getDistance (center1, center2, myangle);
}

double PeriodicBC::getDistance(pGVertex gv1, pGVertex gv2)  // dist bet verts
{
  double pt1[3], pt2[3];

  GV_point(gv1, pt1);
  GV_point(gv2, pt2);

  return getDistance(pt1, pt2, myangle);
}

// dist bet 2 pts
double PeriodicBC::getDistance(double *xyz1, double *xyz2, double theta)
{
 if( theta == 0.0 ) {
   return sqrt(SQ(xyz1[0]-xyz2[0])+SQ(xyz1[1]-xyz2[1])+SQ(xyz1[2]-xyz2[2]));
 } else {
   // rotate the first set to reach the second
   double xyzr[3];
   xyzr[0] = cos(theta) * xyz1[0] - sin(theta) * xyz1[1];
   xyzr[1] = sin(theta) * xyz1[0] + cos(theta) * xyz1[1];
   xyzr[2] = xyz1[2];
   return sqrt(SQ(xyzr[0]-xyz2[0])+SQ(xyzr[1]-xyz2[1])+SQ(xyzr[2]-xyz2[2]));
 }
}

double PeriodicBC::getAngle(pGFace mface, pMesh mesh)
{
  // Ignoring axisym perio for now.. later remove this..
  return 0.0;

  // return the angle between the two axisymmetric model faces
  // since we assume the model faces are planar, the normal is
  // extracted from the first mesh face on each model face

  pGFace sface = (pGFace)this->gf;
  double  parm[2] ={0.5,0.5};
  double x1[3], x2[3];

  GF_normal_flat(mface,mesh,x2);
  GF_normal_flat(sface,mesh,x1);

  double x3[3] = { -x1[1]*x2[2]+x2[1]*x1[2], x1[0]*x2[2]-x2[0]*x1[2],
                   -x1[0]*x2[1]+x2[0]*x1[1] };
  double norm1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]),
         norm2 = sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);

  // we are using asin isntead of acos here to be able to get the sign
  // of the angle correctly.
  double theta = asin(x3[2]/norm1/norm2);
  myangle = theta;
  return theta;
}

double PeriodicBC::getAngle(pGEdge medge)
{
  double *theta=0;

  // return the angle between the two axisymmetric model edges
  // at least one of the bounding faces must be axisymmetric,
  // so just grab the angle from there
  pPList gefaces = GE_faces(medge);
  void* tmp =0;
  pGFace gf1;
  while (gf1 = (pGFace)PList_next(gefaces,&tmp)){
    if (GEN_dataP((pGEntity)gf1,"teta",(void**)&theta)){
      if (fabs(*theta) > 0.0 ) break;
    }
  }
  PList_delete(gefaces);
  // we are getting the theta form the master face here , this theta
  // was negated before we attached it to the masterface. so we have
  // to reverse sign now.

  myangle = (*theta) * -1.0;
  return myangle;
}

double PeriodicBC::getAngle(pGVertex mvert)
{
  double *theta=0;

  // return the angle between the two axisymmetric model vertices
  // at least one of the bounding faces must be axisymmetric,
  // so just grab the angle from there

  pPList gvfaces = GV_faces(mvert);
  void *tmp=0;
  pGFace gf1;
  while ( gf1 = (pGFace)PList_next(gvfaces,&tmp)){
    if (GEN_dataP((pGEntity)gf1,"teta",(void**)&theta)){
      if ( fabs(*theta) > 0.0 ) break;
    }
  }
  PList_delete(gvfaces);
  // comment same as the above function
  myangle = (*theta)* -1.0;
  return myangle;
}

void PeriodicBC::GF_normal_flat(pGFace gface, pMesh mesh, double *xyz)
{
  double xyz0[3];
  double xyz1[3];
//    pPList gf_edges = GF_edges(gface);
//    pPList ge_vertices;
//    void* etmp = 0;
//    void* vtmp = 0;
//    pGEdge gedge;
//    pGVertex gvertex;
  double v[2][3];

  // Get a mesh face on this Model Face
  FIter fIter = M_classifiedFaceIter(mesh,(pGEntity)gface,0);
  pFace face = FIter_next(fIter);
  pEdge edge;
  pVertex pvtx;
  int dir;
  for(int i = 0; i<2; i++) {
    edge = F_edge(face,i);
    dir  = F_dirUsingEdge(face,edge);
    pvtx = E_vertex(edge,0);
    V_coord(pvtx,xyz0);
    pvtx = E_vertex(edge,1);
    V_coord(pvtx,xyz1);
    v[i][0]=xyz1[0]-xyz0[0];
    v[i][1]=xyz1[1]-xyz0[1];
    v[i][2]=xyz1[2]-xyz0[2];
    if(dir==1){
      v[i][0]=-1.0*v[i][0];
      v[i][1]=-1.0*v[i][1];
      v[i][2]=-1.0*v[i][2];
    }
  }
// Calculate two edge vectors
//    for(int i = 0; i<2; i++ ) {
//      vtmp=0;
//      gedge = (pGEdge)PList_next(gf_edges,&etmp);
//      ge_vertices = GE_vertices(gedge);
//      gvertex = (pGVertex)PList_next(ge_vertices,&vtmp);
// another way to do the same thing
//      gvertex = GE_vertex(gedge,0);

//      GV_point(gvertex,xyz0);
//      gvertex = (pGVertex)PList_next(ge_vertices,&vtmp);
// another way to do the same thing
//      gvertex = GE_vertex(gedge,1);
//      GV_point(gvertex,xyz1);
//      v[i][0]=xyz1[0]-xyz0[0];
//      v[i][1]=xyz1[1]-xyz0[1];
//      v[i][2]=xyz1[2]-xyz0[2];
//      PList_delete(ge_vertices);
//    }
//    PList_delete(gf_edges);
  //now take the cross-product to get a normal vector to both edges
  xyz[0]=v[0][1]*v[1][2]-v[0][2]*v[1][1];
  xyz[1]=v[0][2]*v[1][0]-v[0][0]*v[1][2];
  xyz[2]=v[0][0]*v[1][1]-v[0][1]*v[1][0];
  // and normalize
  double mag = xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2];
  mag=sqrt(mag);
  xyz[0]=xyz[0]/mag;
  xyz[1]=xyz[1]/mag;
  xyz[2]=xyz[2]/mag;
}
