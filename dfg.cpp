//C Copyright (C) 2014-2015 Alastair Fyfe
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

/*  dfg : calculate the distance field of the molecular surface from an atomic
          model encoded as a PDB file

 */

// FIXME : when  building a shared object to enable python embedding with cython,
//         the cython build sets NDEBUG and asserts no longer fire. Check for a 
//         cython configuration flag to change stop doing that.
//#undef NDEBUG

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <set>
#include <list>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <limits>
#include <float.h>

#include <clipper/clipper.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/minimol/minimol.h>
#include <clipper/minimol/minimol_utils.h>
#include <clipper/minimol/minimol_io.h>
#include <clipper/clipper-mmdb.h>
#include <iostream>

#include "dfg.h"

using namespace std;
using namespace ::mmdb;
using namespace clipper;

const fmtype DefaultDFRes = 2.;

// debug output can be voluminous (very) so use a vector of flags to control
// which stages should be monitored. each stage goes to a separate file
bool __debugEnabled = false;
// 0 reserved for asserts - always on
bool __debugFlagV[8] = {true, true, true, true, true, true, true, true};
int __debugFileNum;
FILE *__debugCurrFile;

// clipper template utility for non-bonded contacts
// no longer needed ?
//#include "./nonbond.h" 

// Bresenham variant of midpoint circle algorithm adapted for sphere
#include "sphere.h"

#include <ccp4/cmaplib.h>
//#include <ccp4_map_io.h>
void mapToFile(  Xmap<fmtype> xmap,  string filename ); // local export_xmap


#define LOG_ASSERT  do { __debugEnabled = true;  __debugFileNum = 0; }while(0);

void logStr(const char *fmt, ...);

void vlogStr (const char *fmt, va_list ap)
{
  static int currDebugFileNum = -1;
  static char debugFileName[80];
  //printf("changing debug File %d %d\n", currDebugFileNum, __debugFileNum);
  if(currDebugFileNum != __debugFileNum ) {
    if(__debugCurrFile)
      fclose(__debugCurrFile);
    sprintf(debugFileName, "dfg_%d.log", __debugFileNum );
    __debugCurrFile = fopen(debugFileName, "w" );
    printf("changing debug File %d %d\n", currDebugFileNum, __debugFileNum);
    assert (__debugCurrFile);
    setlinebuf(__debugCurrFile);
    currDebugFileNum = __debugFileNum;
  }
  vfprintf(__debugCurrFile, fmt, ap);
}
void logStr(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vlogStr( fmt, ap);
  va_end(ap);
}

#define LOGSTR(x) \
  do { \
    if (__debugEnabled && __debugFlagV[__debugFileNum]) \
      logStr x; \
  } while (0)

// from clipper
class I3 {
  public:
    inline I3() {}
  inline I3( const int& v0, const int& v1, const int &v2 )
    { vec[0] = v0; vec[1] = v1; vec[2] = v2;}
    //! get & set element
    inline const int& operator []( const int& i ) const { return vec[i]; }
    inline int& operator []( const int& i ) { return vec[i]; }

  private:
    int vec[3];
};

inline I3 operator *( const int& s, const I3& v )
{ return I3( s*v[0], s*v[1], s*v[2]); }
inline int operator == ( const I3& v1, const I3& v2 ) 
    { return (v1[0]==v2[0] && v1[1]==v2[1] && v1[2]==v2[2]); }


typedef enum { VState_FAR=0, VState_NARROW, VState_FROZEN, VState_OFFGRID,
	       VState_IGP, VState_INIT }VState;

struct GridIndex : public I3
{
  bool isOffGrid();
  GridIndex(unsigned x, unsigned y, unsigned z);
  const char *toStr() const {
    static char buf[1024];
    sprintf(buf, "(%d,%d,%d)", this->operator[](0),
	    this->operator[](1), this->operator[](2) );
    return buf;
  }
  // new indices offset by "off"
  GridIndex neigh(I3 &off) const {return GridIndex(
				    this->operator[](0) + off[0],
				    this->operator[](1) + off[1],
				    this->operator[](2) + off[2]
				     ); }
  GridIndex(){}
};

struct InitSurfGridPoint {
  //Coord_frac center;
  //fmtype radius;
  fmtype dist;
  class AtomDescr *atmDescr;
  unsigned _flatIx: 30;
  bool useForBoundary:1;
  bool spansZero:1;
  InitSurfGridPoint(unsigned xArg, unsigned yArg, unsigned zArg );
  GridIndex getGridIndex();
};

bool igpCmpFunc(const InitSurfGridPoint &i, const  InitSurfGridPoint&j) {
  return fabs(i.dist)<fabs(j.dist);
}

class MSGrid *_fmGrid;
typedef multimap<fmtype, GridIndex> Heap;
typedef Heap::iterator HeapIt;

struct AtomDescr {
  MAtom *atm;
  Coord_frac cfUnit;
  Coord_orth atmCo;
  fmtype vdw;
  AtomDescr ( MAtom *atmArg,   Coord_frac cfUnitArg,
	      Coord_orth atmCoArg, fmtype vdwArg )
    :atm( atmArg), cfUnit(cfUnitArg), atmCo(atmCoArg),
     vdw(vdwArg) {}
  AtomDescr( const AtomDescr & o ) 
  :atm( o.atm), cfUnit(o.cfUnit), atmCo(o.atmCo),  vdw(o.vdw) {}
};

//FIXME realloc
//const int ALGElem_MAX = 38;
struct ALGElem {
  int nElem;
  AtomDescr *first;
  vector<AtomDescr *> *elemV;
  ALGElem() : nElem(0),elemV(0) {}
  ~ALGElem() {if(elemV) delete elemV; }
};
int maxALGElem;
class AtomListGrid {

  //typedef map<int, vector < AtomDescr*> *>::iterator MapIt_t;
public:
  int nR, nC, nP;
  AtomListGrid( int nR, int nC, int nP ) ;
  ~AtomListGrid() ;
  //vector < AtomDescr*> * elem(int r, int c, int p ) ;
  ALGElem &elem(int r, int c, int p ) ;
  void addAtom(int r, int c, int p, AtomDescr* );
  int maxIndices() { return nR*nC*nP; }
  int usedIndices() { return _map.size(); }
  //vector < AtomDescr*> * elemInit(int i, int j ) ;
  map<int, ALGElem> _map; //grid ix ->AtomDescrs at that ix
private:
  //map<int, vector < AtomDescr*> *> _map; //grid ix ->AtomDescrs at that ix
};
typedef map<int, ALGElem >::iterator MapIt_t;

AtomListGrid::AtomListGrid( int nRArg, int nCArg, int nPArg )
  : nR(nRArg), nC(nCArg), nP(nPArg)
{
  //int nEntries = nR * nC* nP;
  //_map = new  map<int,  vector < AtomDescr*> * >;
  //for( int i = 0; i < nEntries; i++ ) _plane[i] = 0;
}
AtomListGrid::~AtomListGrid() {
  //_map.clear(); not safe?
}

void AtomListGrid::addAtom(int r, int c, int p, AtomDescr *atmDescr) {
  //LOGSTR( ("AtomListPlane try %d %d %d\n", r, c, p ) );
  static ALGElem* newElem;
  if( !newElem ) 
    newElem = new ALGElem;
  //assert(newElem->nElem == 0 );
  int ix = nP*nC*p + nC*c + r;
  pair<MapIt_t, bool> insRC = _map.insert( pair< int, ALGElem  > (ix, *newElem  ));
  ALGElem &el( (*(insRC.first)).second );
  //printf("(%d %d %d ) key %d ins %d nElem  %d\n", r, c, p, (*(insRC.first)).first, insRC.second, el->nElem);
  if(insRC.second == false) {    // key already inserted
    //assert(el->nElem < ALGElem_MAX-1 ) ;
    if(el.nElem == 1 ) {
      el.elemV = new vector<AtomDescr*>;
    }
    if(el.nElem != el.elemV->size()+1 ) {
      printf ("%d %d\n",  el.nElem, el.elemV->size()+1 );
      assert(el.nElem == el.elemV->size()+1 );
    }
    el.elemV->push_back( atmDescr );
  } else {
    newElem = 0;
    assert(el.nElem == 0 );
    el.first = atmDescr;
  }
  el.nElem ++;
  if(  el.nElem > maxALGElem )
    maxALGElem = el.nElem;
  /*
  MapIt_t mit = _map.find(ix);
  if( mit == _map.end()) {
    mit = _map.insert( pair< int, vector < AtomDescr* >*  > (ix,
                       new vector<AtomDescr*>[5] )).first;
  //LOGSTR( ("AtomListPlane populate %d %d %d\n", r, c, p ) );
  }
  (*mit).second->push_back(atmDescr);
  */
}

ALGElem & AtomListGrid::elem(int r, int c, int p ) {
  int ix = nP*nC*p + nC*c + r;
  MapIt_t mit = _map.find(ix);
  static ALGElem nullElem;
  if( mit == _map.end()) 
    return nullElem;
  else
    return (*mit).second;
}


float elemName2vdw(const char * name) {
  if(!strlen(name))
    return -1;
  char buf[100];
  strcpy(buf,name);
  buf[strlen(name)] = 0;
  char *cp = buf;
  while(*cp == ' ') cp++;
  char *tmp=cp;
  while(*tmp ){ 
    char c = *tmp;
    if(isdigit(c) || isspace(c)) {
      *tmp = 0;
      break;
    } 
    *tmp = toupper(c); 
    tmp++; 
  }
  //int vdwTabSz = sizeof(vdwTab)/sizeof(vdwTab[0]);
  for( int i = 0; i < vdwTabSz; i++ ) {
    if(!strcasecmp(cp, vdwTab[i].elemName))
      return vdwTab[i].vdwRadius;
  }
  printf("ERROR: UNMATCHED VDW >%s< tabz %d\n",  cp, vdwTabSz);
  //assert(!"ERROR: UNMATCHED VDW");
  return -1;
}

void prVec(const char* str, vector<AtomDescr*> &v) {
  printf ("%s %u\n", str, v.size() );
  for(int i = 0; i < v.size(); i++ ) {
    printf ("[%d] %p\n", i, v[i] );
  }
}


//
//=========================================
//
fmtype sqrVal(fmtype val) { return val*val; }

const char* stateStr(VState state) {
  switch(state) {
  case VState_FAR: return "FAR"; break;
  case VState_NARROW: return "NARROW"; break;
  case VState_FROZEN: return "FROZEN"; break;
  case VState_OFFGRID: return "OFFGRID"; break;
  case VState_IGP: return "IGP"; break;
  case VState_INIT: return "INIT"; break;
  default: 
    {static char buf[80]; sprintf(buf, "UNKNOWN(%d %x)", state, state); return buf; 
    }break;
  }
}

/*
void MSGrid::checkHeap() {
  fmtype minDist = (*heap.begin()).first;
  //FuzzyComparer fc;
  for (HeapIt it = heap.begin(); it != heap.end();       ++it) {
    fmtype dist = (*it).first;
    if( !(minDist <= dist )) {
      LOGSTR(("min %f dist %f %g %d %d\n", minDist, dist, dist-minDist));
	     // fc(  minDist, dist), fc(3., 4.)));
      assert( minDist <= dist );
    }
  }
}
*/
int solve_quadric(fmtype coeff[3], fmtype sol[2] ) {
  //fmtype a=coeff[0];
  //fmtype b=coeff[1];
  //fmtype c=coeff[2];
  fmtype a=coeff[2];
  fmtype b=coeff[1];
  fmtype c=coeff[0];
  fmtype discrim=b*b-4*a*c;

  if(discrim <0 ) {
    LOGSTR(("ERROR negative discriminant\n"));
    return 0;
  }

  if(a!=0) {
      sol[0]= (-b - sqrt(sqrVal(b)-4.0*a*c)) / (2.0*a);
      sol[1]= (-b + sqrt(sqrVal(b)-4.0*a*c)) / (2.0*a);
  }
  else {
      sol[0]= (2.0*c)/(-b - sqrt(sqrVal(b)-4.0*a*c));
      sol[1]= (2.0*c)/(-b + sqrt(sqrVal(b)-4.0*a*c));
  }
  //LOGSTR(("solve_quadric a %f b %f c %f -> %f %f\n",
  //	 a, b, c, sol[0], sol[1] ));
  return 2;
}

GridIndex::GridIndex(unsigned x, unsigned y, unsigned z) { 
  this->operator[](0) = x;
  this->operator[](1) = y;
  this->operator[](2) = z;
}

struct SurfGridEntry {
  VState state;
  fmtype dist;
  SurfGridEntry() : state(VState_FAR),dist(-1.){}
  SurfGridEntry( const SurfGridEntry& g ) : state(g.state), dist(g.dist) {}
};


class MSGrid {
  I3 N6i[6];
  int nEntries;
  list<InitSurfGridPoint> initList;
  Heap heap;

  MiniMol mMol;
  MMDBfile molMFile;
  Cell origCell, surfCell;
  Spacegroup origSG;
  //Coord_frac surfMaxCf,  surfMinCf;// surfOffsetCf;
  Coord_orth boxMaxCo,  boxMinCo;
  Coord_orth surfMaxCo,  surfMinCo;
  //Coord_frac surfCfRange;
  Coord_orth surfCoRange;
  Grid_sampling surfGridSampling;
  Grid_range surfGridRange;
  Coord_grid surfGridMin,surfGridMax;
  vector<InitSurfGridPoint> igpVec;
  fmtype gridResol[3];
  fmtype idx2_[3];
  fmtype voxelDiag;

  Xmap<fmtype> outMap;

  //SurfGridEntry *grid;
  fmtype *surfGridEntryDistV;
  unsigned char *surfGridEntryStateV;

  fmtype dfResolution;

  void extendSpanZero( InitSurfGridPoint &igp);
  void flagSpanZero( InitSurfGridPoint &igp);
  void heapInsert(fmtype dist, const GridIndex&);
  void heapErase(HeapIt it);
  void checkHeap();
  void neighbors(const GridIndex &gi, GridIndex (&nv)[6]);
  bool recompute(const GridIndex &pi, fmtype &max_sol);
  bool recompute1or2(const GridIndex &pi, fmtype &max_sol, bool use2);
public: 
  int nX, nY, nZ;
  MSGrid(){};
  void init ( const char *pdbFn, fmtype dfResolution);
  Coord_grid frac2SurfGrid(Coord_frac cf);
  Coord_frac surfGrid2CellFrac(Coord_grid cg);

  void allocSurfGrid();
  SurfGridEntry getSurfGridEntry( int i, int j, int k );
  void setSurfGridEntry(SurfGridEntry &ge, int i, int j, int k );
  void showSurfGridEntry(const GridIndex &pi);
  fmtype distance(const GridIndex&);
  VState state(const GridIndex&);
  void freezeVoxel(const GridIndex&, fmtype dist);
  void narrowVoxel(const GridIndex&, fmtype dist);
  void changeVoxelDist(const GridIndex&, fmtype dist);
  void mkInitSurfGridPts();
  void adjustInterior(InitSurfGridPoint&igp, bool checkSpanZero);
  void fmInit();
  void fmLoop();
  void writeMap();
};

inline  InitSurfGridPoint::InitSurfGridPoint(unsigned xArg, unsigned yArg, unsigned zArg )
{ 
  _flatIx = zArg*_fmGrid->nX*_fmGrid->nY + yArg*_fmGrid->nX + xArg;
  dist = 0;
  atmDescr = 0;
  useForBoundary = false;
  spansZero = false;
} 
inline  GridIndex InitSurfGridPoint::getGridIndex() {
  int dZ = _flatIx /  ( _fmGrid->nX* _fmGrid->nY);
  int dY = (_flatIx % ( _fmGrid->nX* _fmGrid->nY)) / _fmGrid->nX;
  int dX = (_flatIx % ( _fmGrid->nX* _fmGrid->nY)) % _fmGrid->nX;
  GridIndex ge(dX, dY, dZ );
  return ge;
}

inline bool GridIndex::isOffGrid() {
  int x, y, z;
  int limX, limY, limZ;
  bool rc = false;
  x = this->operator[](0);
  y = this->operator[](1);
  z = this->operator[](2);
  limX = _fmGrid->nX ;
  limY = _fmGrid->nY ;
  limZ = _fmGrid->nZ ;
  if( x < 0  || x >= (int)limX ||
      y < 0  || y >= (int)limY ||
      z < 0  || z >= (int)limZ )
  {
    rc = true;
  }
  return rc;
}

Coord_frac nearBoundsFrac(Coord_frac cf, Cell &cell) {
  const float fracEpsi = 1e-5;
  /*
  Coord_frac ret = cf.lattice_copy_unit();
  fmtype len = (ret-cf).lengthsq(cell);
  if(len > fracEpsi) {
    LOGSTR(("ERROR nearBounds: %s\n%s %f\n", 
	    cf.format().c_str(),
	    ret.format().c_str(),
	    len));
    assert(len < fracEpsi);
  }
  */
  fmtype u = cf.u(), v= cf.v(), w=cf.w();
  if(fabs(u) < fracEpsi) u= 0;
  if(fabs(v) < fracEpsi) v= 0;
  if(fabs(w) < fracEpsi) w= 0;
  Coord_frac ret(u,v,w);
  assert(ret.u() >= 0 && ret.v() >= 0 && ret.w() >= 0 );
  return ret;
}

Coord_orth nearZeroOrth(Coord_orth co) {
  const float orthEpsi = 1e-4;
  fmtype x = co.x(),  y= co.y(), z=co.z();
  if(fabs(x) < orthEpsi) x= 0;
  if(fabs(y) < orthEpsi) y= 0;
  if(fabs(z) < orthEpsi) z= 0;
  Coord_orth ret(x,y,z);
  return ret;
}

void MSGrid::init( const char *pdbFN, fmtype dfResolutionArg) 
{

  dfResolution = dfResolutionArg;

  N6i[0]=  I3(-1, 0, 0 );
  N6i[1]=  I3( 1, 0, 0 );
  N6i[2]=  I3( 0,-1, 0 );
  N6i[3]=  I3( 0, 1, 0 );
  N6i[4]=  I3( 0, 0,-1 );
  N6i[5]=  I3( 0, 0, 1 );

  MMDBManager& mmdb = static_cast<MMDBManager&>(molMFile);

  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  mmdb.SetFlag( mmdbflags );
  molMFile.read_file( pdbFN );

  assert( mmdb.isSpaceGroup() );
  LOGSTR(("Models %d nAtoms %d sg info :\n", 
	  mmdb.GetNumberOfModels(),mmdb.GetNumberOfAtoms() ));
  molMFile.import_minimol( mMol );

  origSG =  mMol.spacegroup();
  origCell = mMol.cell();
  // FIXME import triggers ObjectCache Leak error
  printf("minimol: read %s isNull %d  cell, sg, asu:\n", pdbFN, mMol.is_null() );
  origCell.debug();
  printf("\n");
  origSG.debug();
  int nAtom = 0;
  boxMaxCo = Coord_orth( -1000., -1000., -1000. );
  boxMinCo = Coord_orth( 1000., 1000., 1000. );
  for ( int p = 0; p < mMol.size(); p++ ) {
    for ( int m = 0; m < mMol[p].size(); m++ ) {
      MMonomer &mono(mMol[p][m]);
      for ( int a = 0; a < mMol[p][m].size(); a++ ) 
	{
	if ( !mMol[p][m][a].is_null() ) {
	  nAtom++;
          MAtom& atom(mMol[p][m][a]);
	  const Coord_orth atmCo =  atom.coord_orth() ;
	  if( atmCo.x() > boxMaxCo.x() )
            boxMaxCo = Coord_orth( atmCo.x(), boxMaxCo.y(), boxMaxCo.z() );
          if( atmCo.x() < boxMinCo.x() )
	    boxMinCo = Coord_orth( atmCo.x(), boxMinCo.y(), boxMinCo.z() );

	  if( atmCo.y() > boxMaxCo.y() )
            boxMaxCo = Coord_orth( boxMaxCo.x(), atmCo.y(), boxMaxCo.z() );
          if( atmCo.y() < boxMinCo.y() )
            boxMinCo = Coord_orth( boxMinCo.x(), atmCo.y(), boxMinCo.z() );

	  if( atmCo.z() > boxMaxCo.z() )
	    boxMaxCo = Coord_orth( boxMaxCo.x(), boxMaxCo.y(), atmCo.z() );
          if( atmCo.z() < boxMinCo.z() )
            boxMinCo = Coord_orth( boxMinCo.x(), boxMinCo.y(),  atmCo.z() );
          
        }
      }
    }
  }
  Coord_orth pad(0.1*origCell.a(),0.1*origCell.b(),0.1*origCell.c());
  //pad=Coord_orth(0.001,0.001,0.001);
  //pad=Coord_orth(1.0,1.0,1.0);
  pad=Coord_orth(0.5,0.5,0.5);
  //pad=Coord_orth(0.,0.,0.);
  //pad=Coord_orth(4.0,4.0,4.0);
  //boxMaxCo = boxMaxCo+pad;
  //boxMinCo = boxMinCo-pad;
  //surfMaxCo = boxMaxCo + pad;
  //surfMinCo = boxMinCo + pad;
  surfMaxCo = boxMaxCo;
  surfMinCo = boxMinCo;
  {
  Grid_sampling outMapGS(origSG, origCell, Resolution(dfResolution), /*rate*/1.5);
  Xmap<float>outMap;
  LOGSTR(("surfMin %s\nsurfMax %s\nboxMin %s\nboxMax %s\n",
	  surfMinCo.format().c_str(),
	  surfMaxCo.format().c_str(),
	  boxMinCo.format().c_str(),
	  boxMaxCo.format().c_str() ));
  outMap.init(origSG, origCell, outMapGS );
  Coord_grid asuMinCg =  outMap.grid_asu().min();
  Coord_grid asuMaxCg =  outMap.grid_asu().max();
  Coord_frac asuMinCf = asuMinCg.coord_frac(outMapGS);
  Coord_frac asuMaxCf = asuMaxCg.coord_frac(outMapGS);
  Coord_orth asuMinCo = asuMinCf.coord_orth(origCell);
  Coord_orth asuMaxCo = asuMaxCf.coord_orth(origCell);
  fmtype  maxX = surfMaxCo.x(), maxY = surfMaxCo.y(), maxZ = surfMaxCo.z();
  if(maxX < asuMaxCo.x()) maxX = asuMaxCo.x();
  if(maxY < asuMaxCo.y()) maxY = asuMaxCo.y();
  if(maxZ < asuMaxCo.z()) maxZ = asuMaxCo.z();
  surfMaxCo=Coord_orth(maxX, maxY, maxZ);

  fmtype  minX = surfMinCo.x(), minY = surfMinCo.y(), minZ = surfMinCo.z();
  if( minX > asuMinCo.x()) minX = asuMinCo.x();
  if( minY > asuMinCo.y()) minY = asuMinCo.y();
  if( minZ > asuMinCo.z()) minZ = asuMinCo.z();
  surfMinCo=Coord_orth(minX, minY, minZ);
  LOGSTR(("asu:\nminCg %s\nmaxCg %s\nminCf %s\nmaxCf %s\nminCo %s\nmaxCo %s\nsurfMin %s\nsurfMax %s\n",
	  asuMinCg.format().c_str(),
	  asuMaxCg.format().c_str(),
	  asuMinCf.format().c_str(),
	  asuMaxCf.format().c_str(),
	  asuMinCo.format().c_str(),
	  asuMaxCo.format().c_str(),
	  surfMinCo.format().c_str(),
	  surfMaxCo.format().c_str() ));
  }
  //surfOffsetCf = -surfMinCf;
  //surfMaxCo = surfMaxCf.coord_orth(origCell);
  //surfMinCo = surfMinCf.coord_orth(origCell);
  //surfCfRange= surfMaxCf -surfMinCf;
  surfCoRange= surfMaxCo -surfMinCo;
  fmtype cellMaxX = surfCoRange.x(), cellMaxY = surfCoRange.y(),cellMaxZ = surfCoRange.z();
  if(cellMaxX < origCell.a() ) cellMaxX = origCell.a();
  if(cellMaxY < origCell.b() ) cellMaxY = origCell.b();
  if(cellMaxZ < origCell.c() ) cellMaxZ = origCell.c();

  // we  will adjust all atom co by surfMinCo to ensure >=0 
  // cell bounds need to include this
  cellMaxX +=  fabs(surfMinCo.x());
  cellMaxY +=  fabs(surfMinCo.y());
  cellMaxZ +=  fabs(surfMinCo.z());

  // lastly round to integral boundary
  cellMaxX = ceil(cellMaxX) ;
  cellMaxY = ceil(cellMaxY) ;
  cellMaxZ = ceil(cellMaxZ) ;
  // surfCell is not original cell, but padded asu in p1
  surfCell=Cell(Cell_descr( cellMaxX, cellMaxY, cellMaxZ,
			    90., 90., 90. )) ;
  LOGSTR(("surfCell: %s\n",  surfCell.format().c_str() ));
  surfGridSampling = Grid_sampling(Spacegroup::p1(), surfCell, Resolution(dfResolution), 
				   /* rate*/ 1.5);
  //surfGridSampling = Grid_sampling(nXtmp, nYtmp, nZtmp);
  surfGridMin = Coord_grid(0, 0, 0);
  surfGridMax = Coord_grid(surfGridSampling.nu()-1, 
                           surfGridSampling.nv()-1, 
                           surfGridSampling.nw()-1);
  surfGridRange=Grid_range(surfGridMin, surfGridMax);
  // NB : in_grid for gridRange tests <= in_grid for gridSampling tests <
  nX = surfGridSampling.nu();
  nY = surfGridSampling.nv();
  nZ = surfGridSampling.nw();
  nEntries = nX*nY*nZ;

  gridResol[0] = surfCell.a()/surfGridSampling.nu();
  gridResol[1] = surfCell.b()/surfGridSampling.nv();
  gridResol[2] = surfCell.c()/surfGridSampling.nw();
  idx2_[0] = 1./gridResol[0]/gridResol[0];
  idx2_[1] = 1./gridResol[1]/gridResol[1];
  idx2_[2] = 1./gridResol[2]/gridResol[2];
  voxelDiag = sqrt(gridResol[0]*gridResol[0] +
                          gridResol[1]*gridResol[1] +
                          gridResol[2]*gridResol[2] );
  /*"maxCf %s\n minCf %s\n "
    " rangeCf %s\n " */
  LOGSTR(("params:\n box max %s\n box min %s\n pad %s\n"
	  "maxCo %s\n minCo %s\n"
	  "rangeCo %s\n"
	  "gs %s grid range min %s max %s\n"
          " gridResol %5.2f %5.2f %5.2f voxelDiag %5.2f\n"
	  " surfGS %d %d %d nX %d nY %d nZ %d nEntries %d",
	  boxMaxCo.format().c_str(),
	  boxMinCo.format().c_str(),
	  pad.format().c_str(),
	  //surfMaxCf.format().c_str(),
	  //surfMinCf.format().c_str(),
	  surfMaxCo.format().c_str(),
	  surfMinCo.format().c_str(),
	  //surfCfRange.format().c_str(),
	  surfCoRange.format().c_str(),
	  surfGridSampling.format().c_str(),
	  surfGridRange.min().format().c_str(),
	  surfGridRange.max().format().c_str(),
	  gridResol[0], gridResol[1], gridResol[2], voxelDiag,
	  surfGridSampling.nu(),
	  surfGridSampling.nv(),
	  surfGridSampling.nw(),
          nX, nY, nZ, nEntries
	  ));
  surfGridEntryDistV = 0;
  surfGridEntryStateV = 0;
  printf("nEntries %d %d x %d  x %d\n", nEntries,  nX, nY, nZ );
}

void MSGrid::allocSurfGrid() {
  assert(nEntries);
  surfGridEntryDistV = new fmtype[nEntries];
  surfGridEntryStateV = new unsigned char[nEntries];
}
inline VState MSGrid::state( const GridIndex &gi ) {
  //Int32Grid::Accessor iAcc = igrid->getAccessor();
  //Coord xyz(gi[0], gi[1], gi[2] );
  //return (VState) iAcc.getValue(xyz);
  unsigned ix = gi[2]*nX*nY + gi[1]*nX + gi[0];
  //return grid[ix].state;
  return  (VState)surfGridEntryStateV[ix];

}

inline fmtype MSGrid::distance( const GridIndex &gi ) {
  //FMTypeGrid::Accessor vAcc = vgrid->getAccessor();
  //Coord xyz(gi[0], gi[1], gi[2] );
  //return vAcc.getValue(xyz);
  unsigned ix = gi[2]*nX*nY + gi[1]*nX + gi[0];
  //return grid[ix].dist;
  return  surfGridEntryDistV[ix];
}

void MSGrid::flagSpanZero( InitSurfGridPoint &igp) {
  GridIndex gi = igp.getGridIndex();
  assert(!gi.isOffGrid() && state(gi)==VState_IGP );
  //SurfGridEntry ge( getSurfGridEntry( igp[0], igp[1], igp[2]));
  SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
  fmtype ourDist = ge.dist;
  assert(! igp.spansZero );
  LOGSTR(("flagSpanZero %d %d %d ourDist %f \n", gi[0], gi[1],gi[2],ourDist));
  for (int dim=0; dim < 3; dim++) { // foreach dim
    fmtype dimDelta = gridResol[dim];
    for( int i = 0; i < 2; i++ ) { // foreach stencil neighbor
      int ix = 2*dim+i;
      GridIndex ngi = gi.neigh(N6i[ix]);
      // neighbors with other states will not have distance set
      if(!ngi.isOffGrid() && 
	 (state(ngi)==VState_IGP || state(ngi) == VState_INIT) ) 
      {
	SurfGridEntry nge( getSurfGridEntry( ngi[0], ngi[1], ngi[2]));
	fmtype neighDist = nge.dist;
	LOGSTR(("adjust neigh %d %d %d ourDist %f \n", ngi[0], ngi[1],ngi[2], neighDist));
	if(neighDist * ourDist < 0 ) {
	  igp.spansZero = true;
	  ge.state=VState_INIT;
	  //setSurfGridEntry( ge,  igp[0], igp[1], igp[2]);
	  setSurfGridEntry( ge,  gi[0], gi[1], gi[2]);
          // not sure why this check : distances of neighbor and this should
          // be less than grid spacing ?

	  //nge.state=VState_INIT;
	  fmtype nDist = fabs(neighDist-ourDist);
	  fmtype distDiff = nDist - dimDelta;
	  /*
	  LOGSTR(("adjust %d %d %d spanZero %d %d %d our %f neigh %f over %f\n",
		  igp[0], igp[1],igp[2],
		  pni[0], pni[1],pni[2],
		  ourDist, neighDist, distDiff));
	  */
	  //assert( distDiff <= 0 );
	  assert( distDiff <= 1e6 );
	  goto done;
	}
      }
    }
  }
 done:
  return;
}

void MSGrid::extendSpanZero( InitSurfGridPoint &igp) {
  GridIndex gi =  igp.getGridIndex();
  assert(!gi.isOffGrid() && state(gi)==VState_IGP );
  //SurfGridEntry ge( getSurfGridEntry( igp[0], igp[1], igp[2]));
  SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
  fmtype ourDist = ge.dist;
  assert(! igp.spansZero && ourDist >= 0 );
  LOGSTR(("extendSpanZero %d %d %d ourDist %f \n", gi[0], gi[1],gi[2], ourDist));
  // if any of igp's  neighbors are spanZero points, igp promoted to INIT
  for (int dim=0; dim < 3; dim++) { // foreach dim
    fmtype dimDelta = gridResol[dim];
    for( int i = 0; i < 2; i++ ) { // foreach stencil neighbor
      int ix = 2*dim+i;
      GridIndex ngi = gi.neigh(N6i[ix]);
      // neighbors with other states will not have distance set
      if(!ngi.isOffGrid() && 
	 (state(ngi)==VState_IGP || state(ngi) == VState_INIT) ) 
      {
	SurfGridEntry nge( getSurfGridEntry( ngi[0], ngi[1], ngi[2]));
	fmtype neighDist = nge.dist;
	if(nge.state == VState_INIT) {
	  ge.state=VState_INIT;
	  //setSurfGridEntry( ge,  igp[0], igp[1], igp[2]);
	  setSurfGridEntry( ge,  gi[0], gi[1], gi[2]);
	}
	/*
	assert(ourDist > neighDist);
	fmtype nDist = ourDist-neighDist;
	fmtype distDiff = nDist - dimDelta;
	assert( distDiff <= 0 );    
	*/
      }
    }
  }
  return;
}

void MSGrid::adjustInterior( InitSurfGridPoint &igp, bool checkSpanZero) {
  GridIndex gi =  igp.getGridIndex();
  assert(!gi.isOffGrid() && state(gi)==VState_IGP );
  //SurfGridEntry ge( getSurfGridEntry( igp[0], igp[1], igp[2]));
  SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
  fmtype lDist[6];
  fmtype nDist[6];
  int nNeigh = 0;
  bool change=false;
  fmtype ourDist = ge.dist;
  LOGSTR(("adjust check %d %d %d span0 %d ourDist %f \n", gi[0], gi[1],gi[2],igp.spansZero, ourDist));
  if(checkSpanZero == false && igp.spansZero == false)
    return;
  do {
    for (int dim=0; dim < 3; dim++) { // foreach dim
      fmtype dimDelta = gridResol[dim];
      for( int i = 0; i < 2; i++ ) { // foreach stencil neighbor
	int ix = 2*dim+i;
	lDist[ix] = 0;
	nDist[ix] = 0;
	GridIndex ngi = gi.neigh(N6i[ix]);
	if(!ngi.isOffGrid() && 
	   (state(ngi)==VState_IGP || state(ngi) == VState_INIT) ) 
        {
	  SurfGridEntry nge( getSurfGridEntry( ngi[0], ngi[1], ngi[2]));
	  fmtype neighDist = nge.dist;
	  if(checkSpanZero) {
	    LOGSTR(("adjust neigh %d %d %d ourDist %f \n", ngi[0], ngi[1], ngi[2], 
		    neighDist));
	    if(neighDist * ourDist < 0 ) {
	      nDist[ix] = neighDist;
	      igp.spansZero = true;
	      ge.state=VState_INIT;
	      setSurfGridEntry( ge,  gi[0], gi[1], gi[2]);
	      //nge.state=VState_INIT;
	      fmtype nDist = fabs(neighDist-ourDist);
	      fmtype distDiff = nDist - dimDelta;
	      LOGSTR(("adjust %d %d %d spanZero %d %d %d our %f neigh %f over %f\n",
		      gi[0], gi[1], gi[2],
		      ngi[0], ngi[1], ngi[2],
		      ourDist, neighDist, distDiff));
	      //assert( distDiff <= 0 );
	      assert( distDiff <= 1e6 );
              /*
	      if( distDiff > 0 )  {
		lDist[ix] = distDiff;
		nNeigh++;
		change = true;
		if(ourDist < 0 ) {
		  ge.dist += distDiff/2;
		  nge.dist -= distDiff/2;
		} else {
		  ge.dist -= distDiff/2;
		  nge.dist += distDiff/2;
		}
	      }
              */
	    }
	  } else {
	    if(nge.state == VState_INIT) {
	      ge.state=VState_INIT;
	      setSurfGridEntry( ge,  gi[0], gi[1], gi[2]);
            }
	    assert(ourDist > neighDist);
	    fmtype nDist = ourDist-neighDist;
	    fmtype distDiff = nDist - dimDelta;
	    assert( distDiff <= 0 );
	    /*
	    if( distDiff > 0 )  {
	      ge.dist -= distDiff;
            }
	    */
	  }
	}
      }
    }
  } while(change);
  if(igp.spansZero && nNeigh ) printf("nNeigh %d (%d %d %d)  oD %f %f %f %f %f %f %f nD %f %f %f %f %f %f\n", 
				     nNeigh, gi[0], gi[1], gi[2], ourDist, 
				     lDist[0],lDist[1],lDist[2],
				     lDist[3],lDist[4],lDist[5],
				     nDist[0],nDist[1],nDist[2],
				     nDist[3],nDist[4],nDist[5]
				     );
}


Coord_grid MSGrid::frac2SurfGrid(Coord_frac cf) {
  assert(cf.u() >= 0 &&
         cf.v() >= 0 &&
         cf.w() >= 0 );
  assert(cf.u() <= 1. &&
         cf.v() <= 1. &&
         cf.w() <= 1. );
  Coord_grid cg  = cf.coord_grid(surfGridSampling);
  return cg;
}
  /*
  assert(cf.u() >= surfMinCf.u() &&
         cf.v() >= surfMinCf.v() &&
         cf.w() >= surfMinCf.w() );
  // FIXME - rework this ; lengthsq  not an option if only
  // one coord close to lim
  if(cf.u() > surfMaxCf.u() || cf.v() > surfMaxCf.v() ||
     cf.w() > surfMaxCf.w() )
  {
    printf("frac2dfGrid out of bounds %f %f %f max %f %f %f \n",
	   cf.u(),cf.v(), cf.w(),
	   surfMaxCf.u(),surfMaxCf.v(),surfMaxCf.w() );
    if( ( surfMaxCf.u() - cf.u() ) < 1e-3)
      cf[0] = surfMaxCf.u();
    if( ( surfMaxCf.v() - cf.v() ) < 1e-3)
      cf[1] = surfMaxCf.v();
    if( ( surfMaxCf.w() - cf.w() ) < 1e-3)
      cf[2] = surfMaxCf.w();
    assert(cf.u() <= surfMaxCf.u() &&
	   cf.v() <= surfMaxCf.v() &&
	   cf.w() <= surfMaxCf.w() );
  }
  Coord_frac adjCf = (cf+surfOffsetCf);
  adjCf = Coord_frac(adjCf.u()/surfCfRange.u(), 
		     adjCf.v()/surfCfRange.v(), 
		     adjCf.w()/surfCfRange.w() );
  Coord_grid cg  = adjCf.coord_grid(surfGridSampling);
  Coord_frac test = surfGrid2CellFrac( cg) ;
  if(0)LOGSTR(("XX frac2SurfGrid cf %6.3f %6.3f %6.3f  adjCf %6.3f %6.3f %6.3f\n cg %d %d %d\n  test %6.3f %6.3f %6.3f \n Min %6.3f %6.3f %6.3f \n",  
	 cf.u(), cf.v(), cf.w(),
	 adjCf.u(), adjCf.v(), adjCf.w(),
	 cg.u(), cg.v(), cg.w(),
	 test.u(), test.v(), test.w(),
	 surfMinCf.u(), surfMinCf.v(), surfMinCf.w() ));

  if(!surfGridRange.in_grid(cg)) {
    printf("ERROR not in grid origCf %s adjCf %s cg %s gs %s\n", 
	   cf.format().c_str(), adjCf.format().c_str(),
	   cg.format().c_str(), surfGridSampling.format().c_str() );
    assert( surfGridRange.in_grid(cg));
  }
  return cg;
  */

Coord_frac MSGrid::surfGrid2CellFrac(Coord_grid cg) {
  if(!surfGridSampling.in_grid(cg)) {
    printf("surfGrid2CellFrac: cg %s\nsurfGS %s range %s-%s in_grid GS %d range %d\n", 
	   cg.format().c_str(),
	   surfGridSampling.format().c_str(),
	   surfGridRange.min().format().c_str(),
	   surfGridRange.max().format().c_str(),
	   surfGridSampling.in_grid(cg),
           surfGridRange.in_grid(cg) );
    assert(surfGridRange.in_grid(cg));
    assert(surfGridSampling.in_grid(cg));
  }
  //assert(surfGridSampling.in_grid(cg));
  Coord_frac cf  = cg.coord_frac(surfGridSampling);
  return cf;
}

void MSGrid::showSurfGridEntry(const GridIndex &pi) {
  SurfGridEntry ge( getSurfGridEntry( pi[0], pi[1], pi[2]));
  LOGSTR(("SurfGridEntry(%d,%d,%d) %s %5.2f\n", pi[0], pi[1], pi[2], 
	  stateStr(ge.state), ge.dist));
  for (int j=0; j < 3; j++) {
    for( int i = 0; i < 2; i++ ) {
      GridIndex ngi = pi.neigh(N6i[2*j+i]);
      if(ngi.isOffGrid()) {
	LOGSTR(("\toffgrid[%d,%d,%d]\n", ngi[0], ngi[1], ngi[2]));
        continue;
      }
      SurfGridEntry nge( getSurfGridEntry(ngi[0], ngi[1], ngi[2]));
      LOGSTR(("\t\t[%d,%d,%d] %s %5.2f\n", ngi[0], ngi[1], ngi[2],
	      stateStr(nge.state), nge.dist));
    }
  }
}



inline SurfGridEntry MSGrid::getSurfGridEntry( int i, int j, int k ) {
  /*
  FMTypeGrid::Accessor vAcc = vgrid->getAccessor();
  Int32Grid::Accessor iAcc = igrid->getAccessor();
  Coord xyz(i, j, k);
  if(isInner)
    xyz.reset(i+nXInner, j+nYInner, k+nZInner);
  GridEntry ge;
  ge.dist = vAcc.getValue(xyz);
  ge.state = (VState)iAcc.getValue(xyz);
  return ge;
  */
  //unsigned ix = k*nX*nY + j*nX + i;
  //assert(ix < nEntries);
  //GridEntry ge2= grid[ix];
  //printf("%f %d %f %d\n", ge.dist, ge.state, ge2.dist, ge2.state);
  //assert( ge.state == ge2.state);
  //assert( fabs(ge.dist - ge2.dist) < 1e6);
  //ge.dist = ge2.dist;
  //return grid[ix];
  assert(i*j*k >= 0 );
  unsigned ix = k*nX*nY + j*nX + i;
  //LOGSTR(("getGridEntry(%d,%d,%d) %d grid %p-%p nX %d nY %d  ix %d -> %p\n",
  //	  i, j, k, nEntries, grid, &grid[nEntries], nX, nY,  ix, &grid[ix] ));
  assert(ix >= 0 && ix < nEntries);
  //return grid[ix];
  SurfGridEntry ret;
  ret.dist=surfGridEntryDistV[ix];
  ret.state=(VState)surfGridEntryStateV[ix];
  return ret;
}

void MSGrid::setSurfGridEntry(SurfGridEntry &ge, int i, int j, int k ) {
  /*
  FMTypeGrid::Accessor vAcc = vgrid->getAccessor();
  Int32Grid::Accessor iAcc = igrid->getAccessor();
  Coord xyz(i, j, k);
  if(isInner)
    xyz.reset(i+nXInner, j+nYInner, k+nZInner);
  vAcc.setValue(xyz, ge.dist);
  iAcc.setValue(xyz, (int)ge.state);
  */
  unsigned ix = k*nX*nY + j*nX + i;
  assert(ix < nEntries);
  surfGridEntryDistV[ix] = ge.dist;
  surfGridEntryStateV[ix] = (unsigned char) ge.state;
}


  /*
  assert(surfGridSampling.in_grid(cg));
  Coord_frac cf  = cg.coord_frac(surfGridSampling);
  Coord_frac adjCf = Coord_frac(cf.u() * surfCfRange.u(), 
				cf.v() * surfCfRange.v(), 
				cf.w() * surfCfRange.w() );
  Coord_frac adj2Cf = (adjCf-surfOffsetCf);
  if(0)LOGSTR(("YY surfGrid2CellFrac cf %6.3f %6.3f %6.3f  adjCf %6.3f %6.3f %6.3f\n cg %d %d %d\n  adj2 %6.3f %6.3f %6.3f \n",
	 cf.u(), cf.v(), cf.w(),
	 adjCf.u(), adjCf.v(), adjCf.w(),
	 cg.u(), cg.v(), cg.w(),
	  adj2Cf.u(), adj2Cf.v(), adj2Cf.w()));

  return adj2Cf;
  */
//struct LStruct {  int i, j, k;};

void MSGrid::mkInitSurfGridPts() {
  //AtomListGrid alg(nX, nY, nZ);
  AtomListGrid alg(surfGridSampling.nu(),
                   surfGridSampling.nv(),
                   surfGridSampling.nw() );
  map<MAtom*, AtomDescr*> atmDescrMap; 
  int nAtom = 0;
  for ( int p = 0; p < mMol.size(); p++ ) {
    for ( int m = 0; m < mMol[p].size(); m++ ) {
      MMonomer &mono(mMol[p][m]);
      LOGSTR(("mono >%s< %s\n", mono.type().c_str(), mono.id().c_str() ));
      for ( int a = 0; a < mMol[p][m].size(); a++ ) 
	{
	if ( !mMol[p][m][a].is_null() ) {
	  nAtom++;
          MAtom& atom(mMol[p][m][a]);
	  // PDB cord
	  const Coord_orth atmCoOrig =  atom.coord_orth() ;
	  Coord_orth atmCo =  atmCoOrig-surfMinCo;
	  atmCo = nearZeroOrth(atmCo);
          LOGSTR(("atmCoOrig %s\natmCo %s\n", atmCoOrig.format().c_str(),
		  atmCo.format().c_str()));
          LOGSTR(("surfMaxCo %s\n fracOrig %s\nfracSurf %s\n", 
		  surfMaxCo.format().c_str(),
		  atmCoOrig.coord_frac(origCell).format().c_str(),  
		  atmCoOrig.coord_frac(surfCell).format().c_str() ));
          /*
	  need <= within tolerance
	  assert( atmCoOrig.x() <= surfMaxCo.x() &&
	          atmCoOrig.y() <= surfMaxCo.y() &&
	          atmCoOrig.z() <= surfMaxCo.z() );
	  */
	  assert( atmCo.x() >= 0 &&
	          atmCo.y() >= 0 &&
	          atmCo.z() >= 0 );
	  assert( atmCo.x() <= surfCell.a() &&
	          atmCo.y() <= surfCell.b() &&
	          atmCo.z() <= surfCell.c() );
          // frac coord relative to orig, "center", cell
	  Coord_frac atmCf = atmCo.coord_frac(surfCell);
	  atmCf = nearBoundsFrac(atmCf, surfCell);
          LOGSTR(("atmCf %s\n", atmCf.format().c_str()));
	  assert( atmCf.u() >= 0 &&
	          atmCf.v() >= 0 &&
	          atmCf.w() >= 0 );
	  assert( atmCf.u() <= 1. &&
	          atmCf.v() <= 1. &&
	          atmCf.w() <= 1. );
          //atmCf = atmCf.lattice_copy_unit();
          //LOGSTR(("atmCf unit %s\n", atmCf.format().c_str()));
	  //Coord_grid atmCg = atmCf.coord_grid(asuGrid);
	  Coord_grid atmCg = frac2SurfGrid(atmCf);
          LOGSTR(("atmCg  %s (%d,%d,%d)-(%d,%d,%d) ix %d\n", 
		  atmCg.format().c_str(),
		  surfGridRange.min().u(), surfGridRange.min().v(), surfGridRange.min().w(),
		  surfGridRange.max().u(), surfGridRange.max().v(), surfGridRange.max().w(),
		  surfGridSampling.index(atmCg) ));
	  LOGSTR(("XX (%f,%f,%f) (%f,%f,%f) (%d %d %d)\n",
		  atmCf.u(),
		  atmCf.v(),
		  atmCf.w(),
		  fmtype(surfGridSampling.nu()),
		  fmtype(surfGridSampling.nv()),
		  fmtype(surfGridSampling.nw()),
		  Util::intr( atmCf.u()*fmtype(surfGridSampling.nu())),
		  Util::intr( atmCf.v()*fmtype(surfGridSampling.nv())),
		  Util::intr( atmCf.w()*fmtype(surfGridSampling.nw()))
		  ));

	  char *atmNam =  strdup(atom.name().c_str());
	  char *elemNam =  strdup( atom.element().c_str() );
          float vdw= elemName2vdw( elemNam );
	  if(vdw<1 ) {
	    // FIXME workaround for hydrogens not recognized as such by clipper
	    LOGSTR(("elem '%s'\n", atom.element().c_str()));
	    if(elemNam[0] == 'H' && strlen(elemNam) > 1)
	      elemNam[1] = 0;
	    vdw= elemName2vdw( elemNam );
          }
	  assert ( vdw > 1. );
	  map<MAtom*, AtomDescr*>::iterator mapItr = atmDescrMap.find(&atom);
	  AtomDescr* atmDescr = 0;
          if(mapItr == atmDescrMap.end() ) {  
	    atmDescr = new AtomDescr( &atom, atmCf, atmCo, vdw );
	    atmDescrMap.insert( pair <MAtom*, AtomDescr*>( &atom, atmDescr ));
          } else {
	    assert(!"duplicate atom");
	    //atmDescr = mapItr->second;
	  }
	  // inner box dims: 2* vdw diameter + voxeldiag should reach
	  // all grid points
          // on both sides of surface
	  fmtype cubeR=2.*vdw+voxelDiag;

	  int nGUx = cubeR/gridResol[0]  + 1;
	  int nGUy = cubeR/gridResol[1]  + 1;
	  int nGUz = cubeR/gridResol[2]  + 1;
	  LOGSTR(("=>'%s' %d '%s' vdw %3.2f  dfcg(%d %d %d) cf(%f %f %f) co(%f %f %f ) nGu (%d %d %d = %d) cubeR %f box(%d-%d,%d-%d,%d-%d) (%5.2f-%5.2f, %5.2f-%5.2f, %5.2f-%5.2f)\n", 
		  atmNam, mono.seqnum(), elemNam, vdw, 
		  atmCg.u(), atmCg.v(), atmCg.w(),
		  atmCf.u(), atmCf.v(), atmCf.w(),
		  atmCo.x(), atmCo.y(), atmCo.z(),
		  nGUx, nGUy, nGUz, nGUx*nGUy*nGUz,  cubeR,
		  atmCg.u()-nGUx, atmCg.u()+nGUx,
		  atmCg.v()-nGUy, atmCg.v()+nGUy,
		  atmCg.w()-nGUz, atmCg.w()+nGUz,
		  atmCo.x()-cubeR, atmCo.x()+cubeR, 
		  atmCo.y()-cubeR, atmCo.y()+cubeR, 
		  atmCo.z()-cubeR, atmCo.z()+cubeR 
		  ));
          free(atmNam);
          free(elemNam);
	  // FIXME - atmCg +nGU should to enclose all interior: is +1 enough?
	  list<LStruct>ixList;
	  /*
          for(int k = atmCg.w()-nGUz; k <= atmCg.w()+nGUz; k++)
            for(int j = atmCg.v()-nGUy; j <= atmCg.v()+nGUy; j++)
              for(int i = atmCg.u()-nGUx; i <= atmCg.u()+nGUx; i++)
		{
		  LStruct lStruct;
		  lStruct.i=i-atmCg.u();
		  lStruct.j=j-atmCg.v();
		  lStruct.k=k-atmCg.w();
		  ixList.push_back(lStruct);
                }
	  */
	  list<LStruct>::iterator lit;
	  //printf("XX %d %d %d %d \n", ixList.size(), nGUx, nGUy, nGUz );
	  int rMax = max(max( nGUx, nGUy), nGUz );
	  ixList.clear();
	  Sphere(rMax, ixList);
	  //printf("XX %d %d %d %d \n", ixList.size(), nGUx, nGUy, nGUz );
	  for(lit=ixList.begin(); lit != ixList.end(); lit++){
	    LStruct lStruct = *lit; 
	    Coord_grid cg( atmCg.u()+lStruct.i, atmCg.v()+lStruct.j, 
			   atmCg.w()+lStruct.k );
	    // FIXME handle in iter range ?
	    if(!surfGridRange.in_grid(cg))
	      continue;
	    Coord_frac cf =  surfGrid2CellFrac(cg);
	    // faster ? Coord_frac cf  = cg.coord_frac(surfGridSampling);
	    Coord_orth co = cf.coord_orth(surfCell);
	    fmtype d = Coord_orth::length(co, atmCo);
	    if( d < cubeR )
	    {  
	      //alg.addAtom(lStruct.i, lStruct.j, lStruct.k, atmDescr);
	      alg.addAtom( atmCg.u()+lStruct.i, atmCg.v()+lStruct.j,  atmCg.w()+lStruct.k, atmDescr );
	      /*LOGSTR(("added %d %d %d co (%5.2f, %5.2f, %5.2f ) d %f %f %f\n", 
		      lStruct.i, lStruct.j, lStruct.k, co.x(), co.y(), co.z(), d, d-vdw, vdw));
	      */
	    } else
	      /*
	      LOGSTR(("skipped %d %d %d co (%5.2f, %5.2f, %5.2f ) d %f %f %f\n", 
		      lStruct.i, lStruct.j, lStruct.k, co.x(), co.y(), co.z(), d, d-vdw, vdw));
	      */
	    continue;
	  }
	}
      }
    }
  }
  LOGSTR(("alg maxIx %d used ix %d  atmDescr %d nAtm %d maxALGElem %d\n", 
	  alg.maxIndices(), alg.usedIndices(), atmDescrMap.size() , nAtom,
          maxALGElem ));
  int nSurfInsideGridPoints = 0;
  int nSurfOutsideGridPoints = 0;

  //for (int xIx = 0; xIx< nX; xIx++ ) 
  //  for (int yIx = 0; yIx< nY; yIx++ )     
  //    for (int zIx = 0; zIx< nZ; zIx++ )     
  for (int xIx = surfGridMin.u(); xIx<= surfGridMax.u(); xIx++ ) {
    for (int yIx = surfGridMin.v(); yIx<= surfGridMax.v(); yIx++ ) {
      for (int zIx = surfGridMin.w(); zIx<= surfGridMax.w(); zIx++ ) {

	ALGElem &ael( alg.elem(xIx,yIx,zIx) );
	if( ael.nElem  ) {
	  LOGSTR(("[%d,%d,%d]atmDescrSet %d\n", xIx, yIx, zIx, 
		  ael.nElem ));
	  int vSz = ael.nElem;
	  assert(vSz);
	  if(vSz > 1 )
	    assert( vSz == ael.elemV->size()+1 );
	  Coord_grid surfCg(xIx, yIx, zIx);
	  //Coord_frac surfCf = surfCg.coord_frac(asuGrid).lattice_copy_unit();
	  Coord_frac surfCf = surfGrid2CellFrac(surfCg);
          Coord_orth surfCo(surfCf.coord_orth(surfCell));
          //Coord_orth surfCo(surfCf.coord_orth(mol.cell()));
	  fmtype minDist = FMTYPE_MAX;
	  AtomDescr *minAtmDescr = 0;
	  // find atom vdw surface we are closest to (and inside)
	  for(int i = 0; i < vSz; i++ )
	  {
	    AtomDescr *atmDescr = 0;
	    if(i == 0 )
	      atmDescr = ael.first;
	    else
	      atmDescr = (*ael.elemV)[i-1];
	      
	    MAtom &atm(*atmDescr->atm);
	    //Coord_grid atmCg = atmDescr->cfUnit.coord_grid(asuGrid);
	    Coord_grid atmCg = frac2SurfGrid(atmDescr->cfUnit);
	    fmtype distAtm = Coord_orth::length(surfCo, atmDescr->atmCo);
	    fmtype dist2vdw = distAtm-atmDescr->vdw;
	    LOGSTR((" dist %f %f atm cg %d %d %d\n", distAtm,
		    dist2vdw, atmCg.u(), atmCg.v(), atmCg.w()    ));
	    // distance must be measured to vdw, which is atom-specific, not to
            // atom center
	    if(dist2vdw < minDist){
	      minDist = dist2vdw;
	      minAtmDescr = atmDescr;
            }
          }
	  assert(minAtmDescr);
	  Coord_grid atmCg = frac2SurfGrid(minAtmDescr->cfUnit);
	  if(minDist < 0) {
	    LOGSTR((" inside %d %d %d %5.2f atm cg %d %d %d\n", 
		    xIx, yIx, zIx, minDist, atmCg.u(), atmCg.v(), atmCg.w() ));
	    nSurfInsideGridPoints ++;
	  } else {
	    LOGSTR((" outside %d %d %d %5.2f atm cg %d %d %d\n", 
		    xIx, yIx, zIx, minDist, atmCg.u(), atmCg.v(), atmCg.w() ));
	    nSurfOutsideGridPoints ++;   
	  }
	  /*
    				xIx+dfGridXOff, 
				yIx+dfGridYOff,
				zIx+dfGridZOff, 
	  */
	  InitSurfGridPoint igp(
				xIx, 
				yIx,
				zIx 
				/*minAtmDescr->cfUnit
				minAtmDescr->vdw*/ );
          //_fmGrid not yet allocated
	  //assert(!igp.isOffGrid());
	  igp.dist = minDist;
	  igp.atmDescr = minAtmDescr;
	  igpVec.push_back(igp);
	}
      }
    }
  }
  LOGSTR(("surface Grid Points inside %d outside %d igpVec %d\n",  
	  nSurfInsideGridPoints, nSurfOutsideGridPoints, igpVec.size()  ));
  printf("surface Grid Points inside %d outside %d igpVec %d\n",  
	 nSurfInsideGridPoints, nSurfOutsideGridPoints, igpVec.size()  );
  /*
  for(MapIt_t mit =alg._map.begin(); mit != alg._map.end(); mit++) {
    ALGElem *el = (*mit).second;
    LOGSTR(("igp %d %d %d\n", el->nElem);
  }
  */
}

void MSGrid::heapInsert(fmtype dist, const GridIndex& newGi) {
  assert(!isnan(dist));
  heap.insert( pair<fmtype, GridIndex>(dist,newGi) );
  /*
  LOGSTR(("heap insert %f at %d %d size %d\n", dist, newGi[0], newGi[1], heap.size() ));
  checkHeap();
  for(HeapIt hi = heap.begin(); hi != heap.end(); hi++) {
    fmtype hd = (*hi).first;
    GridIndex hGi = (*hi).second;
    LOGSTR(("HEAP_I1[%d %d] %f\n", hGi[0], hGi[1], hd ));
  }
  */
  /*
  for(HeapIt hi = heap.begin(); hi != heap.end(); hi++) {
    fmtype hd = (*hi).first;
    GridIndex hGi = (*hi).second;
    LOGSTR(("HEAP_I3[%d %d] %f %d\n", hGi[0], hGi[1], hd, heap.count(hd) ));
  }
  */
}
void MSGrid::heapErase(HeapIt toErase) {
  /*
  fmtype dist = (*toErase).first;
  GridIndex oldGi = (*toErase).second;
  LOGSTR(("heap erase %f at %d %d size %d count %d\n", 
	 dist, oldGi[0], oldGi[1], heap.size(), heap.count(dist) ));
  for(HeapIt hi = heap.begin(); hi != heap.end(); hi++) {
    fmtype hd = (*hi).first;
    GridIndex hGi = (*hi).second;
    LOGSTR(("HEAP_E1[%d %d] %f %d\n", hGi[0], hGi[1], hd, heap.count(hd) ));
  }
  */
  heap.erase(toErase);
  /*
  for(HeapIt hi = heap.begin(); hi != heap.end(); hi++) {
    fmtype hd = (*hi).first;
    GridIndex hGi = (*hi).second;
    LOGSTR(("HEAP_E2[%d %d] %f\n", hGi[0], hGi[1], hd ));
  }
  */
}

void MSGrid::freezeVoxel(const GridIndex &gi, fmtype dist)  {
  SurfGridEntry ge(getSurfGridEntry( gi[0], gi[1], gi[2]) );
  LOGSTR(("freezeVoxel %s with dist %f %s\n", gi.toStr(), dist, stateStr(ge.state) ));
  ge.dist = dist;
  assert(ge.state != VState_FROZEN);
  ge.state = VState_FROZEN;
  setSurfGridEntry( ge, gi[0], gi[1], gi[2] );
}

void MSGrid::narrowVoxel(const GridIndex &gi, fmtype dist)  {
  SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]) );
  //LOGSTR(("narrowVoxel: %s with dist %f %s\n", gi.toStr(), dist, stateStr(ge.state) ));
  assert(ge.state != VState_NARROW);
  assert(ge.state != VState_FROZEN);
  ge.dist = dist;
  //LOGSTR(("%p %p %s\n", &ge.state, &grid[nEntries], stateStr(ge.state)));
  ge.state = VState_NARROW;
  //LOGSTR(("%p %p %s\n", &ge.state, &grid[nEntries], stateStr(ge.state)));
  //assert ((char*)&ge.state < (char*)&grid[nEntries]);
  //checkHeap();
  setSurfGridEntry( ge, gi[0], gi[1], gi[2] );
  heapInsert(dist ,gi);
}

void MSGrid::changeVoxelDist(const GridIndex &gi, fmtype dist)  {
  SurfGridEntry ge (getSurfGridEntry( gi[0], gi[1], gi[2]));
  assert(ge.state == VState_NARROW);
  /*LOGSTR(("changeVoxelDist (%d,%d,%d) %f -> %f %f\n", 
  	 gi[0], gi[1], gi[2], ge.dist, dist, ge.dist-dist ));
  */
  fmtype oldDist = ge.dist;
  ge.dist = dist;
  setSurfGridEntry(ge, gi[0], gi[1], gi[2] );
  // need to also update heap
  pair<HeapIt,HeapIt> range = heap.equal_range( oldDist );
  HeapIt toErase = heap.end();
  for( HeapIt it = range.first; it != range.second; it++) {
    GridIndex &rgi((*it).second);
    if( rgi == gi ) {
      //LOGSTR(("heap match %s\n", gi.toStr()));
      assert(toErase == heap.end() );     
      toErase = it;     
    }
  }
  assert(toErase != heap.end() );     
  if ( oldDist != (*toErase).first) {
    LOGSTR(("DIFF %f %f %f\n",  oldDist, (*toErase).first, 
	    oldDist- (*toErase).first ));
    assert( oldDist == (*toErase).first);
  }
  GridIndex oldGi = (*toErase).second;
  GridIndex newGi(oldGi[0], oldGi[1], oldGi[2]);
  heapErase( toErase);
  heapInsert( dist, newGi );
}
bool MSGrid::recompute(const GridIndex &pi, fmtype &max_sol)
{
  if(!recompute1or2(pi, max_sol, true)) {
    bool ret = recompute1or2(pi, max_sol, false);
    if(!ret) {
      LOG_ASSERT;
      showSurfGridEntry(pi);
      ret = recompute1or2(pi, max_sol, false);
      assert(ret);
    }
  }
  return true;
}
bool MSGrid::recompute1or2(const GridIndex &pi, fmtype &max_sol, bool use2)
{
  //vector<GridIndex> contribVec;
  //                  c b a
  fmtype coeff[3] = { -1, 0, 0 };
  //GridIndex N6i[6];
  //neighbors(pi, N6i);
  LOGSTR(("recompute for %s start \n", pi.toStr() ));
  //if(pi[0] == 0 && pi[1] == 5 && pi[2] == 5 )   assert(!"BINGO");
  int selDir[3];
  fmtype selVal[3];
  for (int j=0; j < 3; j++) {
    fmtype val1 = FMTYPE_MAX;
    fmtype val2 = FMTYPE_MAX;
    selVal[j] = -100;
    selDir[j] = -1;
    //GridIndex *contrib = 0;
    for( int i = 0; i < 2; i++ ) {

      GridIndex ngi = pi.neigh(N6i[2*j+i]);
      fmtype _val1;
      if( !ngi.isOffGrid() && state(ngi) == VState_FROZEN ) {
	_val1 = distance(ngi);
	//
	LOGSTR(("pi[%d,%d,%d] j %d i %d val %f ngi %d %d %d \n", 
	       pi[0], pi[1], pi[2], j, i, 
		_val1, ngi[0], ngi[1], ngi[2] ));
	//
	//assert(_val1 != -1);
	//assert(_val1 > 0);
	//if( (_val1 > 0 && val1 > 0 && _val1 < val1 ) || (_val1 < 0 && _val1 > val1 ) )
	if( _val1 < val1 ) 
	  {
	  LOGSTR(("recompute: val1 neighbor %s state %s dist %f \n", 
		  ngi.toStr(), stateStr(state(ngi)), _val1 ));
	  val1 = _val1;
	  selDir[j] = 2*j+i;
	  selVal[j] = val1;
	  I3 tmp = 2*N6i[2*j+i];
	  GridIndex ngi2 = pi.neigh(tmp);
	  VState vState2;
	  fmtype _val2;
	  if(ngi2.isOffGrid()) {
	    //LOGSTR(("recompute neighbor %s offgrid\n", ngi2.toStr().c_str() ));
	    _val2 = FMTYPE_MAX;
	    vState2  = VState_OFFGRID;
	    val2 = FMTYPE_MAX;
          } else {
	    _val2 = distance(ngi2);
	    vState2  = state(ngi2);
	    if( vState2 == VState_FROZEN && 
		( ( _val2 < val1  && val1 >= 0 ) ||
		  ( _val2 >= val1 && val1 <= 0 ) ) ) 
           {
	      LOGSTR(("recompute: val2 neighbor %s state %s dist %f \n", 
		      ngi2.toStr(), stateStr(state(ngi2)), _val2 ));
	      val2 = _val2;
	      selDir[j] = (2*j+i)+200;
	      selVal[j] = val2;
	    } else
	      val2 = FMTYPE_MAX;
          }
        }
      }
    }
    //if(contrib) {      contribVec.push_back( *contrib );    }
//#if HIGH_ACCURACY
#if 1
    if( use2 && val2 != FMTYPE_MAX) {
      fmtype tp = (1.0/3) * (4 *val1-val2);
      fmtype a = 9.0/4;
      coeff[2] += idx2_[j]*a;
      coeff[1] -= idx2_[j]*2 * a * tp;
      coeff[0] += idx2_[j]*a*sqrVal( tp ) ;
    } else
#endif
    if(val1 != FMTYPE_MAX) {
      coeff[2] += idx2_[j] ;
      coeff[1] -= idx2_[j]*2*val1;
      coeff[0] += idx2_[j]*sqrVal(val1);
    }
  }
  //showGridEntry(pi);
  max_sol = FMTYPE_MAX;
  fmtype sol[2] = { FMTYPE_MAX, FMTYPE_MAX};
  //
  fmtype a= coeff[2], b = coeff[1], c = coeff[0];
  LOGSTR(("pi[%d,%d,%d] a %f b %f c %f b^2 %f 4ac %f b^2-4ac %f\n selDir %d %d %d selVal %f %f %f\n", 
	  pi[0], pi[1], pi[2], a, b, c, b*b, 4*a*c, b*b-4*a*c, 
	  selDir[0], selDir[1], selDir[2], 
	  selVal[0], selVal[1], selVal[2] ));
  //
  if( int no_solutions = solve_quadric(coeff, sol) ) {
    if( no_solutions == 2 ) 
      max_sol = max( sol[1], sol[0] );
    else
      max_sol = sol[0];
    //assert( !isnan(max_sol));
    return true;
  } 
  //assert(!"complex sol");
  return false;
}


void MSGrid::neighbors(const GridIndex&gi, GridIndex (&nv)[6])
{
  for(int i = 0; i < 6; i++ )
    nv[i] = gi.neigh(N6i[i]);
}
void MSGrid::fmInit() {
  //igpVec.reserve(50000000);
  __debugFileNum = 2;
  mkInitSurfGridPts();
  allocSurfGrid();
  //grid = new SurfGridEntry[nEntries];
  __debugFileNum = 3;
  vector<InitSurfGridPoint> &iv(igpVec);
  // sort by abs value of distance from vdw surface
  sort(iv.begin(), iv.end(), igpCmpFunc);

  //FIXME leave as vector
  for(unsigned i = 0; i < iv.size(); i++) {
    InitSurfGridPoint igp = iv[i];
    //SurfGridEntry ge(getSurfGridEntry(igp[0], igp[1], igp[2]) );
    GridIndex gi =  igp.getGridIndex();
    SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
    assert(ge.state== VState_FAR);
    ge.state=VState_IGP;
    //ge.state=VState_FAR;
    ge.dist = igp.dist;
    setSurfGridEntry( ge,  gi[0], gi[1], gi[2]);
    //LOGSTR(("X1(%d %d %d) %s %f\n", igp[0], igp[1], igp[2], stateStr(ge.state), ge.dist));
    initList.push_back(igp);
  }

  list<InitSurfGridPoint>::iterator lit;
  
  fmtype minVoxelWidth = gridResol[0];
  minVoxelWidth = min(minVoxelWidth, gridResol[1]);
  minVoxelWidth = min(minVoxelWidth, gridResol[2]);
  
  LOGSTR(("fmInit %d voxelDiag %f minVoxelWidth %f\n", initList.size(), 
	  voxelDiag, minVoxelWidth ));
  printf("fmInit %u voxelDiag %f minVoxelWidth %f\n", 
	 initList.size(), voxelDiag, minVoxelWidth );
  // traverse from largest to smallest, looking for points that "contact" vdw
  for ( list<InitSurfGridPoint>::reverse_iterator rit = initList.rbegin(); 
	rit != initList.rend();	rit++) 
  {
    InitSurfGridPoint &(igp) = *rit;
    //SurfGridEntry ge(getSurfGridEntry( igp[0], igp[1], igp[2]) );
    GridIndex gi =  igp.getGridIndex();
    SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
    //LOGSTR(("X2(%d %d %d) %s %f\n", igp[0], igp[1], igp[2], stateStr(ge.state), ge.dist));
    if(ge.state == VState_IGP && fabs(igp.dist) < 1.5*voxelDiag) {
      flagSpanZero( igp);
    }
  }
  //
  // smallest to largest : points that contact those that span boundary
  for ( lit = initList.begin(); lit != initList.end();	lit++) {
    InitSurfGridPoint &(igp) = *lit;
    //SurfGridEntry ge(getSurfGridEntry( igp[0], igp[1], igp[2]) );
    GridIndex gi =  igp.getGridIndex();
    SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
    if(ge.state == VState_IGP && igp.dist >= 0 && igp.dist < 2*voxelDiag) {
      // igp may be a neighbor of a spanZero point
      extendSpanZero(igp);
    }
  }
  //  
  __debugFileNum = 4;
  int nDiscard = 0;
  for ( lit = initList.begin(); lit != initList.end();	lit++) {
    InitSurfGridPoint &(igp) = *lit;
    //SurfGridEntry ge(getSurfGridEntry( igp[0], igp[1], igp[2]) );
    GridIndex gi =  igp.getGridIndex();
    SurfGridEntry ge( getSurfGridEntry( gi[0], gi[1], gi[2]));
    if(ge.state == VState_INIT) {
      freezeVoxel(gi, igp.dist);
      igp.useForBoundary = true;
    } else if(ge.state == VState_IGP) {
      if( igp.dist < 0)
	freezeVoxel(gi, igp.dist);
      else {
	ge.state= VState_FAR;
	ge.dist = -2;
	setSurfGridEntry( ge,  gi[0], gi[1], gi[2]);
      }
      nDiscard++;
    }
  }
  // points flagged during initial sphere search but discarded
  printf("discard %d + used %u = %u\n", nDiscard, 
	 initList.size()-nDiscard, initList.size()  );

  for ( lit = initList.begin(); lit != initList.end();	lit++) {
    InitSurfGridPoint &(igp) = *lit;
    GridIndex gi = igp.getGridIndex();
    if(!igp.useForBoundary)  
      continue;
    //continue;
    //assert( state(igp) !=VState_IGP && state(igp) != VState_FAR );
    assert( state(gi) ==VState_FROZEN);
    LOGSTR(("narrow neighbors of (%d,%d,%d)\n", gi[0], gi[1], gi[2])); 
    showSurfGridEntry(gi);
    GridIndex nVec[6];
    neighbors(gi, nVec);
    for( int i = 0; i< 6; i++ ) {
      GridIndex &ngi(nVec[i]);
      if( nVec[i].isOffGrid() )
	continue;
      if( state(nVec[i])==VState_FROZEN ) {
	continue;  
      } 
      if( state(nVec[i])==VState_NARROW ) {
	continue;
      }
      fmtype dist;
      recompute(nVec[i], dist);
      narrowVoxel(nVec[i], dist);
    }
  }
  igpVec.clear();
}
void MSGrid::fmLoop() {
  int loopCount = 0;
  while(heap.size() > 0 ) {
    loopCount++;
    HeapIt hit = heap.begin();
    fmtype dist = (*hit).first;
    GridIndex gi = (*hit).second;
    LOGSTR(("loopCount %d (%d,%d,%d) dist %f\n", loopCount, gi[0], gi[1], gi[2], dist ));
    //checkHeap();
    freezeVoxel(gi, dist);
    //LOGSTR(("top heap erase %d %d\n", gi[0], gi[1] ));
    heapErase( hit );

    GridIndex nVec[6];
    neighbors(gi, nVec);
    for( int i = 0; i< 6; i++ ) {
      //checkHeap();
      if( !nVec[i].isOffGrid() && state(nVec[i]) != VState_FROZEN ) {
	fmtype maxSol;
	recompute( nVec[i], maxSol );
	if( state(nVec[i]) == VState_NARROW ) {
          //checkHeap();
    	  changeVoxelDist(nVec[i], maxSol);
          //checkHeap();
        } else {
	  narrowVoxel(nVec[i], maxSol);
          //checkHeap();
        }
      }
      //checkHeap();
    }
  }
  //LOGSTR(("heap size %d\n", heap.size()));
}

void MSGrid::writeMap() {
  Grid_sampling outMapGS(origSG, origCell, Resolution(dfResolution), /*rate*/1.5);
  outMap.init(origSG, origCell, outMapGS );
  printf("gs map %s surf %s\n", 
	 outMapGS.format().c_str(),surfGridSampling.format().c_str() );
  for( Xmap<fmtype>::Map_reference_coord outRef = outMap.first_coord(); 
       !outRef.last(); outRef = outRef.next() ) 
  {
    Coord_orth outCo = outRef.coord_orth();
    Coord_orth outCo2 = outCo-surfMinCo;
    Coord_frac surfCf = outCo2.coord_frac(surfCell);
    Coord_grid surfCg = surfCf.coord_grid(surfGridSampling);
    fmtype val = distance( GridIndex(surfCg.u(), surfCg.v(), surfCg.w() ) ); 
    outMap[outRef] = val;
    Coord_grid mapCg =  outRef.coord();
    /*
    printf("%f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f (%d %d %d) (%d %d %d) (%5.2f %5.2f %5.2f )\n", 
	   val, 
	   outCo.x(), outCo.y(), outCo.z(), 
	   outCo2.x(), outCo2.y(), outCo2.z(), 
           surfCg.u(), surfCg.v(), surfCg.w(),
           mapCg.u(),  mapCg.v(),  mapCg.w(),
	   surfCf.u(), surfCf.v(), surfCf.w() );
    */
    if(!surfGridRange.in_grid(surfCg) || !surfGridSampling.in_grid(surfCg)) {
      printf("range %s-%s\n", 
	     surfGridRange.min().format().c_str(),
	     surfGridRange.max().format().c_str() );
      assert(surfGridRange.in_grid(surfCg));
      assert(surfGridSampling.in_grid(surfCg));
    }
    /*
    */
  }
  mapToFile(  outMap,  string("dfclp.map" ) );
  
  //NXmap<fmtype> outMap( mol.cell(), gridSam, Grid_range(clipper::Grid(nX,nY,nZ), Coord_frac(0,0,0), Coord_frac(1.,1.,1.)));
#if 0
  Coord_grid asuMin =  outMap.grid_asu().min();
  Coord_grid asuMax =  outMap.grid_asu().max();
  LOGSTR(("writemap outMapGrid: min %d %d %d max %d %d %d Nx 0..%d, Ny 0..%d, Nz 0..%d\n", 
	  asuMin.u(), asuMin.v(), asuMin.w(),
	  asuMax.u(), asuMax.v(), asuMax.w(),
          nX-1, nY-1, nZ-1));
  int i, j, k, ix = 0;
#if 0
  for(i=0;i<nX;i++) {
    for(j=0;j<nY;j++) {
      for(k=0;k<nZ;k++) {
      //LOGSTR(("[%d] ", i));
	Coord_grid cg(i, j, k);
        if(mapAsuGrid.in_grid(cg)) {
	  GridEntry ge (getGridEntry( i, j, k ) );
	  fmtype dist = ge.dist;
	  //Coord_frac cf = cg.coord_frac( asuGrid );
	  Coord_frac cf = dfGrid2CellFrac(cg);
	  Coord_orth co = cf.coord_orth( cell );
	  //printf("%d %d %d %f %6.3f %6.3f %6.3f \n", i, j, k, dist,
	  //     co.x(), co.y(), co.z() );
	//if( dist == FMTYPE_MAX)
	//  dist = 0;
	  outMap.set_data(cg, dist);
	//wGrid[ix] = ge.dist;
	//dObj[ix] = 1;
	//if( nCol > 30 )
	//  distM[ix] = 0;
	  ix++;
	}
      }
    }
    //LOGSTR(("\n"));
  }
#endif
  LOGSTR(("%d map values written\n", ix));
#if 0
  clipper::CCP4MAPfile ccp4Map;
  ccp4Map.open_write( "dfclp.map" );      // write map
  ccp4Map.export_xmap( outMap );
  ccp4Map.close_write();
#endif
  mapToFile(  outMap,  string("dfclp.map" ) );
#endif
  //testMap("dfclp.map");
}


int main(int argc, char **argv) { 
  __debugFileNum = 1; // 0 reserved for asserts

  char *pdbFN = NULL;
  float dfResolution=DefaultDFRes;
  --argc,++argv;
  while( argc ) {
    if(argc && !strcmp(*argv,"-c")) {
      --argc,++argv;
      pdbFN = *argv;
    } else if(argc && !strcmp(*argv,"-r")) {
      --argc,++argv;
      dfResolution = atof(*argv);
    } else if(argc && !strcmp(*argv,"-d")) {
      __debugEnabled = true;
    } else {
      assert(!"invalid argument");
    }
    --argc,++argv;
  }
  if(!pdbFN) {
    fprintf(stderr,"usage : dfg -c PDB_FILE [-r resolution] [-d]\n ");
    return -1;
  }
  MSGrid msGrid;
  _fmGrid = &msGrid;

  msGrid.init( pdbFN, dfResolution );
  //__debugFileNum  2 - 4
  msGrid.fmInit();
  __debugFileNum = 5;
  printf("start loop\n");
  msGrid.fmLoop();
  __debugFileNum = 4;
  msGrid.writeMap();
  return 0;
}

void mapToFile(  Xmap<fmtype> xmap,  string filename ) {

  const char* title = "From clipper Xmap                                                               ";
  char symop[80];
  int grid[3], orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  float cp[6];

  int spg = xmap.spacegroup().descr().spacegroup_number();  // spacegroup
  // get axis order
  switch ( spg ) {
  case 1:  case 2:  case 3:  case 4:
  case 10: case 16: case 17: case 18:
  case 20: case 21: case 23:
    orderfms[0] = 2; orderfms[1] = 1; orderfms[2] = 3; break;
  default:
    orderfms[0] = 3; orderfms[1] = 1; orderfms[2] = 2; break;
  }
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  // grids
  for ( int i = 0; i < 3; i++ ) {
    grid[i] = xmap.grid_sampling()[i];
    gfms0[orderxyz[i]] = xmap.grid_asu().min()[i];
    gfms1[orderxyz[i]] = xmap.grid_asu().max()[i];
  }
  LOGSTR(("mapToFile asu min %d %d %d\n", 
	 xmap.grid_asu().min()[0],
	 xmap.grid_asu().min()[1],
	  xmap.grid_asu().min()[2]));
  LOGSTR(("mapToFile asu max %d %d %d\n", 
	 xmap.grid_asu().max()[0],
	 xmap.grid_asu().max()[1],
	  xmap.grid_asu().max()[2]));
  LOGSTR(("XXX gfms0 %d %d %d\n", gfms0[0], gfms0[1], gfms0[2]));
  LOGSTR(("XXX gfms1 %d %d %d\n", gfms1[0], gfms1[1], gfms1[2]));
  Cell_descr cd = xmap.cell().descr();
  cp[0] = cd.a(); cp[3] = cd.alpha_deg();
  cp[1] = cd.b(); cp[4] = cd.beta_deg ();
  cp[2] = cd.c(); cp[5] = cd.gamma_deg();
  for ( int i = 0; i < 3; i++ ) dim[i] = gfms1[i] - gfms0[i] + 1;

  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_WRONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: export_xmap - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_set_cell( file, cp );
  CMap_io::ccp4_cmap_set_grid( file, grid );
  LOGSTR(("XXX grid %d %d %d\n", grid[0], grid[1], grid[2]));
  CMap_io::ccp4_cmap_set_order( file, orderfms );
  LOGSTR(("XXX dim %d %d %d\n", dim[0], dim[1], dim[2]));
  CMap_io::ccp4_cmap_set_dim( file, dim );
  CMap_io::ccp4_cmap_set_origin( file, gfms0 );

  CMap_io::ccp4_cmap_set_spacegroup( file, spg );
  CMap_io::ccp4_cmap_set_title( file, title );
  CMap_io::ccp4_cmap_set_datamode( file, 2 );

  // write symops
  for ( int i = 0; i < xmap.spacegroup().num_symops(); i++ ) {
    String strop = xmap.spacegroup().symop(i).format();
    for ( int j = 0; j < 80; j++ ) symop[j] = ' ';
    for ( int j = 0; j < strop.length(); j++ ) symop[j] = strop[j];
    CMap_io::ccp4_cmap_set_symop( file, symop );
  }

  // write the map data
  int n0 = (gfms1[0]-gfms0[0]+1);
  int n1 = n0 * (gfms1[1]-gfms0[1]+1);
  std::vector<float> section( n1 );
  int index, g[3];
  Xmap_base::Map_reference_coord x( xmap );
  for ( g[2] = gfms0[2]; g[2] <= gfms1[2]; g[2]++ ) {
    index = 0;
    for ( g[1] = gfms0[1]; g[1] <= gfms1[1]; g[1]++ ) {
      for ( g[0] = gfms0[0]; g[0] <= gfms1[0]; g[0]++ ) {
	x.set_coord( Coord_grid( g[orderxyz[0]], g[orderxyz[1]],
				 g[orderxyz[2]] ) );
        section[ index++ ] = float( xmap[x] );
      }
    }
    CMap_io::ccp4_cmap_write_section( file, &section[0] );
  }

  // done
  CMap_io::ccp4_cmap_close( file );
}

