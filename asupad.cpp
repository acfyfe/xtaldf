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

/*  asupad : preprocessing utility for distance field computation

   Given an input PDB file, refer all atoms to the ASU then generate an output
   PDB ("dfgen.pdb") consisting of either (1) the asu expanded to a specified
   width, given by the "asuPad" argument as a fractional coordinate, or (2) 
   the symmetry-expanded  atoms that  occupy a given quadrilateral box 
   given by "boxSize" and padded by  "borderMult".

   For example "asupad -ap 1. -c foo.pdb" gives the 26-neighbor covering 
   of the starting asu.

   For DF processing, this is usually run with "-ap 0.5"

 */
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <set>
#include <list>
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



using namespace std;
using namespace clipper;
using namespace ::mmdb;
typedef float fmtype; // aka floating type for "fast marching" computation

class ASUAtom {  // wrapper for item in clipper  MiniMol hierarchy 
public:
  MAtom &atom;
  MMonomer &mono;
  MPolymer &poly;
  Residue *res;
  MAtomIndex atmIx;
  ASUAtom(MAtom &aArg, MMonomer &mArg, MPolymer &pArg, Residue *resArg,
	    MAtomIndex atmIxArg) 
    : atom(aArg), mono(mArg), poly(pArg), res(resArg), atmIx(atmIxArg)
  {
  }
};

class CoordOrthComp { // put an ordering on 3D coordinates
   public:
      bool  operator()  ( const Coord_orth l, const Coord_orth r) const { 
	const float epsi = 1e-4;
        if(fabs(l.x()-r.x()) > epsi) {
	  return l.x() < r.x();
        } else  if(fabs(l.y()-r.y()) > epsi) {
	  return l.y() < r.y();
        } else  if(fabs(l.z()-r.z()) > epsi) {
	  return l.z() < r.z();
        } 
        return false;
      }
};

list<Coord_frac> &genCellTrans(int max) {
  list <Coord_frac> &ret(* new list<Coord_frac> );
  for( int i = -max; i<=max; i++) {
    for( int j = -max; j<=max; j++) {
      for( int k = -max; k<=max; k++) {
	ret.push_back(Coord_frac( i, j, k));
      }
    }
  }
  return ret;
}

// clipper ouput methods don't provide sufficient access to  file details, hence mmdb
void writeAtom(MAtom &atom, const Coord_orth &aCo, Residue *res) {
  PAtom cAtom = new ::mmdb::Atom;
  //cAtom->SetAtomName ( atom.name().c_str() );// it has to be a PDB name!
  cAtom->SetAtomName (atom.name().c_str() );// it has to be a PDB name!
  if(atom.id().find(":") != string::npos) {
    AltLoc altLoc;
    string tmp = atom.id().substr(atom.id().find(":")+1);
    strncpy(altLoc, tmp.c_str(), sizeof(altLoc)-1 );
    strncpy(cAtom->altLoc, altLoc, sizeof(cAtom->altLoc)-1);
  }  
  cAtom->SetCoordinates(aCo.x(), aCo.y(), aCo.z(), 
			atom.occupancy(), Util::u2b(atom.u_iso()));
  cAtom->SetElementName(atom.element().c_str());
  cAtom->SetResidue(res); 
  res->AddAtom(cAtom);
}

// along with main output file, a couple of "guide" pdbs are emitted to help
// visualize dimensions of the asu and padded box 
void writeBox(Coord_frac minCf, Coord_frac maxCf, 
	      const Spacegroup &sg, Cell &cell, const char *fName) 
{
  Coord_frac cf000 = Coord_frac(minCf.u(), minCf.v(), minCf.w() );
  Coord_frac cf010 = Coord_frac(minCf.u(), maxCf.v(), minCf.w() );
  Coord_frac cf100 = Coord_frac(maxCf.u(), minCf.v(), minCf.w() );
  Coord_frac cf110 = Coord_frac(maxCf.u(), maxCf.v(), minCf.w() );
  Coord_frac cf001 = Coord_frac(minCf.u(), minCf.v(), maxCf.w() );
  Coord_frac cf011 = Coord_frac(minCf.u(), maxCf.v(), maxCf.w() );
  Coord_frac cf101 = Coord_frac(maxCf.u(), minCf.v(), maxCf.w() );
  Coord_frac cf111 = Coord_frac(maxCf.u(), maxCf.v(), maxCf.w() );

  Coord_orth c000 = cf000.coord_orth( cell );
  Coord_orth c010 = cf010.coord_orth( cell );
  Coord_orth c100 = cf100.coord_orth( cell );
  Coord_orth c110 = cf110.coord_orth( cell );
  Coord_orth c001 = cf001.coord_orth( cell );
  Coord_orth c011 = cf011.coord_orth( cell );
  Coord_orth c101 = cf101.coord_orth( cell );
  Coord_orth c111 = cf111.coord_orth( cell );
  Coord_frac cfV[]={ cf000, cf010, cf100, cf110,
                     cf001, cf011, cf101, cf111 };
  Coord_orth coV[]={ c000, c010, c100, c110,
                     c001, c011, c101, c111 };
  int nCo=sizeof(coV)/sizeof(coV[0]);
  MiniMol outMM=MiniMol(sg,cell);
  MModel outMMod;
  int newMonoSeqNum = 0;
  MPolymer newPol;
  newPol.set_id(string(1, 'B'));
  string id(" SR ");
  for(int i = 0; i < nCo; i++) {
    MAtom &atom(*new MAtom());
    atom.set_id(id);
    atom.set_name(id);
    printf("XX [%d] %5.2f %5.2f %5.2f\n", i, cfV[i].u(), cfV[i].v(), cfV[i].w() );
    atom.set_coord_orth(coV[i]);
    atom.set_occupancy(Util::u2b(30.));
    atom.set_u_iso(1);
 
    MMonomer &nMono( *new MMonomer);
    nMono.set_id(id);
    nMono.set_type(id);
    nMono.set_seqnum(++newMonoSeqNum);
    nMono.insert(atom);
    newPol.insert(nMono);
  }
  outMMod.insert(newPol);
  outMM.model() = outMMod;
  MMDBManager cmmdbOF;
  MMDBfile *mmdbOF = static_cast<MMDBfile*>(&cmmdbOF);
  mmdbOF->export_minimol(outMM);
  mmdbOF->write_file(fName);
}

inline bool approxEq( float x, float y ) {
  const float epsi = 1e-4;
  return (fabs(x-y) <= epsi);
}

inline bool approxLe( float x, float y ) {
  if( x <= y  )
    return true;
  else  return  approxEq( x, y );
}

inline bool approxGe( float x, float y ) {
  if( x >= y  )
    return true;
  else  return  approxEq( x, y );
}


/* use either boxSize and borderMult (eg for octamer extension) or asupad (for df)
 */
void genPaddedAsuPdb( const char *molFileName,   fmtype borderMult,
		      fmtype boxSize, fmtype  asuPad, bool debug=true )
{

  MiniMol mmMol;
  MMDBfile molMFile;
  MMDBManager& mmdb = static_cast<MMDBManager&>(molMFile);
  // FIXME  check how permissive these need to be to process entire PDB archive
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  mmdb.SetFlag( mmdbflags );
  molMFile.read_file( molFileName );
  if(debug ) {
    printf("mmdb read: isSpaceGroup %d %s crystready %d \n", mmdb.isSpaceGroup(),  
	   mmdb.GetSpaceGroup(), mmdb.CrystReady() );
    molMFile.spacegroup().debug();
  }
  assert( mmdb.isSpaceGroup() );
  molMFile.import_minimol( mmMol );
  const Spacegroup &sg( mmMol.spacegroup() );
  Cell cell = mmMol.cell();

  Coord_frac asuMax = sg.asu_max();
  Coord_frac asuMin = sg.asu_min();
  printf("asuMax        %s\n", asuMax.format().c_str() );
  printf("asuMin        %s\n", asuMin.format().c_str() );
  Coord_orth boxMin, boxMax;
  if(asuPad == 0 ) {
    boxMin = Coord_orth(-boxSize,-boxSize,-boxSize);
    boxMax = Coord_orth(boxSize,boxSize,boxSize);
    for ( int p = 0; p < mmMol.size(); p++ ) { // polymers/chains
      MPolymer &poly(mmMol[p]);
      for ( int m = 0; m < mmMol[p].size(); m++ ) { // monomers
	MMonomer &mono(mmMol[p][m]);
	for ( int a = 0; a < mmMol[p][m].size(); a++ ) { //atoms
	  MAtom& atom(mmMol[p][m][a]);
	  const Coord_orth &co(atom.coord_orth());
	  if(co.x() < boxMin.x()) boxMin= Coord_orth( co.x(), boxMin.y(), boxMin.z());
	  if(co.y() < boxMin.y()) boxMin= Coord_orth( boxMin.x(), co.y(), boxMin.z());
	  if(co.z() < boxMin.z()) boxMin= Coord_orth( boxMin.x(), boxMin.y(), co.z());

	  if(co.x() > boxMax.x()) boxMax= Coord_orth( co.x(), boxMax.y(), boxMax.z());
	  if(co.y() > boxMax.y()) boxMax= Coord_orth( boxMax.x(), co.y(), boxMax.z());
	  if(co.z() > boxMax.z()) boxMax= Coord_orth( boxMax.x(), boxMax.y(), co.z());
	}
      }
    }
    Coord_orth boxBorder = 0.1*(boxMax-boxMin);
    boxMax = boxMax + borderMult*boxBorder;
    boxMin = boxMin - borderMult*boxBorder;
  } else {
    boxMin = Coord_frac(asuMin.u()-asuPad, 
			asuMin.v()-asuPad, 
			asuMin.w()-asuPad ).coord_orth(cell);
    boxMax = Coord_frac(asuMax.u()+asuPad, 
			asuMax.v()+asuPad, 
			asuMax.w()+asuPad ).coord_orth(cell);
  }
  Coord_frac boxMinCf=boxMin.coord_frac(cell);
  Coord_frac boxMaxCf=boxMax.coord_frac(cell);
  if(debug ) {
    printf("box %s %s\n%s %s\n", 
	   boxMin.format().c_str(),
	   boxMinCf.format().c_str(),
	   boxMax.format().c_str(),
	   boxMaxCf.format().c_str() );
  }
  writeBox(boxMinCf, boxMaxCf, sg, cell, "box.pdb"); 
  writeBox(asuMin, asuMax, sg, cell, "asubox.pdb");
  
  MMDBManager cMMDB;
  cMMDB.set_cell(mmMol.cell());
  cMMDB.set_spacegroup( mmMol.spacegroup() );
  Model* model = new Model;
  cMMDB.AddModel(model);

  Chain* chain = 0;
  Residue* res = 0;
  fmtype rad_ = 1.0;
  // use nonbond grid spacing to 1A resol
  //Grid_sampling grid = Grid_sampling( Util::intc(1.0/(rad_*cell.a_star())), Util::intc(1.0/(rad_*cell.b_star())), Util::intc(1.0/(rad_*cell.c_star())));
  //printf("grid %d %d %d = %d\n", grid.nu(), grid.nv(), grid.nw(),
  //	 grid.nu()* grid.nv()* grid.nw() );
  //map<int, list<MAtomIndexSymmetry>* > gridMap;
  list<pair<ASUAtom*,Coord_frac> > asuCoL;
  for ( int p = 0; p < mmMol.size(); p++ ) { // polymers/chains
    MPolymer &poly(mmMol[p]);
    printf("XX minimol chain '%s'\n", mmMol[p].id().c_str() );  
    chain = new Chain;
    chain->SetChainID ( mmMol[p].id().c_str() );  
    model->AddChain(chain);
    for ( int m = 0; m < mmMol[p].size(); m++ ) { // monomers
      MMonomer &mono(mmMol[p][m]);
      // FIXME monomers can be encoded as 0 or negative, see 3G8T - legit?
      //assert(mono.seqnum());
      res = new Residue;
      res->seqNum = mono.seqnum();
      res->SetResName(mono.type().c_str());
      chain->AddResidue ( res );
      for ( int a = 0; a < mmMol[p][m].size(); a++ ) { //atoms
        MAtom& atom(mmMol[p][m][a]);
	ASUAtom *ad = new ASUAtom(atom, mono,poly, res, MAtomIndex(p, m, a));
	const Coord_orth &aCo( atom.coord_orth());
	const Coord_frac &aCf( aCo.coord_frac(cell));
	Coord_frac cfz(0, 0, 0);
	//Coord_frac aCf0 = aCf.symmetry_copy_near( msg, cell, cfz );
	Coord_frac aCfU = aCf.lattice_copy_unit();
	Coord_orth aCo_U  = aCfU.coord_orth(cell);
        double d1 = ( aCf - cfz ).lengthsq( cell );
        double d2 = ( aCfU - cfz ).lengthsq( cell );
	if( debug) 
        {
          Coord_frac aCf2= Coord_frac(cell.matrix_frac()* aCo);
	  printf("atom[%d] %s(%d) co (%4.1f,%4.1f,%4.1f) cf (%4.2f,%4.2f,%4.2f) d %4.2f "
		 "near-0 (%4.2f,%4.2f,%4.2f) d %4.2f  unit (%4.2f,%4.2f,%4.2f) \n"
		 "aCf2 (%4.2f,%4.2f, %4.2f) \n",
		 a, atom.id().c_str(), mono.seqnum(),
		 aCo.x(), aCo.y(), aCo.z(), 
		 aCf.u(), aCf.v(), aCf.w(),  d1,
		 aCfU.u(),aCfU.v(),aCfU.w(), d2,
		 aCo_U.x(),aCo_U.y(),aCo_U.z(),
		 aCf2.u(), aCf2.v(), aCf2.w() );
	}
        int s;
	Coord_frac tCfU;
	for ( s = 0; s < sg.num_symops(); s++ ) {
	  tCfU = sg.symop(s) * aCfU;
	  tCfU = tCfU.lattice_copy_unit();
	  bool inASU=false;
	  if(  approxLe(asuMin.u(),tCfU.u()) && 
               approxGe(asuMax.u(),tCfU.u()) &&
	       approxLe(asuMin.v(),tCfU.v()) && 
	       approxGe(asuMax.v(),tCfU.v()) &&
	       approxLe(asuMin.w(),tCfU.w()) && 
	       approxGe(asuMax.w(),tCfU.w())  ) 
	    {
	      inASU=true;
	    }
	    if( debug) {
	      const Coord_orth tCo = tCfU.coord_orth(cell);
              if(debug )
	      printf("\t%s inASU %d  (%4.2f,%4.2f,%4.2f) co (%4.1f,%4.1f,%4.1f)\n",
	           sg.symop(s).format().c_str(), inASU,
	           tCfU.u(),tCfU.v(),tCfU.w(),
	           tCo.x(), tCo.y(), tCo.z() );
	    }
	    if(inASU)
	      break;
	}
	if( s == sg.num_symops() ) {
	  printf("ERROR : %s(%d) cell %s co %s\n"
		 " cf %s cf-near %s\n", 
		 atom.id().c_str(), mono.seqnum(), cell.format().c_str(),
		 aCo.format().c_str(), 
		 aCf.format().c_str(), 
		 aCfU.format().c_str() );
	  for ( s = 0; s < sg.num_symops(); s++ ) {
	    tCfU = sg.symop(s) * aCfU;
	    printf("[%d]%s  tcfNear %s \n", 
		   s, sg.symop(s).format().c_str(), 
		   tCfU.format().c_str() );
	  }
	  assert( s < sg.num_symops() );
	}
	aCo_U = tCfU.coord_orth(cell);
	asuCoL.push_back(pair<ASUAtom*,Coord_frac>(ad, tCfU));
      }
    }
  }
  Coord_frac cellTrans[] =  { 
    Coord_frac( 0, 0, 0), 
    Coord_frac( 1, 0, 0), 
    Coord_frac(-1, 0, 0), 
    Coord_frac( 0, 1, 0), 
    Coord_frac( 1, 1, 0), 
    Coord_frac(-1, 1, 0), 
    Coord_frac( 0,-1, 0), 
    Coord_frac( 1,-1, 0), 
    Coord_frac(-1,-1, 0), 

    Coord_frac( 0, 0, 1), 
    Coord_frac( 1, 0, 1), 
    Coord_frac(-1, 0, 1), 
    Coord_frac( 0, 1, 1), 
    Coord_frac( 1, 1, 1), 
    Coord_frac(-1, 1, 1), 
    Coord_frac( 0,-1, 1), 
    Coord_frac( 1,-1, 1), 
    Coord_frac(-1,-1, 1), 

    Coord_frac( 0, 0,-1), 
    Coord_frac( 1, 0,-1), 
    Coord_frac(-1, 0,-1), 
    Coord_frac( 0, 1,-1), 
    Coord_frac( 1, 1,-1), 
    Coord_frac(-1, 1,-1), 
    Coord_frac( 0,-1,-1), 
    Coord_frac( 1,-1,-1), 
    Coord_frac(-1,-1,-1) 
  };
  int nCellTrans = sizeof(cellTrans)/sizeof(cellTrans[0]);
  
  map< Coord_orth, Coord_orth, CoordOrthComp > uniqMap;
  int uniqCount = 0;
  
  list<Coord_frac> &cellTransL(genCellTrans(2));
  for( list<pair<ASUAtom*,Coord_frac> >::iterator asuIt=asuCoL.begin(); 
       asuIt != asuCoL.end(); asuIt++)
  {
    ASUAtom &ad(*(asuIt->first));
    Coord_frac &asuCf(asuIt->second);
    bool inASU=(approxLe(asuMin.u(),asuCf.u()) && 
                approxGe(asuMax.u(),asuCf.u()) &&
	        approxLe(asuMin.v(),asuCf.v()) && 
	        approxGe(asuMax.v(),asuCf.v()) &&
	        approxLe(asuMin.w(),asuCf.w()) && 
	        approxGe(asuMax.w(),asuCf.w())  ) ;
    /*
    Coord_frac vDiff1 = asuMax-asuCf;
    Coord_frac vDiff2 = asuCf-asuMin;
    printf("%s %s %s %s %s\n", 
	   asuMin.format().c_str(),
	   asuMax.format().c_str(),
	   asuCf.format().c_str(),
	   vDiff1.format().c_str(), 
	   vDiff2.format().c_str()
	   );
    */
    assert( inASU);
    for( list<Coord_frac>::iterator ctIt = cellTransL.begin();
	 ctIt != cellTransL.end(); ctIt++ )
    {
      Coord_frac &ct(*ctIt);
      for ( int s = 0; s < sg.num_symops(); s++ ) {
	Coord_frac tCfU = sg.symop(s) * asuCf;
	tCfU = tCfU.lattice_copy_unit();
	//tCfU = tCfU+cellMult*cellTrans[transN];
	tCfU = tCfU+ct;
        /*
	bool inASU=
	  ( tCfU.u()   >= asuMin.u()  && tCfU.u() <= asuMax.u() ) &&
	  ( tCfU.v()   >= asuMin.v()  && tCfU.v() <= asuMax.v() ) &&
	  ( tCfU.w()   >= asuMin.w()  && tCfU.w() <= asuMax.w() );
	assert( inASU);
        */
	Coord_orth tCo_U = tCfU.coord_orth(cell);
	bool useBoxLimit = asuPad == 0;
        bool useTheseCoord = false;
	if( useBoxLimit ) {
          if ( ( tCo_U.x()   >= boxMin.x() && tCo_U.x() <= boxMax.x() ) &&
	       ( tCo_U.y()   >= boxMin.y() && tCo_U.y() <= boxMax.y() ) &&
	       ( tCo_U.z()   >= boxMin.z() && tCo_U.z() <= boxMax.z() ) )
          {
	    useTheseCoord = true;
	  }
	} else {
          if ( ( tCfU.u()   >= boxMinCf.u() && tCfU.u() <= boxMaxCf.u() ) &&
	       ( tCfU.v()   >= boxMinCf.v() && tCfU.v() <= boxMaxCf.v() ) &&
	       ( tCfU.w()   >= boxMinCf.w() && tCfU.w() <= boxMaxCf.w() ) )
          {
	    useTheseCoord = true;
	  }
        }
        /* in symmetry-expansion, multiple operations may map an atom to a given coord
           yet we only want to write it once : identify duplicates
	*/
        if( useTheseCoord ) {
	  if( uniqMap.count( tCo_U ) == 0) {
	    uniqMap[tCo_U] = tCo_U;
	    writeAtom(ad.atom, tCo_U, ad.res);
	  }
        }
      }
    }
  }
  asuCoL.clear(); 
  
  cMMDB.PDBCleanup ( PDBCLEAN_SERIAL | PDBCLEAN_INDEX );
  string pdbName = "dfgen.pdb";
  ERROR_CODE writeRC = cMMDB.WritePDBASCII( const_cast<char*>(pdbName.c_str()) );
  if(writeRC != 0) {
    printf("chsolv PDB write err %s\n", GetErrorDescription(writeRC) );
  }
}

int main(int argc, char**argv) {
  --argc,++argv;
  char *pdbFN = NULL;
  // use either boxsize borderMult (for making symm equiv beehives) or asupad
  // for dfgen
  float boxSize = 1000;
  float borderMult=7.5;
  float asuPad=0;
  while( argc ) {
    if(argc && !strcmp(*argv,"-c")) {
      --argc,++argv;
      pdbFN = *argv;
    } else if(argc && !strcmp(*argv,"-bm")) {
      --argc,++argv;
      borderMult = atof(*argv);
    } else if(argc && !strcmp(*argv,"-bs")) {
      --argc,++argv;
      boxSize = atof(*argv);
    } else if(argc && !strcmp(*argv,"-ap")) {
      --argc,++argv;
      asuPad = atof(*argv);
    }
    --argc,++argv;
  }
  assert(pdbFN);
  genPaddedAsuPdb( pdbFN, boxSize,  borderMult, asuPad );
}
