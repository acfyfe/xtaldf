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

/*  dfsel : given a map encoding  the distance field from the molecular surface
            of an atomic model and either (1) a pdb file or  (2) an additional map
            emit either a PDB or map which includes the subset of (1) or (2)
            that lies between the distance field range specified by "d1" and 
	    "d2" arguments
*/
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <set>
#include <list>
#include <vector>
#include <string>
#include <map>
#include <iostream>

#include <clipper/clipper.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/minimol/minimol.h>
#include <clipper/minimol/minimol_io.h>
#include <clipper/clipper-mmdb.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>

const double  MAP_NO_VAL = 99999. ;

using namespace std;
using namespace ::mmdb;
using namespace clipper;
typedef float ftype_t;

void readMap(Xmap<ftype_t>& xmap, char *mapFN) {
  clipper::CCP4MAPfile file;
  file.open_read( mapFN );
  file.import_xmap( xmap );
  file.close_read();
  Grid_range gr = xmap.grid_asu();
#if 0
  printf("min (%d,%d,%d) max (%d, %d, %d)\n", 
	 gr.min().u(),gr.min().v(),gr.min().w(), 
	 gr.max().u(),gr.max().v(),gr.max().w());
  gr.debug();
#endif
}

void checkPDB( Xmap<ftype_t> &dfMap, const char *pdbFN ) {
  MiniMol mmMol;
  MMDBfile molMFile;
  molMFile.read_file( pdbFN);
  MMDBManager& mmdb = static_cast<MMDBManager&>(molMFile);
  assert( mmdb.isSpaceGroup() );
  const Spacegroup &msg( mmMol.spacegroup() );
  molMFile.import_minimol( mmMol );
  Cell cell = mmMol.cell();
  Grid_sampling dfCG = dfMap.grid_sampling();

  for ( int p = 0; p < mmMol.size(); p++ ) { // polymers/chains
    for ( int m = 0; m < mmMol[p].size(); m++ ) { // monomers
      MMonomer mono(mmMol[p][m]);
      for ( int a = 0; a < mmMol[p][m].size(); a++ ) { //atoms
        MAtom& atom(mmMol[p][m][a]);
	const Coord_orth &aCo( atom.coord_orth());
	const Coord_frac &aCf( aCo.coord_frac(cell));
	const Coord_grid &aCg( aCf.coord_grid(dfCG));
	Xmap_base::Map_reference_coord refDf(dfMap);
	refDf.set_coord( aCg );
	ftype_t dVal = dfMap[ refDf ];
	printf("XX %d %s %s %s %s dVal %f\n", mono.seqnum(), atom.id().c_str(),
	       aCo.format().c_str(), 
	       aCf.format().c_str(), 
	       aCg.format().c_str(), 
	       dVal);
      }
    }
  }
}

void doMapSel( Xmap<ftype_t> &dfMap, char *rhoMapFN, char *fcMapFN, char * diffMapFN, 
	       float d1Cut, float d2Cut ) {
  Xmap<ftype_t> rhoMap, fcMap, diffMap;
  Xmap<ftype_t>  rhoSelMap, fcSelMap, diffSelMap;
  readMap(rhoMap, rhoMapFN);
  const Grid_sampling rhoGS(rhoMap.grid_sampling());
  const Grid_sampling dfGS(dfMap.grid_sampling());
  const Cell rhoCell( rhoMap.cell() );
  const Cell dfCell( dfMap.cell() );

  /*
  assert(rhoCell.a() == dfCell.a() &&
	 rhoCell.b() == dfCell.b() &&
	 rhoCell.c() == dfCell.c() &&
	 rhoCell.alpha() == dfCell.alpha() &&
	 rhoCell.beta() == dfCell.beta() &&
	 rhoCell.gamma() == dfCell.gamma() );
  */
  rhoMap.spacegroup().debug();
  dfMap.spacegroup().debug();
  printf("GS rho %s df %s\n", rhoGS.format().c_str(), dfGS.format().c_str() );
  printf("cell rho %s df %s\n", 
	 rhoCell.descr().format().c_str(), 
	 dfCell.descr().format().c_str() );
  rhoSelMap.init(rhoMap.spacegroup(), rhoCell, rhoGS );
  if(fcMapFN) {
    readMap(fcMap, fcMapFN);
    fcSelMap.init(rhoMap.spacegroup(), rhoCell, rhoGS );
  }
  if(diffMapFN) {
    readMap(diffMap, diffMapFN);
    diffSelMap.init(rhoMap.spacegroup(), rhoCell, rhoGS );
  }
  Xmap_base::Map_reference_index ix;
  int mapIxCount = 0, mapIxCopyCount = 0, mapIxZeroCount = 0;
  //for ( ix = rhoMap.first(); !ix.last(); ix.next() ) 
  //for ( ix = dfMap.first(); !ix.last(); ix.next() ) 
 
  vector<double> foV, fcV, diffV;
  for ( ix = rhoSelMap.first(); !ix.last(); ix.next() ) 
  {
    mapIxCount++;
    Coord_grid cg = ix.coord();
    Coord_frac cf  = cg.coord_frac(rhoGS);
    Coord_grid dfCg = cf.coord_grid(dfGS);
    Xmap_base::Map_reference_coord refDf(dfMap, dfCg);
    Xmap_base::Map_reference_coord refRho(rhoMap, cg);
    //refDf.set_coord(dfCg);
    //refOut.set_coord( rhoCg );
    //ftype_t dVal = dfMap[ refDf ];
    ftype_t dVal = dfMap[ refDf ];
    if( (dVal > d1Cut) && (dVal <= d2Cut) ) {
      //printf("%s %6.4f %6.4f\n", ix.coord().format().c_str(), dVal, rhoMap[ ix ]);
      //outMap[ refOut ] = rhoMap[ix] + 50.;
      ftype_t val = rhoMap[ refRho];
      foV.push_back ( val);
      rhoSelMap[ ix ] = val;
      //outMap[ refOut ] = dfMap[refDf];
      if(fcMapFN) {
	ftype_t fcVal = fcMap[ refRho];
	fcSelMap[ ix ] = fcVal;
	fcV.push_back ( fcVal);
      }
      if(diffMapFN) {
	ftype_t diffVal = diffMap[ refRho];
	diffSelMap[ ix ] = diffVal;
	diffV.push_back ( diffVal );
      }
      mapIxCopyCount++;
    } else {
      //outMap[ refOut ] = -(50 + dfMap[refDf]);
      rhoSelMap[ ix ] = MAP_NO_VAL;
      mapIxZeroCount++;
      if(fcMapFN)
	fcSelMap[ ix ] = MAP_NO_VAL;
      if(diffMapFN) 
	diffSelMap[ ix ] = MAP_NO_VAL;
    }
  }
  printf("=====\nstatsDict={\n");
  printf("  'mapVoxCount' : %d,\n  'mapZeroCount' : %d,\n  'mapCopiedCount' : %d,\n",
	 mapIxCount, mapIxZeroCount, mapIxCopyCount );

  double aVoxW = rhoCell.a()/rhoGS.nu();
  double bVoxW = rhoCell.b()/rhoGS.nv();
  double cVoxW = rhoCell.c()/rhoGS.nw();
  double voxVol =  aVoxW*bVoxW*cVoxW *
    sqrt( 2.0*cos(rhoCell.alpha())*cos(rhoCell.beta())*cos(rhoCell.gamma())
	- cos(rhoCell.alpha())*cos(rhoCell.alpha())
	- cos( rhoCell.beta())*cos( rhoCell.beta())
	- cos(rhoCell.gamma())*cos(rhoCell.gamma()) + 1.0 );

  Grid_range asuRange = rhoMap.grid_asu();
  Coord_grid gMin =asuRange.min();
  Coord_grid gMax =asuRange.max();
  assert(gMin.u() == 0 &&  gMin.v() == 0 && gMin.w() == 0 );
  Coord_frac fMax =gMax.coord_frac(rhoGS);
  double calcVol =  rhoCell.a() * rhoCell.b() * rhoCell.c() * 
    sqrt( 2.0*cos(rhoCell.alpha())*cos(rhoCell.beta())*cos(rhoCell.gamma())
	- cos(rhoCell.alpha())*cos(rhoCell.alpha())
	- cos( rhoCell.beta())*cos( rhoCell.beta())
	- cos(rhoCell.gamma())*cos(rhoCell.gamma()) + 1.0 );
  /*0
  rhoCell.a()*rhoCell.b()*rhoCell.c()*
      (1. - (cos(rhoCell.alpha())*cos(rhoCell.alpha())) -
	    (cos(rhoCell.beta())*cos(rhoCell.beta())) -
	    (cos(rhoCell.gamma())*cos(rhoCell.gamma()))) +
      2*sqrt(cos(rhoCell.alpha()) * cos(rhoCell.beta()) * cos(rhoCell.gamma()));
  */

  double mapVolAdjust = rhoCell.volume() /
    (voxVol*mapIxCount*rhoMap.spacegroup().num_symops());

    printf("  'cellVol' : %.3f,\n  'calcVolClp' : %.3f,\n  'voxVol' : %.3f,\n  'cellVolFromVox' : %.3f , #voxVol*nMapVox*nOps\n  'cellVolFromGrid' : %.3f, #voxVol*nGridVox\n  'mapVolAdjust' : %f,\n  'fracCopied' : %.3f,\n  'totMapPts' : %d,\n  'setMapPts' : %d,\n", 
	   rhoCell.volume(), calcVol, 
	   voxVol, voxVol*mapIxCount*rhoMap.spacegroup().num_symops(),
	   voxVol*rhoGS.nu()*rhoGS.nv()*rhoGS.nw(),
	   mapVolAdjust, 
	   ((double(mapIxCopyCount)/mapIxCount)*mapVolAdjust)*100.,
	   mapIxCount,
	   mapIxCopyCount
	   //float(mapIxCount)/float((rhoGS.nu()-1)*(rhoGS.nv()-1)*(rhoGS.nw()-1)) 
           );
  if( fcMapFN ) {
    double pearson = 0., fcMean = 0., fcSD = 0.;
    if( foV.size() >  0 && fcV.size() > 0) {
      assert( foV.size() ==  fcV.size() );
      gsl_vector_const_view fo = gsl_vector_const_view_array( &foV[0], foV.size());
      gsl_vector_const_view fc = gsl_vector_const_view_array( &fcV[0], fcV.size());

      pearson = gsl_stats_correlation( (double*) fo.vector.data, 1,
					    (double*) fc.vector.data, 1,
				       foV.size()  );
      fcSD = gsl_stats_sd( (double*)     fc.vector.data, 1, fcV.size()  );
      fcMean = gsl_stats_mean( (double*) fc.vector.data, 1, fcV.size()  );
    }
    printf ("  'fcSD' : %f,\n",    fcSD);
    printf ("  'fcMean' : %f,\n",  fcMean);
    printf ("  'pearson' : %f,\n", pearson);
  }
  if( diffMapFN && diffV.size() > 0) {
    double diffSD = 0., diffMean = 9999.;
    gsl_vector_const_view fd = gsl_vector_const_view_array( &diffV[0], diffV.size());
    diffSD = gsl_stats_sd( (double*) fd.vector.data, 1, diffV.size()  );
    diffMean = gsl_stats_mean( (double*) fd.vector.data, 1, diffV.size()  );
    printf ("  'diffSD' : %f,\n", diffSD);
    printf ("  'diffMean' : %f,\n", diffMean);
  }
  printf("}\n"); // close statsDict
  CCP4MAPfile mapWrite1;
  mapWrite1.open_write( "2fofc_sel.map" );
  mapWrite1.export_xmap( rhoSelMap );
  mapWrite1.close_write();	
  if( fcMapFN && fcV.size() ) {
    CCP4MAPfile mapWrite2;
    mapWrite2.open_write( "Dfc_sel.map" );
    mapWrite2.export_xmap( fcSelMap );
    mapWrite2.close_write();	
  }
  if( diffMapFN && diffV.size() > 0) {
    CCP4MAPfile mapWrite3;
    mapWrite3.open_write( "fofc_sel.map" );
    mapWrite3.export_xmap( diffSelMap );
    mapWrite3.close_write();	
  }
}

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
  cAtom->SetResidue(res); 
  res->AddAtom(cAtom);
}

void doPDBSel( Xmap<ftype_t> &dfMap, char *pdbFN ) {
  MiniMol mmMol;
  MMDBfile molMFile;
  molMFile.read_file( pdbFN);
  MMDBManager& mmdb = static_cast<MMDBManager&>(molMFile);
  assert( mmdb.isSpaceGroup() );
  const Spacegroup &msg( mmMol.spacegroup() );
  molMFile.import_minimol( mmMol );
  Cell cell = mmMol.cell();

  MMDBManager cMMDB;
  cMMDB.set_cell(mmMol.cell());
  cMMDB.set_spacegroup( mmMol.spacegroup() );
  Model* model = new Model;
  cMMDB.AddModel(model);

  Grid_sampling dfCG = dfMap.grid_sampling();
  Chain* chain = 0;
  Residue* res = 0;
  for ( int p = 0; p < mmMol.size(); p++ ) { // polymers/chains
    printf("XX minimol chain '%s'\n", mmMol[p].id().c_str() );  
    chain = new Chain;
    chain->SetChainID ( mmMol[p].id().c_str() );  
    model->AddChain(chain);
    for ( int m = 0; m < mmMol[p].size(); m++ ) { // monomers
      MMonomer mono(mmMol[p][m]);
      res = new Residue;
      res->seqNum = mono.seqnum();
      res->SetResName(mono.type().c_str());
      chain->AddResidue ( res );
      for ( int a = 0; a < mmMol[p][m].size(); a++ ) { //atoms
        MAtom& atom(mmMol[p][m][a]);
	const Coord_orth &aCo( atom.coord_orth());
	const Coord_frac &aCf( aCo.coord_frac(cell));
	const Coord_grid &aCg( aCf.coord_grid(dfCG));
	Xmap_base::Map_reference_coord refDf(dfMap);
	refDf.set_coord( aCg );
	ftype_t dVal = dfMap[ refDf ];
	printf("XX %d dVal %f", mono.seqnum(), dVal);
        if(dVal >0 && dVal < 2. )  {
	  writeAtom(atom, aCo, res);
	  printf("\n");
	} else
	  printf("excluded\n");
      }
    }
  }

  cMMDB.PDBCleanup ( PDBCLEAN_SERIAL | PDBCLEAN_INDEX );
  string pdbName = "dfsel.pdb";
  ERROR_CODE writeRC = cMMDB.WritePDBASCII( const_cast<char*>(pdbName.c_str()) );
  if(writeRC != 0) {
    printf("dfsel: PDB write err %s\n", GetErrorDescription(writeRC) );
  }
}

int main(int argc, char **argv) { 
  --argc, ++argv;
  char *dfMapFN= 0, *rhoMapFN = 0, *fcMapFN = 0, *pdbFN = 0, *diffMapFN = 0;
  float d1Cut = 0.25, d2Cut= FLT_MAX;

  if(argc ) do {
    if(!strcmp(*argv,"-df" )) {
      --argc, ++argv;
      dfMapFN= *argv;
    } else if(!strcmp(*argv,"-m" )) {
      --argc, ++argv;
      rhoMapFN = *argv;
    } else if(!strcmp(*argv,"-fc" )) {
      --argc, ++argv;
      fcMapFN = *argv;
    } else if(!strcmp(*argv,"-diff" )) {
      --argc, ++argv;
      diffMapFN = *argv;
    } else if(!strcmp(*argv,"-c" )) {
      --argc, ++argv;
      pdbFN = *argv;
    } else if(!strcmp(*argv,"-d1" )) {
      --argc, ++argv;
      d1Cut = atof(*argv);
    } else if(!strcmp(*argv,"-d2" )) {
      --argc, ++argv;
      d2Cut = atof(*argv);
    } else {
      assert(!"bad args");
    }
    --argc, ++argv;
  } while (argc);
  assert( dfMapFN && (rhoMapFN || pdbFN) );
  Xmap<ftype_t> dfMap;
  readMap(dfMap, dfMapFN);
  const Grid_sampling dfGS(dfMap.grid_sampling());
  const Cell dfCell( dfMap.cell() );
  if( rhoMapFN) {
    //checkPDB( dfMap, "./p1a-cell.pdb" );
    doMapSel( dfMap, rhoMapFN, fcMapFN, diffMapFN, d1Cut, d2Cut );
  }
  if( pdbFN) {
    doPDBSel( dfMap, pdbFN );
  }
}
