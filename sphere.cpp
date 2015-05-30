// Bresenham variant of midpoint circle algorithm adapted for sphere
// from Graphics Gems "Spheres-to-voxels conversion" chapter by Claudio Montani 
// and Roberto Scopigno

#include <stdio.h>
#include <stdlib.h>
#include <list>
using namespace std;
#include "sphere.h"

void Sphere(int R, list<LStruct> &voxList);
void Slice(int r, int y);
void Track(int x, int y, int z);
void Voxel(int x, int y, int z);
#if 0
int main(int argc, char **argv) {
  int R = 3;
  --argc, ++argv;
  if(argc == 1)
    R = atol(*argv);
  Sphere(R);
}
#endif
static  list<LStruct> *voxListPtr;
void Sphere(int R, list<LStruct> &voxListArg)
{
  voxListPtr = &voxListArg;
  int x = 0, y = R;
  int delta=2*(1-R), limit = 0;
  while (y >= limit) {
    if (delta < 0 ) {
      int deltaLC = 2*delta + 2*y - 1;
      if(deltaLC > 0 ) {
	Slice(x,y);
	x = x + 1; y = y - 1;
        delta = delta + 2*x - 2*y + 2;
      } else {
	x = x+1;
	delta = delta + 2*x + 1;
      }
    } else if (delta > 0 ) {
      int deltaLC = 2*delta - 2*x - 1;
      if(deltaLC > 0 ) {
	Slice(x,y);
	y = y - 1; 
        delta = delta - 2*y + 1;
      } else {
	Slice(x,y);
	x = x + 1; y = y - 1;
        delta = delta + 2*x - 2*y + 2;
      }
    } else {
      Slice(x,y);
      x = x+1; y = y-1; 
      delta = delta + 2*x - 2*y + 2;
    }
  }
}

void Slice(int r, int y) {
  int x = 0, z = r;
  int delta = 2*(1-r), limit = 0;
  Track( x, y, z );
  while( z >= limit ) {
    if( delta < 0 ) {
      int deltaLC = 2*delta + 2*z - 1;
      if( deltaLC > 0 ) {
	x = x + 1; z = z - 1;
	delta = delta + 2*x - 2*z +2;
	Track( x, y, z );
	Track( -x, y, z );
      } else {
	x = x + 1; delta = delta + 2*x  +1;
	Track( x, y, z );
	Track( -x, y, z );
      }
    } else if( delta > 0 ) {
      int deltaLC = 2*delta - 2*x - 1;
      if( deltaLC > 0 ) {
	z = z - 1;
	delta = delta - 2*z + 1;
      } else {
	x = x + 1; z = z - 1;
	delta = delta + 2*x -2*z +2;
	Track( x, y, z );
	Track( -x, y, z );
      }
    } else {
 	x = x + 1; z = z - 1;
	delta = delta + 2*x - 2*z  + 2;
	Track( x, y, z );
	Track( -x, y, z );
   }
  }
}

void Track(int x, int y, int z) {
  for(int k = -z; k <= z;  k= k+ 1 ){
    Voxel(x, y, k);
  }
  if( y != 0 ) {
    for(int k = -z; k <= z;  k= k+ 1 ){    
      Voxel(x, -y, k);
    }
  }
}
void Voxel(int x, int y, int z) {
  LStruct lStruct;
  lStruct.i = x;
  lStruct.j = y;
  lStruct.k = z;
  voxListPtr->push_back(lStruct);
  //printf("%d %d %d\n", x, y, z);
}
