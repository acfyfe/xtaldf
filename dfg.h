//typedef double fmtype;
//#define FMTYPE_MAX DBL_MAX

typedef float fmtype;
#define FMTYPE_MAX FLT_MAX
#define FMTYPE_MIN FLT_MIN
void mkGrid(char* pdbFn, fmtype *grid,  int *dObj, int nRow, int nCol, int nPlanes);
struct ElemDescr {
  char * elemName;
  //int elemNum;
  float vdwRadius;
};
extern ElemDescr vdwTab[];
extern int vdwTabSz;
