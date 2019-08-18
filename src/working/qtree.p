main(int argc, char **argv)
;
void MapSpheres(PDB *pdb, SPHERE *spheres, int NSphere)
;
SPHERE *CreateSphereList(PDB *pdb, int NAtom)
;
BOOL SpaceFill(SPHERE *AllSpheres, int NSphere)
;
void SplitPic(int x0, int y0, int x1, int y1, SPHERE **spheres,
              int NSphere)
;
SPHERE **UpdateSphereList(REAL   x0, 
                          REAL   y0, 
                          REAL   x1, 
                          REAL   y1, 
                          SPHERE **spheres,
                          int    NSphere,
                          int    *NSphOut)
;
SPHERE **SortSpheresOnX(SPHERE *AllSpheres, int NSphere)
;
SPHERE **SortSpheresOnZ(SPHERE **Spheres, int NSphere)
;
void ColourPixel(REAL x, REAL y, SPHERE **spheres, int NSphere)
;
int FarLeftSearch(SPHERE **spheres, int NSphere, REAL x)
;
int FarRightSearch(SPHERE **spheres, int NSphere, REAL x)
;
int CtrlCExit(void)
;
int CtrlCNoExit(void)
;
void ShadePixel(REAL x, REAL y, REAL z, SPHERE *sphere)
;
void UsageExit(BOOL ShowHelp)
;
