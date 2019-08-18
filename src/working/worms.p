int main(int argc, char **argv)
;
PDB *SplinePDB(PDB *pdb, int SplineIter, int *nspline)
;
PDB *Spline(PDB *InWorm, int NIn)
;
void WriteWorm(FILE *fp, PDB *worm, int nworm)
;
PDB *WormPDB(PDB *spline, int nspline, int NDivide, int *nworm)
;
PDB *FindChainPDB(PDB *pdb)
;
