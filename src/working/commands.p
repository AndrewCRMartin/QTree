BOOL SetupParser(void)
;
void HandleControl(char *file, PDB *pdb, SPHERE *spheres, int NSphere, 
                   BOOL ReportError)
;
void DoDefault(SPHERE *spheres, int NSphere, REAL RGB[3])
;
void DoPhong(SPHERE *spheres, int NSphere, REAL shine, REAL metallic)
;
void DoZone(SPHERE *spheres, PDB *pdb, int NSphere, char *start, 
            char *end, char *red, char *green, char *blue)
;
void ParseResSpec(char *spec, char *chain, int *resnum, char *insert)
;
void DoAtom(SPHERE *spheres, PDB *pdb, int NSphere, char *atom, 
            char *red, char *green, char *blue)
;
void DoResidue(SPHERE *spheres, PDB *pdb, int NSphere, char *resnam, 
               char *red, char *green, char *blue)
;
void DoRotate(PDB *pdb, char *direction, char *amount)
;
void DoCentre(PDB *pdb, char *resspec, char *atomspec)
;
void DoBackground(REAL r1, REAL g1, REAL b1, REAL r2, REAL g2, REAL b2)
;