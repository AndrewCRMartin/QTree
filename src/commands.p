BOOL SetupParser(void)
;
void HandleControl(char *file, PDB *pdb, SPHERE *spheres, int NSphere,
                   BOOL ReportError)
;
void DoDefault(SPHERE *spheres, int NSphere, REAL RGB[3],
               PDB *pdb, BOOL ColourByTemp)
;
void DoPhong(SPHERE *spheres, int NSphere, REAL shine, REAL metallic)
;
void DoZone(SPHERE *spheres, PDB *pdb, int NSphere, char *start, 
            char *end, char *red, char *green, char *blue, int type)
;
void DoResidue(SPHERE *spheres, PDB *pdb, int NSphere, char *resnam, 
               char *red, char *green, char *blue)
;
void DoRotate(PDB *pdb, char *direction, char *amount)
;
void DoCentre(PDB *pdb, char *resspec, char *atom)
;
void DoBackground(REAL r1, REAL g1, REAL b1, REAL r2, REAL g2, REAL b2)
;
void DoSlab(PDB *pdb, char *resspec, char *atom)
;
void HSL2RGB(REAL hue, REAL saturation, REAL luminance,
             REAL *red, REAL *green, REAL *blue)
;
void DoChain(SPHERE *spheres, PDB *pdb, int NSphere, char *chain, 
             char *red, char *green, char *blue)
;
void DoRadius(char *atomspec, char *radius_str)
;
void DoAtom(SPHERE *spheres, PDB *pdb, int NSphere, char *atom, 
            char *red, char *green, char *blue)
;
