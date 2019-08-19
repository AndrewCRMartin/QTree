BOOL InitGraphics(void)
;
void EndGraphics(char *outFile, int outFormat)
;
void SetPixel(int x0, int y0, REAL r, REAL g, REAL b)
;
void SetAbsPixel(int x0, int y0, REAL r, REAL g, REAL b)
;
BOOL WriteMTVFile(char *FileName, int xsize, int ysize);
BOOL WritePNGFile(char *FileName, int xsize, int ysize);
