/* Example code to read a .mtv file
*/




/************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines
*/
#define XSIZE 1280
#define YSIZE 1024

/************************************************************************/
/* Global arrays to contain the picture once read. Each position contains
   a value between 0 and 255 for R, G and B.

   These will need to be accessed by whatever software is used to display
   the image. N.B. If a workstation with only 256 colours is in use,
   some conversion via a colour lookup table (CLUT) will be necessary
   before display.
*/
unsigned char  gRed[XSIZE][YSIZE],
               gGreen[XSIZE][YSIZE],
               gBlue[XSIZE][YSIZE];

/************************************************************************/
/*>BOOL ReadMTV(FILE *fp, int *xres, int *yres)
   --------------------------------------------
   Read an MTV file. Takes a file pointer as input and outputs the x- and
   y-resolutions of the picture which is read into the global arrays,
   gRed[][], gGreen[][] and gBlue[][]
   
   03.08.93 Original    By: ACRM
*/
BOOL ReadMTV(FILE *fp, int *xres, int *yres)
{
   int   x,
         y;
         
   if(fscanf(fp, "%d %d",xres,yres) != 2) return(FALSE);
   
   if(*xres > XSIZE || *yres > YSIZE) return(FALSE);
   
   if(getc(fp) != '\n') return(FALSE);
   
   for(y=0; y < (*yres); y++)
   {
      for(x=0; x < (*xres); x++)
      {
         gRed[x][y]   = getc(fp);
         gGreen[x][y] = getc(fp);
         gBlue[x][y]  = getc(fp);
      }
   }
   
   return(TRUE);
}

