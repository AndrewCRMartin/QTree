/*************************************************************************

   Program:    QTree
   File:       graphics.c
   
   Version:    V3.0
   Date:       19.08.19
   Function:   Display routines for QTree
   
   Copyright:  (c) SciTech Software 1993-2019
   Author:     Prof. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain.

   It may not be copied or made available to third parties, but may be
   freely used by non-profit-making organisations who have obtained it
   directly from the author or by FTP.
   
   You are requested to send EMail to the author to say that you are
   using this code so that you may be informed of future updates.
   
   The code may not be made available on other FTP sites without express
   permission from the author.
   
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If
   someone else breaks this code, the author doesn't want to be blamed
   for code that does not work! You may not distribute any
   modifications, but are encouraged to send them to the author so
   that they may be incorporated into future versions of the code.
   
   The code may not be sold commercially or used for commercial purposes
   without prior permission from the author.
                                                
**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  19.07.93 Original
   V1.1  28.07.93 Added support for ball and stick
   V1.2  29.07.93 Added MTV file output support and background colouring
   V1.3  04.08.93 Corrected usage message
   V1.4  12.08.93 Added force option to create output file of dimensions
                  equal to image rather than XSIZE/YSIZE
                  Changed to require output file. Removed direct graphic
                  display
   V1.5  14.09.93 Skipped
   V1.6  04.01.94 Added casts for GCC
   V1.7  24.03.94 Minimal tidying up
   V1.8  09.05.94 Skipped
   V1.9  13.05.94 Skipped
   V1.10 24.06.94 Handles TERMPERATURE (changes in commands.c)
   V1.11 04.10.94 Skipped
   V1.12 21.12.94 Skipped
   V2.0  28.03.95 Modified for output on stdout
   V2.1  23.10.95 Skipped
   V2.2  14.10.03 Skipped
   V2.3  18.10.03 Skipped
   V2.4  27.01.15 Skipped
   V2.5  18.08.19 General cleanup and moved into GitHub
   V3.0  19.08.19 Added PNG support

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/array.h"

#include "qtree.h"
#include "graphics.p"

/************************************************************************/
/* Arrays to contain the image
*/
static unsigned char 
#ifdef _AMIGA
                     __far 
#endif
                           **sRed   = NULL,
                           **sGreen = NULL,
                           **sBlue  = NULL;
               
/************************************************************************/
/*>BOOL InitGraphics(void)
   -----------------------
   Perform any initialisation for graphics.
   19.07.93 Original    By: ACRM
   12.08.93 Modified for file output only.
   04.01.94 Added casts on blFreeArray2D
*/
BOOL InitGraphics(void)
{
   int   x, y;
   
   sRed   = (unsigned char **)blArray2D(sizeof(unsigned char), 
                                      gScreen[0], gScreen[1]);
   sGreen = (unsigned char **)blArray2D(sizeof(unsigned char), 
                                      gScreen[0], gScreen[1]);
   sBlue  = (unsigned char **)blArray2D(sizeof(unsigned char), 
                                      gScreen[0], gScreen[1]);

   if(sRed == NULL || sGreen == NULL || sBlue == NULL)
   {
      if(sRed   == NULL) blFreeArray2D((char **)sRed,  
                                     gScreen[0], gScreen[1]);
      if(sGreen == NULL) blFreeArray2D((char **)sGreen, 
                                     gScreen[0], gScreen[1]);
      if(sBlue  == NULL) blFreeArray2D((char **)sBlue,  
                                     gScreen[0], gScreen[1]);
      
      sRed   = NULL;
      sGreen = NULL;
      sBlue  = NULL;
      
      return(FALSE);
   }
   
   /* Clear the screen                                                  */
   for(y=0; y<gScreen[1]; y++)
      for(x=0; x<gScreen[0]; x++)
         sRed[x][y] = sGreen[x][y] = sBlue[x][y] = 0;

   return(TRUE);
}


/************************************************************************/
/*>void EndGraphics(void)
   ----------------------
   Performs any remaining graphics display routines, waits for return to
   be pressed, and cleans up the graphics. Also writes the graphics to
   a file if specified.
   19.07.93 Original    By: ACRM
   29.07.93 Added file output
   12.08.93 Modified for screen size specification
   04.01.94 Added casts on blFreeArray2D
   28.03.95 No longer checks file specified
   19.08.19 Now takes filename and type as a parameter
*/
void EndGraphics(char *outFile, int outFormat)
{
   /* Write the file in MTV format                                      */
   switch(outFormat)
   {
   case OUTPUT_MTV:
      WriteMTVFile(outFile,gScreen[0],gScreen[1]);
      break;
#ifdef SUPPORT_PNG
   case OUTPUT_PNG:
      WritePNGFile(outFile,gScreen[0],gScreen[1]);
      break;
#endif
   default:
      fprintf(stderr, "Error: Unknown output format! Using MTV");
      WriteMTVFile(outFile,gScreen[0],gScreen[1]);
   }

   /* Free memory for RGB arrays                                        */
   if(sRed   != NULL) blFreeArray2D((char **)sRed,   
                                  gScreen[0], gScreen[1]);
   if(sGreen != NULL) blFreeArray2D((char **)sGreen, 
                                  gScreen[0], gScreen[1]);
   if(sBlue  != NULL) blFreeArray2D((char **)sBlue,  
                                  gScreen[0], gScreen[1]);
   
   sRed   = NULL;
   sGreen = NULL;
   sBlue  = NULL;
}


/************************************************************************/
/*>void SetPixel(int x0, int y0, REAL r, REAL g, REAL b)
   -----------------------------------------------------
   Sets a pixel applying offsets to centre image on screen. This is just
   done into the colour arrays. Currently nothing is actually placed on 
   the screen.
   19.07.93 Original    By: ACRM
   21.07.93 Corrected y calculation
   29.07.93 Added centering
   12.08.93 Modified for screen size spec
   19.10.07 Checks that pixels are in range
*/
void SetPixel(int x0, int y0, REAL r, REAL g, REAL b)
{
   int temp;

   if((x0 >= 0) && (x0 < gScreen[0]) &&
      (y0 >= 0) && (y0 < gScreen[1]))
   {
      x0 += (gScreen[0]-gSize)/2;
      y0 += (gScreen[1]-gSize)/2;
      
      y0 = gScreen[1] - y0 - 1;
      
      temp = (int)(256.0 * r + 0.5);
      sRed[x0][y0]   = (temp > 255) ? 255 : temp;
      
      temp = (int)(256.0 * g + 0.5);
      sGreen[x0][y0] = (temp > 255) ? 255 : temp;
      
      temp = (int)(256.0 * b + 0.5);
      sBlue[x0][y0]  = (temp > 255) ? 255 : temp;
   }
}


/************************************************************************/
/*>void SetAbsPixel(int x0, int y0, REAL r, REAL g, REAL b)
   --------------------------------------------------------
   Sets a pixel applying no offsets. Used for colouring background. This 
   is just one into the colour arrays. Currently nothing is actually 
   placed on the screen.

   29.07.93 Original based on SetPixel()     By: ACRM
   12.08.93 Modified for screen size spec
*/
void SetAbsPixel(int x0, int y0, REAL r, REAL g, REAL b)
{
   int temp;


   temp = (int)(256.0 * r + 0.5);
   sRed[x0][y0]   = (temp > 255) ? 255 : temp;
   
   temp = (int)(256.0 * g + 0.5);
   sGreen[x0][y0] = (temp > 255) ? 255 : temp;
   
   temp = (int)(256.0 * b + 0.5);
   sBlue[x0][y0]  = (temp > 255) ? 255 : temp;
}

/************************************************************************/
/*>BOOL WriteMTVFile(char *FileName, int xsize, int ysize)
   -------------------------------------------------------
   Write a graphics file in MTV format
   29.07.93 Original    By: ACRM
   12.08.93 Added check on arrays
   28.03.95 Modified for output on stdout if blank filename
*/
BOOL WriteMTVFile(char *FileName, int xsize, int ysize)
{
   FILE  *fp;
   int   x, y;


   if(sRed==NULL || sGreen == NULL || sBlue==NULL) 
      return(FALSE);

   if(FileName[0])
      fp = fopen(FileName,"w");
   else
      fp = stdout;

   if(fp != NULL)
   {
      /* Header                                                         */
      fprintf(fp,"%d %d\n",xsize,ysize);
      
      /* And the graphics                                               */
      for(y=0; y<ysize; y++)
      {
         for(x=0; x<xsize; x++)
         {
            fputc((int)sRed[x][y],   fp);
            fputc((int)sGreen[x][y], fp);
            fputc((int)sBlue[x][y],  fp);
         }
      }
      
      fclose(fp);
      return(TRUE);
   }
   
   return(FALSE);
}


#ifdef SUPPORT_PNG
#include "writepng.h"
/************************************************************************/
/*>BOOL WritePNGFile(char *FileName, int xsize, int ysize)
   -------------------------------------------------------
   Write a graphics file in PNG format

   19.08.19 Original    By: ACRM
*/
BOOL WritePNGFile(char *FileName, int xsize, int ysize)
{
   int   x, y;
   blPNGIMAGE image;
   BOOL retval = TRUE;
   
   if(sRed==NULL || sGreen == NULL || sBlue==NULL) 
      return(FALSE);

   image.width  = xsize;
   image.height = ysize;
   if((image.pixels = (blPNGPIXEL *)calloc(xsize*ysize,
                                           sizeof(blPNGPIXEL)))==NULL)
   {
      return(FALSE);
   }
   
   /* Copy the image over into the blPNGIMAGE bitmap                     */
   for(y=0; y<ysize; y++)
   {
      for(x=0; x<xsize; x++)
      {
         blPNGPIXEL *pixel = blPNGPixelAt(&image, x, y);

         pixel->red   = (int)sRed[x][y];
         pixel->green = (int)sGreen[x][y];
         pixel->blue  = (int)sBlue[x][y];
      }
   }

   /* Write the image to the file or stdout if FileName is a NULL string */
   if(!blSavePNGToFile(&image, FileName)) 
   {
      fprintf(stderr, "Error writing file.\n");
      retval = FALSE;
   }
   
   free(image.pixels);
   
   return(retval);
}

#endif


