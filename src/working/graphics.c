/*************************************************************************

   Program:    QTree
   File:       graphics.c
   
   Version:    V1.2
   Date:       29.07.93
   Function:   Display routines for QTree
   
   Copyright:  (c) SciTech Software 1993
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0372) 275775
   EMail:      UUCP:  cbmehq!cbmuk!scitec!amartin
                      amartin@scitec.adsp.sub.org
               JANET: andrew@uk.ac.ox.biop
               
**************************************************************************

   This program is not in the public domain, but it may be freely copied
   and distributed for no charge providing this header is included.
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! The code may not be sold commercially without prior permission 
   from the author, although it may be given away free with commercial 
   products, providing it is made clear that this program is free and that 
   the source code is provided with the program.

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

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include <intuition/intuition.h>
#include <proto/all.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/array.h"

#include "qtree.h"
#include "graphics.p"

/************************************************************************/
/* Arrays to contain the image
*/
static unsigned char __far **sRed   = NULL,
                           **sGreen = NULL,
                           **sBlue  = NULL;
               
/************************************************************************/
/*>BOOL InitGraphics(void)
   -----------------------
   Perform any initialisation for graphics.
   19.07.93 Original    By: ACRM
   12.08.93 Modified for file output only.
*/
BOOL InitGraphics(void)
{
   int   x, y;
   
   sRed   = (unsigned char **)Array2D(sizeof(unsigned char), 
                                      gScreen[0], gScreen[1]);
   sGreen = (unsigned char **)Array2D(sizeof(unsigned char), 
                                      gScreen[0], gScreen[1]);
   sBlue  = (unsigned char **)Array2D(sizeof(unsigned char), 
                                      gScreen[0], gScreen[1]);

   if(sRed == NULL || sGreen == NULL || sBlue == NULL)
   {
      if(sRed   == NULL) FreeArray2D(sRed,   gScreen[0], gScreen[1]);
      if(sGreen == NULL) FreeArray2D(sGreen, gScreen[0], gScreen[1]);
      if(sBlue  == NULL) FreeArray2D(sBlue,  gScreen[0], gScreen[1]);
      
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
*/
void EndGraphics(void)
{
   /* If a file name was specified, write the file in MTV format        */
   if(gOutFile[0])
      WriteMTVFile(gOutFile,gScreen[0],gScreen[1]);

   /* Free memory for RGB arrays                                        */
   if(sRed   != NULL) FreeArray2D(sRed,   gScreen[0], gScreen[1]);
   if(sGreen != NULL) FreeArray2D(sGreen, gScreen[0], gScreen[1]);
   if(sBlue  != NULL) FreeArray2D(sBlue,  gScreen[0], gScreen[1]);
   
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
*/
void SetPixel(int x0, int y0, REAL r, REAL g, REAL b)
{
   int temp;

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

//   y0 = gScreen[1] - y0 - 1;
   
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
*/
BOOL WriteMTVFile(char *FileName, int xsize, int ysize)
{
   FILE  *fp;
   int   x, y;


   if(sRed==NULL || sGreen == NULL || sBlue==NULL) 
      return(FALSE);

   if((fp=fopen(FileName,"w")) != NULL)
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



