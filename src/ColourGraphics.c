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
   This is an example file for direct graphics support from within QTree.
   This has now been abandoned in favour of writing a file.

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

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include <intuition/intuition.h>
#include <proto/all.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

#include "qtree.h"
#include "graphics.p"

/************************************************************************/
/* Defines
*/
#define THRESHHOLD 4

/************************************************************************/
/* Amiga graphics setup
*/
#ifdef _AMIGA  /*-------------------- AMIGA ONLY -----------------------*/
struct IntuitionBase *IntuitionBase = NULL;
struct GfxBase       *GfxBase       = NULL;
struct Screen        *scrn          = NULL;
struct Window        *wind          = NULL;
struct RastPort      *rp            = NULL;
struct ViewPort      *viewport      = NULL;
int                  nallocr        = 0;
#endif         /*------------------ END OF AMIGA CODE ------------------*/

/************************************************************************/
/* Arrays to contain the image
*/
static unsigned char __far sRed[XSIZE][YSIZE],
                           sGreen[XSIZE][YSIZE],
                           sBlue[XSIZE][YSIZE];
               
/************************************************************************/
/*>BOOL InitGraphics(void)
   -----------------------
   Perform any initialisation for graphics.
   Supports Amiga HAM mode
   19.07.93 Original    By: ACRM
*/
BOOL InitGraphics(void)
{
#ifdef _AMIGA  /*-------------------- AMIGA ONLY -----------------------*/
   int      i;

   struct NewScreen newscreen = 
   {
      0,0,              /* Leftedge, topedge                            */
      XSIZE,YSIZE,      /* Width, height                                */
      6,                /* Depth                                        */
      0,0,              /* pens                                         */
      HAM,              /* ViewModes                                    */
      CUSTOMSCREEN,     /* Flags                                        */
      NULL,             /* Font                                         */
      NULL,             /* title                                        */
      NULL,             /* gadgets                                      */
      NULL              /* bitmap                                       */
   };
   struct NewWindow newwindow = 
   {
      0,0,              /* LeftEdge, TopEdge                            */
      XSIZE,YSIZE,      /* Width, Height                                */
      0,1,              /* DetailPen, BlockPen                          */
      MENUPICK,         /* IDCMPFlags                                   */
      BORDERLESS,       /* Flags                                        */
      NULL,             /* FirstGadget                                  */
      NULL,             /* CheckMark                                    */
      NULL,             /* Title                                        */
      NULL,             /* Screen                                       */
      NULL,             /* BitMap                                       */
      0,0,              /* MinWidth, MinHeight                          */
      0,0,              /* MaxWidth, MaxHeight                          */
      CUSTOMSCREEN      /* Type                                         */
   };
 
 
   /* Open intuition library                                            */
   if((IntuitionBase = (struct IntuitionBase *)
                       OpenLibrary("intuition.library",36)) == NULL)
   {
      return(FALSE);
   }
      
   /* Open graphics library                                             */
   if((GfxBase = (struct GfxBase *)
                 OpenLibrary("graphics.library",36)) == NULL)
   {
      return(FALSE);
   }
      
                   
   newscreen.Width   = XSIZE;
   newscreen.Height  = YSIZE;
   newscreen.TopEdge = 0;

   if(newscreen.Width < 362)
      newscreen.LeftEdge = (362 - newscreen.Width)/2;

   if((scrn = OpenScreen(&newscreen)) == NULL)
   {
      return(FALSE);
   }

   newwindow.Screen    = scrn;
   newwindow.Width     = XSIZE;
   newwindow.MinWidth  = XSIZE;
   newwindow.MaxWidth  = XSIZE;
   newwindow.Height    = YSIZE;
   newwindow.MinHeight = YSIZE;
   newwindow.MaxHeight = YSIZE;

   if((wind = OpenWindow(&newwindow)) == NULL)
   {
      return(FALSE);
   }

   rp = wind->RPort;
 
   viewport = &(scrn->ViewPort);

   for(i=0; i<32; i++)
      SetRGB4(viewport,i,0,0,0);

   SetRGB4(viewport,1,15,15,15);
   
   nallocr = 2;      /* Reset number of allocated colour registers      */
#endif         /*------------------ END OF AMIGA CODE ------------------*/

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
*/
void EndGraphics(void)
{
   char  buffer[80];

#ifdef _AMIGA  /*-------------------- AMIGA ONLY -----------------------*/
   int   x,y,pix[3];

   printf(" Displaying...");
      
   for(y=0;y<YSIZE;y++)
   {
      for(x=0;x<XSIZE;x++)
      {
         pix[0] = sRed[x][y]   / 17;
         pix[1] = sGreen[x][y] / 17;
         pix[2] = sBlue[x][y]  / 17;

         if(pix[0] > 15) pix[0] = 15;
         if(pix[1] > 15) pix[1] = 15;
         if(pix[2] > 15) pix[2] = 15;
         HAM2(x,y,pix);
      }
   }
   
   printf("Complete.\n");
#endif         /*------------------ END OF AMIGA CODE ------------------*/

   /* Wait for return to be pressed unless gQuickExit is set            */
   if(!gQuickExit)
   {
      printf("Press <return> to exit...");
      gets(buffer);
   }

#ifdef _AMIGA  /*-------------------- AMIGA ONLY -----------------------*/
   if(wind != NULL)           CloseWindow(wind);
   if(scrn != NULL)           CloseScreen(scrn);
   if(GfxBase != NULL)        CloseLibrary((struct Library *)
                                           GfxBase);
   if(IntuitionBase != NULL)  CloseLibrary((struct Library *)
                                           IntuitionBase);
#endif         /*------------------ END OF AMIGA CODE ------------------*/

   /* If a file name was specified, write the file in MTV format        */
   if(gOutFile[0])
   {
      WriteMTVFile(gOutFile,XSIZE,YSIZE);
   }
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
*/
void SetPixel(int x0, int y0, REAL r, REAL g, REAL b)
{
   int temp;

   x0 += (XSIZE-gSize)/2;
   y0 += (YSIZE-gSize)/2;
   
   y0 = YSIZE - y0 - 1;
   
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
*/
void SetAbsPixel(int x0, int y0, REAL r, REAL g, REAL b)
{
   int temp;

   y0 = YSIZE - y0 - 1;
   
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
*/
BOOL WriteMTVFile(char *FileName, int xsize, int ysize)
{
   FILE  *fp;
   int   x, y;
   
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


#ifdef _AMIGA  /*-------------------- AMIGA ONLY -----------------------*/
/************************************************************************/
/*                                                                      */
/*                    AMIGA SPECIFIC GRAPHICS ROUTINES                  */
/*                                                                      */
/************************************************************************/
/* Colour register array
*/
static int  ColourRegister[32][3]=
{  {0,0,0},
   {15,15,15},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0},
   {0,0,0}
};

/************************************************************************/
/*>void HAM2(int i, int j, int pix[3])
   -----------------------------------
   Takes a coordinate and a colour specified as 3 integers of 0--15 and
   plots a HAM pixel.
   19.07.93 Original taken from CPK    By: ACRM
*/
void HAM2(int i, int j, int pix[3])
{
   int         k;
   static int  prevpix[3],
               map[3] = {0x20,0x30,0x10};
   int         dif,  dif2,
               id,
               maxdif,
               pen,
               regist = FALSE;
               
   /* If it's the first pixel on a line, use a register                 */
   if(!i)
   {
      pen = NearestPen(pix,&dif2);
reg:
      for(k=0; k<3; ++k)
         prevpix[k] = ColourRegister[pen][k];

      regist = TRUE;
   }
   else
   {
      /* What change from the last pixel?                               */
      dif = ColourDistance(pix,prevpix);
      if(dif)
      {
         /* Which register is the nearest?                              */
         pen = NearestPen(pix,&dif2);
         if(dif2 < dif)
            goto reg;
         else
            regist = FALSE;
      }
      else
      {
         /* See if last pixel was a register                            */
         if(regist)
         {
            pen = NearestPen(pix,&dif2);
            goto reg;
         }
      }

      id = maxdif = 0;
      
      for(k=0; k<3; k++)
      {
         dif = pix[k] - prevpix[k];
         if(dif < 0) dif = (-dif);
         if(dif > maxdif)
         {
            maxdif=dif;
            id=k;
         }
      }

      pen = map[id]+pix[id];
      prevpix[id] = pix[id];        /* Use a HAM pixel                  */
   }
   
   SetAPen(rp,pen);
   WritePixel(rp,i,j);
}
 
/************************************************************************/
/*>int NearestPen(int *c, int *dist)
   ---------------------------------
   Return the pen nearest to the color c.
   Allocate a new register if necessary and possible.
   19.07.93 Original taken from CPK    By: ACRM
*/
int NearestPen(int *c, int *dist)
{
   int   i,
         mindist,
         nearest,
         d;
 
   mindist=32000;
 
   for(i=0; i<nallocr; ++i)
   {
      d=ColourDistance2(c,ColourRegister[i]);
      if(d < mindist)
      {
         mindist=d;
         nearest=i;
      }
   }
   
   if(mindist > THRESHHOLD && nallocr < 16)
   {
      for(i=0; i<3; ++i)
         ColourRegister[nallocr][i]=c[i];
      SetRGB4(viewport,nallocr, c[0],c[1],c[2]);
      nearest = nallocr++;
      mindist = 0;
   }
   
   *dist = mindist;

   return(nearest);
}
 
/************************************************************************/
/*>int ColourDistance(int *a, int *b)
   ----------------------------------
   The 'distance' between two colours.
   After we correct the worst colour.
   19.07.93 Original taken from CPK    By: ACRM
*/
int ColourDistance(int *a, int *b)
{
   int   k,
         r,
         d,
         m;
 
   for(r=k=m=0; k<3; ++k)
   {
      d = a[k]-b[k];
      if(d < 0) d = -d;

      r += d;
      if(d > m) m = d;
   }
   
   return(r-m);
}
 
 
/************************************************************************/
/*>int ColourDistance2(int *a, int *b)
   -----------------------------------
   Calculate the the 'distance' between two colors with no correction
   19.07.93 Original taken from CPK    By: ACRM
*/
int ColourDistance2(int *a, int *b)
{
   int   k,
         r,
         d;
 
   for(r=k=0; k<3; ++k)
   {
      d = a[k]-b[k];
      if(d < 0) d= -d;

      r += d;
   }

   return(r);
}

#endif         /*------------------ END OF AMIGA CODE ------------------*/

