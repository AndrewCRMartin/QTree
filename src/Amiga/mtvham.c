/*************************************************************************

   Program:    MTVHam
   File:       MTVHAM.c
   
   Version:    V1.0
   Date:       17.08.93
   Function:   Display an MTV file on an Amiga HAM screen
   
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
   Reads an MTV file and displays it as a HAM image. Somewhat memory 
   hungry since the image is read into arrays first. This is not really 
   necessary, the image could just be placed straight on the screen, but 
   this was knocked up in a hurry!

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  17.08.93 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>

#include <intuition/intuition.h>
#include <proto/all.h>

#include "bioplib/SysDefs.h"
#include "bioplib/array.h"

/************************************************************************/
/* Prototypes
*/
int  main(int, char **);
BOOL InitGraphics(void);
void EndGraphics(void);
void DispMTV(int, int);
void HAM2(int x, int y, int pix[3]);
int  NearestPen(int *, int *);
int  ColourDistance(int *, int *);
int  ColourDistance2(int *, int *);
BOOL ReadMTV(FILE *, int *, int *);
BOOL AllocArrays(int, int);
void CleanExit(int, int);

/************************************************************************/
/* Defines
*/
#define THRESHHOLD   4        /* Register allocation threshold          */
#define XSIZE        320      /* Screen dimensions                      */
#define YSIZE        256

/************************************************************************/
/* Amiga graphics setup
*/
struct IntuitionBase *IntuitionBase = NULL;
struct GfxBase       *GfxBase       = NULL;
struct Screen        *scrn          = NULL;
struct Window        *wind          = NULL;
struct RastPort      *rp            = NULL;
struct ViewPort      *viewport      = NULL;
int                  nallocr        = 0;

/************************************************************************/
/* Arrays to contain the image
*/
static unsigned char __far **sRed   = NULL,
                           **sGreen = NULL,
                           **sBlue  = NULL;
               
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
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for mtvham.
   17.08.93 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE  *fp   = NULL;
   int   xres  = 0,
         yres  = 0;

   /* Check command line                                                */
   if(argc != 2)
   {
      printf("Usage: mtvham <file.mtv>\n");
      CleanExit(0,0);
   }
   
   /* Open file                                                         */
   if((fp=fopen(argv[1],"r"))==NULL)
   {
      printf("Unable to open MTV file: %s\n",argv[1]);
      CleanExit(0,0);
   }
   
   /* Read file                                                         */
   if(ReadMTV(fp, &xres, &yres))
   {
      /* Check screen dimensions are OK                                 */
      if(xres > XSIZE || yres > YSIZE)
      {
         printf("MTV image dimensions (%d x %d) too large for screen.\n",
                xres, yres);
         CleanExit(xres, yres);
      }
      
      /* Open HAM screen                                                */
      if(InitGraphics())
      {
         /* Display image                                               */
         printf(" Displaying...");
         DispMTV(xres, yres);
         printf("Complete\n");
         
         /* End graphics                                                */
         EndGraphics();
      }
   }
   else
   {
      printf("Unable to read file %s as MTV file\n",argv[1]);
   }
   
   CleanExit(xres, yres);
}

/************************************************************************/
/*>BOOL InitGraphics(void)
   -----------------------
   Perform any initialisation for graphics.
   Supports Amiga HAM mode
   19.07.93 Original    By: ACRM
*/
BOOL InitGraphics(void)
{
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

   return(TRUE);
}

/************************************************************************/
/*>void EndGraphics(void)
   ----------------------
   Performs any remaining graphics display routines, waits for return to
   be pressed, and cleans up the graphics. Also writes the graphics to
   a file if specified.
   17.08.93 Original    By: ACRM
*/
void EndGraphics(void)
{
   char  buffer[80];

   /* Wait for return to be pressed                                     */
   printf("Press <return> to exit...");
   gets(buffer);

   if(wind != NULL)           CloseWindow(wind);
   if(scrn != NULL)           CloseScreen(scrn);
   if(GfxBase != NULL)        CloseLibrary((struct Library *)
                                           GfxBase);
   if(IntuitionBase != NULL)  CloseLibrary((struct Library *)
                                           IntuitionBase);
}

/************************************************************************/
/*>void DispMTV(int xsize, int ysize)
   ----------------------------------
   Display an MTV image stored in 3 2D arrays on a HAM screen.
   17.08.93 Original    By: ACRM
*/
void DispMTV(int xsize, int ysize)
{
   int   x,y,pix[3];

   for(y=0;y<ysize;y++)
   {
      for(x=0;x<xsize;x++)
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
}


/************************************************************************/
/*>BOOL ReadMTV(FILE *fp, int *xres, int *yres)
   --------------------------------------------
   Read an MTV file
   
   03.08.93 Original    By: ACRM
   12.08.93 Allocates arrays
*/
BOOL ReadMTV(FILE *fp, int *xres, int *yres)
{
   int   x,
         y;
         
   if(fscanf(fp, "%d %d",xres,yres) != 2) return(FALSE);
   
   if(!AllocArrays(*xres, *yres))
   {
      printf("Unable to allocate memory\n");
      CleanExit(*xres, *yres);
   }
   
   if(getc(fp) != '\n') return(FALSE);
   
   for(y=0; y < (*yres); y++)
   {
      for(x=0; x < (*xres); x++)
      {
         sRed[x][y]   = getc(fp);
         sGreen[x][y] = getc(fp);
         sBlue[x][y]  = getc(fp);
      }
   }
   
   return(TRUE);
}

/************************************************************************/
/*>BOOL AllocArrays(int xsize, int ysize)
   --------------------------------------
   Allocates memory for the arrays filled by ReadMTV()
   12.08.93 Original    By: ACRM
*/
BOOL AllocArrays(int xsize, int ysize)
{
   sRed   = (unsigned char **)Array2D(sizeof(unsigned char),xsize,ysize);
   sGreen = (unsigned char **)Array2D(sizeof(unsigned char),xsize,ysize);
   sBlue  = (unsigned char **)Array2D(sizeof(unsigned char),xsize,ysize);
   
   if(sRed == NULL || sGreen == NULL || sBlue == NULL)
      return(FALSE);
      
   return(TRUE);
}

/************************************************************************/
/*>void CleanExit(int xsize, int ysize)
   ------------------------------------
   Free memory allocated for image arrays and exit.
   
   17.08.93 Original    By: ACRM
*/
void CleanExit(int xsize, int ysize)
{
   if(xsize != 0 || ysize != 0)
   {
      if(sRed   != NULL) FreeArray2D(sRed,   xsize, ysize);
      if(sGreen != NULL) FreeArray2D(sGreen, xsize, ysize);
      if(sBlue  != NULL) FreeArray2D(sBlue,  xsize, ysize);
      
      sRed = sGreen = sBlue = 0;
   }
   
   exit(0);
}

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

