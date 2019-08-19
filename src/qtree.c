/*************************************************************************

   Program:    QTree
   File:       qtree.c
   
   Version:    V3.0
   Date:       19.08.19
   Function:   Use quad-tree algorithm to display a molecule
   
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
   QTree implements the quad tree algorithm for generating shaded 
   pictures.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======
   FarLeftSearch() and FarRightSearch() simply step along the sorted array
   to find the required points. This is actually *more* efficient than
   using a binary search since, at the deeper recursion levels (which 
   take the most time) there are few atoms in the list.
   
   Because of the division into squares, the fastest check for the 
   presence of atoms is to assume they are square (rather than circular).
   Attempts to optimize the search for the front pixel by sorting on z
   thus fail since atoms not really in this pixel get included.
   
   Conditional compilation flags are defined in qtree.h: 
   
      SPEC        - Include specular reflection
      DEPTHCUE    - Perform depth cueing
      SHOW_INFO   - Show run statistics

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
   V1.5  14.09.93 Added sphere scaling option
   V1.6  04.01.94 Fixed bug in argument parsing
   V1.7  28.03.94 Handles SLAB. Applies sphere scaling when radius comes
                  from B-value
   V1.8  09.05.94 Skipped
   V1.9  13.05.94 Skipped
   V1.10 24.06.94 Handles TEMPERATURE - changes in commands.c
   V1.11 04.10.94 With -b, reads radii from occ rather than bval
   V1.12 21.12.94 Improved Usage message
   V1.12a21.12.94 Fixed Usage text for -b
   V2.0  28.03.95 Modified to allow I/O through pipes
                  Ctrl-C handled properly for non-AmigaDOS systems
   V2.1  23.10.95 Changes in commands.c to send warnings and errors to
                  stderr
   V2.1a 06.12.95 Added CHAIN command
   V2.1b 08.02.96 Fixed bug when handling inserts within zones
   V2.1c 18.06.96 Moved InZone() into bioplib as InPDBZone()
   V2.2  14.10.03 Added BOUNDS and RADIUS
   V2.2a 18.10.07 Some cleanup for ANSI C
   V2.3  18.10.07 Added HIGHLIGHT
   V2.4  27.01.15 Various changes for new bioplib to ensure chain and
                  insert are handled as strings
   V2.5  18.08.19 General cleanup and moved into GitHub
   V3.0  19.08.19 Started to add different output formats

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <setjmp.h>
#include <time.h>

#ifdef _AMIGA
#include <dos.h>
#else
#include <signal.h>   /* If you don't have signal.h, it's no great loss */
typedef void (*sighandler_t)(int);
#endif

#include "bioplib/macros.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/help.h"

#define MAIN
#include "qtree.h"

/************************************************************************/
/* Defines
*/
#define DEF_CONTROL  "qtree.def"    /* Default control file             */
#define HELPFILE     "qtree.hlp"    /* Help file                        */

/************************************************************************/
/* Prototypes
*/
#include "qtree.p"
#include "graphics.p"
#include "commands.p"

/************************************************************************/
/* Variables global to this file only
*/
static jmp_buf sSaveUnwind;         /* For setjmp() / longjmp()         */
static BOOL    sOKFlag    = TRUE,   /* Set FALSE on calling longjmp()   */
               sBallStick = FALSE;  /* Default to CPK images            */

#ifdef SHOW_INFO
static int     sNPixels = 0;        /* Number of pixels coloured        */
#endif

#ifdef _AMIGA
/* Version string                                                       */
static unsigned char 
   *sVers="\0$VER: QTree V3.0 - SciTech Software, 1993-2019";
#endif


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main routine for the quad tree PDB space filling program
   
   19.07.93 Original    By: ACRM
   21.07.93 Added depth cue handling
   22.07.93 Calls MapSpheres() separately from CreateSphereList()
   23.07.93 Added default control file
            Added check for -h
            Added onbreak()
            Modified argument parsing method
   26.07.93 Corrected check for -h with new parsing method
   28.07.93 Added Ball & Stick support
            Improved exit if allocation failed.
   29.07.93 Added File output & resolution setting
   12.08.93 Added handling of force option.
            Requires output file to be specified
   14.09.93 Added gSphScale initialisation
   07.10.93 Argument while() loop checks argc
   28.03.94 Initialises gSlab.flag, etc.
   04.10.94 Calls blReadPDBAll() if it's ball & stick, so all occupancy
            atoms are read (radius is now stored in occ rather then
            in bval).
   28.03.95 Moved screen text to stderr.
            Added Quiet handling and stdio I/O
   23.10.95 V2.1
   06.12.95 V2.1a
   08.02.96 V2.1b
   18.02.96 V2.1c
   14.10.03 V2.2
   18.10.07 V2.2a Cast added to onbreak()
   19.08.19 Added PNG output
*/
int main(int argc, char **argv)
{
   PDB      *pdb           = NULL;
   FILE     *fp            = NULL;
   SPHERE   *spheres       = NULL;
   BOOL     DoControl      = FALSE,
            OK             = TRUE,
            DoResolution   = FALSE,
            Quiet          = FALSE;
   int      NAtom          = 0,
            resolution     = 0,
            outFormat      = OUTPUT_MTV;
   char     ControlFile[160],
            InFile[160],
            outFile[160];
            
#ifdef SHOW_INFO
   clock_t  StartTime,
            StopTime;
            
   StartTime = clock();
#endif

   /* Establish a NULL CTRL-C trap                                      */
   onbreak((void *)&CtrlCNoExit);
   
#ifdef DEPTHCUE
   /* Set default depth cueing parameter                                */
   gDepthCue.contrast = 0.75;
#endif

   /* Default scaling                                                   */
   gScale    = 0.9;
   gSphScale = 1.0;

   /* Default to work out boundaries internally                         */
   gBounds.flag = FALSE;
   
   /* Default image size                                                */
   gSize = SIZE;
   
   /* Set flag for non-ball & stick images                              */
   sBallStick = FALSE;
   
   /* Default to asking for automatic calculation of centre             */
   gMidPoint.x = gMidPoint.y = gMidPoint.z = -9999.0;

   /* Set default screen size                                           */
   gScreen[0] = XSIZE;
   gScreen[1] = YSIZE;

   /* Initialise for slabbing                                           */
   gSlab.flag  = FALSE;
   gSlab.z     = 0.0;
   gSlab.depth = 1000.0;

   if(ParseCmdLine(argc, argv, InFile, outFile, &DoControl, ControlFile,
                   &sBallStick, &DoResolution, &resolution, &Quiet,
                   &(gScreen[0]), &(gScreen[1]), &outFormat))
   {
      /* If the resolution flag has been set, calculate the resolution  */
      if(DoResolution)
      {
         int size;
         
         size = MIN(gScreen[0],gScreen[1]);
         if(resolution > size) resolution = size;
         
         /* Ensure size a power of 2                                    */
         for(gSize = 2; gSize<resolution; gSize*=2);
         if(gSize > resolution) gSize /= 2;
      }
      
      /* Make sure gSize is OK (if screen size specified, but not res   */
      while(gSize > MIN(gScreen[0],gScreen[1]))
         gSize /= 2;
      
      /* Set up default lighting condition                              */
      gLight.x    = (REAL)gSize*2;
      gLight.y    = (REAL)gSize*2;
      gLight.z    = (REAL)gSize*5;
      gLight.amb  = 0.3;
      gLight.spec = FALSE;
      
      /* Banner message                                                 */
      if(!Quiet)
      {
         fprintf(stderr,"\nQTree V3.0\n");
         fprintf(stderr,"========== \n");
         fprintf(stderr,"CPK program for PDB files. SciTech Software\n");
         fprintf(stderr,"Copyright (C) 1993-2019 SciTech Software. All \
Rights Reserved.\n");
         fprintf(stderr,"This program is freely distributable providing \
no profit is made in so doing.\n\n");
         fprintf(stderr,"Rendering...");
      }
      
      /* Open file for reading (Modified for V2.0)                      */
      if(InFile[0])
      {
         if((fp=fopen(InFile,"r")) == NULL)
         {
            fprintf(stderr,"Unable to read file: %s\n",InFile);
            exit(1);
         }
      }
      else
      {
         fp = stdin;
      }
      
      if(InitGraphics())
      {
         /* Read the PDB file                                           */
         if(sBallStick)
            pdb = blReadPDBAll(fp, &NAtom);
         else
            pdb = blReadPDB(fp, &NAtom);
         
         if(pdb != NULL)
         {
            /* Convert to sphere list                                   */
            if((spheres = CreateSphereList(pdb, NAtom)) != NULL)
            {
               /* Handle control file if specified                      */
               if(DoControl) 
                  HandleControl(ControlFile, pdb, spheres, NAtom, TRUE);
               else
                  HandleControl(DEF_CONTROL, pdb, spheres, NAtom, FALSE);
               
               /* Set and scale coords in sphere list                   */
               MapSpheres(pdb, spheres, NAtom);
               
               /* Remove spheres outside slab range                     */
               if(gSlab.flag)
                  spheres = SlabSphereList(spheres, &NAtom);
               
               /* Free memory of PDB linked list                        */
               FREELIST(pdb, PDB);
               pdb = NULL;
               
               /* Run the space fill                                    */
               if(!SpaceFill(spheres, NAtom))
               {
                  fprintf(stderr,"Memory allocation failed or Ctrl-C \
pressed.\n");
                  OK = FALSE;
               }
               
               /* Free the allocated space                              */
               free(spheres);
            }
            else
            {
               fprintf(stderr,"Unable to allocate memory for initial \
sphere list.\n");
               OK = FALSE;
            }
            
            /* Free PDB linked list if CreateSphereList() failed        */
            if(pdb != NULL)
               FREELIST(pdb, PDB);
         }
         else
         {
            fprintf(stderr,"Unable to read atoms from PDB file.\n");
            OK = FALSE;
         }
      }
      
      if(OK && !Quiet) 
         fprintf(stderr,"Complete.\n");
      
#ifdef SHOW_INFO
      StopTime = clock();
#endif
      
      EndGraphics(outFile, outFormat);
      
#ifdef SHOW_INFO
      if(OK && !Quiet)
      {
         fprintf(stderr,"Pixel coverage: %.3f\n",
                 (double)sNPixels/(double)(gSize*gSize));
         fprintf(stderr,"CPU Time:       %.3f seconds\n",
                 (double)(StopTime-StartTime)/CLOCKS_PER_SEC);
      }
#endif
   }
   else
   {
      UsageExit(FALSE);
   }

   return(0);  
}


/************************************************************************/
/*>void MapSpheres(PDB *pdb, SPHERE *spheres, int NSphere)
   -------------------------------------------------------
   Maps the pdb linked list into the spheres array.

   22.07.93 Original (Extracted from CreateSphereList())    By: ACRM
   23.07.93 Changed to place centre of molecule on the 0 z-plane.
            Added worm support
   28.07.93 Added ball & stick support
   14.09.93 Multiply radii by sphere scaling factor
   29.03.94 Applies z scaling to the slab data as well.
            Radius multiplied by gSphScale when using bval
   04.10.94 Sphere radius taken from oxx rather than bval
   14.10.03 Added BOUNDS and RADII stuff
*/
void MapSpheres(PDB *pdb, SPHERE *spheres, int NSphere)
{
   PDB   *p;
   int   i;
   REAL  xmin, xmax,
         ymin, ymax,
         zmin, zmax,
         size;
   BOOL  found;
   RADII *r;
   
    
   /* If we have specified boundaries, use those, otherwise initialize
      boundaries so we can find them from the coordinates
   */
   if(gBounds.flag)
   {
      xmin = gBounds.xmin;
      xmax = gBounds.xmax;
      ymin = gBounds.ymin;
      ymax = gBounds.ymax;
      zmin = gBounds.zmin;
      zmax = gBounds.zmax;
   }
   else
   {
      xmin = xmax = pdb->x;
      ymin = ymax = pdb->y;
      zmin = zmax = pdb->z;
   }

   /* Transfer coordinates and find limits of size                      */
   for(p=pdb,i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      spheres[i].x  = p->x;
      spheres[i].y  = p->y;
      spheres[i].z  = p->z;
   
      if(sBallStick)
      {
         spheres[i].rad = p->occ * gSphScale;
      }
      else
      {
         /* If we have specified atom radii, see if this atom is in the
            list
         */
         found = FALSE;
         if(gRadii)
         {
            for(r=gRadii; r!=NULL; NEXT(r))
            {
               if(((r->atnam[0] == '\0') ||
                   !strncmp(r->atnam, p->atnam, 4)) &&
                  ((r->resnam[0] == '\0') || 
                   !strncmp(r->resnam, p->resnam, 4)))
               {
                  spheres[i].rad = r->radius;
                  found = TRUE;
                  break;
               }
            }
         }

         /* If not, use default values                                  */
         if(!found)
         {
            switch(p->atnam[0])
            {
            case 'C':                        /* Carbon:      r = 1.7    */
               spheres[i].rad    = 1.7  * gSphScale;
               break;
            case 'N':                        /* Nitrogen:    r = 1.7    */
               spheres[i].rad    = 1.7  * gSphScale;
               break;
            case 'O':                        /* Oxygen:      r = 1.35   */
               spheres[i].rad    = 1.35 * gSphScale;
               break;
            case 'S':                        /* Sulphur:     r = 1.7    */
               spheres[i].rad    = 1.7  * gSphScale;
               break;
            case 'H':                        /* Hydrogen:    r = 1.0    */
               spheres[i].rad    = 1.0  * gSphScale;
               break;
            case 'W':                        /* Worm sphere: r = 0.75   */
               if(!strncmp(p->atnam,"WRM ",4))
                  spheres[i].rad = 0.75 * gSphScale;
               break;
            default:                         /* Default:     r = 1.7    */
               spheres[i].rad    = 1.7  * gSphScale;
               break;
            }
         }
      }
      
      if(!gBounds.flag)
      {
         if(spheres[i].x > xmax) xmax = spheres[i].x + spheres[i].rad;
         if(spheres[i].x < xmin) xmin = spheres[i].x - spheres[i].rad;
         if(spheres[i].y > ymax) ymax = spheres[i].y + spheres[i].rad;
         if(spheres[i].y < ymin) ymin = spheres[i].y - spheres[i].rad;
         if(spheres[i].z > zmax) zmax = spheres[i].z + spheres[i].rad;
         if(spheres[i].z < zmin) zmin = spheres[i].z - spheres[i].rad;
      }
   }

      
   /* If midpoint undefined, calculate it                               */
   if(gMidPoint.x == -9999.0 && 
      gMidPoint.y == -9999.0 &&
      gMidPoint.z == -9999.0)
   {
      /* Calculate centre of molecule                                   */
      gMidPoint.x = xmin + (xmax - xmin) / 2.0;
      gMidPoint.y = ymin + (ymax - ymin) / 2.0;
      gMidPoint.z = zmin + (zmax - zmin) / 2.0;
   }

   size = MAX((xmax-xmin),(ymax-ymin));
      
   /* Move atoms to centre picture and scale                            */
   for(i=0; i<NSphere; i++)
   {
      spheres[i].x    -= gMidPoint.x;
      spheres[i].x    *= gSize * gScale / size;
      spheres[i].x    += gSize / 2;
      
      spheres[i].y    -= gMidPoint.y;
      spheres[i].y    *= gSize * gScale / size;
      spheres[i].y    += gSize / 2;
      
      spheres[i].z    -= gMidPoint.z;
      spheres[i].z    *= gSize * gScale / size;

      spheres[i].rad  *= gSize * gScale / size;

      spheres[i].xmax = spheres[i].x + spheres[i].rad;
      spheres[i].xmin = spheres[i].x - spheres[i].rad;
      spheres[i].ymax = spheres[i].y + spheres[i].rad;
      spheres[i].ymin = spheres[i].y - spheres[i].rad;
   }

   /* Apply z scaling to the Slab information                           */
   gSlab.z     -= gMidPoint.z;
   gSlab.z     *= gSize * gScale / size;
   gSlab.depth *= gSize * gScale / size;
   
#ifdef DEPTHCUE
   /* Calculate values for depth cueing                                 */
   zmin             -= gMidPoint.z;
   gDepthCue.ZMin    = zmin * gSize * gScale / size;
   
   zmax             -= gMidPoint.z;
   zmax             *= gSize * gScale / size;
   
   gDepthCue.ZRange  = zmax - gDepthCue.ZMin;
#endif
}


/************************************************************************/
/*>SPHERE *CreateSphereList(PDB *pdb, int NAtom)
   ---------------------------------------------
   Convert a PDB linked list to an array of sphere structures. Simply
   creates the array and fills in default colour and specular shading 
   data.

   19.07.93 Original    By: ACRM
   20.07.93 Added shine and metallic setting
   21.07.93 Pre-calc x and y bounds for each sphere
            Corrected centering. Proper radii.
            Added depth cue handling
   22.07.93 Separated out MapSpheres()
   17.10.07 Sets .highlight
*/
SPHERE *CreateSphereList(PDB *pdb, int NAtom)
{
   PDB      *p    = NULL;
   SPHERE   *sp   = NULL;
   int      NSphere;
   
   if((sp = (SPHERE *)malloc(NAtom * sizeof(SPHERE))) != NULL)
   {
      /* Set colour info in sphere list                                 */
      for(p=pdb,NSphere=0; p!=NULL; NEXT(p), NSphere++)
      {
         sp[NSphere].set         = FALSE;
         sp[NSphere].highlight   = 0;
         
         switch(p->atnam[0])
         {
         case 'C':                        /* Carbon: white              */
            sp[NSphere].r        = 1.0;
            sp[NSphere].g        = 1.0;
            sp[NSphere].b        = 1.0;
            break;
         case 'N':                        /* Nitrogen: blue             */
            sp[NSphere].r        = 0.0;
            sp[NSphere].g        = 0.0;
            sp[NSphere].b        = 1.0;
            break;
         case 'O':                        /* Oxygen: red                */
            sp[NSphere].r        = 1.0;
            sp[NSphere].g        = 0.0;
            sp[NSphere].b        = 0.0;
            break;
         case 'S':                        /* Sulphur: yellow            */
            sp[NSphere].r        = 1.0;
            sp[NSphere].g        = 1.0;
            sp[NSphere].b        = 0.0;
            break;
         case 'H':                        /* Hydrogen: grey             */
            sp[NSphere].r        = 0.6;
            sp[NSphere].g        = 0.6;
            sp[NSphere].b        = 0.6;
            break;
         default:                         /* Default: green             */
            sp[NSphere].r        = 0.0;
            sp[NSphere].g        = 1.0;
            sp[NSphere].b        = 0.0;
            break;
         }

         sp[NSphere].shine    = 0.2; /* Larger value makes brighter spot*/
         sp[NSphere].metallic = 8;   /* Larger value makes smaller spot */
      }
   }
   
   return(sp);
}


/************************************************************************/
/*>BOOL SpaceFill(SPHERE *AllSpheres, int NSphere)
   -----------------------------------------------
   Main entry point for quad tree space filling given an array of SPHERE
   structures. Takes the array and the number of spheres as parameters.
   Returns TRUE or FALSE to indicate success or failure.

   19.07.93 Original    By: ACRM
   20.07.93 Added onbreak() for Amiga. Corrected call to SplitPic() to
            use NSphOut, not NSphere
   18.10.07 Added casts to onbreak()
*/
BOOL SpaceFill(SPHERE *AllSpheres, int NSphere)
{
   int      NSphOut     = 0;
   SPHERE   **SrtSph    = NULL,
            **spheres   = NULL;

   
   /* Assume all OK (no error has occurred)                             */
   sOKFlag = TRUE;
   
   /* Do a setjmp() so we can unwind from errors                        */
   setjmp(sSaveUnwind);
   
   /* Establish a CTRL-C trap                                           */
   onbreak((void *)&CtrlCExit);
   
   /* sOKFlag is cleared if an error occurs during the recursive 
      splitting process and will only be set if we have returned here 
      via longjmp()
   */
   if(sOKFlag)
   {
      /* Create an index into AllSpheres sorted on x                    */
      if((SrtSph = SortSpheresOnX(AllSpheres, NSphere)) != NULL)
      {
         /* Extract list which is within the bounds of the screen       */
         if((spheres = UpdateSphereList((REAL)0,    (REAL)0,
                                        (REAL)gSize, (REAL)gSize,
                                        SrtSph,     NSphere,
                                        &NSphOut)) != NULL)
         {
            /* Call recursive quad-tree routine.                        */
            SplitPic(0,0,gSize,gSize,spheres,NSphOut);
         }
      }
   }
   
   /* Establish a NULL CTRL-C trap                                      */
   onbreak((void *)&CtrlCNoExit);
   
   /* Free memory                                                       */
   if(spheres != NULL)  free(spheres);
   if(SrtSph  != NULL)  free(SrtSph);
   
   return(sOKFlag);
}


/************************************************************************/
/*>void SplitPic(int x0, int y0, int x1, int y1, SPHERE **spheres,
                 int NSphere)
   ---------------------------------------------------------------
   This is the recursive quad-tree algorithm. Takes the top left and
   bottom right of the pixel block and the current array of sphere 
   structure pointers. Splits the block into 4 quadrants; for each, 
   updates the sphere list and if any spheres are present, recurses. 
   If the picture contains only 1 pixel, the shading algorithm is called 
   instead, ending the recursion
   
   19.07.93 Original    By: ACRM
   20.07.93 Added chkabort() for Amiga
   21.07.93 Cast x0 and y0 to REAL in call to ColourPixel
*/
void SplitPic(int x0, int y0, int x1, int y1, SPHERE **spheres,
              int NSphere)
{
   SPHERE   **SplitSpheres = NULL;
   int      xm,
            ym,
            NSphOut;

#ifdef _AMIGA
   chkabort();
#endif

   /* Check for remaining pixel to be coloured                          */
   if(x1-x0 == 1 && y1-y0 == 1)
   {
      ColourPixel(x0, y0, spheres, NSphere);
   }
   else
   {
      /* Find mid point of current square                               */
      xm = x0 + (x1-x0)/2;
      ym = y0 + (y1-y0)/2;
      
      /* Recurse for the top left quadrants                             */
      if((SplitSpheres = UpdateSphereList((REAL)x0, (REAL)y0,
                                          (REAL)xm, (REAL)ym,
                                          spheres,  NSphere,
                                          &NSphOut)) != NULL)
      {
         SplitPic(x0,y0,xm,ym,SplitSpheres,NSphOut);
         free(SplitSpheres);
      }
      
      /* Recurse for the top right quadrants                            */
      if((SplitSpheres = UpdateSphereList((REAL)xm, (REAL)y0,
                                          (REAL)x1, (REAL)ym,
                                          spheres,  NSphere,
                                          &NSphOut)) != NULL)
      {
         SplitPic(xm,y0,x1,ym,SplitSpheres,NSphOut);
         free(SplitSpheres);
      }
      
      /* Recurse for the bottom left quadrants                          */
      if((SplitSpheres = UpdateSphereList((REAL)x0, (REAL)ym,
                                          (REAL)xm, (REAL)y1,
                                          spheres,  NSphere,
                                          &NSphOut)) != NULL)
      {
         SplitPic(x0,ym,xm,y1,SplitSpheres,NSphOut);
         free(SplitSpheres);
      }
      
      /* Recurse for the bottom right quadrants                         */
      if((SplitSpheres = UpdateSphereList((REAL)xm, (REAL)ym,
                                          (REAL)x1, (REAL)y1,
                                          spheres,  NSphere,
                                          &NSphOut)) != NULL)
      {
         SplitPic(xm,ym,x1,y1,SplitSpheres,NSphOut);
         free(SplitSpheres);
      }
   }
}


/************************************************************************/
/*>SPHERE **UpdateSphereList(REAL x0, REAL y0, REAL x1, REAL y1, 
                             SPHERE **spheres, int NSphere, int *NSphOut)
   ----------------------------------------------------------------------
   Take screen coordinates (x0,y0)--(x1,y1) and an array of sphere 
   pointers. Returns a new array of sphere pointers for those spheres
   which are in the bounds of the screen coordinates. If no atoms are in 
   range, returns NULL instead.
   
   19.07.93 Original    By: ACRM
   20.07.93 Made `in' a register int
   23.07.93 Added swap of offsets if needed.
*/
SPHERE **UpdateSphereList(REAL   x0, 
                          REAL   y0, 
                          REAL   x1, 
                          REAL   y1, 
                          SPHERE **spheres,
                          int    NSphere,
                          int    *NSphOut)
{
   int            LOffset,
                  ROffset;
   register int   in;
   SPHERE         **sp = NULL;


   /* Binary search for far left sphere                                 */
   LOffset = FarLeftSearch(spheres,  NSphere, x0);
   ROffset = FarRightSearch(spheres, NSphere, x1);

   /* Check to see if either is out of range                            */
   if(LOffset == (-1) || ROffset == (-1))
      return(NULL);

   /* Swap them if the sphere diameters has resulted in positions being
      reversed
   */
   if(ROffset < LOffset)
   {
      int temp;
      
      temp    = ROffset;
      ROffset = LOffset;
      LOffset = temp;
   }
      
   /* Find number of spheres in range                                   */
   NSphere = ROffset - LOffset + 1;
   
   /* Allocate memory for all spheres in range                          */
   if((sp = (SPHERE **)malloc(NSphere * sizeof(SPHERE *))) == NULL)
   {
      /* Jump out through the recursion if we ran out of memory         */
      sOKFlag = FALSE;
      longjmp(sSaveUnwind,1);
   }
   
   /* Copy in the pointers if y's are in range                          */
   *NSphOut = 0;
   for(in = LOffset; in<=ROffset; in++)
   {
      if((spheres[in])->ymax >= y0 && 
         (spheres[in])->ymin <= y1)
         sp[(*NSphOut)++] = spheres[in];
   }
   
   if(*NSphOut)
   {
      return(sp);
   }
   else
   {
      free(sp);
      return(NULL);
   }
}


/************************************************************************/
/*>SPHERE **SortSpheresOnX(SPHERE *AllSpheres, int NSphere)
   --------------------------------------------------------
   Performs a heapsort returning an array of pointers which represent the 
   input array AllSpheres sorted on x.
   
   19.07.93 Original    By: ACRM
*/
SPHERE **SortSpheresOnX(SPHERE *AllSpheres, int NSphere)
{
   SPHERE   **sp = NULL,
            *temp;
   int      i, j,
            l, ir;
   REAL     q;
   
   /* Return NULL if no atoms                                           */
   if(NSphere == 0) return(NULL);
   
   /* Allocate memory for index                                         */
   if((sp = (SPHERE **)malloc(NSphere * sizeof(SPHERE *))) == NULL)
      return(NULL);

   /* Set all pointers to input ordering                                */
   for(j=0; j<NSphere; j++) 
      sp[j] = (AllSpheres + j);

   if(NSphere > 1)
   {
      l  = NSphere/2 + 1;
      ir = NSphere;
   
      for(;;)
      {
         if(l>1)
         {
            temp = sp[--l - 1];
            q    = temp->x;
         }
         else
         {
            temp = sp[ir-1];
            q    = temp->x;
            
            sp[ir-1]=sp[0];
            if(--ir == 1)
            {
               sp[0]=temp;
               return(sp);
            }
         }
         
         i = l;
         j = l+l;
   
         while(j<=ir)
         {
            if(j<ir)
            {
               if((sp[j-1])->x < (sp[j])->x) j++;
            }
            
            if(q<(sp[j-1])->x)
            {
               sp[i-1] = sp[j-1];
               i       = j;
               j      += j;
            }
            else
            {
               j       = ir+1;
            }
         }
         sp[i-1]       = temp;
      }
   }
   
   /* Return the sorted array if there was only one atom                */
   return(sp);
}


/************************************************************************/
/*>void ColourPixel(int x, int y, SPHERE **spheres, int NSphere)
   ---------------------------------------------------------------
   Search through the sphere list for this pixel to identify the 
   front-most sphere. When found, call the shading routine.

   19.07.93 Original    By: ACRM
   20.07.93 Made q and z register; Fixed Z search just to look at nearest
            5 centre points.
   21.07.93 Made x and y real
   23.07.93 Removed the Z sorting; always run through the whole list.
   18.10.07 Made x and y ints the cast them inside here
*/
void ColourPixel(int xi, int yi, SPHERE **spheres, int NSphere)
{
   REAL           x, y,
                  MaxZ,
                  junk;
   int            FrontSphere = (-1),
                  sph,
                  xx, yy;
   BOOL           border;

   /* Cast x and y as REALs                                             */
   x = (REAL)xi;
   y = (REAL)yi;

   FrontSphere = FindSphere(x, y, spheres, NSphere, &MaxZ);
   
   /* Shade the pixel                                                   */
   if(FrontSphere != (-1))
   {
      if(spheres[FrontSphere]->highlight)
      {
         border = FALSE;
         for(xx = xi-BORDER_NEIGHBOUR; xx <= xi+BORDER_NEIGHBOUR; xx++)
         {
            for(yy = yi-BORDER_NEIGHBOUR; yy <= yi+BORDER_NEIGHBOUR; yy++)
            {
               sph = FindSphere((REAL)xx, (REAL)yy,
                                   spheres, NSphere, &junk);
               if((sph==(-1)) ||
                  (spheres[sph]->highlight != 
                   spheres[FrontSphere]->highlight))
               {
                  border = TRUE;
                  xx = xi+BORDER_NEIGHBOUR+BORDER_NEIGHBOUR;
                  break;
               }
            }
         }
         
         if(border)
         {
            for(xx = xi-gBorderWidth; xx <= xi; xx++)
            {
               for(yy = yi-gBorderWidth; yy <= yi; yy++)
               {
                  SetPixel(xx, yy, 
                           spheres[FrontSphere]->hr,
                           spheres[FrontSphere]->hg,
                           spheres[FrontSphere]->hb);
               }
            }
         }
         else
         {
            ShadePixel(x, y, MaxZ, spheres[FrontSphere]);
         }
      }
      else
      {
         ShadePixel(x, y, MaxZ, spheres[FrontSphere]);
      }
   }
}

/************************************************************************/
int FindSphere(REAL x, REAL y, SPHERE **spheres, int NSphere, REAL *MaxZ)
{
   REAL           XOff,
                  YOff;
   register REAL  q, z;
   int            i,
                  FrontSphere = (-1);

   /* Search back through the spheres for the nearest z position at this
      pixel.
   */
   for(i=NSphere-1; i>=0; i--)
   {
      XOff = x - spheres[i]->x;
      YOff = y - spheres[i]->y;
      
      q = (spheres[i]->rad * spheres[i]->rad) - 
          (XOff * XOff) -
          (YOff * YOff);
      
      if(q >= 0.0)
      {
         /* Find z for this sphere on this pixel                        */
         z = sqrt(q) + spheres[i]->z;
         
         if(FrontSphere == (-1))
         {
            *MaxZ = z;
            FrontSphere = i;
         }
         else
         {
            if(z > *MaxZ)
            {
               *MaxZ = z;
               FrontSphere = i;
            }
         }
      }
   }
   return(FrontSphere);
}   

/************************************************************************/
/*>int FarLeftSearch(SPHERE **spheres, int NSphere, REAL x)
   --------------------------------------------------------
   Searches along the spheres array (sorted on X) for the first sphere
   which is at least partially to the right of x. i.e. the first sphere
   whose right boundary is >= x

   20.07.93 Original    By: ACRM
   22.07.93 Changed not to make out of bounds test first since this can
            give wrong results with spheres of different sizes
*/
int FarLeftSearch(SPHERE **spheres, int NSphere, REAL x)
{
   register int i;
   
   for(i=0; i<NSphere; i++)
   {
      if(spheres[i]->xmax >= x)
         return(i);
   }

   return(-1);
}


/************************************************************************/
/*>int FarRightSearch(SPHERE **spheres, int NSphere, REAL x)
   ---------------------------------------------------------
   Searches backwards along the spheres array (sorted on X) for the first
   sphere which is at least partially to the left of x. i.e. the first 
   sphere whose left boundary is <= x

   20.07.93 Original    By: ACRM
   22.07.93 Changed not to make out of bounds test first since this can
            give wrong results with spheres of different sizes
*/
int FarRightSearch(SPHERE **spheres, int NSphere, REAL x)
{
   register int i;
   
   for(i=NSphere-1; i>=0; i--)
   {
      if(spheres[i]->xmin <= x)
         return(i);
   }

   return(-1);
}


/************************************************************************/
/*>int CtrlCExit(void)
   -------------------
   Control-C was hit so we longjmp() back out of the recursion

   20.07.93 Original    By: ACRM
*/
int CtrlCExit(void)
{
   sOKFlag = FALSE;
   longjmp(sSaveUnwind,1);

   return(0);
}


/************************************************************************/
/*>int CtrlCNoExit(void)
   ---------------------
   Control-C was hit after completing. Do nothing.

   22.07.93 Original    By: ACRM
*/
int CtrlCNoExit(void)
{
   return(0);
}


#ifndef _AMIGA
/************************************************************************/
/*>void onbreak(void *func)
   ------------------------
   Establish a trap for Ctrl-C under OS other than AmigaDOS.

   If you're system doesn't support signal() properly, just comment out
   the signal() call so this is a NULL routine

   30.03.95 Original    By: ACRM
   18.10.07 Changed (void *) case to (sighander_t)
*/
void onbreak(void *func)
{
   signal((int)SIGINT, (sighandler_t)func);
}
#endif


/************************************************************************/
/*>void ShadePixel(REAL x, REAL y, REAL z, SPHERE *sphere)
   -------------------------------------------------------
   This routine performs the actual work of calculating the colour of a
   pixel and calls the SetPixel() routine to colour the pixel.
   If SPEC is defined, specular reflections will be considered.
   If DEPTHCUE is defined, handle depth cueing.
   
   19.07.93 Original    By: ACRM
   20.07.93 Modified and corrected specular reflection code
   21.07.93 Added depth cue handling
   23.07.93 Added pixel count
*/
void ShadePixel(REAL x, REAL y, REAL z, SPHERE *sphere)
{
   VEC3F L,             /* Vector from surface point to light           */
         N;             /* Surface normal vector                        */
   REAL  dot,
         rr, gg, bb,
         cosval,
         NLen,
         LLen;

#ifdef SPEC
   VEC3F dotN,          /* Vector N scaled by L.N dot produce = N(L.N)  */
         Temp,
         V,             /* Vector from surface point to observer        */
         R;             /* Reflected light ray vector                   */
   REAL  VLen,
         RLen,
         spec = 0.0;
#endif


#ifdef SHOW_INFO
   sNPixels++;
#endif


   /* Calculate surface normal                                          */
   N.x        = x - sphere->x;
   N.y        = y - sphere->y;
   N.z        = z - sphere->z;
   NLen       = sqrt(N.x * N.x +
                     N.y * N.y +
                     N.z * N.z);
   
   /* Calculate vector from surface point to light                      */
   L.x     = gLight.x - x;
   L.y     = gLight.y - y;
   L.z     = gLight.z - z;
   LLen    = sqrt(L.x * L.x +
                  L.y * L.y +
                  L.z * L.z);
   
   /* Calculate angle between surface normal and this vector            */
   dot = N.x * L.x +
         N.y * L.y +
         N.z * L.z;
   cosval = dot/(NLen * LLen);
   
#ifdef DEPTHCUE
   cosval *= (1 - gDepthCue.contrast + 
             gDepthCue.contrast * (z - gDepthCue.ZMin) / gDepthCue.ZRange);
#endif

   if(cosval < 0.0) cosval = 0.0;
   
   /* Calculate diffuse reflection colour components                    */
   rr = sphere->r * ((1.0-gLight.amb)*cosval + gLight.amb);
   gg = sphere->g * ((1.0-gLight.amb)*cosval + gLight.amb);
   bb = sphere->b * ((1.0-gLight.amb)*cosval + gLight.amb);


#ifdef SPEC
   if(gLight.spec)
   {
      /* Now specular reflection                                        */
      dotN.x = dot * N.x;
      dotN.y = dot * N.y;
      dotN.z = dot * N.z;
      Temp.x = L.x - dotN.x;
      Temp.y = L.y - dotN.y;
      Temp.z = L.z - dotN.z;
      R.x    = N.x - Temp.x;
      R.y    = N.y - Temp.y;
      R.z    = N.z - Temp.z;
      RLen   = sqrt(R.x * R.x + 
                    R.y * R.y + 
                    R.z * R.z);
      
      /* Observer                                                       */
      V.x = gSize/2.0 - x;
      V.y = gSize/2.0 - y;
      V.z = gSize*5.0 - z;
      VLen   = sqrt(V.x * V.x + 
                    V.y * V.y + 
                    V.z * V.z);
      
   
      /* Calc dot product of reflected and observer                     */
      dot = R.x * V.x +
            R.y * V.y +
            R.z * V.z;
      cosval = dot/(RLen * VLen);
   
      cosval = 2.0 * cosval * cosval - 1.0;
      if(cosval < 0.00001)
         spec = 0.0;
      else
         spec = sphere->shine * pow(cosval,sphere->metallic);
   
      rr += spec;
      gg += spec;
      bb += spec;
   }
#endif

   if(rr > 1.0) rr = 1.0;
   if(gg > 1.0) gg = 1.0;
   if(bb > 1.0) bb = 1.0;
   
   SetPixel((int)x, (int)y, rr, gg, bb);
}


/************************************************************************/
/*>SPHERE *SlabSphereList(SPHERE *spheres, int *Natom)
   ---------------------------------------------------
   Trim the sphere list by removing atoms outside the slabbing range
   Natom is required both for input and output. The old sphere list is
   freed before the routine returns.
   The routine allocates a new list of the same length as the old list
   which is slightly wasteful of memory, but the old list is freed, so
   it's not a memory loss compared with not calling the routine.
   If memory allocation fails returns the input list unmodified.

   Input:    SPHERE  *spheres   Array of spheres
   I/O:      int     *Natom     Length of input and output sphere lists
   Returns:  SPHERE  *          Updated sphere list

   28.03.94 Original    By: ACRM
   29.03.94 Modified such that any atom which overlaps the slab will
            be included when OVERLAP_SLAB is defined
*/
SPHERE *SlabSphereList(SPHERE *spheres, int *Natom)
{
   SPHERE *spl;
   int    i,
          NOut;
   REAL   SlabMin,
          SlabMax;

   /* Allocate memory for new sphere list                               */
   if((spl = (SPHERE *)malloc(*Natom * sizeof(SPHERE)))==NULL)
      return(spheres);

   /* Calculate bounds of the slab                                      */
   SlabMin = gSlab.z - gSlab.depth/(REAL)2.0;
   SlabMax = gSlab.z + gSlab.depth/(REAL)2.0;

   /* Copy spheres within slab                                          */
   for(i=0, NOut=0; i<(*Natom); i++)
   {
#ifdef OVERLAP_SLAB
      REAL zmin = (spheres[i].z - spheres[i].rad),
           zmax = (spheres[i].z + spheres[i].rad);

      if((zmax >= SlabMin && zmax <= SlabMax) ||
         (zmin >= SlabMin && zmin <= SlabMax) ||
         (zmin <= SlabMin && zmax >= SlabMax))
#else
      if(spheres[i].z >= SlabMin && spheres[i].z <= SlabMax)
#endif
      {
         spl[NOut] = spheres[i];
         NOut++;
      }
   }

   *Natom = NOut;
   return(spl);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *DoControl, char *ControlFile, 
                     BOOL *DoBallStick, 
                     BOOL *DoResolution, int *resolution, BOOL *quiet,
                     int *screenx, int *screeny, int *outFormat)
   ---------------------------------------------------------------------
   Input:   int    argc               Argument count
            char   **argv             Argument array
   Output:  char   *infile            Input file (or blank string)
            char   *outfile           Output file (or blank string)
            BOOL   *DoControl         A control file has been given
            char   *ControlFile       Name of control file
            BOOL   *DoBallStick       Read raidus from occ column
            BOOL   *DoResolution      A resolution has been given
            int    *resolution        The resolution
            BOOL   *Quiet             Operate quietly
            int    *screenx           X image size
            int    *screeny           Y image size
            int    *outFormat         Output format
   Returns: BOOL                      Success?

   Parse the command line
   
   28.03.95 Original    By: ACRM
   19.08.19 Added outFormat
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *DoControl, char *ControlFile, BOOL *DoBallStick, 
                  BOOL *DoResolution, int *resolution, BOOL *quiet,
                  int *screenx, int *screeny, int *outFormat)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(!strncmp(argv[0],"-help",5)  ||
         !strncmp(argv[0],"-HELP",5)  ||
         !strncmp(argv[0],"--help",5) ||
         !strncmp(argv[0],"--HELP",5))
      {
         UsageExit(TRUE);
         return(FALSE);
      }
      
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
         case 'C':
            argc--;  argv++;
            *DoControl = TRUE;
            strcpy(ControlFile,argv[0]);
            break;
         case 'h':
         case 'H':
            return(FALSE);
            break;
         case 'b':
         case 'B':
            *DoBallStick = TRUE;
            break;
         case 'r':
         case 'R':
            argc--;  argv++;
            sscanf(argv[0],"%d",resolution);
            *DoResolution = TRUE;
            break;
         case 'q':
         case 'Q':
            *quiet = TRUE;
            break;
         case 's':
         case 'S':
            argc--;  argv++;
            sscanf(argv[0],"%d",screenx);
            argc--;  argv++;
            sscanf(argv[0],"%d",screeny);
            break;
         case 'f':
         case 'F':
            argc--;  argv++;
            LOWER(argv[0]);
            if(!strncmp(argv[0], "mtv", 3))
            {
               *outFormat = OUTPUT_MTV;
            }
#ifdef SUPPORT_PNG
            else if(!strncmp(argv[0], "png", 3))
            {
               *outFormat = OUTPUT_PNG;
            }
#endif
            else
            {
               fprintf(stderr, "Unknown output format: %s\n", argv[0]);
               exit(1);
            }
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void UsageExit(BOOL ShowHelp)
   -----------------------------
   Print usage info and exit. If ShowHelp set, then display help file
   if present.

   21.07.93 Original (split from main())     By: ACRM
   23.07.93 Added help file support
   28.07.93 Added ball & stick support
   12.08.93 Added -s option
   21.12.94 Added Copyright/version. Corrected text for -b
   28.03.95 Changed to stderr and added -q
   23.10.95 V2.1
   06.12.95 V2.1a
   08.02.96 V2.1b
   18.06.96 V2.1c
   14.10.03 V2.2
   18.10.07 V2.2a
   18.10.07 V2.3
   27.01.15 V2.4
   18.08.19 V2.5
   19.08.19 V3.0
*/
void UsageExit(BOOL ShowHelp)
{
   if(ShowHelp)
   {
      blDoHelp("HELP",HELPFILE);
      blHelp("Dummy","CLOSE");
   }
   else
   {
      fprintf(stderr,"\nQTree V3.0 (c) 1993-2019 Prof. Andrew C.R. \
Martin, SciTech Software\n\n");
      
      fprintf(stderr,"Usage: qtree [-q] [-b] [-c <control.dat>] [-r <n>] \
[-f fmt] [-s <x> <y>] [<file.pdb> [<file.mtv>]]\n");
      fprintf(stderr,"       qtree [--help]\n\n");
      fprintf(stderr,"       -q Operate quietly\n");
      fprintf(stderr,"       -b Interpret occupancy as radius for ball \
& stick\n");
      fprintf(stderr,"       -c Specify control file\n");
      fprintf(stderr,"       -r Specify pixel resolution (power of 2) \
[%d]\n", SIZE);
      fprintf(stderr,"       -s Specify screen size (%d %d)\n",
             XSIZE,YSIZE);
      fprintf(stderr,"       -h Enter help utility\n");
      fprintf(stderr,"       -f Specify output format (mtv|png)\n");
      fprintf(stderr,"          Default output is in MTV raytracer \
format\n\n");
      fprintf(stderr,"       Render a space filling picture of a PDB \
file\n\n");
      fprintf(stderr,"       Enter 'qtree --help' to enter the help \
program\n\n");
   }
   exit(0);
}


