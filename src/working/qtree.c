/*************************************************************************

   Program:    QTree
   File:       qtree.c
   
   Version:    V1.4
   Date:       12.08.93
   Function:   Use quad-tree algorithm to display a molecule
   
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
static unsigned char *sVers="\0$VER: QTree V1.4 - SciTech Software, 1993";
#endif

/************************************************************************/
/*>main(int argc, char **argv)
   ---------------------------
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
*/
main(int argc, char **argv)
{
   PDB      *pdb           = NULL;
   FILE     *fp            = NULL;
   SPHERE   *spheres       = NULL;
   BOOL     DoControl      = FALSE,
            OK             = TRUE,
            DoResolution   = FALSE;
   int      NAtom          = 0,
            resolution     = 0;
   char     ControlFile[160];
            
#ifdef SHOW_INFO
   clock_t  StartTime,
            StopTime;
            
   StartTime = clock();
#endif

#ifdef _AMIGA
   /* Establish a NULL CTRL-C trap                                      */
   onbreak(&CtrlCNoExit);
#endif
   
#ifdef DEPTHCUE
   /* Set default depth cueing parameter                                */
   gDepthCue.contrast = 0.75;
#endif

   /* Default scaling                                                   */
   gScale = 0.9;
   
   /* Default image size                                                */
   gSize = SIZE;
   
   /* Set flag for non-ball & stick images                              */
   sBallStick = FALSE;
   
   /* Default to asking for automatic calculation of centre             */
   gMidPoint.x = gMidPoint.y = gMidPoint.z = -9999.0;

   /* Set default screen size                                           */
   gScreen[0] = XSIZE;
   gScreen[1] = YSIZE;


   /* Parse the command line                                            */
   argc--;  argv++;
   
   while(argv[0][0] == '-')
   {
      switch(argv[0][1])
      {
      case 'c':
      case 'C':
         argc--;  argv++;
         DoControl = TRUE;
         strcpy(ControlFile,argv[0]);
         break;
      case 'h':
      case 'H':
         UsageExit(TRUE);
         break;
      case 'b':
      case 'B':
         sBallStick = TRUE;
         break;
      case 'r':
      case 'R':
         argc--;  argv++;
         sscanf(argv[0],"%d",&resolution);
         DoResolution = TRUE;
         break;
      case 's':
      case 'S':
         argc--;  argv++;
         sscanf(argv[0],"%d",&(gScreen[0]));
         argc--;  argv++;
         sscanf(argv[0],"%d",&(gScreen[1]));
         break;
      default:
         break;
      }
      
      argc--;  argv++;
   }

   if(argc != 2)
      UsageExit(FALSE);

   /* Copy in name of output filename                                   */
   strcpy(gOutFile,argv[1]);

   /* If the resolution flag has been set, calculate the resolution     */
   if(DoResolution)
   {
      int size;
      
      size = MIN(gScreen[0],gScreen[1]);
      if(resolution > size) resolution = size;

      /* Ensure size a power of 2                                       */
      for(gSize = 2; gSize<resolution; gSize*=2);
      if(gSize > resolution) gSize /= 2;
      
   }

   /* Make sure gSize is OK (if screen size specified, but not res      */
   while(gSize > MIN(gScreen[0],gScreen[1]))
      gSize /= 2;

   /* Set up default lighting condition                                 */
   gLight.x    = (REAL)gSize*2;
   gLight.y    = (REAL)gSize*2;
   gLight.z    = (REAL)gSize*5;
   gLight.amb  = 0.3;
   gLight.spec = FALSE;
   
   /* Banner message                                                    */
   printf("\nQTree V1.4\n");
   printf("==========\n");
   printf("CPK program for PDB files. SciTech Software\n");
   printf("Copyright (C) 1993 SciTech Software. All Rights Reserved.\n");
   printf("This program is freely distributable providing no profit is \
made in so doing.\n\n");
   printf("Rendering...");
      
   /* Open file for reading                                             */
   if((fp=fopen(argv[0],"r")) == NULL)
   {
      printf("Unable to read file: %s\n",argv[0]);
      exit(0);
   }

   if(InitGraphics())
   {
      /* Read the PDB file                                              */
      if((pdb = ReadPDB(fp, &NAtom)) != NULL)
      {
         /* Convert to sphere list                                      */
         if((spheres = CreateSphereList(pdb, NAtom)) != NULL)
         {
            /* Handle control file if specified                         */
            if(DoControl) 
               HandleControl(ControlFile, pdb, spheres, NAtom, TRUE);
            else
               HandleControl(DEF_CONTROL, pdb, spheres, NAtom, FALSE);
      
            /* Set and scale coords in sphere list                      */
            MapSpheres(pdb, spheres, NAtom);
            
            /* Free memory of PDB linked list                           */
            FREELIST(pdb, PDB);
            pdb = NULL;
            
            /* Run the space fill                                       */
            if(!SpaceFill(spheres, NAtom))
            {
               printf("Memory allocation failed or Ctrl-C pressed.\n");
               OK = FALSE;
            }
            
            /* Free the allocated space                                 */
            free(spheres);
         }
         else
         {
            printf("Unable to allocate memory for initial sphere \
list.\n");
            OK = FALSE;
         }
         
         /* Free PDB linked list if CreateSphereList() failed           */
         if(pdb != NULL)
            FREELIST(pdb, PDB);
      }
      else
      {
         printf("Unable to read atoms from PDB file.\n");
         OK = FALSE;
      }
   }
   
   if(OK) printf("Complete.\n");
   
#ifdef SHOW_INFO
   StopTime = clock();
#endif
   
   EndGraphics();
   
#ifdef SHOW_INFO
   if(OK)
   {
      printf("Pixel coverage: %.3lf\n",
             (double)sNPixels/(double)(gSize*gSize));
      printf("CPU Time:       %.3lf seconds\n",
             (double)(StopTime-StartTime)/CLOCKS_PER_SEC);
   }
#endif
}


/************************************************************************/
/*>void MapSpheres(PDB *pdb, SPHERE *spheres, int NSphere)
   -------------------------------------------------------
   Maps the pdb linked list into the spheres array.

   22.07.93 Original (Extracted from CreateSphereList())    By: ACRM
   23.07.93 Changed to place centre of molecule on the 0 z-plane.
            Added worm support
   28.07.93 Added ball & stick support
*/
void MapSpheres(PDB *pdb, SPHERE *spheres, int NSphere)
{
   PDB   *p;
   int   i;
   REAL  xmin, xmax,
         ymin, ymax,
         zmin, zmax,
         size;
         
   xmin = xmax = pdb->x;
   ymin = ymax = pdb->y;
   zmin = zmax = pdb->z;

   /* Transfer coordinates and find limits of size                      */
   for(p=pdb,i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      spheres[i].x  = p->x;
      spheres[i].y  = p->y;
      spheres[i].z  = p->z;
   
      if(sBallStick)
      {
         spheres[i].rad = p->bval;
      }
      else
      {
         switch(p->atnam[0])
         {
         case 'C':                           /* Carbon: r = 1.7         */
            spheres[i].rad      = 1.7;
            break;
         case 'N':                           /* Nitrogen: r = 1.7       */
            spheres[i].rad      = 1.7;
            break;
         case 'O':                           /* Oxygen: r = 1.35        */
            spheres[i].rad      = 1.35;
            break;
         case 'S':                           /* Sulphur: r = 1.7        */
            spheres[i].rad      = 1.7;
            break;
         case 'H':                           /* Hydrogen: r = 1.0       */
            spheres[i].rad      = 1.0;
            break;
         case 'W':                           /* Worm sphere: r = 0.75   */
            if(!strncmp(p->atnam,"WRM ",4))
               spheres[i].rad   = 0.75;
            break;
         default:                            /* Default: r = 1.7        */
            spheres[i].rad      = 1.7;
            break;
         }
      }
      
      if(spheres[i].x > xmax) xmax = spheres[i].x + spheres[i].rad;
      if(spheres[i].x < xmin) xmin = spheres[i].x - spheres[i].rad;
      if(spheres[i].y > ymax) ymax = spheres[i].y + spheres[i].rad;
      if(spheres[i].y < ymin) ymin = spheres[i].y - spheres[i].rad;
      if(spheres[i].z > zmax) zmax = spheres[i].z + spheres[i].rad;
      if(spheres[i].z < zmin) zmin = spheres[i].z - spheres[i].rad;
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
   
#ifdef _AMIGA
   /* Establish a CTRL-C trap                                           */
   onbreak(&CtrlCExit);
#endif
   
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
   
#ifdef _AMIGA
   /* Establish a NULL CTRL-C trap                                      */
   onbreak(&CtrlCNoExit);
#endif
   
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
      ColourPixel((REAL)x0, (REAL)y0, spheres, NSphere);
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
/*>void ColourPixel(REAL x, REAL y, SPHERE **spheres, int NSphere)
   ---------------------------------------------------------------
   Search through the sphere list for this pixel to identify the 
   front-most sphere. When found, call the shading routine.

   19.07.93 Original    By: ACRM
   20.07.93 Made q and z register; Fixed Z search just to look at nearest
            5 centre points.
   21.07.93 Made x and y real
   23.07.93 Removed the Z sorting; always run through the whole list.
*/
void ColourPixel(REAL x, REAL y, SPHERE **spheres, int NSphere)
{
   REAL           XOff,
                  YOff,
                  MaxZ;
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
            MaxZ = z;
            FrontSphere = i;
         }
         else
         {
            if(z > MaxZ)
            {
               MaxZ = z;
               FrontSphere = i;
            }
         }
      }
          
   }
   
   /* Shade the pixel                                                   */
   if(FrontSphere != (-1))
      ShadePixel(x, y, (REAL)MaxZ, spheres[FrontSphere]);
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


#ifdef _AMIGA
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
/*>void UsageExit(BOOL ShowHelp)
   -----------------------------
   Print usage info and exit. If ShowHelp set, then display help file
   if present.

   21.07.93 Original (split from main())     By: ACRM
   23.07.93 Added help file support
   28.07.93 Added ball & stick support
   12.08.93 Added -s option
*/
void UsageExit(BOOL ShowHelp)
{
   if(ShowHelp)
   {
      DoHelp("HELP",HELPFILE);
      Help("Dummy","CLOSE");
   }
   else
   {
      printf("Usage: qtree [-b] [-c <control.dat>] [-r <n>] [-s <x> <y>] \
<file.pdb> <file.mtv>\n");
      printf("       qtree [-h]\n\n");
      printf("       -b Interpret b-val as radius for ball & stick\n");
      printf("       -c Specify control file\n");
      printf("       -r Specify pixel resolution (power of 2) [%d]\n",
             SIZE);
      printf("       -s Specify screen size [%d %d]\n",XSIZE,YSIZE);
      printf("       -h Enter help utility\n\n");
      printf("       Output is in MTV raytracer format\n\n");
      printf("       Render a space filling picture of a PDB file\n");
      printf("       Enter qtree -h to enter the help program\n");
   }
   exit(0);
}

