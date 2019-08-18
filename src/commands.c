/*************************************************************************

   Program:    QTree
   File:       commands.c
   
   Version:    V2.5
   Date:       18.08.19
   Function:   Handle command files for QTree program
   
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
   V1.3           Skipped
   V1.4           Skipped
   V1.5  14.09.93 Added sphere scaling option
   V1.6  04.01.94 Skipped
   V1.7  24.03.94 Added SLAB option and warning messages.
                  Removed blParseResSpec() as this is now in the library
   V1.8  09.05.94 Skipped
   V1.9  13.05.94 Skipped
   V1.10 24.06.94 Default colouring may now be done on temperature
                  factor using command TEMPERATURE
   V1.11 04.10.94 Skipped
   V1.12 21.12.94 Skipped
   V2.0  28.03.95 Fixed DoZone() which wasn't working properly for
                  Ball and Stick images which split the PDB data into
                  two sections.
   V2.1  23.10.95 Errors and warnings got to stderr
   V2.1a 06.12.95 Added CHAIN command for setting colour
   V2.1b 08.02.96 Fixed bug in specifying zones. InZone() didn't
                  correctly handle residues within a zone which had
                  insertion codes
   V2.1c 18.06.96 Changed InZone() to blInPDBZone() and moved to bioplib
   V2.2  14.10.03 Added BOUNDS and RADIUS commands
   V2.3  18.10.07 Added HIGHLIGHT command
   V2.4  27.01.15 Modifications for new version of BiopLib
   V2.5  18.08.19 General cleanup and moved into GitHub

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/parse.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/matrix.h"

#include "qtree.h"

/************************************************************************/
/* Prototypes
*/
#include "commands.p"
#include "graphics.p"
#include "qtree.p"

/************************************************************************/
/* Parser setup
*/
#define PARSER_MAXSTRPARAM    5
#define PARSER_MAXREALPARAM   9
#define PARSER_MAXSTRLEN      80 /* Must be >= 8; gets padded           */

#define COM_ZONE              0
#define COM_ATOM              1
#define COM_RESIDUE           2
#define COM_DEFAULT           3
#define COM_AMBIENT           4
#define COM_SPEC              5
#define COM_CONTRAST          6
#define COM_PHONG             7
#define COM_SCALE             8
#define COM_ROTATE            9
#define COM_MATRIX            10
#define COM_CENTRE            11
#define COM_CENTER            12
#define COM_XMATRIX           13
#define COM_LIGHT             14
#define COM_BACKGRND          15
#define COM_SPHSCALE          16
#define COM_SLAB              17
#define COM_TEMP              18
#define COM_CHAIN             19
#define COM_BOUNDS            20
#define COM_RADIUS            21
#define COM_HIGHLIGHT         22
#define COM_BORDERWIDTH       23
#define PARSER_NCOMM          24

/************************************************************************/
KeyWd sKeyWords[PARSER_NCOMM];         /* Parser keywords               */
char  *sStrParam[PARSER_MAXSTRPARAM];  /* Parser string parameters      */
REAL  sRealParam[PARSER_MAXREALPARAM]; /* Parser real parameters        */


/************************************************************************/
/*>BOOL SetupParser(void)
   ----------------------
   Set up the command parser
   21.07.93 Original    By: ACRM
   22.07.93 Added Rotate, matrix, and centre/center
   23.07.93 Added light & background
   14.09.93 Added spherescale
   24.03.94 Added SLAB
   24.06.94 Added TEMPERATURE
   06.12.95 Added CHAIN
   14.10.03 Added BOUNDS and RADIUS
   18.10.07 Added HIGHLIGHT
*/
BOOL SetupParser(void)
{
   int i;
   
   /* Allocate memory for the string parameters                         */
   for(i=0; i<PARSER_MAXSTRPARAM; i++)
   {
      if((sStrParam[i] = 
         (char *)malloc(PARSER_MAXSTRLEN * sizeof(char)))==NULL)
      {
         int j;
         
         for(j=0;j<i;j++) 
            free(sStrParam[j]);
         return(FALSE);
      }
   }
   
   /* Set up the keywords                                               */
   MAKEKEY(sKeyWords[COM_ZONE],       "ZONE",        STRING,5);
   MAKEKEY(sKeyWords[COM_ATOM],       "ATOM",        STRING,4);
   MAKEKEY(sKeyWords[COM_RESIDUE],    "RESIDUE",     STRING,4);
   MAKEKEY(sKeyWords[COM_DEFAULT],    "DEFAULT",     NUMBER,3);
   MAKEKEY(sKeyWords[COM_AMBIENT],    "AMBIENT",     NUMBER,1);
   MAKEKEY(sKeyWords[COM_SPEC],       "SPECULAR",    NUMBER,0);
   MAKEKEY(sKeyWords[COM_CONTRAST],   "CONTRAST",    NUMBER,1);
   MAKEKEY(sKeyWords[COM_PHONG],      "PHONG",       NUMBER,2);
   MAKEKEY(sKeyWords[COM_SCALE],      "SCALE",       NUMBER,1);
   MAKEKEY(sKeyWords[COM_ROTATE],     "ROTATE",      STRING,2);
   MAKEKEY(sKeyWords[COM_MATRIX],     "MATRIX",      NUMBER,9);
   MAKEKEY(sKeyWords[COM_CENTRE],     "CENTRE",      STRING,2);
   MAKEKEY(sKeyWords[COM_CENTER],     "CENTER",      STRING,2);
   MAKEKEY(sKeyWords[COM_XMATRIX],    "XMATRIX",     NUMBER,9);
   MAKEKEY(sKeyWords[COM_LIGHT],      "LIGHT",       NUMBER,3);
   MAKEKEY(sKeyWords[COM_BACKGRND],   "BACKGROUND",  NUMBER,6);
   MAKEKEY(sKeyWords[COM_SPHSCALE],   "SPHERESCALE", NUMBER,1);
   MAKEKEY(sKeyWords[COM_SLAB],       "SLAB",        STRING,3);
   MAKEKEY(sKeyWords[COM_TEMP],       "TEMPERATURE", NUMBER,0);
   MAKEKEY(sKeyWords[COM_CHAIN],      "CHAIN",       STRING,4);
   MAKEKEY(sKeyWords[COM_BOUNDS],     "BOUNDS",      NUMBER,6);
   MAKEKEY(sKeyWords[COM_RADIUS],     "RADIUS",      STRING,2);
   MAKEKEY(sKeyWords[COM_HIGHLIGHT],  "HIGHLIGHT",   STRING,5);
   MAKEKEY(sKeyWords[COM_BORDERWIDTH],"BORDERWIDTH", NUMBER,1);
   
   /* Check all allocations OK                                          */
   for(i=0; i<PARSER_NCOMM; i++)
   {
      if(sKeyWords[i].name == NULL)
         return(FALSE);
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void HandleControl(char *file, PDB *pdb, SPHERE *spheres, int NSphere,
                      BOOL ReportError)
   ----------------------------------------------------------------------
   Modify the sphere list based on the control file
   21.07.93 Original    By: ACRM
   22.07.93 Removed Remap flag since mapping is now done after file
            parsing anyway
            Added Rotate, matrix, and centre/center
   23.07.93 Added ReportError; added light & exit
   29.07.93 Added background support
   03.08.93 Corrected call to ApplyMatrixPDB to blRotatePDB
   14.09.93 Added spherescale
   24.03.94 Added SLAB
   24.06.94 Added TEMPERATURE colouring. SetDefault() now has a different
            set of parameters
   23.10.95 Changed error messages to go to stderr
   06.12.95 Added CHAIN
   14.10.03 Added BOUNDS and RADIUS
   18.10.07 Added HIGHLIGHT
*/
void HandleControl(char *file, PDB *pdb, SPHERE *spheres, int NSphere,
                   BOOL ReportError)
{
   FILE  *fp = NULL;
   char  buffer[160],
         CentreRes[16],
         CentreAtom[8],
         SlabRes[16],
         SlabAtom[8];
   int   key,
         i, j,
         nhighlight = 0;
   REAL  DefaultRGB[3],
         matrix[3][3];
   BOOL  SetDefault   = FALSE,
         SetCentre    = FALSE,
         ColourByTemp = FALSE;
   
   /* Default spheres to white                                          */
   for(i=0; i<3; i++)
      DefaultRGB[i] = 1.0;

   if(SetupParser())
   {
      if((fp=fopen(file,"r")) == NULL)
      {
         if(ReportError)
            fprintf(stderr,"Unable to open control file: %s\n",file);

         return;
      }
      
      while(fgets(buffer,159,fp))
      {
         TERMINATE(buffer);
         
         key = blParse(buffer, PARSER_NCOMM, sKeyWords, sRealParam,
                       sStrParam);
         
         switch(key)
         {
         case PARSE_ERRC:
         case PARSE_ERRP:
            fprintf(stderr,"Error in command file line:\n%s\n",buffer);
            break;
         case COM_ZONE:
            DoZone(spheres,pdb,NSphere,sStrParam[0],sStrParam[1],
                           sStrParam[2],sStrParam[3],sStrParam[4],
                           0);
            SetDefault = TRUE;
            break;
         case COM_ATOM:
            DoAtom(spheres,pdb,NSphere,sStrParam[0],sStrParam[1],
                           sStrParam[2],sStrParam[3]);
            SetDefault = TRUE;
            break;
         case COM_RESIDUE:
            DoResidue(spheres,pdb,NSphere,sStrParam[0],sStrParam[1],
                              sStrParam[2],sStrParam[3]);
            SetDefault = TRUE;
            break;
         case COM_DEFAULT:
            for(i=0; i<3; i++)
               DefaultRGB[i] = sRealParam[i];
            SetDefault = TRUE;
            break;
         case COM_AMBIENT:
            gLight.amb = sRealParam[0];
            break;
         case COM_SPEC:
            gLight.spec = TRUE;
            break;
         case COM_CONTRAST:
#ifdef DEPTHCUE
            gDepthCue.contrast = sRealParam[0];
#endif
            break;
         case COM_PHONG:
            DoPhong(spheres,NSphere,sRealParam[0],sRealParam[1]);
            break;
         case COM_SCALE:
            gScale = sRealParam[0];
            break;
         case COM_ROTATE:
            DoRotate(pdb,sStrParam[0],sStrParam[1]);
            break;
         case COM_MATRIX:
            for(i=0; i<3; i++)
            {
               for(j=0; j<3; j++)
               {
                  matrix[i][j] = sRealParam[i*3 + j];
               }
            }
            blRotatePDB(pdb,matrix);
            break;
         case COM_XMATRIX:
            for(i=0; i<3; i++)
            {
               for(j=0; j<3; j++)
               {
                  matrix[j][i] = sRealParam[i*3 + j];
               }
            }
            blRotatePDB(pdb,matrix);
            break;
         case COM_CENTRE:
         case COM_CENTER:
            SetCentre = TRUE;
            strcpy(CentreRes,  sStrParam[0]);
            strcpy(CentreAtom, sStrParam[1]);
            break;
         case COM_LIGHT:
            gLight.x = (REAL)gSize * sRealParam[0];
            gLight.y = (REAL)gSize * sRealParam[1];
            gLight.z = (REAL)gSize * sRealParam[2];
            break;
         case COM_BACKGRND:
            DoBackground(sRealParam[0],sRealParam[1],sRealParam[2],
                         sRealParam[3],sRealParam[4],sRealParam[5]);
            break;
         case COM_SPHSCALE:
            gSphScale = sRealParam[0];
            break;
         case COM_SLAB:
            strcpy(SlabRes,sStrParam[0]);
            strcpy(SlabAtom,sStrParam[1]);
            sscanf(sStrParam[2],"%lf",&(gSlab.depth));
            gSlab.flag = TRUE;
            break;
         case COM_TEMP:
            SetDefault   = TRUE;
            ColourByTemp = TRUE;
            break;
         case COM_CHAIN:
            DoChain(spheres,pdb,NSphere,sStrParam[0],sStrParam[1],
                    sStrParam[2],sStrParam[3]);
            SetDefault = TRUE;
            break;
         case COM_BOUNDS:
            gBounds.flag = TRUE;
            gBounds.xmin = sRealParam[0];
            gBounds.xmax = sRealParam[1];
            gBounds.ymin = sRealParam[2];
            gBounds.ymax = sRealParam[3];
            gBounds.zmin = sRealParam[4];
            gBounds.zmax = sRealParam[5];
            break;
         case COM_RADIUS:
            DoRadius(sStrParam[0], sStrParam[1]);
            break;
         case COM_HIGHLIGHT:
            nhighlight++;
            DoZone(spheres, pdb, NSphere, sStrParam[0],sStrParam[1], 
                   sStrParam[2], sStrParam[3], sStrParam[4], 
                   nhighlight);
            break;
         case COM_BORDERWIDTH:
            gBorderWidth = (int)sRealParam[0];
            gBorderWidth--;
            if(gBorderWidth<0)
            {
               gBorderWidth = 0;
            }
            break;
         default:
            break;
         }
      }
      
      if(SetDefault)
         DoDefault(spheres,NSphere,DefaultRGB,pdb,ColourByTemp);
      if(SetCentre)
         DoCentre(pdb,CentreRes,CentreAtom);
      if(gSlab.flag)
         DoSlab(pdb,SlabRes,SlabAtom);
   }
   else
   {
      fprintf(stderr,"Unable to set up parser\n");
   }
}


/************************************************************************/
/*>void DoDefault(SPHERE *spheres, int NSphere, REAL RGB[3], 
                  PDB *pdb, BOOL ColourByTemp)
   ---------------------------------------------------------
   Set the default colour for atoms which have not otherwise been 
   coloured.
   21.07.93 Original    By: ACRM
   24.06.94 Modified to handle colour by temperature factor
            pdb and ColourByTemp parameters added
*/
void DoDefault(SPHERE *spheres, int NSphere, REAL RGB[3],
               PDB *pdb, BOOL ColourByTemp)
{
   int  i;
   PDB  *p;
   REAL colour, 
        MinBVal, 
        MaxBVal;
   
   if(ColourByTemp)
   {
      /* Find min and max BVal                                          */
      MinBVal = MaxBVal = pdb->bval;
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(p->bval > MaxBVal) MaxBVal = p->bval;
         if(p->bval < MinBVal) MinBVal = p->bval;
      }
         
      for(p=pdb, i=0; p!=NULL && i<NSphere; NEXT(p), i++)
      {
         if(!spheres[i].set)
         {
            /* Convert the BVal to a value between 0.0 and 1.0          */
            colour = (p->bval - MinBVal)/(MaxBVal - MinBVal);

            /* Since HSL has 0.0 == red, we convert this to a sensible
               hue to give colours from blue...green...red
            */
            colour = (REAL)2.0*((REAL)1.0 - colour)/(REAL)3.0;

            /* Calculate equivalent RGB colours                         */
            HSL2RGB(colour,1.0,1.0,&(RGB[0]),&(RGB[1]),&(RGB[2]));
            
            /* Set the colour of the sphere                             */
            spheres[i].r = RGB[0];
            spheres[i].g = RGB[1];
            spheres[i].b = RGB[2];
         }
      }
   }
   else
   {
      for(i=0; i<NSphere; i++)
      {
         if(!spheres[i].set)
         {
            spheres[i].r = RGB[0];
            spheres[i].g = RGB[1];
            spheres[i].b = RGB[2];
         }
      }
   }
}


/************************************************************************/
/*>void DoPhong(SPHERE *spheres, int NSphere, REAL shine, REAL metallic)
   ---------------------------------------------------------------------
   Set Phong shading parameters for all spheres
   21.07.93 Original    By: ACRM
*/
void DoPhong(SPHERE *spheres, int NSphere, REAL shine, REAL metallic)
{
   int i;
   
   for(i=0; i<NSphere; i++)
   {
      spheres[i].shine    = shine;
      spheres[i].metallic = metallic;
   }
}


/************************************************************************/
/*>void DoZone(SPHERE *spheres, PDB *pdb, int NSphere, char *start, 
               char *end, char *red, char *green, char *blue, 
               int highlight)
   ----------------------------------------------------------------
   Set colours of spheres based on zone information.
   
   21.07.93 Original    By: ACRM
   29.03.95 Modified to work correctly with b&s images where the 
            appropriate zone occurs in two parts
   18.06.96 Changed InZone() to blInPDBZone()
   18.10.07 Added type parameter
   27.01.15 Updated for new bioplib
*/
void DoZone(SPHERE *spheres, PDB *pdb, int NSphere, char *start, 
            char *end, char *red, char *green, char *blue, 
            int highlight)
{
   char     chain1[8],  chain2[8],
            insert1[8], insert2[8];
   int      resnum1,    resnum2,
            i;
   PDB      *p;
   double   r, g, b;
   
   sscanf(red,   "%lf",&r);
   sscanf(green, "%lf",&g);
   sscanf(blue,  "%lf",&b);
   
   blParseResSpec(start, chain1, &resnum1, insert1);
   blParseResSpec(end,   chain2, &resnum2, insert2);
   
   if(!CHAINMATCH(chain1,chain2)) return;
   
   for(p=pdb, i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      if(blInPDBZone(p, chain1, resnum1, insert1, resnum2, insert2))
      {
         if(highlight)
         {
            spheres[i].hr   = r;
            spheres[i].hg   = g;
            spheres[i].hb   = b;
            spheres[i].highlight = highlight;
         }
         else
         {
            spheres[i].r   = r;
            spheres[i].g   = g;
            spheres[i].b   = b;
            spheres[i].set = TRUE;
         }
      }
   }
}


/************************************************************************/
/*>void DoResidue(SPHERE *spheres, PDB *pdb, int NSphere, char *resnam, 
                  char *red, char *green, char *blue)
   --------------------------------------------------------------------
   Set colours of spheres based on residue information. 
   
   21.07.93 Original    By: ACRM
   22.07.93 Padded residue name to 4 chars
   06.12.95 Corrected comments
*/
void DoResidue(SPHERE *spheres, PDB *pdb, int NSphere, char *resnam, 
               char *red, char *green, char *blue)
{
   int      i;
   PDB      *p;
   double   r, g, b;
   
   sscanf(red,   "%lf",&r);
   sscanf(green, "%lf",&g);
   sscanf(blue,  "%lf",&b);

   UPPER(resnam);

   /* Pad residue name to 4 chars. N.B. We assume resnam is large enough
      to handle this!
   */
   blPadterm(resnam,4);
   
   for(p=pdb, i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      if(!strncmp(p->resnam,resnam,4))
      {
         spheres[i].r   = r;
         spheres[i].g   = g;
         spheres[i].b   = b;
         spheres[i].set = TRUE;
      }
   }
}


/************************************************************************/
/*>void DoRotate(PDB *pdb, char *direction, char *amount)
   ------------------------------------------------------
   Create a rotation matrix and apply to the pdb linked list.
   22.07.93 Original    By: ACRM
*/
void DoRotate(PDB *pdb, char *direction, char *amount)
{
   REAL  matrix[3][3],
         angle;
   
   angle  = (REAL)atof(amount);
   angle *= PI/180.0;
   
   blCreateRotMat(*direction, angle, matrix);
   
   blRotatePDB(pdb, matrix);
}


/************************************************************************/
/*>void DoCentre(PDB *pdb, char *resspec, char *atom)
   --------------------------------------------------
   Parse a residue spec and centre the display on the specified atom in
   this residue.
   23.07.93 Original    By: ACRM
   28.03.94 Added warning if atom not found
   23.10.95 Warnings go to stderr
   27.01.15 Updated for new bioplib
*/
void DoCentre(PDB *pdb, char *resspec, char *atom)
{
   char  chain[8],
         insert[8];
   int   resnum;
   PDB   *p;
   BOOL  found = FALSE;
   
   /* Parse the residue specification                                   */
   blParseResSpec(resspec, chain, &resnum, insert);
   
   /* Tidy up atom specification                                        */
   UPPER(atom);
   blPadterm(atom,4);

   /* Walk the pdb linked list                                          */
   for(p=pdb;p!=NULL;NEXT(p))
   {
      /* Check we are in the correct residue                            */
      if(p->resnum    == resnum &&
         CHAINMATCH(p->insert, insert) &&
         CHAINMATCH(p->chain, chain))
      {
         /* Residue found; set flag                                     */
         found = TRUE;

         /* If we've got the correct atom, set midpoint and return      */
         if(!strcmp(p->atnam,atom))
         {
            gMidPoint.x = p->x;
            gMidPoint.y = p->y;
            gMidPoint.z = p->z;
            
            return;
         }
         
         /* If it's the CA, record coords; this will be over-ridden by
            the requested atom if found, since we don't return from this
            part
         */
         if(!strcmp(p->atnam,"CA  "))
         {
            gMidPoint.x = p->x;
            gMidPoint.y = p->y;
            gMidPoint.z = p->z;
         }
      }
   }

   /* If residue not found, issue warning                               */
   if(!found)
      fprintf(stderr,"Warning: Residue for centre of display not \
found.\n");
}


/************************************************************************/
/*>void DoBackground(REAL r1, REAL g1, REAL b1, REAL r2, REAL g2, REAL b2)
   -----------------------------------------------------------------------
   Colour in the background shading from r1,g1,b1 at the top to r2,g2,b2 
   at the bottom.
   29.07.93 Original    By: ACRM
*/
void DoBackground(REAL r1, REAL g1, REAL b1, REAL r2, REAL g2, REAL b2)
{
   int   x, y;
   REAL  r, g, b;
   
   for(y=0; y<gScreen[1]; y++)
   {
      r = r1 + y * (r2 - r1)/gScreen[1];
      g = g1 + y * (g2 - g1)/gScreen[1];
      b = b1 + y * (b2 - b1)/gScreen[1];
      
      for(x=0; x<gScreen[0]; x++)
         SetAbsPixel(x,y,r,g,b);
   }
}


/************************************************************************/
/*>void DoSlab(PDB *pdb, char *resspec, char *atom)
   ------------------------------------------------
   Parse a residue spec and set coords to centre the slab wedge
   24.03.94 Original (based on DoCentre())    By: ACRM
   28.03.94 Changed check on atom name to use strncmp()
            Added warning if not found.
   23.10.95 Warnings go to stderr
   27.01.15 Updated for new bioplib
*/
void DoSlab(PDB *pdb, char *resspec, char *atom)
{
   char  chain[8],
         insert[8];
   int   resnum;
   PDB   *p;
   BOOL  found = FALSE;
   
   /* Parse the residue specification                                   */
   blParseResSpec(resspec, chain, &resnum, insert);
   
   /* Tidy up atom specification                                        */
   UPPER(atom);
   blPadterm(atom,4);

   /* Walk the pdb linked list                                          */
   for(p=pdb;p!=NULL;NEXT(p))
   {
      /* Check we are in the correct residue                            */
      if(p->resnum    == resnum &&
         CHAINMATCH(p->insert, insert) &&
         CHAINMATCH(p->chain, chain))
      {
         /* Residue found; set flag                                     */
         found = TRUE;

         /* If we've got the correct atom, set midpoint and return      */
         if(!strncmp(p->atnam,atom,4))
         {
            gSlab.z = p->z;
            
            return;
         }
         
         /* If it's the CA, record coords; this will be over-ridden by
            the requested atom if found, since we don't return from this
            part
         */
         if(!strcmp(p->atnam,"CA  "))
         {
            gSlab.z = p->z;
         }
      }
   }

   /* If not found residue, print warning                               */
   if(!found)
      fprintf(stderr,"Warning: Residue for centre of slab not found.\n");
}


/************************************************************************/
/*>void HSL2RGB(REAL hue, REAL saturation, REAL luminance,
                REAL *red, REAL *green, REAL *blue)
   -------------------------------------------------------
   Input:   REAL  hue            HSL hue value (0.0--1.0)
            REAL  saturation     HSL saturation value (0.0--1.0)
            REAL  luminance      HSL luminance value (0.0--1.0)
   Output:  REAL  *red           RGB red value (0.0--1.0)
            REAL  *green         RGB green value (0.0--1.0)
            REAL  *blue          RGB blue value (0.0--1.0)

   Converts an HSL colour value to an RGB colour value
   Loosely based on code by R.J. Mical from Book 1 of the Amiga 
   Programmers' Suite.

   Although the code has been completely re-written here is his
   copyright notice:

   Copyright (C) 1986, 1987, Robert J. Mical All Rights Reserved.  

   Any or all of this code can be used in any program as long as this
   entire copyright notice is retained, ok?  Thanks.

   The Amiga Programmer's Suite Book 1 is copyrighted but freely
   distributable.  All copyright notices and all file headers must be
   retained intact.

   The Amiga Programmer's Suite Book 1 may be compiled and assembled,
   and the resultant object code may be included in any software
   product.  However, no portion of the source listings or
   documentation of the Amiga Programmer's Suite Book 1 may be
   distributed or sold for profit or in a for-profit product without
   the written authorization of the author, RJ Mical.


   24.06.94 Original    By: ACRM
*/
void HSL2RGB(REAL hue, REAL saturation, REAL luminance,
             REAL *red, REAL *green, REAL *blue)
{
   REAL rising, falling, InvSat;
   int  sixth;
   
   /* Find which sixth of the hue spectrum we are in                    */
   sixth   = (int)((REAL)6.0 * hue);

   rising  = (hue - ((REAL)sixth / (REAL)6.0)) * (REAL)6.0;
   falling = (REAL)1.0 - rising;

   InvSat  = (REAL)1.0 - saturation;

   switch(sixth)
   {
   case 0:
   case 6:
      *red   = (REAL)1.0;
      *green = rising;
      *blue  = (REAL)0.0;
      break;
   case 1:
      *red   = falling;
      *green = (REAL)1.0;
      *blue  = (REAL)0.0;
      break;
   case 2:
      *red   = (REAL)0.0;
      *green = (REAL)1.0;
      *blue  = rising;
      break;
   case 3:
      *red   = (REAL)0.0;
      *green = falling;
      *blue  = (REAL)1.0;
      break;
   case 4:
      *red   = rising;
      *green = (REAL)0.0;
      *blue  = (REAL)1.0;
      break;
   case 5:
      *red   = (REAL)1.0;
      *green = (REAL)0.0;
      *blue  = falling;
      break;
   }

   *red   *= luminance;
   *green *= luminance;
   *blue  *= luminance;

   *red   += ((luminance-(*red))   * InvSat);
   *green += ((luminance-(*green)) * InvSat);
   *blue  += ((luminance-(*blue))  * InvSat);
}


/************************************************************************/
/*>void DoChain(SPHERE *spheres, PDB *pdb, int NSphere, char *chain,
                char *red, char *green, char *blue)
   ------------------------------------------------------------------
   Set colours of spheres based on chain name.
   
   06.12.95 Original    By: ACRM
   27.01.15 Updated for new bioplib
*/
void DoChain(SPHERE *spheres, PDB *pdb, int NSphere, char *chain, 
             char *red, char *green, char *blue)
{
   int      i;
   PDB      *p;
   double   r, g, b;
   
   sscanf(red,   "%lf",&r);
   sscanf(green, "%lf",&g);
   sscanf(blue,  "%lf",&b);

#ifdef DO_UPPER
   UPPER(chain);
#endif 

   for(p=pdb, i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      if(CHAINMATCH(p->chain, chain))
      {
         spheres[i].r   = r;
         spheres[i].g   = g;
         spheres[i].b   = b;
         spheres[i].set = TRUE;
      }
   }
}

/************************************************************************/
/*>void DoRadius(char *atomspec, char *radius_str)
   -----------------------------------------------
   Input:   char    *atomspec      Atom specification
            char    *radius_str    Radius (as a string)
   
*/
void DoRadius(char *atomspec, char *radius_str)
{
   static RADII *r   = NULL;
   char         *chp = NULL;
   REAL         radius;

   if(!sscanf(radius_str, "%lf", &radius))
   {
      fprintf(stderr,"Warning: RADIUS command has invalid \
radius: %s (ignored)\n", radius_str);
   }
   else
   {
      /* Make space in linked list to store atom spec and radius        */
      if(gRadii == NULL)
      {
         INIT(gRadii, RADII);
         r = gRadii;
      }
      else
      {
         ALLOCNEXT(r, RADII);
      }
      if(r==NULL)
      {
         fprintf(stderr,"Error: No memory for radius list\n");
         exit(1);
      }
      
      /* Store the radius                                               */
      r->radius = radius;

      /* Now parse out atom name and residue name                       */
      if((chp=strchr(atomspec, (int)'.'))==NULL)
      {
         /* No '.' in atomspec so it's just an atom name                */
         if(atomspec[0] == '*')
         {
            r->atnam[0] = '\0';
         }
         else
         {
            strncpy(r->atnam, atomspec, MAXATNAM);
            PADMINTERM(r->atnam, 4);
         }
         r->resnam[0] = '\0';
      }
      else
      {
         /* We found a '.', so terminate string there, copy last part to
            atomspec and grab residue name
         */
         *chp = '\0';
         if(*(chp+1) == '*')
         {
            r->atnam[0] = '\0';
         }
         else
         {
            strncpy(r->atnam,  chp+1, MAXATNAM);
            PADMINTERM(r->atnam, 4);
         }

         strncpy(r->resnam, atomspec, MAXATNAM);
         PADMINTERM(r->resnam, 4);
      }
   }
}

/************************************************************************/
/*>void DoAtom(SPHERE *spheres, PDB *pdb, int NSphere, char *atom, 
               char *red, char *green, char *blue)
   ---------------------------------------------------------------
   Set colours of spheres based on atom information. Handles * as wild
   card and ' is translated to a * for use in PDB files
   
   21.07.93 Original    By: ACRM
   22.07.93 Padded atom name to 4 chars
*/
void DoAtom(SPHERE *spheres, PDB *pdb, int NSphere, char *atom, 
            char *red, char *green, char *blue)
{
   int      i,
            NComp = 4;  /* Normally compare all 4 chars of atom name    */
   PDB      *p;
   char     *ptr;
   double   r, g, b;
   
   sscanf(red,   "%lf",&r);
   sscanf(green, "%lf",&g);
   sscanf(blue,  "%lf",&b);
   
   UPPER(atom);
   
   /* See if there is a wild card in the atom spec                      */
   if((ptr = strchr(atom,'*')) != NULL)
   {
      *ptr = '\0';
      NComp = strlen(atom);   /* Compare fewer characters               */
      
      if(NComp == 0) NComp = 4;
   }
   else
   {
      /* No *, pad atom to 4 chars. N.B. We assume atom is large enough
         to handle this!
      */
      blPadterm(atom,4);
   }
   
   /* Change a ' to a * for comparison                                  */
   if((ptr = strchr(atom,'\'')) != NULL)
   {
      *ptr = '*';
   }
   
   for(p=pdb, i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      if(!strncmp(p->atnam,atom,NComp))
      {
         spheres[i].r   = r;
         spheres[i].g   = g;
         spheres[i].b   = b;
         spheres[i].set = TRUE;
      }
   }
}


