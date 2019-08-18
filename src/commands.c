/*************************************************************************

   Program:    QTree
   File:       commands.c
   
   Version:    V1.7
   Date:       24.03.94
   Function:   Handle command files for QTree program
   
   Copyright:  (c) SciTech Software 1993-4
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
   V1.3           Skipped
   V1.4           Skipped
   V1.5  14.09.93 Added sphere scaling option
   V1.6  04.01.94 Skipped
   V1.7  24.03.94 Added SLAB option and warning messages.
                  Removed ParseResSpec() as this is now in the library

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
#define PARSER_NCOMM          18

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
   MAKEKEY(sKeyWords[COM_ZONE],     "ZONE",        STRING,5);
   MAKEKEY(sKeyWords[COM_ATOM],     "ATOM",        STRING,4);
   MAKEKEY(sKeyWords[COM_RESIDUE],  "RESIDUE",     STRING,4);
   MAKEKEY(sKeyWords[COM_DEFAULT],  "DEFAULT",     NUMBER,3);
   MAKEKEY(sKeyWords[COM_AMBIENT],  "AMBIENT",     NUMBER,1);
   MAKEKEY(sKeyWords[COM_SPEC],     "SPECULAR",    NUMBER,0);
   MAKEKEY(sKeyWords[COM_CONTRAST], "CONTRAST",    NUMBER,1);
   MAKEKEY(sKeyWords[COM_PHONG],    "PHONG",       NUMBER,2);
   MAKEKEY(sKeyWords[COM_SCALE],    "SCALE",       NUMBER,1);
   MAKEKEY(sKeyWords[COM_ROTATE],   "ROTATE",      STRING,2);
   MAKEKEY(sKeyWords[COM_MATRIX],   "MATRIX",      NUMBER,9);
   MAKEKEY(sKeyWords[COM_CENTRE],   "CENTRE",      STRING,2);
   MAKEKEY(sKeyWords[COM_CENTER],   "CENTER",      STRING,2);
   MAKEKEY(sKeyWords[COM_XMATRIX],  "XMATRIX",     NUMBER,9);
   MAKEKEY(sKeyWords[COM_LIGHT],    "LIGHT",       NUMBER,3);
   MAKEKEY(sKeyWords[COM_BACKGRND], "BACKGROUND",  NUMBER,6);
   MAKEKEY(sKeyWords[COM_SPHSCALE], "SPHERESCALE", NUMBER,1);
   MAKEKEY(sKeyWords[COM_SLAB],     "SLAB",        STRING,3);
   
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
   03.08.93 Corrected call to ApplyMatrixPDB to RotatePDB
   14.09.93 Added spherescale
   24.03.94 Added SLAB
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
         i, j;
   REAL  DefaultRGB[3],
         matrix[3][3];
   BOOL  SetDefault   = FALSE,
         SetCentre    = FALSE;
   
   /* Default spheres to white                                          */
   for(i=0; i<3; i++)
      DefaultRGB[i] = 1.0;
   
   if(SetupParser())
   {
      if((fp=fopen(file,"r")) == NULL)
      {
         if(ReportError)
            printf("Unable to open control file: %s\n",file);

         return;
      }
      
      while(fgets(buffer,159,fp))
      {
         TERMINATE(buffer);
         
         key = parse(buffer,PARSER_NCOMM,sKeyWords,sRealParam,sStrParam);
         
         switch(key)
         {
         case PARSE_ERRC:
         case PARSE_ERRP:
            printf("Error in command file line:\n%s\n",buffer);
            break;
         case COM_ZONE:
            DoZone(spheres,pdb,NSphere,sStrParam[0],sStrParam[1],
                           sStrParam[2],sStrParam[3],sStrParam[4]);
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
            RotatePDB(pdb,matrix);
            break;
         case COM_XMATRIX:
            for(i=0; i<3; i++)
            {
               for(j=0; j<3; j++)
               {
                  matrix[j][i] = sRealParam[i*3 + j];
               }
            }
            RotatePDB(pdb,matrix);
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
         default:
            break;
         }
      }
      
      if(SetDefault)
         DoDefault(spheres,NSphere,DefaultRGB);
      if(SetCentre)
         DoCentre(pdb,CentreRes,CentreAtom);
      if(gSlab.flag)
         DoSlab(pdb,SlabRes,SlabAtom);
   }
   else
   {
      printf("Unable to set up parser\n");
   }
}

/************************************************************************/
/*>void DoDefault(SPHERE *spheres, int NSphere, REAL RGB[3])
   ---------------------------------------------------------
   Set the default colour for atoms which have not otherwise been 
   coloured.
   21.07.93 Original    By: ACRM
*/
void DoDefault(SPHERE *spheres, int NSphere, REAL RGB[3])
{
   int i;
   
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
               char *end, char *red, char *green, char *blue)
   ----------------------------------------------------------------
   Set colours of spheres based on zone information.
   
   21.07.93 Original    By: ACRM
*/
void DoZone(SPHERE *spheres, PDB *pdb, int NSphere, char *start, 
            char *end, char *red, char *green, char *blue)
{
   char     chain1,  chain2,
            insert1, insert2;
   int      resnum1, resnum2,
            i;
   PDB      *p;
   double   r, g, b;
   
   sscanf(red,   "%lf",&r);
   sscanf(green, "%lf",&g);
   sscanf(blue,  "%lf",&b);
   
   ParseResSpec(start, &chain1, &resnum1, &insert1);
   ParseResSpec(end,   &chain2, &resnum2, &insert2);
   
   if(chain1 != chain2) return;
   
   for(p=pdb, i=0; p!=NULL && i<NSphere; NEXT(p), i++)
   {
      if(p->chain[0]  == chain1   &&
         p->resnum    == resnum1  &&
         p->insert[0] == insert1)
      {
         /* Start found, step through to end                            */
         for(; p!=NULL && i<NSphere; NEXT(p), i++)
         {
            if(p->chain[0]  == chain2   &&
               p->resnum    == resnum2  &&
               p->insert[0] == insert2)
               break;
            spheres[i].r   = r;
            spheres[i].g   = g;
            spheres[i].b   = b;
            spheres[i].set = TRUE;
         }
         
         /* Check for abnormal ending                                   */
         if(p==NULL || i>=NSphere) break;
         
         /* Step through the end residue                                */
         for(; p!=NULL && i<NSphere; NEXT(p), i++)
         {
            if(p->chain[0]  != chain2   ||
               p->resnum    != resnum2  ||
               p->insert[0] != insert2)
               break;
            spheres[i].r   = r;
            spheres[i].g   = g;
            spheres[i].b   = b;
            spheres[i].set = TRUE;
         }
         
         /* This zone finished so break out                             */
         break;
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
      padterm(atom,4);
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

/************************************************************************/
/*>void DoResidue(SPHERE *spheres, PDB *pdb, int NSphere, char *type, 
                  char *red, char *green, char *blue)
   ------------------------------------------------------------------
   Set colours of spheres based on residue information. 
   
   21.07.93 Original    By: ACRM
   22.07.93 Padded residue name to 4 chars
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

   /* Pad residue name to 4 chars. N.B. We assume atom is large enough
      to handle this!
   */
   padterm(resnam,4);
   
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
   
   CreateRotMat(*direction, angle, matrix);
   
   RotatePDB(pdb, matrix);
}

/************************************************************************/
/*>void DoCentre(PDB *pdb, char *resspec, char *atom)
   --------------------------------------------------
   Parse a residue spec and centre the display on the specified atom in
   this residue.
   23.07.93 Original    By: ACRM
   28.03.94 Added warning if atom not found
*/
void DoCentre(PDB *pdb, char *resspec, char *atom)
{
   char  chain,
         insert;
   int   resnum;
   PDB   *p;
   BOOL  found = FALSE;
   
   /* Parse the residue specification                                   */
   ParseResSpec(resspec, &chain, &resnum, &insert);
   
   /* Tidy up atom specification                                        */
   UPPER(atom);
   padterm(atom,4);

   /* Walk the pdb linked list                                          */
   for(p=pdb;p!=NULL;NEXT(p))
   {
      /* Check we are in the correct residue                            */
      if(p->resnum    == resnum &&
         p->insert[0] == insert &&
         p->chain[0]  == chain)
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
      printf("Warning: Residue for centre of display not found.\n");
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
*/
void DoSlab(PDB *pdb, char *resspec, char *atom)
{
   char  chain,
         insert;
   int   resnum;
   PDB   *p;
   BOOL  found = FALSE;
   
   /* Parse the residue specification                                   */
   ParseResSpec(resspec, &chain, &resnum, &insert);
   
   /* Tidy up atom specification                                        */
   UPPER(atom);
   padterm(atom,4);

   /* Walk the pdb linked list                                          */
   for(p=pdb;p!=NULL;NEXT(p))
   {
      /* Check we are in the correct residue                            */
      if(p->resnum    == resnum &&
         p->insert[0] == insert &&
         p->chain[0]  == chain)
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
      printf("Warning: Residue for centre of slab not found.\n");
}

