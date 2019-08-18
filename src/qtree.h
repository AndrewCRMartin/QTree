/*************************************************************************

   Program:    QTree
   File:       qtree.h
   
   Version:    V2.2
   Date:       14.10.03
   Function:   Include file for QTree
   
   Copyright:  (c) SciTech Software 1993-2003
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
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
   SPEC, DEPTHCUE, OVERLAP_SLAB and SHOW_INFO may be defined for 
   conditional compilation of additional features.

**************************************************************************

   Revision History:
   =================
   V1.0  19.07.93 Original
   V1.1           Skipped
   V1.2           Skipped
   V1.3           Skipped
   V1.4           Skipped
   V1.5  14.09.93 Added gSphScale
   V1.6           Skipped
   V1.7  28.03.94 Added Slab
   V1.8  09.05.94 Skipped
   V1.9  13.05.94 Skipped
   V1.10 24.06.94 Handles TEMPERATURE (Changes in commands.c)
   V1.11 04.10.94 Skipped
   V1.12 21.12.94 Skipped
   V2.0  30.03.95 Skipped
   V2.1  23.10.95 Skipped
   V2.2  14.10.03 Added BOUNDS and RADIUS stuff

*************************************************************************/

#ifndef _QTREE_H
#define _QTREE_H

/************************************************************************/
/* Includes for types used here
*/
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Conditional compilation flags
*/
#define DEPTHCUE        /* Handle depth cueing                          */
#define SPEC            /* Handle specular reflections                  */
#define SHOW_INFO       /* Show program statistics                      */
#define OVERLAP_SLAB    /* Slabs will include any atom which overlaps the
                           slab region                                  */

/************************************************************************/
/* Defines
*/
#define SIZE  512       /* Display size (square) 256                    */
#define XSIZE 800       /* Screen width          320                    */
#define YSIZE 600       /* Screen height         256                    */
#define MAXATNAM 8      /* Max atom name array size                     */

/************************************************************************/
/* Structure type definitions
*/

typedef struct
{
   REAL  x, y, z,
         xmin, xmax,
         ymin, ymax,
         rad,
         r, g, b,
         shine,
         metallic;
   BOOL  set;
}  SPHERE;

typedef struct
{
   REAL  x, y, z,
         amb;
   BOOL  spec;
}  LIGHT;

typedef struct
{
   REAL  ZMin,
         ZRange,
         contrast;
}  DCUE;

typedef struct
{
   REAL z,
        depth;
   BOOL flag;
}  SLAB;

typedef struct
{
   REAL xmin, xmax,
        ymin, ymax,
        zmin, zmax;
   BOOL flag;
}  BOUNDS;

typedef struct _radii
{
   struct _radii *next;
   REAL radius;
   char atnam[MAXATNAM],
        resnam[MAXATNAM];
}  RADII;

/************************************************************************/
/* Global variables
*/
#ifdef MAIN    /*---------------- Main program definitions -------------*/
LIGHT  gLight;             /* The light                                 */
DCUE   gDepthCue;          /* Depthcue info                             */
REAL   gScale     = 0.9,   /* Scaling factor                            */
       gSphScale  = 1.0;   /* Sphere scaling factor                     */
VEC3F  gMidPoint;          /* Mid point of structure for centering      */
char   gOutFile[160];      /* Output file name                          */
int    gSize = SIZE,       /* Display size                              */
       gScreen[2];         /* Screen size                               */
SLAB   gSlab;              /* Slabbing                                  */
BOUNDS gBounds;            /* User specified boundary of image          */
RADII  *gRadii = NULL;     /* Linked list of atom radii                 */
#else          /*----------------------- Externals ---------------------*/
extern LIGHT  gLight;
extern DCUE   gDepthCue;
extern REAL   gScale,
              gSphScale;
extern VEC3F  gMidPoint;
extern char   gOutFile[160];
extern int    gSize,
              gScreen[2];
extern SLAB   gSlab;
extern BOUNDS gBounds;
extern RADII  *gRadii;
#endif         /*-------------------------------------------------------*/


#endif
