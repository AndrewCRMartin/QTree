/*************************************************************************

   Program:    QTree
   File:       qtree.h
   
   Version:    V1.0
   Date:       19.07.93
   Function:   Include file for QTree
   
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
   SPEC, DEPTHCUE and SHOW_INFO may be defined for conditional compilation
   of additional features.

**************************************************************************

   Revision History:
   =================

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

/************************************************************************/
/* Defines
*/
#define SIZE  512       /* Display size (square) 256                    */
#define XSIZE 800       /* Screen width          320                    */
#define YSIZE 600       /* Screen height         256                    */

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

/************************************************************************/
/* Global variables
*/
#ifdef MAIN    /*---------------- Main program definitions -------------*/
LIGHT gLight;              /* The light                                 */
DCUE  gDepthCue;           /* Depthcue info                             */
REAL  gScale     = 0.9;    /* Scaling factor                            */
VEC3F gMidPoint;           /* Mid point of structure for centering      */
char  gOutFile[160];       /* Output file name                          */
int   gSize = SIZE,        /* Display size                              */
      gScreen[2];          /* Screen size                               */
#else          /*----------------------- Externals ---------------------*/
extern LIGHT gLight;
extern DCUE  gDepthCue;
extern REAL  gScale;
extern VEC3F gMidPoint;
extern char  gOutFile[160];
extern int   gSize,
             gScreen[2];
#endif         /*-------------------------------------------------------*/


#endif