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
   V3.0  19.08.19 Added PNG support

*************************************************************************/
/* Includes
*/
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines
*/
typedef struct  /* A coloured pixel                                      */
{
   uint8_t red;
   uint8_t green;
   uint8_t blue;
}  blPNGPIXEL;

typedef struct  /* An image                                              */
{
   blPNGPIXEL *pixels;
   size_t width;
   size_t height;
}  blPNGIMAGE;
    
/************************************************************************/
/* Prototypes
*/
blPNGPIXEL *blPNGPixelAt(blPNGIMAGE *bitmap, int x, int y);
BOOL blSavePNGToFile(blPNGIMAGE *bitmap, const char *path);
