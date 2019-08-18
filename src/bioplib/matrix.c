/*************************************************************************

   Program:    
   File:       matrix.c
   
   Version:    V1.6R
   Date:       27.09.95
   Function:   Simple matrix and vector operations
   
   Copyright:  (c) SciTech Software 1991-5
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  06.09.91 Original
   V1.0a 01.06.92 Documented
   V1.1  30.09.92 Matrix multiplication added
   V1.2  10.06.93 void return from matrix multiplication
   V1.3  22.07.93 Added CreateRotMat()
   V1.4  03.08.93 Changed matrix multiplication to standard direction
   V1.5  28.07.95 Added VecDist()
   V1.6  27.09.95 Added MatMult33_33()

*************************************************************************/
/* Includes
*/
#include <math.h>
#include "MathType.h"

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/

/************************************************************************/
/*>void MatMult3_33(VEC3F vecin, REAL matin[3][3], VEC3F *vecout)
   -------------------------------------------------------------
   Input:   VEC3F vecin        Vector to be multiplied
            REAL  matin[3][3]  Rotation matrix
   Output:  VEC3F *vecout      Output multiplied vector

   Multiply a 3-vector by a 3x3 matrix

   30.09.92 Original
   03.08.93 Changed multiplication to standard direction
*/
void MatMult3_33(VEC3F vecin, 
                 REAL  matin[3][3], 
                 VEC3F *vecout)
{
   vecout->x = vecin.x * matin[0][0] +
               vecin.y * matin[1][0] +
               vecin.z * matin[2][0];
   vecout->y = vecin.x * matin[0][1] +
               vecin.y * matin[1][1] +
               vecin.z * matin[2][1];
   vecout->z = vecin.x * matin[0][2] +
               vecin.y * matin[1][2] +
               vecin.z * matin[2][2];
}

/************************************************************************/
/*>void invert33(REAL s[3][3], REAL ss[3][3])
   ------------------------------------------
   Input:   REAL    s[3][3]  Input matrix
   Output:  REAL    ss[3][3] Ouput inverted matrix

   Invert a 3x3 matrix

   06.09.91 Original
   01.06.92 Documented
   10.06.93 void return
*/
void invert33(REAL s[3][3],
              REAL ss[3][3])
{
   int   i,  j,
         i1, j1,
         i2, j2;
   REAL  det;
       
   for(i=0;i<3;i++)
      for(j=0;j<3;j++)
         ss[i][j] = 0.0;
         
   det = 0.0;
   for(j=0;j<3;j++)
   {
      for(i=0;i<3;i++)
      {
         switch(i)
         {
         case 0:
            i1 = 2;
            i2 = 3;
            break;
         case 1:
            i1 = 1;
            i2 = 3;
            break;
         case 3:
            i1 = 1;
            i2 = 2;
            break;
         }
         switch(j)
         {
         case 0:
            j1 = 2;
            j2 = 3;
            break;
         case 1:
            j1 = 1;
            j2 = 3;
            break;
         case 3:
            j1 = 1;
            j2 = 2;
            break;
         }
         ss[i][j] = (REAL)pow((double)-1.0,(double)(i+j)) * 
                    (s[j1][i1] * s[j2][i2] - s[j2][i1] * s[j1][i2]);
      }
   }
   det = s[0][0]*ss[0][0] + s[0][1]*ss[1][0] + s[0][2]*ss[2][0];
   det = 1.0/det;
   
   for(i=0;i<3;i++)
      for(j=0;j<3;j++)
         ss[i][j] = det*ss[i][j];
}


/************************************************************************/
/*>void CreateRotMat(char direction, REAL angle, REAL matrix[3][3])
   ----------------------------------------------------------------
   Input:   char direction    Axis about which to rotate
            REAL angle        Angle (in rads) to rotate
   Output:  REAL matrix[3][3] Rotation matrix

   Create a 3x3 rotation matrix. Takes a direction as a single character
   ('x', 'y', or 'z'), an angle (in rads) and outputs a rotation matrix
   
   22.07.93 Original    By: ACRM
*/
void CreateRotMat(char direction, REAL angle, REAL matrix[3][3])
{
   int   i, j,
         m, 
         m1,
         m2;
   REAL  CosTheta,
         SinTheta;
         
   /* Initialise matrix to all 0.0                                      */
   for(i=0; i<3; i++)
      for(j=0; j<3; j++)
         matrix[i][j] = 0.0;
   
   /* Select the items that need to be filled in                        */
   switch(direction)
   {
   case 'x':   case 'X':
      m = 0;
      break;
   case 'y':   case 'Y':
      m = 1;
      break;
   case 'z':   case 'Z':
      m = 2;
      break;
   default:                   /* Just return the unit matrix            */
      for(i=0; i<3; i++)
         matrix[i][i] = 1.0;
      return;
   }
   
   /* Find which items these relate to                                  */
   m1 = (m+1)  % 3;
   m2 = (m1+1) % 3;
   
   /* Fill in the values                                                */
   matrix[m][m]   = 1.0;
   CosTheta       = (REAL)cos((double)angle);
   SinTheta       = (REAL)sin((double)angle);
   matrix[m1][m1] = CosTheta;
   matrix[m2][m2] = CosTheta;
   matrix[m1][m2] = SinTheta;
   matrix[m2][m1] = -SinTheta;
}


/************************************************************************/
/*>REAL VecDist(REAL *a, REAL *b, int len)
   ---------------------------------------
   Input:   REAL    *a     An arbitrary length vector (as an array)
            REAL    *b     An arbitrary length vector (as an array)
            int     len    The dimensionality of the vectors (array
                           length)
   Returns: REAL           The distance between the points described by
                           the two vectors

   Finds the distance between two vectors of arbitrary length

   28.07.95 Original    By: ACRM
*/
REAL VecDist(REAL *a, REAL *b, int len)
{
   REAL sumsq = 0.0;
   int  i;
   
   for(i=0; i<len; i++)
      sumsq += (a[i] - b[i]) * (a[i] - b[i]);

   return((REAL)sqrt((double)sumsq));
}


/************************************************************************/
/*>void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
   -------------------------------------------------------------
   Input:   REAL  a[3][3]      Matrix to be multiplied
            REAL  b[3][3]      Matrix to be multiplied
   Output:  REAL  out[3][3]    Output matrix

   Multiply two 3x3 matrices

   27.09.95 Original
*/
void MatMult33_33(REAL a[3][3], REAL b[3][3], REAL out[3][3])
{
   int  i, j, k;
   REAL ab;
   
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         ab = (REAL)0.0;
         for(k=0; k<3; k++)
         {
            ab += a[i][k]*b[k][j];
         }
         out[i][j]=ab;
      }
   }
}

         
   
      
   
   
