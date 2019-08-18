/*************************************************************************

   Program:    Worms
   File:       worms.c
   
   Version:    V1.4
   Date:       07.10.93
   Function:   Preprocessor for QTree to create a worms image
   
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
   Reads a PDB file extracting the C-alpha coordinates. These are then
   smoothed using a B-spline and coordinates are interpolated between the
   smoothed coordinates. A new PDB file is then written with these
   coordinates and atom type WRM. QTree will then read this file and
   generate a worms type image.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  23.07.93 Original
   V1.1  26.07.93 Added B-spline smoothing
   V1.2  28.07.93 Corrected bug in failure to open file
   V1.3  11.08.93 Frees other allocated memory
   V1.4  07.10.93 Fixed bug in FindChainPDB()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Variables global to this file only
*/
int   sTotalCAlpha = 0;

/************************************************************************/
#ifdef _AMIGA
/* Version string                                                       */
static unsigned char *sVers="\0$VER: Worms V1.4 - SciTech Software, 1993";
#endif

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
PDB *DivideSmoothPDB(PDB *pdb, int DivideSmoothIter, int *nsmooth);
PDB *DivideSmooth(PDB *InWorm, int NIn);
void WriteWorm(FILE *fp, PDB *worm, int nworm);
PDB *InterpPDB(PDB *spline, int nsmooth, int NDivide, int *nworm);
PDB *FindChainPDB(PDB *pdb);
PDB *BSplineSmoothPDB(PDB *pdb, int NStep, int *nspline);
void BSplineCoef(REAL x0,  REAL x1,  REAL x2,  REAL x3, 
                 REAL *a0, REAL *a1, REAL *a2, REAL *a3);
void BuildCA(REAL x1, REAL y1, REAL z1,
             REAL x2, REAL y2, REAL z2,
             REAL x3, REAL y3, REAL z3,
             REAL *outx, REAL *outy, REAL *outz);
PDB *BuildCAList(PDB *pdb);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for creating worms type space filling pictures
   
   23.07.93 Original    By: ACRM
   28.07.93 Corrected bug in unable to open file
   11.08.93 Free's memory allocated as divide and worm
*/
int main(int argc, char **argv)
{
   FILE  *in      = NULL,
         *out     = NULL;
   PDB   *pdb     = NULL,
         *divide  = NULL,
         *worm    = NULL;
   int   natom,
         nsmooth,
         nworm,
         DivideSmoothIter  = 4,
         NDivide           = 30,
         TotalOut          = 0;
   BOOL  DivSmooth         = FALSE;

   /* Parse command line arguments                                      */
   argc--;  argv++;
   
   while(argc > 2)
   {
      if(argv[0][0] != '-') break;
      
      switch(argv[0][1])
      {
      case 's':
      case 'S':
         argc--;  argv++;
         sscanf(argv[0],"%d",&DivideSmoothIter);
         break;
      case 'n':
      case 'N':
         argc--;  argv++;
         sscanf(argv[0],"%d",&NDivide);
         break;
      case 'd':
      case 'D':
         DivSmooth = TRUE;
      default:
         break;
      }
      
      argc--;  argv++;
   }
         
   /* Error in command line                                             */
   if(argc != 2)
   {
      printf("Usage: worms [-n <n>] [-d] [-s <n>] <in.pdb> <out.pdb>\n");
      printf("       -n Specify the number of spheres to be placed \
between splined atoms [30]\n");
      printf("       -d Use division smoothing rather than B-spline\n");
      printf("       -s Specify the division smoothing factor [4]\n");
      printf("\nCreates a set of spheres along the smoothed C-alpha \
chain.\n");
      exit(0);
   }

   /* Open input and output files                                       */
   if((in = fopen(*argv,"r")) == NULL)
   {
      printf("Unable to open input file %s\n",*argv);
      exit(0);
   }
   argv++;

   if((out = fopen(*argv,"w")) == NULL)
   {
      printf("Unable to open output file %s\n",*argv);
      exit(0);
   }
   
   /* Banner message                                                    */
   printf("\nWorms V1.4\n");
   printf("==========\n");
   printf("Worms program for use with QTree. SciTech Software\n");
   printf("Copyright (C) 1993 SciTech Software. All Rights Reserved.\n");
   printf("This program is freely distributable providing no profit is \
made in so doing.\n\n");

   /* Read PDB file                                                     */
   pdb = ReadPDB(in, &natom);
   
   if(pdb != NULL)
   {
      PDB   *start,
            *end,
            *p;
   
      /* Handle each chain in turn                                      */
      start = pdb;
      
      while(start != NULL)
      {
         end = FindChainPDB(start);
         for(p=start; p->next != end; NEXT(p)) ;
         p->next = NULL;

         if(DivSmooth)     /* Use division smoothing                    */
         {
            if((divide = DivideSmoothPDB(start,DivideSmoothIter,&nsmooth))
                         != NULL)
            {
               /* Interpolate positions between smoothed points         */
               if((worm = InterpPDB(divide,nsmooth,NDivide,&nworm)) 
                          != NULL)
               {
                  TotalOut += nworm;
                  WriteWorm(out, worm, nworm);
               }
            }
         }
         else              /* Use B-spline smoothing                    */
         {
            if((worm = BSplineSmoothPDB(start,NDivide,&nworm))
                       != NULL)
            {
               TotalOut += nworm;
               WriteWorm(out, worm, nworm);
            }
         }

         if(divide != NULL)   free(divide);
         if(worm   != NULL)   free(worm);
         divide = NULL;
         worm   = NULL;

         FREELIST(start, PDB);

         start = end;
      }
      
      /* Print information                                              */
      printf("Smoothing method    = %s\n",
             (DivSmooth?"Division":"B-spline"));
      if(DivSmooth)
         printf("DivideSmooth level  = %d\n",DivideSmoothIter);
      printf("Divisions           = %d\n",NDivide);
      printf("Input C-alpha atoms = %d\n",sTotalCAlpha);
      printf("Output worm atoms   = %d\n",TotalOut);
   }
}

/************************************************************************/
/*>PDB *DivideSmoothPDB(PDB *pdb, int DivideSmoothIter, int *nsmooth)
   ------------------------------------------------------------------
   Takes a PDB linked list, extracts the C-alphas and calls a division
   smothing routine DivideSmoothIter times over. Returns an array of PDB 
   structures containing the smoothed C-alpha points and outputs the 
   number of points.
   
   23.07.93 Original    By: ACRM
*/
PDB *DivideSmoothPDB(PDB *pdb, int DivideSmoothIter, int *nsmooth)
{
   int   NCalpha = 0,
         i;
   PDB   *p,
         *spline,
         *spline2;
   
   /* Count the C-alphas                                                */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ",4))
         NCalpha++;
   }
   
   /* Update total C-alphas                                             */
   sTotalCAlpha += NCalpha;
   
   /* Allocate sufficient memory                                        */
   if((spline = (PDB *)malloc(NCalpha * sizeof(PDB))) == NULL)
      return(NULL);
      
   /* Copy in the C-alpha data                                          */
   for(p=pdb, i=0; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ",4))
         CopyPDB(&(spline[i++]), p);
   }
   
   /* Do the division smoothing                                         */
   for(i=0; i<DivideSmoothIter; i++)
   {
      spline2 = DivideSmooth(spline, NCalpha++);
         if(spline2==NULL) return(spline);
      free(spline);
      spline  = spline2;
   }
   
   *nsmooth = NCalpha;

   return(spline);   
}

/************************************************************************/
/*>PDB *DivideSmooth(PDB *InWorm, int NIn)
   ---------------------------------------
   Takes an array of PDB structures and the number of items in the array 
   and creates a division smoothed version (1 division level only). 
   Returns a pointer to the output array which will be one longer than 
   the input array.
   
   23.07.93 Original    By: ACRM
*/
PDB *DivideSmooth(PDB *InWorm, int NIn)
{
   PDB   *OutWorm;
   int   i;
   
   if((OutWorm = (PDB *)malloc((NIn + 1) * sizeof(PDB))) == NULL)
      return(NULL);
      
   CopyPDB(&(OutWorm[0]),   &(InWorm[0]));
   CopyPDB(&(OutWorm[NIn]), &(InWorm[NIn-1]));
   
   for(i=1; i<NIn; i++)
   {
      CopyPDB(&(OutWorm[i]), &(InWorm[i]));
      
      OutWorm[i].x = InWorm[i-1].x + (InWorm[i].x - InWorm[i-1].x)/2.0;
      OutWorm[i].y = InWorm[i-1].y + (InWorm[i].y - InWorm[i-1].y)/2.0;
      OutWorm[i].z = InWorm[i-1].z + (InWorm[i].z - InWorm[i-1].z)/2.0;
   }
   
   return(OutWorm);
}

/************************************************************************/
/*>void WriteWorm(FILE *fp, PDB *worm, int nworm)
   ----------------------------------------------
   Writes the array of smoothed and interpolated atom positions setting 
   the atom name to WRM. Form is as a PDB file.
   
   23.07.93 Original    By: ACRM
*/
void WriteWorm(FILE *fp, PDB *worm, int nworm)
{
   int i;
   
   for(i=0; i<nworm; i++)
   {
      fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   %8.3lf%8.3lf%8.3lf\
%6.2lf%6.2lf\n",
              worm[i].junk,
              worm[i].atnum,
              "WRM ",
              worm[i].resnam,
              worm[i].chain,
              worm[i].resnum,
              worm[i].insert,
              worm[i].x,
              worm[i].y,
              worm[i].z,
              worm[i].occ,
              worm[i].bval);
   }
}

/************************************************************************/
/*>PDB *InterpPDB(PDB *spline, int nsmooth, int NDivide, int *nworm)
   -----------------------------------------------------------------
   Takes a divide smoothed set of C-alpha positions and divides the set 
   up introducing additional interpolated positions.
   Input:   spline   Array of smoothed C-alpha positions
            nsmooth  Number of positions
            NDivide  Number of divisions to make
   Output:  nworm    Resulting number of positions created
   Returns:          Array of interpolated positions
   
   23.07.93 Original    By: ACRM
*/
PDB *InterpPDB(PDB *spline, int nsmooth, int NDivide, int *nworm)
{
   PDB   *worm = NULL;
   int   i, j;
   
   *nworm = 1 + (nsmooth-1) * NDivide;
   
   if((worm = (PDB *)malloc(*nworm * sizeof(PDB))) == NULL)
      return(NULL);
      
   for(i=0; i<nsmooth-1; i++)
   {
      for(j=0;j<NDivide; j++)
      {
         int OutPos = i * NDivide + j;
         
         CopyPDB(&(worm[OutPos]), &(spline[i]));
         worm[OutPos].x = spline[i].x + 
                          j * (spline[i+1].x - spline[i].x)/NDivide;
         worm[OutPos].y = spline[i].y + 
                          j * (spline[i+1].y - spline[i].y)/NDivide;
         worm[OutPos].z = spline[i].z + 
                          j * (spline[i+1].z - spline[i].z)/NDivide;
      }
   }
   CopyPDB(&(worm[*nworm - 1]), &(spline[nsmooth - 1]));

   return(worm);   
}

/************************************************************************/
/*>PDB *FindChainPDB(PDB *pdb)
   ---------------------------
   Find the end of a chain using an inter C-alpha distance criterion. 
   Returns a pointer to the start of the next chain or NULL.
   23.07.93 Original    By: ACRM
   07.10.93 Now returns NULL if no new chain found (as it should
            have all along)
*/
PDB *FindChainPDB(PDB *pdb)
{
   PDB   *p,
         *q;
   COOR  CAlpha;
   
   CAlpha.x = -9999.0;
   CAlpha.y = -9999.0;
   CAlpha.z = -9999.0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If it's a C-alpha                                              */
      if(!strncmp(p->atnam,"CA  ",4))
      {
         /* If we've had a C-alpha already                              */
         if(CAlpha.x != -9999.0)
         {
            /* If the distance between C-alphas is >5A we've got a new
               chain
            */
            if(DISTSQ(&CAlpha, p) > 25.0)
            {
               /* Find start of this residue and return it              */
               for(q=pdb; q!=NULL; NEXT(q))
               {
                  if(p->resnum    == q->resnum  &&
                     p->chain[0]  == q->chain[0] &&
                     p->insert[0] == q->insert[0])
                     return(q);
               }
               return(p);
            }
         }
         
         /* Update CAlpha                                               */
         CAlpha.x = p->x;
         CAlpha.y = p->y;
         CAlpha.z = p->z;
      }
   }
   return((PDB *)NULL);
}

/************************************************************************/
/*>PDB *BSplineSmoothPDB(PDB *pdb, int NStep, int *nspline)
   --------------------------------------------------------
   Takes a PDB linked list, extracts the C-alphas and calls a B-spline
   smoothing routine. 
   Takes a PDB linked list and number of division steps as input; outputs 
   the number of spline-dmoothed points created and returns an array of 
   PDB structures containing the smoothed C-alpha points and outputs the
   number of points.
   
   23.07.93 Original    By: ACRM
*/
PDB *BSplineSmoothPDB(PDB *pdb, int NStep, int *nspline)
{
   int   NCalpha = 0,
         i, j, k;
   PDB   *p,
         *p0,
         *p1,
         *p2,
         *p3,
         *spline,
         *ca_pdb;
   REAL  a0,
         a1,
         a2,
         a3,
         t, tstep;
   
   /* Count the C-alphas                                                */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ",4))
         NCalpha++;
   }
   
   /* Update total C-alphas                                             */
   sTotalCAlpha += NCalpha;
   
   /* Make NStep a multiple of 2                                        */
   if(NStep%2) NStep++;
   
   /* Calculate total number of points and spline step (tstep)          */
   *nspline = NStep * (NCalpha - 1);
   tstep    = 1.0/NStep;
   
   /* Allocate sufficient memory for returned PDB array                 */
   if((spline = (PDB *)malloc(*nspline * sizeof(PDB))) == NULL)
      return(NULL);

   /* Build CA linked list with additional extrapolated CAs at start and
      finish
   */
   if((ca_pdb = BuildCAList(pdb))==NULL)
      return(NULL);

   /* Copy residue data into the output array                           */
   for(p=ca_pdb->next, i=0; p->next!=NULL; NEXT(p))
   {
      for(j=(-NStep)/2; j<NStep/2; j++)
      {
         k = i*NStep + j;
         
         if(k>=0 && k<*nspline)
            CopyPDB(&(spline[k]), p);
      }
      i++;
   }
   
   /* Do the B-spline smoothing                                         */
   for(p=ca_pdb, i=0; p->next->next->next!=NULL; NEXT(p))
   {
      /* Find pointers                                                  */
      p0 = p;
      p1 = p->next;
      p2 = p1->next;
      p3 = p2->next;
      
      /* Do B-spline for x                                              */
      BSplineCoef(p0->x, p1->x, p2->x, p3->x, &a0, &a1, &a2, &a3);
      for(t=0.0, j=i; t<1.0; t += tstep, j++)
         spline[j].x = ((a3 * t + a2)*t + a1)*t + a0;

      /* Do B-spline for y                                              */
      BSplineCoef(p0->y, p1->y, p2->y, p3->y, &a0, &a1, &a2, &a3);
      for(t=0.0, j=i; t<1.0; t += tstep, j++)
         spline[j].y = ((a3 * t + a2)*t + a1)*t + a0;

      /* Do B-spline for z                                              */
      BSplineCoef(p0->z, p1->z, p2->z, p3->z, &a0, &a1, &a2, &a3);
      for(t=0.0, j=i; t<1.0; t += tstep, j++)
         spline[j].z = ((a3 * t + a2)*t + a1)*t + a0;
         
      /* Reset i for next group                                         */
      i = j;
   }
   
   /* Free the Calpha linked list                                       */
   FREELIST(ca_pdb, PDB);
   
   *nspline = i;
   return(spline);   
}

/************************************************************************/
/*>void BSplineCoef(REAL x0,  REAL x1,  REAL x2,  REAL x3, 
                    REAL *a0, REAL *a1, REAL *a2, REAL *a3)
   --------------------------------------------------------
   Calculate B-spline coefficients.
   
   26.07.93 Original    By: ACRM
*/
void BSplineCoef(REAL x0,  REAL x1,  REAL x2,  REAL x3, 
                 REAL *a0, REAL *a1, REAL *a2, REAL *a3)
{
   *a0 = (x0 + 4.0*x1 + x2)/6.0;
   *a1 = (-x0 + x2)/2.0;
   *a2 = (x0 - 2.0*x1 + x2)/2.0;
   *a3 = (-x0 - 3.0*(x2 - x1) + x3)/6.0;
}



/************************************************************************/
/*>void BuildCA(REAL x1, REAL y1, REAL z1, 
                REAL x2, REAL y2, REAL z2,
                REAL x3, REAL y3, REAL z3, 
                REAL *outx, REAL *outy, REAL *outz)
   ------------------------------------------------
   Calculate the coordinates for a C-alpha given the coordinates for
   3 others.
   
   26.07.93 Original    By: ACRM
*/
void BuildCA(REAL x1, REAL y1, REAL z1,
             REAL x2, REAL y2, REAL z2,
             REAL x3, REAL y3, REAL z3,
             REAL *outx, REAL *outy, REAL *outz)
{

   REAL  alpha       = 123.0,    /* CA-CA-CA angle                      */
         beta        = 180.0,    /* Torsion angle                       */
         length      = 3.8,      /* CA-CA distance                      */
         sina,cosa,
         sinb,cosb,
         cosax,cosay,cosaz,
         x21,  y21,  z21,
         x23,  y23,  z23,  r23,
         x32,  y32,  z32,
         xh,   yh,   zh,
         xp,   yp,   zp,
         x21p, y21p, z21p, r21p,
         xs,   ys,   zs,
         xv,   yv,   zv,
         scalpr;
         
   cosa = (REAL)cos((double)alpha);
   sina = (REAL)sin((double)alpha);

   x21=x2-x1;
   y21=y2-y1;
   z21=z2-z1;

   x23=x2-x3;
   y23=y2-y3;
   z23=z2-z3;
   r23=(REAL)sqrt((double)(x23*x23 + y23*y23 + z23*z23));

   x32=x3-x2;
   y32=y3-y2;
   z32=z3-z2;
   
   xh=x32/r23;
   yh=y32/r23;
   zh=z32/r23;
   
   scalpr=(x21*x32+y21*y32+z21*z32)/r23;
   xp=scalpr*xh;
   yp=scalpr*yh;
   zp=scalpr*zh;
   
   x21p=x21-xp;
   y21p=y21-yp;
   z21p=z21-zp;
   r21p=(REAL)sqrt((double)(x21p*x21p + y21p*y21p + z21p*z21p));
   
   xv=x21p/r21p;
   yv=y21p/r21p;
   zv=z21p/r21p;
   
   xs=yh*zv-zh*yv;
   ys=zh*xv-xh*zv;
   zs=xh*yv-yh*xv;
   
   cosax=cosa*xh;
   cosay=cosa*yh;
   cosaz=cosa*zh;

   cosb=(REAL)cos((double)beta);
   sinb=(REAL)sin((double)beta);
   *outx=x3+length*(cosax+sina*(cosb*xv+sinb*xs));
   *outy=y3+length*(cosay+sina*(cosb*yv+sinb*ys));
   *outz=z3+length*(cosaz+sina*(cosb*zv+sinb*zs));
}

/************************************************************************/
/*>PDB *BuildCAList(PDB *pdb)
   --------------------------
   Builds a new linked list containing only the CAs. If possible, we 
   create extra records at the start and end and calculate coordinates 
   for these. This is because, to calculate the spline curve for atoms i 
   and i+1, we need to know about i-1 and i+2. The extrapolated 
   coordinates for i-1 and i+2 at the ends of the chain allow us to 
   include all the C-alphas (otherwise we would loose the termini).
   
   Takes a PDB linked list as input and returns another.

   26.07.93 Original    By: ACRM
*/
PDB *BuildCAList(PDB *pdb)
{
   PDB   *q,
         *p,
         *p1,
         *p2,
         *p3,
         *ca_pdb;
   int   i;
         
   q = ca_pdb = NULL;
   
   /* Create initial extra record                                       */
   INIT(ca_pdb, PDB);
   q = ca_pdb;
   if(q==NULL) return(NULL);
   CopyPDB(q,pdb);   /* We're only interested in the residue info       */

   /* Copy C-alphas into new linked list                                */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam, "CA  ",4))
      {
         /* Allocate memory                                             */
         ALLOCNEXT(q, PDB);
         
         /* Check allocation                                            */
         if(q==NULL)
         {
            if(ca_pdb != NULL) FREELIST(ca_pdb, PDB);
            return(NULL);
         }
         
         /* Copy PDB record                                             */
         CopyPDB(q, p);
      }
   }
      
   /* Create extra ending record                                        */
   ALLOCNEXT(q, PDB);
   if(q==NULL)
   {
      FREELIST(ca_pdb, PDB);
      return(NULL);
   }
   p=pdb;
   LAST(p);
   CopyPDB(q,p);     /* Again, only the residue info is relevant        */
   
   /* Find pointers to adjacent CAs for N terminus                      */
   p1 = p2 = p3 = NULL;
   for(p=ca_pdb,i=0; p!=NULL && i<4; NEXT(p),i++)
   {
      switch(i)
      {
      case 1:
         p3 = p;
         break;
      case 2:
         p2 = p;
         break;
      case 3:
         p1 = p;
         break;
      default:
         break;
      }
   }
   
   /* If any of these was NULL, we don't have enough CAs                */
   if(p1==NULL || p2==NULL || p3==NULL)
   {
      p = ca_pdb->next;
      free(ca_pdb);
      ca_pdb = p;
      
      for(p=ca_pdb; p->next != NULL; NEXT(p)) ;
      
      free(p->next);
      p->next = NULL;
   }
   else  /* Calculate the coords for the additional Nter position       */
   {
      BuildCA(p1->x, p1->y, p1->z, 
              p2->x, p2->y, p2->z, 
              p3->x, p3->y, p3->z, 
              &(ca_pdb->x), &(ca_pdb->y), &(ca_pdb->z));
   }
   

   /* Find pointers to adjacent CAs for C terminus                      */
   p1 = p2 = p3 = NULL;
   for(p=ca_pdb; p->next!=NULL; NEXT(p))
   {
      p1 = p2;
      p2 = p3;
      p3 = p;
   }
   
   /* If any of these was NULL, we don't have enough CAs                */
   if(p1==NULL || p2==NULL || p3==NULL)
   {
      p = ca_pdb->next;
      free(ca_pdb);
      ca_pdb = p;
      
      for(p=ca_pdb; p->next != NULL; NEXT(p)) ;
      
      free(p->next);
      p->next = NULL;
   }
   else  /* Calculate the coords for the additional Nter position       */
   {
      q = ca_pdb;
      LAST(q);
      
      BuildCA(p1->x, p1->y, p1->z, 
              p2->x, p2->y, p2->z, 
              p3->x, p3->y, p3->z, 
              &(q->x), &(q->y), &(q->z));
   }
   
   return(ca_pdb);
}
