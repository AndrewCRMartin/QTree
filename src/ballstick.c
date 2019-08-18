/*************************************************************************

   Program:    BallStick
   File:       BallStick.c
   
   Version:    V1.12
   Date:       21.12.94
   Function:   Preprocessor for QTree to create a Ball & Stick image
   
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
   Reads a PDB file and writes a new pdb file containing atoms at linearly
   interpolated positions along each bond to create ball and stick type 
   information. QTree will then read this file and generate a ball and 
   stick image.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  28.07.93 Original
   V1.1           Skipped
   V1.2  29.07.93 Added support for disulphide bonds
   V1.3           Skipped
   V1.4           Skipped
   V1.5           Skipped
   V1.6  04.01.94 Made externals static
   V1.7  29.03.94 Skipped
   V1.8  09.05.94 Skipped
   V1.9  13.05.94 Skipped
   V1.10 24.06.94 Skipped
   V1.11 04.10.94 Sphere size now in occ rather than BVal
   V1.12 21.12.94 Improved Usage message

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
static int sTotalCAlpha = 0;
static PDB *sStickArray = NULL;

/************************************************************************/
#ifdef _AMIGA
/* Version string                                                       */
static unsigned char *sVers="\0$VER: BallStick V1.12 - SciTech Software, \
1993-4";
#endif

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
long int WriteSticks(FILE *out, PDB *pdb, int NDivide, REAL StickRad,
                     BOOL Disulphides);
void WriteBalls(FILE *out, PDB *sStickArray, int NDivide);
void UsageExit(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for creating Ball & Stick type space filling pictures
   
   28.07.93 Original    By: ACRM
   29.07.93 Added disulphide flag
   04.10.94 Writes to occ rather than bval
*/
int main(int argc, char **argv)
{
   FILE     *in      = NULL,
            *out     = NULL;
   PDB      *pdb,
            *p;
   int      natom,
            NDivide  = 30;
   long int TotalOut = 0;
   REAL     BallRad  = 0.4,
            StickRad = 0.2;
   BOOL     Disulphides = TRUE;

   /* Parse command line arguments                                      */
   argc--;  argv++;
   
   while(argc > 2)
   {
      if(argv[0][0] != '-') break;
      
      switch(argv[0][1])
      {
      case 'n':
      case 'N':
         argc--;  argv++;
         sscanf(argv[0],"%d",&NDivide);
         if(NDivide%2) NDivide++;         /* Round to multiple of 2     */
         break;
      case 'b':
      case 'B':
         argc--; argv++;
         sscanf(argv[0],"%lf",&BallRad);
         break;
      case 's':
      case 'S':
         argc--; argv++;
         sscanf(argv[0],"%lf",&StickRad);
         break;
      case 'd':
      case 'D':
         Disulphides = FALSE;
         break;
      default:
         UsageExit();
         break;
      }
      
      argc--;  argv++;
   }
         
   /* Error in command line                                             */
   if(argc != 2)
   {
      UsageExit();
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
   printf("\nBallStick V1.12\n");
   printf("===============\n");
   printf("Ball and Stick program for use with QTree. \
SciTech Software\n");
   printf("Copyright (C) 1993-4 SciTech Software. All Rights \
Reserved.\n");
   printf("This program is freely distributable providing no profit is \
made in so doing.\n\n");

   /* Read PDB file                                                     */
   pdb = ReadPDB(in, &natom);
   TotalOut = natom;
   
   if(pdb != NULL)
   {
      /* Allocate memory for stick array                                */
      if((sStickArray = (PDB *)malloc(NDivide * sizeof(PDB))) == NULL)
      {
         printf("No memory for stick array\n");
         exit(0);
      }
      
      /* V1.11 Sets occup rather than B-val                             */
      /* Set the occup of all atoms to the ball radius                  */
      for(p=pdb; p!=NULL; NEXT(p))
         p->occ = BallRad;
         
      /* Re-write the current atom information                          */
      WritePDB(out,pdb);
      
      /* Write the stick information                                    */
      TotalOut += WriteSticks(out,pdb,NDivide,StickRad,Disulphides);
      
      /* Print information                                              */
      printf("Input atoms  = %d\n",natom);
      printf("Output atoms = %d\n",TotalOut);
   }
}

/************************************************************************/
/*>long WriteSticks(FILE *out, PDB *pdb, int NDivide, REAL StickRad,
                    BOOL Disulphides)
   -----------------------------------------------------------------
   Write sets of small spheres for the sticks
   28.07.93 Original    By: ACRM
   29.07.93 Added disulphide support
   04.10.94 Stores in occ rather than Bval
*/
long int WriteSticks(FILE *out, PDB *pdb, int NDivide, REAL StickRad,
                     BOOL Disulphides)
{
   PDB      *start,
            *end,
            *p,
            *q;
   int      i,
            HalfNDiv;
   long int NBall = 0;
   REAL     xstep,
            ystep,
            zstep;
         
   HalfNDiv = NDivide / 2;
         
   /* In residue links                                                  */
   for(start=pdb; start!=NULL; start=end)
   {
      end = FindEndPDB(start);
      
      for(p=start; p!= end; NEXT(p))
      {
         for(q=p->next; q!=end; NEXT(q))
         {
            if(DIST(p,q) < 2.0)
            {
               xstep = (q->x - p->x) / NDivide;
               ystep = (q->y - p->y) / NDivide;
               zstep = (q->z - p->z) / NDivide;
               
               for(i=0;i<NDivide;i++)
               {
                  /* Copy PDB information from parent atom              */
                  CopyPDB(&(sStickArray[i]),(i<HalfNDiv?p:q));
                  
                  /* V1.11 occ rather than bval                         */
                  /* Set radius                                         */
                  sStickArray[i].occ = StickRad;
                  
                  /* Set coords                                         */
                  sStickArray[i].x = p->x + i * xstep;
                  sStickArray[i].y = p->y + i * ystep;
                  sStickArray[i].z = p->z + i * zstep;
               }
               
               NBall += NDivide;
               
               WriteBalls(out, sStickArray, NDivide);
            }
         }
      }
   }
   
   /* Between residue links (C--N and O3*--P)                           */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"C   ",4) || !strncmp(p->atnam,"O3* ",4))
      {
         for(q=p->next; q!=NULL; NEXT(q))
         {
            if(p->chain[0] != q->chain[0]) break;
            
            if(!strncmp(q->atnam,"N   ",4) || 
               !strncmp(q->atnam,"P   ",4))
            {
               xstep = (q->x - p->x) / NDivide;
               ystep = (q->y - p->y) / NDivide;
               zstep = (q->z - p->z) / NDivide;
               
               for(i=0;i<NDivide;i++)
               {
                  /* Copy PDB information from parent atom              */
                  CopyPDB(&(sStickArray[i]),(i<HalfNDiv?p:q));
                  
                  /* V1.11 occ rather than bval                         */
                  /* Set radius                                         */
                  sStickArray[i].occ = StickRad;
                  
                  /* Set coords                                         */
                  sStickArray[i].x = p->x + i * xstep;
                  sStickArray[i].y = p->y + i * ystep;
                  sStickArray[i].z = p->z + i * zstep;
               }
                           
               NBall += NDivide;

               WriteBalls(out, sStickArray, NDivide);

               p = q;
               break;
            }
         }
      }
   }
   
   if(Disulphides)
   {
      /* Disulphide bonds                                               */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->resnam, "CYS ", 4) && 
            !strncmp(p->atnam,  "SG  ", 4))
         {
            for(q=p->next; q!=NULL; NEXT(q))
            {
               if(!strncmp(p->resnam, "CYS ", 4) && 
                  !strncmp(p->atnam,  "SG  ", 4))
               {
                  if(DIST(p,q) < 2.5)
                  {
                     xstep = (q->x - p->x) / NDivide;
                     ystep = (q->y - p->y) / NDivide;
                     zstep = (q->z - p->z) / NDivide;
                     
                     for(i=0;i<NDivide;i++)
                     {
                        /* Copy PDB information from parent atom        */
                        CopyPDB(&(sStickArray[i]),(i<HalfNDiv?p:q));
                        
                        /* V1.11 occ rather than bval                   */
                        /* Set radius                                   */
                        sStickArray[i].occ = StickRad;
                        
                        /* Set coords                                   */
                        sStickArray[i].x = p->x + i * xstep;
                        sStickArray[i].y = p->y + i * ystep;
                        sStickArray[i].z = p->z + i * zstep;
                     }
                                 
                     NBall += NDivide;
      
                     WriteBalls(out, sStickArray, NDivide);
                     break;
                  }
               }
            }
         }
      }
   }
   
   return(NBall);
}   

/************************************************************************/
/*>void WriteBalls(FILE *fp, PDB *balls, int nballs)
   -------------------------------------------------
   Write the balls which form a stick.
   28.07.93 Original    By: ACRM
*/
void WriteBalls(FILE *fp, PDB *balls, int nballs)
{
   int i;
   
   for(i=0; i<nballs; i++)
   {
      fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   %8.3lf%8.3lf%8.3lf\
%6.2lf%6.2lf\n",
              balls[i].junk,
              balls[i].atnum,
              balls[i].atnam,
              balls[i].resnam,
              balls[i].chain,
              balls[i].resnum,
              balls[i].insert,
              balls[i].x,
              balls[i].y,
              balls[i].z,
              balls[i].occ,
              balls[i].bval);
   }
}

/************************************************************************/
/*>void UsageExit(void)
   --------------------
   Give usage message and exit
   29.07.93 Original extracted from main()   By: ACRM
   21.12.94 Added Copyright/version message
*/
void UsageExit(void)
{
   printf("\nBallStick V1.12 (c) 1993-4 Dr. Andrew C.R. Martin, \
SciTech Software\n\n");
   
   printf("Usage: BallStick [-n <n>] [-b <r>] [-s <r>] [-d] <in.pdb> \
<out.pdb>\n");
   printf("       -n Specify the number of spheres to be placed \
between atoms [30]\n");
   printf("       -b Specify ball radius [0.4]\n");
   printf("       -s Specify stick radius [0.2]\n");
   printf("       -d Don't do disulphides\n");
   printf("\nCreates a set of spheres between each atom pair.\n");

   exit(0);
}
