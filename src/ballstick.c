/*************************************************************************

   Program:    BallStick
   File:       BallStick.c
   
   Version:    V2.5
   Date:       28.08.19
   Function:   Preprocessor for QTree to create a Ball & Stick image
   
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
   V2.0  28.03.95 Supports stdio
   V2.1  23.10.95 Skipped
   V2.1a 18.09.97 Added command line option to set the maximum distance
                  between atoms to count as a bond
   V2.1b 30.09.97 Tidied for clean compile under gcc
   V2.1c 30.06.98 If -m is specified this now also applies to between-
                  residue links
   V2.2  14.10.03 Skipped
   V2.2  14.10.03 Changed for new PDB structure
   V2.3  18.10.07 Skipped
   V2.4  27.01.15 Updated for new bioplib
   V2.5  18.08.19 General cleanup and moved into GitHub

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
#include "bioplib/general.h"

/************************************************************************/
/* Variables global to this file only
*/
static PDB *sStickArray   = NULL;
static REAL gMaxBondLenSq = (REAL)4.0;
static BOOL gMaxSpecified = FALSE;

/************************************************************************/
/* Defines
*/
#define MAXBUFF 160

#ifdef _AMIGA
/* Version string                                                       */
static unsigned char *sVers="\0$VER: BallStick V2.5 - SciTech Software, \
1993-2007";
#endif

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
long int WriteSticks(FILE *out, PDB *pdb, int NDivide, REAL StickRad,
                     BOOL Disulphides);
void WriteBalls(FILE *out, PDB *sStickArray, int NDivide);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *NDivide, REAL *BallRad, REAL *StickRad, 
                  BOOL *Disulphides, BOOL *Quiet);
void UsageExit(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for creating Ball & Stick type space filling pictures
   
   28.07.93 Original    By: ACRM
   29.07.93 Added disulphide flag
   04.10.94 Writes to occ rather than bval
   28.03.95 Modified to support stdio. Now calls ParseCmdLine() and
            blOpenStdFiles(). Text goes to stderr
   23.10.95 V2.1
   30.09.97 V2.1b
   30.06.98 V2.1c
   14.10.03 V2.2
   18.10.07 V2.3
   27.01.15 V2.4
   18.08.19 V2.5
*/
int main(int argc, char **argv)
{
   FILE     *in      = stdin,
            *out     = stdout;
   PDB      *pdb,
            *p;
   int      natom,
            NDivide  = 30;
   long int TotalOut = 0;
   REAL     BallRad  = 0.4,
            StickRad = 0.2;
   BOOL     Disulphides = TRUE,
            Quiet       = FALSE;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];

   if(ParseCmdLine(argc, argv, infile, outfile, &NDivide, &BallRad,
                   &StickRad, &Disulphides, &Quiet))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         /* Banner message                                              */
         if(!Quiet)
         {
            fprintf(stderr,"\nBallStick V2.5\n");
            fprintf(stderr,"==============\n");
            fprintf(stderr,"Ball and Stick program for use with QTree. \
SciTech Software\n");
            fprintf(stderr,"Copyright (C) 1993-2019 SciTech Software. \
All Rights Reserved.\n");
            fprintf(stderr,"This program is freely distributable \
providing no profit is made in so doing.\n\n");
         }

         /* Read PDB file                                               */
         pdb = blReadPDB(in, &natom);
         TotalOut = natom;
         
         if(pdb != NULL)
         {
            /* Allocate memory for stick array                          */
            if((sStickArray = (PDB *)malloc(NDivide * sizeof(PDB))) == 
               NULL)
            {
               fprintf(stderr,"No memory for stick array\n");
               exit(0);
            }
            
            /* V1.11 Sets occup rather than B-val                       */
            /* Set the occup of all atoms to the ball radius            */
            for(p=pdb; p!=NULL; NEXT(p))
               p->occ = BallRad;
            
            /* Re-write the current atom information                    */
            blWritePDB(out,pdb);
            
            /* Write the stick information                              */
            TotalOut += WriteSticks(out,pdb,NDivide,StickRad,Disulphides);
            
            /* Print information                                        */
            if(!Quiet)
            {
               fprintf(stderr,"Input atoms  = %d\n",natom);
               fprintf(stderr,"Output atoms = %ld\n",TotalOut);
            }
         }
      }
   }
   else
   {
      UsageExit();
   }

   return(0);
}

/************************************************************************/
/*>long WriteSticks(FILE *out, PDB *pdb, int NDivide, REAL StickRad,
                    BOOL Disulphides)
   -----------------------------------------------------------------
   Write sets of small spheres for the sticks
   28.07.93 Original    By: ACRM
   29.07.93 Added disulphide support
   04.10.94 Stores in occ rather than Bval
   18.09.97 Max bond length in gMaxBondSq rather than hard coded
   27.01.15 Uses CHAINMATCH()
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
      end = blFindNextResidue(start);
      
      for(p=start; p!= end; NEXT(p))
      {
         for(q=p->next; q!=end; NEXT(q))
         {
            if(DISTSQ(p,q) < gMaxBondLenSq)
            {
               xstep = (q->x - p->x) / NDivide;
               ystep = (q->y - p->y) / NDivide;
               zstep = (q->z - p->z) / NDivide;
               
               for(i=0;i<NDivide;i++)
               {
                  /* Copy PDB information from parent atom              */
                  blCopyPDB(&(sStickArray[i]),(i<HalfNDiv?p:q));
                  
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
            if(!CHAINMATCH(p->chain, q->chain)) break;
            
            if(!strncmp(q->atnam,"N   ",4) || 
               !strncmp(q->atnam,"P   ",4))
            {
               if(!gMaxSpecified || (DISTSQ(p,q) < gMaxBondLenSq))
               {
                  xstep = (q->x - p->x) / NDivide;
                  ystep = (q->y - p->y) / NDivide;
                  zstep = (q->z - p->z) / NDivide;
                  
                  for(i=0;i<NDivide;i++)
                  {
                     /* Copy PDB information from parent atom           */
                     blCopyPDB(&(sStickArray[i]),(i<HalfNDiv?p:q));
                     
                     /* V1.11 occ rather than bval                      */
                     /* Set radius                                      */
                     sStickArray[i].occ = StickRad;
                     
                     /* Set coords                                      */
                     sStickArray[i].x = p->x + i * xstep;
                     sStickArray[i].y = p->y + i * ystep;
                     sStickArray[i].z = p->z + i * zstep;
                  }
                  
                  NBall += NDivide;
                  
                  WriteBalls(out, sStickArray, NDivide);
               }
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
                        blCopyPDB(&(sStickArray[i]),(i<HalfNDiv?p:q));
                        
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
   14.10.03 changed 'junk' to 'record_type' for new PDB structure
*/
void WriteBalls(FILE *fp, PDB *balls, int nballs)
{
   int i;
   
   for(i=0; i<nballs; i++)
   {
      fprintf(fp,"%-6s%5d  %-4s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f\
%6.2f%6.2f\n",
              balls[i].record_type,
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     int *NDivide, REAL *BallRad, REAL *StickRad, 
                     BOOL *Disulphides, BOOL *Quiet)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            int    *NDivide     Number of divisions
            REAL   *BallRad     Ball radius
            REAL   *StickRad    Stick radius
            BOOL   *Disulphides Do disulphides
            BOOL   *Quiet       Operate quietly
   Returns: BOOL                Success?

   Parse the command line
   
   28.03.95 Original    By: ACRM
   18.09.97 Added -m to specify max bond length. Added checks on sscanf()
            returns
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *NDivide, REAL *BallRad, REAL *StickRad, 
                  BOOL *Disulphides, BOOL *Quiet)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
         case 'N':
            argc--;  argv++;
            if(!sscanf(argv[0],"%d",NDivide))
               return(FALSE);
            if((*NDivide)%2) (*NDivide)++;  /* Round to multiple of 2   */
            break;
         case 'b':
         case 'B':
            argc--; argv++;
            if(!sscanf(argv[0],"%lf",BallRad))
               return(FALSE);
            break;
         case 's':
         case 'S':
            argc--; argv++;
            if(!sscanf(argv[0],"%lf",StickRad))
               return(FALSE);
            break;
         case 'd':
         case 'D':
            *Disulphides = FALSE;
            break;
         case 'q':
         case 'Q':
            *Quiet = TRUE;
            break;
         case 'm':
         case 'M':
            argc--; argv++;
            if(!sscanf(argv[0],"%lf",&gMaxBondLenSq))
               return(FALSE);
            gMaxBondLenSq *= gMaxBondLenSq;
            gMaxSpecified = TRUE;
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
/*>void UsageExit(void)
   --------------------
   Give usage message and exit
   29.07.93 Original extracted from main()   By: ACRM
   21.12.94 Added Copyright/version message
   28.03.95 V2.0 Uses stderr
   23.10.95 V2.1
   18.09.97 V2.1a Added -m
   30.09.97 V2.1b
   30.06.98 V2.1c -m now also applies to between residue links
   14.10.03 V2.2
   18.10.07 V2.3
   27.01.15 V2.4
   18.08.19 V2.5
*/
void UsageExit(void)
{
   fprintf(stderr,"\nBallStick V2.5 (c) 1993-2019 Prof. Andrew C.R. \
Martin, SciTech Software\n\n");
   
   fprintf(stderr,"Usage: BallStick [-q] [-n <n>] [-b <r>] [-s <r>] [-d] \
[-m <d>] [<in.pdb> [<out.pdb>]]\n");
   fprintf(stderr,"       -q Operate quietly\n");
   fprintf(stderr,"       -n Specify the number of spheres to be placed \
between atoms [30]\n");
   fprintf(stderr,"       -b Specify ball radius [0.4]\n");
   fprintf(stderr,"       -s Specify stick radius [0.2]\n");
   fprintf(stderr,"       -d Don't do disulphides\n");
   fprintf(stderr,"       -m Specify maximum distance between bonded \
atoms\n");
   fprintf(stderr,"          By default this is 2.0A within a residue \
and unlimited\n");
   fprintf(stderr,"          between residues. -m sets a limit for \
both\n");

   fprintf(stderr,"\nCreates a set of spheres between each atom pair.\n");

   exit(0);
}
