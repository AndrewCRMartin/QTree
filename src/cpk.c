/*************************************************************************

   Program:    CPK
   File:       cpk.c
   
   Version:    V2.4
   Date:       27.01.15
   Function:   Size pre-processor for QTree to create CPK images
   
   Copyright:  (c) SciTech Software 1995-2015
   Author:     Dr. Andrew C. R. Martin
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
   CPK simply patches the atom radius into the occupancy column. This
   gives a file which should generate a normal CPK image when used with
   qtree -b (i.e. the results will be the same as qtree without the -b).
   This allows one to combine files from this program with B&S images
   generated with the ballstick preprocessor.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V2.0  28.03.95 Original (numbered to match other programs)
   V2.1  23.10.95 Cosmetic
   V2.2  14.10.03 Skipped
   V2.3  18.10.07 Skipped
   V2.4  27.01.15 Skipped

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "bioplib/macros.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines
*/
#define MAXBUFF 160

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Banner(void);
void PatchSizes(PDB *pdb);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *quiet);
void Usage(void);


/************************************************************************/
/* Variables global to this file only
*/
#ifdef _AMIGA
/* Version string                                                       */
static unsigned char 
   *sVers="\0$VER: CPK V2.4 - SciTech Software, 1995-2015";
#endif

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main routine for the CPK QTree preprocessor
   
   28.03.95 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   PDB  *pdb = NULL;
   FILE *in  = stdin,
        *out = stdout;
   BOOL quiet = FALSE;
   int  natom;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
            
   if(ParseCmdLine(argc, argv, infile, outfile, &quiet))
   {
      if(!quiet)
         Banner();
      
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in, &natom))!=NULL)
         {
            PatchSizes(pdb);
            WritePDB(out, pdb);
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>void Banner(void)
   ----------------
   Prints a banner message

   28.03.95 Original    By: ACRM
   23.10.95 V2.1
   14.10.03 V2.2
   18.10.07 V2.3
   27.01.15 V2.4
*/
void Banner(void)
{
   fprintf(stderr,"\nCPK V2.4\n");
   fprintf(stderr,"========\n");
   fprintf(stderr,"CPK preprocessor for QTree. SciTech Software\n");
   fprintf(stderr,"Copyright (C) 1995-2015 SciTech Software. All Rights \
Reserved.\n");
   fprintf(stderr,"This program is freely distributable providing no \
profit is made in so doing.\n\n");
}

      
/************************************************************************/
/*>void PatchSizes(PDB *pdb)
   ------------------------
   Places default sphere radii into occ column of a PDB linked list

   28.03.95 Original    By: ACRM
*/
void PatchSizes(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      switch(p->atnam[0])
      {
      case 'C':                           /* Carbon:      r = 1.7    */
         p->occ  = 1.7;
         break;
      case 'N':                           /* Nitrogen:    r = 1.7    */
         p->occ  = 1.7;
         break;
      case 'O':                           /* Oxygen:      r = 1.35   */
         p->occ  = 1.35;
         break;
      case 'S':                           /* Sulphur:     r = 1.7    */
         p->occ  = 1.7;
         break;
      case 'H':                           /* Hydrogen:    r = 1.0    */
         p->occ  = 1.0;
         break;
      default:                            /* Default:     r = 1.7    */
         p->occ  = 1.7;
         break;
      }
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *quiet)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *quiet       Should text be displayed
   Returns: BOOL                Success?

   Parse the command line
   
   28.03.95 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *quiet)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *quiet    = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'q':
            *quiet = TRUE;
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
/*>void Usage(void)
   ----------------
   Prints a usage message

   28.03.95 Original    By: ACRM
   23.10.95 V2.1; Changed `Work' to `Operate'
   14.10.03 V2.2
   18.10.07 V2.3
   27.01.15 V2.4
*/
void Usage(void)
{
   fprintf(stderr,"\nCPK V2.4 (c) 1995-2015, SciTech Software\n\n");

   fprintf(stderr,"Usage: cpk [-q] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -q Operate quietly\n\n");

   fprintf(stderr,"CPK takes a PDB file as input and writes a PDB file \
with radii in the\n");
   fprintf(stderr,"occupancy column. The output of this program may be \
merged with the\n");
   fprintf(stderr,"output of the BallStick preprocessor to obtain mixed \
CPK---ball and\n");
   fprintf(stderr,"stick images. I/O through standard input/output if \
files not specified\n\n");
}
