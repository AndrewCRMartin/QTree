/*************************************************************************

   Program:    
   File:       ZonePDB.c
   
   Version:    V1.3R
   Date:       19.09.96
   Function:   Routines for handling zones in PDB linked lists
   
   Copyright:  (c) SciTech Software 1993-6
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
   V1.0  30.09.92 Original
   V1.1  16.06.93 Tidied for book. Mode now a char.
   V1.2  18.06.96 Added InPDBZone() from QTree program
   V1.3  19.09.96 Added InPDBZoneSpec()

*************************************************************************/

/* Includes
*/
#include "SysDefs.h"
#include "MathType.h"

#include "pdb.h"
#include "macros.h"

/************************************************************************/
/*>BOOL FindZonePDB(PDB *pdb, int start, char startinsert, int stop, 
                    char stopinsert, char chain, int mode, 
                    PDB **pdb_start, PDB **pdb_stop)
   -------------------------------------------------------------
   Input:   PDB   *pdb        PDB linked list
            int   start       Resnum of start of zone
            char  startinsert Insert code for start of zone
            int   stop        Resnum of end of zone
            char  stopinsert  Insert code for end of zone
            char  chain       Chain name
            int   mode        ZONE_MODE_RESNUM:     Use PDB residue 
                                                    numbers/chain
                              ZONE_MODE_SEQUENTIAL: Use sequential 
                                                    numbering
   Output:  PDB   **pdb_start Start of zone
            PDB   **pdb_stop  End of zone
   Returns: BOOL              OK?

   Finds pointers to the start and end of a zone in a PDB linked list. The
   end is the atom *after* the specified zone

   30.09.92 Original
   17.07.95 Chain name was being ignored in specs like L* (for whole
            of light chain)
   18.08.95 Now handles inserts
   31.07.95 Fixed bug when zone end==chain end
*/
BOOL FindZonePDB(PDB   *pdb,
                 int   start,
                 char  startinsert,
                 int   stop,
                 char  stopinsert,
                 char  chain,
                 int   mode,
                 PDB   **pdb_start,
                 PDB   **pdb_stop)
{
   PDB   *p;
   int   rescount,
         resnum,
         InStop = FALSE;
   char  insert;
   
   /* To start, we don't know where either are                          */
   *pdb_start = NULL;
   *pdb_stop  = NULL;
   
   /* If both start and stop are -1, then the whole structure (or a whole
      chain) is being specified
   */
   if((start == (-1)) && (stop == (-1)))
   {
      if(chain == ' ')                    /* Whole structure            */
      {
         *pdb_start = pdb;
         *pdb_stop  = NULL;
         return(TRUE);
      }
      else                                /* An individual chain        */
      {
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(p->chain[0] == chain)
            {
               if(*pdb_start == NULL)
               {
                  *pdb_start = p;
               }
            }
            else if(*pdb_start != NULL)   /* We've aleady got the start */
            {
               *pdb_stop = p;
               return(TRUE);
            }
         }
         if(*pdb_start==NULL) 
            return(FALSE);                /* Chain not found            */
         else
            return(TRUE);
      }
   }
   
   /* Handle one end of a zone being set to -1                          */
   if(start == -1) *pdb_start = pdb;
   if(stop  == -1) *pdb_stop  = NULL;
   
   /* If either end is still undefined                                  */
   if(*pdb_start == NULL || *pdb_stop == NULL)
   {
      /* Search reference structure for start and end of zone           */
      rescount = 1;
      resnum   = pdb->resnum;
      insert   = pdb->insert[0];
      InStop   = FALSE;
   
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(mode == ZONE_MODE_RESNUM)
         {
            if(chain == ' ' || chain == p->chain[0])
            {  /* We are in the correct chain                           */

               /* If start undefined, see if residue matches            */
               if(*pdb_start == NULL)
               {
                  if((p->resnum == start) && 
                     (p->insert[0] == startinsert))
                     *pdb_start = p;
               }
               
               /* If stop undefined, then find the following residue    */
               if(*pdb_stop == NULL)
               {
                  /* See if we have just moved out of the stop residue.
                     If so, set the stop position and return
                  */
                  if(InStop && 
                     (p->resnum != stop || p->insert[0] != stopinsert))
                  {
                     *pdb_stop = p;
                     return((*pdb_start==NULL)?FALSE:TRUE);
                  }

                  /* Residue matches, so set flag to say we're in the
                     last residue of the zone.
                  */
                  if((p->resnum == stop) &&
                     (p->insert[0] == stopinsert))
                     InStop = TRUE;
               }
               if(*pdb_start != NULL && *pdb_stop != NULL)  /*Found both*/
                  break;
            }
            else if(InStop)
            {
               /* We will get here if InStop has been set without having 
                  found the start of the next residue. This will occur
                  if the last residue of a zone was also the last 
                  residue of a chain, since the chain name will now have
                  changed.
                  We just set *pdb_stop to this pointer and return.
               */
               *pdb_stop = p;
               return((*pdb_start==NULL)?FALSE:TRUE);
            }
         }  /* End of ZONE_MODE_RESNUM                                  */
         else if(mode == ZONE_MODE_SEQUENTIAL)
         {
            /* Correct the residue count                                */
            if(p->resnum != resnum || p->insert[0] != insert)
            {
               rescount++;
               resnum = p->resnum;
               insert = p->insert[0];
            }
            
            if(*pdb_start == NULL)           /* Identify zone start     */
               if(rescount == start) *pdb_start = p;
            if(*pdb_stop == NULL)            /* Identify zone stop      */
            {
               if(InStop && rescount != stop) *pdb_stop = p;
               
               if(rescount == stop) InStop = TRUE;
            }
            if(*pdb_start != NULL && *pdb_stop != NULL)   /* Found both */
               break;
         }
         else
         {
            printf("   Error==> CreateFitAtoms(): Internal confusion!\n");
         }
      }  /* End of loop through PDB linked list                         */
   }  /* End of if() one pointer undefined                              */

   /* Check start of range has been found and return                    */
   return((*pdb_start==NULL)?FALSE:TRUE);
}


/************************************************************************/
/*>BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1, 
                  int resnum2, char insert2)
   ----------------------------------------------------------
   Input:   PDB    *p         Pointer to a PDB record
            char   chain      Chain name
            int    resnum1    First residue
            char   insert1    First insert code
            int    resnum2    Second residue
            char   insert2    Second insert code
   Returns: BOOL              Is p in the range specified?

   Checks that atom stored in PDB pointer p is within the specified 
   residue range.

   N.B. This assumes ASCII coding.

   29.03.95 Original    By: ACRM
   08.02.96 Insert residues inside a zone were not handled correctly!
   18.06.96 Added to bioplib from QTree (was called InZone())
*/
BOOL InPDBZone(PDB *p, char chain, int resnum1, char insert1, 
               int resnum2, char insert2)
{
   if(p->chain[0] == chain)
   {
      
      /* If residue number is *within* the range, return TRUE           */
      if((p->resnum > resnum1) && (p->resnum < resnum2))
         return(TRUE);
      
      /* If the range has a single residue number, check both inserts   */
      if((p->resnum == resnum1) && (p->resnum == resnum2))
      {
         if(((int)p->insert[0] >= (int)insert1) &&
            ((int)p->insert[0] <= (int)insert2))
            return(TRUE);
      }
      
      /* If residue number matches ends of range check insert           */
      if(((p->resnum == resnum1) &&
          ((int)p->insert[0] >= (int)insert1)) ||
         ((p->resnum == resnum2) &&
          ((int)p->insert[0] <= (int)insert2)))
         return(TRUE);
   }
   
   return(FALSE);
}


/************************************************************************/
/*>BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2)
   ----------------------------------------------------------
   Input:   PDB    *p         Pointer to a PDB record
            char   *resspec1  Res spec for first residue
            char   *resspec2  Res spec for last residue
   Returns: BOOL              Is p in the range specified?

   Determines whether a PDB pointer is within a residue range specified
   using standard format: [c]nnn[i]

   Also handles the residue spec of c* (i.e. chain name and a * to
   indicate all residues in a given chain). This must be given as
   resspec1 (resspec2 is then ignored).

   Calls InPDBZone() to do the actual work

   19.09.96 Original  By: ACRM
*/
BOOL InPDBZoneSpec(PDB *p, char *resspec1, char *resspec2)
{
   char chain1,  chain2,
        insert1, insert2;
   int  res1,    res2;

   /* Check for wildcard specification of whole chain                   */
   if(resspec1[1] == '*')
   {
      UPPER(resspec1);
      if(p->chain[0] == resspec1[0])
      {
         return(TRUE);
      }
      else
      {
         return(FALSE);
      }
   }
   
   ParseResSpec(resspec1, &chain1, &res1, &insert1);
   ParseResSpec(resspec2, &chain2, &res2, &insert2);

   if(chain1 != chain2)
      return(FALSE);
   
   return(InPDBZone(p, chain1, res1, insert1, res2, insert2));
}
