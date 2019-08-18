/*************************************************************************

   Program:    
   File:       PDBList.c
   
   Version:    V1.10R
   Date:       08.10.99
   Function:   PDB linked list manipulation
   
   Copyright:  (c) SciTech Software 1992-6
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0) 1372 275775
   EMail:      andrew@stagleys.demon.co.uk
               
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
   V1.0  22.02.94 Original release
   V1.1  23.05.94 Added FindNextChainPDB()
   V1.2  05.10.94 KillSidechain() uses BOOL rather than int
   V1.3  24.07.95 Added TermPDB()
   V1.4  25.07.95 Added GetPDBChainLabels()
   V1.5  26.09.95 Fixed bug in TermPDB()
   V1.6  12.10.95 Added DupePDB(), CopyPDBCoords()
   V1.7  23.10.95 Moved FindResidueSpec() to ParseRes.c
   V1.8  10.01.96 Added ExtractZonePDB()
   V1.9  14.03.96 Added FindAtomInRes()
   V1.10 08.10.99 Initialised some variables

*************************************************************************/
/* Includes
*/
#include <math.h>
#include <stdlib.h>

#include "MathType.h"
#include "SysDefs.h"
#include "pdb.h"
#include "macros.h"
#include "general.h"

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
/*>PDB *KillPDB(PDB *pdb, PDB *prev)
   ---------------------------------
   Input:   PDB  *pdb    Pointer to item in PDB linked list to be removed
            PDB  *prev   Pointer to previous item in linked list
   Returns: PDB  *       Next item in PDB linked list

   Kill an item in the PDB linked list and re-link correctly. Returns the
   next item in the list, so will be NULL when the last item in the list
   is killed.

   12.05.92 Original
   11.03.94 Now handles prev==NULL to delete first item in a list
*/
PDB *KillPDB(PDB *pdb,              /* Pointer to record to kill        */
             PDB *prev)             /* Pointer to previous record       */
{
   PDB *p;

/* Old action was just to return if prev==NULL
   if(prev == NULL) return(NULL);
*/
   p = pdb->next;

   if(prev!=NULL)
      prev->next = pdb->next;       /* Relink the list                  */
   free(pdb);                       /* Free the item                    */

   return(p);
}

/************************************************************************/
/*>void CopyPDB(PDB *out, PDB *in)
   -------------------------------
   Input:   PDB  *in     Input PDB record pointer
   Output:  PDB  *out    Output PDB record pointer

   Copy a PDB record, except that the ->next is set to NULL;

   12.05.92 Original    By: ACRM
*/
void CopyPDB(PDB *out,
             PDB *in)
{
   strcpy(out->junk,   in->junk);
   strcpy(out->atnam,  in->atnam);
   strcpy(out->resnam, in->resnam);
   strcpy(out->insert, in->insert);
   strcpy(out->chain,  in->chain);
   out->atnum  = in->atnum;
   out->resnum = in->resnum;
   out->x      = in->x;
   out->y      = in->y;
   out->z      = in->z;
   out->occ    = in->occ;
   out->bval   = in->bval;
   out->next   = NULL;
}

/************************************************************************/
/*>BOOL MovePDB(PDB *move, PDB **from, PDB **to)
   ---------------------------------------------
   Input:   PDB    *move     PDB record to be moved
   I/O:     PDB    **from    Start of PDB linked list containing record
            PDB    **to      Start of output linked list
   Returns: BOOL             Success?

   Moves a PDB record from one linked list to another. from and ret should
   point to the start of the 2 lists. If the ret list hasn't been started,
   ret should be NULL. Returns TRUE if moved, FALSE otherwise.

   13.05.92 Original
   19.06.92 Changed p=*to, etc. for crappy compilers
*/
BOOL MovePDB(PDB *move, PDB **from, PDB **to)
{
   PDB *p;
   BOOL ret = FALSE;
   
   if(move != NULL && *from != NULL)
   {
      /* Find the item before move in the *from list                    */
      if(move == (*from))           /* Start of list                    */
      {
         p = NULL;
      }
      else                          /* Middle of list                   */
      {
         /* Move p to item before move                                  */
         for(p = (*from); p->next && p->next != move; NEXT(p)) ;
      }
      
      /* Unlink move from *from                                         */
      if(p)          /* We're moving something in the middle of the list*/
      {
         /* Unlink move                                                 */
         p->next = move->next;
      }
      else           /* We're moving the first one in the list          */
      {
         /* If first in *from list, reset *from list                    */
         *from = move->next;
      }

      /* Add move onto the end of *to                                   */
      move->next = NULL;
      if(*to)
      {
         /* Move p to end of *to list                                   */
         for(p=(*to); p->next; NEXT(p)) ;
         /* Link in move                                                */
         p->next = move;
      }
      else
      {
         /* Initialise *to list                                         */
         *to = move;
      }
      ret = TRUE;
   }
   return(ret);
}

/************************************************************************/
/*>PDB *AppendPDB(PDB *first, PDB *second)
   ---------------------------------------
   Input:   PDB   *first    First linked list (may be NULL)
            PDB   *second   Second linked list
   Returns: PDB   *         Start of list 
   
   Appends list second onto first. Returns start of new list (useful if 
   first was NULL).
   13.05.92 Original
   09.07.93 Changed to use LAST()
*/
PDB *AppendPDB(PDB *first,
               PDB *second)
{
   PDB *p;

   if(first == NULL)
      return(second);
   
   p = first;
   LAST(p);

   p->next = second;
   return(first);
}

/************************************************************************/
/*>PDB *FindEndPDB(PDB *start)
   ---------------------------
   Input:   PDB   *start    PDB linked list
   Returns: PDB   *         pointer to next residue

   Step along a PDB linked list from start until we find a different
   residue. Return a pointer to this PDB item.
   
   08.07.93 Original    By: ACRM
   09.08.95 Now simply calls FindNextResidue() which is a rather more
            sensible name. Retained for backwards compatibility
*/
PDB *FindEndPDB(PDB *start)
{
   return(FindNextResidue(start));
}

/************************************************************************/
/*>BOOL KillSidechain(PDB *ResStart, PDB *NextRes, BOOL doCB)
   ----------------------------------------------------------
   Input:   PDB   *ResStart     Start of a residue in linked list
            PDB   *NextRes      Start of next residue
            BOOL  doCB          Flag to kill CB as part of s/c
   Returns: BOOL                Success?
   
   Kill a sidechain, by calls to KillPDB(). If doCB is set, will kill 
   the CB.
   N.B. At least 1 backbone atom must occur in the linked list before the
   sidechain.
   
   12.05.92 Original
   05.10.94 doCB is now a BOOL as is the return
*/
BOOL KillSidechain(PDB *ResStart,   /* Pointer to start of residue      */
                   PDB *NextRes,    /* Pointer to start if next residue */
                   BOOL doCB)       /* Flag to kill the CB              */
{
   PDB *p,
       *prev = NULL;
   
   for(p=ResStart; p && p!=NextRes; NEXT(p))
   {
      if(strcmp(p->atnam, "N   ") &&
         strcmp(p->atnam, "CA  ") &&
         strcmp(p->atnam, "C   ") &&
         strcmp(p->atnam, "O   ") &&
         strcmp(p->atnam, "CB  "))
      {
         if(prev == NULL) return(FALSE); /* No b/b atom before s/c      */

         /* KillPDB() returns the next in list, so exit if list ended   */
         if(KillPDB(p, prev) == NULL) break;
         p = prev;
      }
      
      /* Kill the CB if required                                        */
      if(doCB && !strcmp(p->atnam,"CB  "))
      {
         if(prev == NULL) return(FALSE);  /* No b/b atom before s/c     */

         /* KillPDB() returns the next in list, so exit if list ended   */
         if(KillPDB(p, prev) == NULL) break;
         p = prev;
      }
      
      prev = p;
   }
   return(TRUE);
}

/************************************************************************/
/*>PDB *GetPDBByN(PDB *pdb, int n)
   -------------------------------
   Input:   PDB   *pdb    PDB linked list
            int   n       Offset into linked list
   Returns: PDB   *       Pointer to n'th item in linked list

   Gets a pointer to a pdb item by taking a PDB linked list and an 
   integer.
   The pointer returned is the n'th item in the list
   
   13.05.92 Original
*/
PDB *GetPDBByN(PDB *pdb,
               int n)
{
   PDB *p;
   int i;
   
   for(i=0, p=pdb; p && i<n ; NEXT(p), i++) ;

   return(p);
}

/************************************************************************/
/*>PDB *FindNextChainPDB(PDB *pdb)
   -------------------------------
   I/O:     PDB   *pdb      PDB linked list
   Returns: PDB   *         Pointer to start of next chain in linked list

   Terminates the linked list at the end of the current chain and 
   returns a pointer to the start of the next chain.

   23.05.94 Original    By: ACRM
*/
PDB *FindNextChainPDB(PDB *pdb)
{
   PDB  *p, *ret = NULL;
   char chain;
   
   chain = pdb->chain[0];
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->next == NULL) || (p->next->chain[0] != chain))
      {
         ret = p->next;
         p->next = NULL;
         break;
      }
   }

   return(ret);
}


/************************************************************************/
/*>PDB *TermPDB(PDB *pdb, int length)
   ----------------------------------
   Input:   PDB   *pdb         PDB linked list
            int   length       Number of residues after which to terminate
   Returns: PDB   *            Pointer to next residue after terminated
                               list. NULL if not enough residues in linked
                               list.

   Terminate a PDB linked list after length residues, returning a pointer
   to the next residue. 

   Note that the number of residues may cross chain boundaries.

   06.07.95 Original    By: ACRM
   26.09.95 Corrected update of resnum etc to use p-> not pdb-> (!!)
*/
PDB *TermPDB(PDB *pdb, int length)
{
   int  resnum,
        count;
   char insert,
        chain;
   PDB  *p,
        *prev = NULL;
   
   resnum = pdb->resnum;
   insert = pdb->insert[0];
   chain  = pdb->chain[0];

   for(p=pdb, count=1; p!=NULL; NEXT(p))
   {
      if((p->resnum    != resnum) ||
         (p->chain[0]  != chain)  ||
         (p->insert[0] != insert))
      {
         if(++count > length)
         {
            prev->next = NULL;
            return(p);
         }

         resnum = p->resnum;
         insert = p->insert[0];
         chain  = p->chain[0];
      }
      prev = p;
   }

   return(NULL);
}


/************************************************************************/
/*>char *GetPDBChainLabels(PDB *pdb)
   ---------------------------------
   Input:   PDB    *pdb      PDB linked list
   Returns: char   *         Allocated string containing chain labels
                             NULL if unable to allocate memory

   Scans a PDB linked list for chain names. Allocates memory for a 
   string containing these labels which is returned.

   N.B. You must free the allocated memory when you've finished with it!

   25.07.95 Original    By: ACRM
*/
char *GetPDBChainLabels(PDB *pdb)
{
   char *chains;
   int  nchains   = 0,
        maxchains = 16;
   PDB  *p;
   
   /* Just return if linked list is NULL                                */
   if(pdb==NULL)
      return(NULL);

   /* Allocate a chunk for storing the chains                           */
   if((chains = (char *)malloc(maxchains * sizeof(char)))==NULL)
      return(NULL);

   /* Set up first chain label                                          */
   chains[nchains] = pdb->chain[0];

   /* Run through the linked list                                       */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If chain label has changed                                     */
      if(p->chain[0] != chains[nchains])
      {
         /* Increment chain count and reallocate memory if needed       */
         if(++nchains == maxchains)
         {
            maxchains += 16;
            if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
               return(NULL);
         }
         /* Store this new chain label                                  */
         chains[nchains] = p->chain[0];
      }
   }

   /* Increment chain count and reallocate memory if needed             */
   if(++nchains == maxchains)
   {
      maxchains += 16;
      if((chains = realloc(chains, maxchains * sizeof(char)))==NULL)
         return(NULL);
   }

   /* Terminate the chain list with a NUL character                     */
   chains[nchains] = '\0';

   return(chains);
}


/************************************************************************/
/*>PDB *FindNextResidue(PDB *pdb)
   ------------------------------
   Input:   PDB   *pdb      PDB linked list
   Returns: PDB   *         Next residue in PDB linked list or NULL if
                            there is none.

   Finds the next residue in a PDB linked list.

   08.08.95 Original    By: ACRM
*/
PDB *FindNextResidue(PDB *pdb)
{
   PDB  *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    != pdb->resnum) ||
         (p->insert[0] != pdb->insert[0]) ||
         (p->chain[0]  != pdb->chain[0]))
         return(p);
   }

   return(NULL);
}

/************************************************************************/
/*>PDB *DupePDB(PDB *in)
   ---------------------
   Input:   PDB    *in     Input PDB linked list
   Returns: PDB    *       Duplicated PDB linked list
                           (NULL on allocation failure)

   Duplicates a PDB linked list. Allocates new linked list with identical
   data.

   11.10.95 Original   By: ACRM
   08.10.99 Initialise q to NULL
*/
PDB *DupePDB(PDB *in)
{
   PDB *out = NULL,
       *p, *q = NULL;

   for(p=in; p!=NULL; NEXT(p))
   {
      if(out==NULL)
      {
         INIT(out, PDB);
         q=out;
      }
      else
      {
         ALLOCNEXT(q, PDB);
      }
      if(q==NULL)
      {
         FREELIST(out, PDB);
         return(NULL);
      }
      
      CopyPDB(q, p);
   }
   
   return(out);
}


/************************************************************************/
/*>BOOL CopyPDBCoords(PDB *out, PDB *in)
   -------------------------------------
   Input:   PDB  *in      Input PDB linked list
   Output:  PDB  *out     Output PDB linked list
   Returns: BOOL          Success?

   Applies the coordinates of `in' to `out'. Assumes that the structures
   are equivalent with identical atom ordering. Makes a simple check on
   resnam and atnam at each position.

   11.10.95 Original   By: ACRM
*/
BOOL CopyPDBCoords(PDB *out, PDB *in)
{
   PDB *p, *q;
   
   for(p=in, q=out; p!=NULL && q!=NULL; NEXT(p), NEXT(q))
   {
      if(strncmp(p->atnam,  q->atnam,  4) ||
         strncmp(p->resnam, q->resnam, 4))
         return(FALSE);
      
      q->x = p->x;
      q->y = p->y;
      q->z = p->z;
   }

   if(p!=NULL || q!=NULL)
      return(FALSE);
   
   return(TRUE);
}


/************************************************************************/
/*>PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, 
                       char *insert1, char *chain2, int resnum2, 
                       char *insert2)
   -----------------------------------------------------------------------
   Input:   PDB    *inpdb   Input PDB linked list
            char   *chain1  Start residue chain name
            int    resnum1  Start residue number
            char   *insert1 Start residue insert code
            char   *chain2  End residue chain name
            int    resnum2  End residue number
            char   *insert2 End residue insert code
   Returns: PDB    *        PDB linked list of the region of interest.

   Reduces a PDB linked list to those residues within a specified zone.
   Note that the PDB linked list is duplicated before extraction so
   pointers do not match those in the input PDB linked list. Excess
   records in the new PDB linked list are freed.

   10.01.96 Original   By: ACRM
*/
PDB *ExtractZonePDB(PDB *inpdb, char *chain1, int resnum1, char *insert1,
                    char *chain2, int resnum2, char *insert2)
{
   PDB *pdb, *p, 
       *start = NULL, 
       *last  = NULL, 
       *prev  = NULL;

   /* Duplicate the PDB linked list                                     */
   if((pdb = DupePDB(inpdb))==NULL)
      return(NULL);

   /* Find the first residue in the PDB linked list                     */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    == resnum1)   &&
         (p->chain[0]  == chain1[0]) &&
         (p->insert[0] == insert1[0]))
      {
         start = p;
         break;
      }
      prev = p;
   }

   if(start==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Find the last residue in the PDB linked list                      */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    == resnum2)   &&
         (p->chain[0]  == chain2[0]) &&
         (p->insert[0] == insert2[0]))
      {
         last = p;
         break;
      }
   }

   if(last==NULL)
   {
      FREELIST(pdb, PDB);
      return(NULL);
   }

   /* Step last onto the final atom in that residue                     */
   for(; last->next!=NULL; NEXT(last))
   {
      if((last->next->resnum    != resnum2)   ||
         (last->next->chain[0]  != chain2[0]) ||
         (last->next->insert[0] != insert2[0]))
      {
         break;
      }
   }

   /* Free linked list after 'last'                                     */
   if(last->next != NULL)
   {
      FREELIST(last->next, PDB);
      last->next = NULL;
   }
   
   /* Unlink 'start' from rest of linked list and free memory before 
      'start'
   */
   if(prev != NULL)
   {
      prev->next = NULL;
      FREELIST(pdb, PDB);
   }

   return(start);
}


/************************************************************************/
/*>PDB *FindAtomInRes(PDB *pdb, char *atnam_in)
   --------------------------------------------
   Input:   PDB    *pdb         The beginning of a residue in a PDB 
                                linked list
            char   *atnam_in    An atom name to search for (doesn't need
                                to be space-padded)
   Returns: PDB    *            Pointer to required atom, NULL if not
                                found

   14.03.96 Original   By: ACRM
*/
PDB *FindAtomInRes(PDB *pdb, char *atnam_in)
{
   PDB *end,
       *p;

   char atnam[8];
   
   /* First copy the specified atom name and pad to 4 chars             */
   strcpy(atnam,atnam_in);
   padterm(atnam,4);
   
   /* Find the end of this residue                                      */
   end = FindNextResidue(pdb);
   
   /* Search for the required atom                                      */
   for(p=pdb; p!=end; NEXT(p))
   {
      if(!strncmp(p->atnam,atnam,4))
         return(p);
   }
   
   return(NULL);
}
