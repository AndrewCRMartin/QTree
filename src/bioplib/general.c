/*************************************************************************

   Program:    
   File:       general.c
   
   Version:    V1.20R
   Date:       18.09.96
   Function:   General purpose routines
   
   Copyright:  (c) SciTech Software 1991-6
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
   These are generally useful C routines, mostly string handling.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.1  08.02.91 Added KillLine()
   V1.2  10.02.91 Added setextn() and index()
   V1.3  20.03.91 Added Word()
   V1.4  28.05.92 ANSIed
   V1.5  22.06.92 Added tab check to Word(). Improved setextn().
                  Added WordN(). Documented other routines.
   V1.6  27.07.93 Corrected fsscanf() for double precision
   V1.7  07.10.93 Checks made on case before toupper()/tolower()
                  for SysV compatibility. Also index() becomes
                  chindex()
   V1.8  18.03.94 getc() -> fgetc()
   V1.9  11.05.94 Added GetFilestem(), upstrcmp(), upstrncmp() &
                  GetWord()
   V1.10 24.08.94 Added OpenStdFiles()
   V1.11 08.03.95 Corrected OpenFile() for non-UNIX
   V1.12 09.03.95 Added check on non-NULL filename in OpenFile()
   V1.13 17.07.95 Added countchar()
   V1.14 18.10.95 Moved YorN() to WindIO.c
   V1.15 06.11.95 Added StoreString(), InStringList() and FreeStringList()
   V1.16 22.11.95 Moved ftostr() to generam.c
   V1.17 15.12.95 Added QueryStrStr()
   V1.18 18.12.95 OpenStdFiles() treats filename of - as stdin/stdout
   V1.19 05.02.96 OpenStdFiles() allows NULL pointers instead if filenames
   V1.20 18.09.96 Added padchar()

*************************************************************************/
/* Includes
*/
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MathType.h"
#include "SysDefs.h"
#include "general.h"
#include "macros.h"

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
/*>void StringToLower(char *string1, char *string2)
   ------------------------------------------------
   Input:      char  *string1       A character string
   Output:     char  *string2       Lower case version of string1
   
   This routine converts a lower or mixed case string to lower case.
   
   06.02.91 Original
   28.05.92 ANSIed
   07.01.93 Checks case before converting for SysV
*/
void StringToLower(char *string1,
                   char *string2)
{
   int i;
   
   for(i=0;i<strlen(string1);i++)
      if(isupper(string1[i])) 
         string2[i]=tolower(string1[i]);
      else
         string2[i]=string1[i];
   string2[i]='\0';
}


/************************************************************************/
/*>void StringToUpper(char *string1, char *string2)
   ------------------------------------------------
   Input:      char  *string1       A character string
   Output:     char  *string2       Upper case version of string1
   
   This routine converts a lower or mixed case string to upper case.
   
   06.02.91 Original
   28.05.92 ANSIed
   07.01.93 Checks case before converting for SysV
*/
void StringToUpper(char *string1,
                   char *string2)
{
   int i;
   
   for(i=0;i<strlen(string1);i++)
      if(islower(string1[i])) 
         string2[i]=toupper(string1[i]);
      else
         string2[i]=string1[i];
   string2[i]='\0';
}


/************************************************************************/
/*>char *KillLeadSpaces(char *string)
   ----------------------------------
   Input:      char  *string        A character string
   Returns:   *char                 A pointer to the string with the 
                                    leading spaces removed

   This routine strips leading spaces and tabs from a string returning
   a pointer to the first non-whitespace character.

   N.B. THE MACRO KILLLEADSPACES() MAY NOW BE USED INSTEAD
   
   06.02.91 Original
   28.05.92 ANSIed
   06.07.93 Added tab skipping
*/
char *KillLeadSpaces(char *string)
{
   while (*string == ' ' || *string == '\t') string++;
   return(string);
}


/************************************************************************/
/*>void KillLine(FILE *fp)
   -----------------------
   Input:      FILE  *fp            A file pointer
   
   This routine reads characters from a file till it reaches a '\n' or 
   the end of file.

   08.02.91 Original
   28.05.92 ANSIed
   18.03.94 getc() -> fgetc()
*/
void KillLine(FILE *fp)
{
   int ch;
   
   ch = fgetc(fp);
   while(ch!='\n' && !feof(fp)) ch=fgetc(fp);
}


/************************************************************************/
/*>void SetExtn(char *File, char *Ext)
   -----------------------------------
   I/O:    char  *File     Filename to be modified
   Input:  char  *Ext      New extension

   Force a filename extension. Modifies the input filename to have the
   specified extension. Note that the string File should be large enough
   to cope with longer extensions, if required. Searches back through the
   string to find a `.' Anything after this is replaced with the 
   specified extension. The search for a `.' stops as soon as a `/',`\',
   or a `:' is found as these indicate a directory. If no `.' is found,
   one is appended to the string and the extension is added.

   10.02.91 Original
   28.05.92 ANSIed
   22.06.92 Improved to work from end of string for Unix filenames, etc.
   11.03.94 Added check on '\' for MS-DOS
*/
void SetExtn(char *File, char *Ext)
{
   int   pos,
         InDir    = FALSE;
   
   pos = strlen(File);
   
   while(--pos >= 0 && File[pos] != '.')
   {
      if(File[pos] == ':' || File[pos] == '/' || File[pos] == '\\')
      {
         InDir = TRUE;
         break;
      }
   }
   
   if(pos >= 0 && !InDir)           /* Dot found, change extension      */
   {
      strcpy(File+pos+1, Ext);
      File[pos+1+strlen(Ext)] = '\0';
   }
   else                             /* No dot found, just append ext.   */
   {
      strcat(File,".");
      strcat(File,Ext);
   }
}


/************************************************************************/
/*>int chindex(char *string, char ch)
   ----------------------------------
   Input:      char  *string        A string
                     ch             A character for which to search
   Returns:    int                  The offset of ch in string.
   
   Returns the offset of a character in a string. -1 if not found. This is
   used in a similar manner to strchr(), but gives an offset in the string
   rather than a pointer to the character.

   10.02.91 Original
   28.05.92 ANSIed
   06.10.93 Changed name to chindex() to avoid UNIX name clash
*/
int chindex(char  *string,
            char  ch)
{
   int count;
   
   for(count=0;count<strlen(string);count++)
      if(string[count] == ch) break;
      
   if(count >= strlen(string)) count = -1;
   
   return(count);
}


/************************************************************************/
/*>void Word(char *string1, char *string2)
   ---------------------------------------
   Input:      char  *string1       A string
   Output:     char  *string2       A new string

   Removes leading spaces and extracts a space/tab delimited word. string2
   must have the same amount of space as string1.

   20.03.91 Original
   28.05.92 ANSIed
   22.06.92 Added tab check.
*/
void Word(char  *string1,
          char  *string2)
{
   int j;
   char *str;
   
   str = KillLeadSpaces(string1);
   strcpy(string2,str);
   
   for(j=0;j<strlen(string2);j++)
   {
      if(string2[j] == ' ' || string2[j] == '\t')
      {
         string2[j] = '\0';
         break;
      }
   }
}


/************************************************************************/
/*>void WordN(char *string1, char *string2, int MaxChar)
   -----------------------------------------------------
   Input:      char  *string1       A string
               int   MaxChar        Max characters to copy
   Output:     char  *string2       A new string

   Removes leading spaces and extracts a space/tab delimited word. A 
   maximum of MaxChar characters will be copied into the word. If 
   necessary a terminating NULL will be added *after* MaxChar other 
   characters. Thus MaxChar should always be at least 1 smaller than the 
   size of string2.

   22.06.92 Original based on Word()
*/
void WordN(char  *string1,
           char  *string2,
           int   MaxChar)
{
   int j;
   char *str;
   
   str = KillLeadSpaces(string1);
   strncpy(string2,str,MaxChar);
   string2[MaxChar] = '\0';
   
   for(j=0; j<strlen(string2); j++)
   {
      if(string2[j] == ' ' || string2[j] == '\t')
      {
         string2[j] = '\0';
         break;
      }
   }
}


/************************************************************************/
/*>void padterm(char *string, int  length)
   ---------------------------------------
   I/O:     char  *string   String to be padded with spaces
   Input:   int   length    Required size for string

   Pads a string with spaces to length characters, then terminates it.

   06.09.91 Original    By: ACRM   
*/
void padterm(char *string,
             int  length)
{
   int i;
   
   for(i=strlen(string); i<length; i++)
      string[i] = ' ';
   string[length] = '\0';
}


/************************************************************************/
/*>void padchar(char *string, int length, char ch)
   -----------------------------------------------
   I/O:     char  *string   String to be padded with spaces
   Input:   int   length    Required size for string
            char  ch        Character with which to pad

   Pads a string with a specified character to length characters, then 
   terminates it.

   18.09.96 Original based on padterm()    By: ACRM   
*/
void padchar(char *string,
             int  length,
             char ch)
{
   int i;
   
   for(i=strlen(string); i<length; i++)
      string[i] = ch;
   string[length] = '\0';
}


/************************************************************************/
/*>BOOL CheckExtn(char *string, char *ext)
   ---------------------------------------
   Input:   char *string    String to be checked for given extension
            char *extn      Extension to check for
   Returns: BOOL            Found?

   Check the extension of a filename. For use on machines like VAXes,
   MS-DOS and Amigas, everything is converted to upper case first.
   18.06.93 Original    By: ACRM
*/
BOOL CheckExtn(char  *string,
               char  *ext)
{
   int   extl     = strlen(ext),
         strl     = strlen(string);
   char  *buff1,
         *buff2;
   BOOL  RetVal   = TRUE;
         
   buff1 = (char *)malloc(strl * sizeof(char));
   buff2 = (char *)malloc(extl * sizeof(char));
   
   if(buff1==NULL || buff2==NULL)
   {
      if(buff1) free(buff1);
      if(buff2) free(buff2);
      return(FALSE);
   }
         
   StringToUpper(string,buff1);
   StringToUpper(ext,   buff2);
         
   if(strncmp(buff2,buff1+(strl-extl),extl))
      RetVal = FALSE;
      
   free(buff1);
   free(buff2);
   return(RetVal);
}


/************************************************************************/
/*>void GetFilestem(char *filename, char *stem)
   --------------------------------------------
   Input:   char  *filename      Complete filename
   Output:  char  *stem          The filestem

   Extracts the filestem from a complete filename. Should work under
   Unix, VMS, MS-DOS, AmigaDOS, etc.

   14.04.94 Original    By: ACRM
*/
void GetFilestem(char *filename, char *stem)
{
   char *p, 
        *q;

   q = filename;

   /* First step past a ] if found (VMS)                                */
   if((p=strchr(q,']'))!=NULL)      q=p+1;

   /* Step past any colons/double colons (VMS, AmigaDOS, Unix, MS-DOS)  */
   while((p=strchr(q,':'))!=NULL)   q=p+1;

   /* Step past any / (Unix, AmigaDOS)                                  */
   while((p=strchr(q,'/'))!=NULL)   q=p+1;
   
   /* Step past any \ (MS-DOS)                                          */
   while((p=strchr(q,'\\'))!=NULL)  q=p+1;
   
   /* We should now have the actual filename, with the path removed.
      Copy it into our output array
   */
   strcpy(stem, q);
   
   /* Terminate at the last . to remove the extension.                  */
   q = stem;
   for(p = q + strlen(q) - 1; p!=q; p--)
   {
      if(*p == '.')
      {
         *p = '\0';
         break;
      }
   }
}


/************************************************************************/
/*>int upstrcmp(char *word1, char *word2)
   --------------------------------------
   Input:   char *word1     First word
            char *word2     Second word
   Returns: int             0 if strings match or offset of first 
                            mismatched character

   Like strcmp(), but upcases each character before comparison

   20.04.94 Original   By: ACRM
*/
int upstrcmp(char *word1, char *word2)
{
   int i;

   for(i=0; word1[i] && word2[i]; i++)
      if((islower(word1[i])?toupper(word1[i]):word1[i]) != 
         (islower(word2[i])?toupper(word1[i]):word2[i])) return(i+1);

   if(word1[i] || word2[i]) return(i+1);
   
   return(0);
}


/************************************************************************/
/*>int upstrncmp(char *word1, char *word2, int ncomp)
   --------------------------------------------------
   Input:   char *word1     First word
            char *word2     Second word
            int  ncomp      Number of characters to compare
   Returns: int             0 if strings match or offset of first 
                            mismatched character

   Like strncmp(), but upcases each character before comparison

   20.04.94 Original   By: ACRM
*/
int upstrncmp(char *word1, char *word2, int ncomp)
{
   int i;

   for(i=0; i<ncomp; i++)
   {
      if(!word1[i] || !word2[i]) return(i+1);
      
      if((islower(word1[i])?toupper(word1[i]):word1[i]) != 
         (islower(word2[i])?toupper(word1[i]):word2[i])) return(i+1);
   }
   
   return(0);
}


/************************************************************************/
/*>char *GetWord(char *buffer, char *word)
   ---------------------------------------
   Input:   char  *buffer        String buffer
   Output:  char  *word          Word extracted from buffer
   Returns: char  *              Pointer into buffer after word

   Extracts a comma or white space delimited word from a buffer and
   returns a pointer into the buffer after the word has been removed.
   ' or " may be used to define a word containing white space.

   By default, the bounding ' or " are not returned as part of the word.
   Calling the routine with 2 NULL pointers will switch the mode such
   that the bounding inverted commas are returned.
   Calling again will switch back to the default mode.

   12.04.94 Original    By: ACRM
   26.04.94 Now handles groups of characters in inverted commas as a
            single word.
   11.05.94 Added special call; if called with 2 NULL pointers, switches
            the treatment of ' and " to be included or not in the output
*/
char *GetWord(char *buffer, char *word)
{
   char        *p;
   int         i         = 0,
               j         = 0,
               Commas    = 0;     /* 1: In singles, 2: In doubles       */
   BOOL        GotComma  = FALSE;
   static BOOL CommaMode = FALSE;
   
   /* Check for special calls                                           */
   if(buffer == NULL && word == NULL)
   {
      CommaMode = !CommaMode;
      return(NULL);
   }

   /* Return a blank string if the input buffer is NULL                 */
   word[0] = '\0';
   if(buffer==NULL) return(NULL);

   /* Remove leading spaces                                             */
   KILLLEADSPACES(p, buffer);

   /* Copy up to next comma or white space                              */
   for(i=0; p[i]; i++)
   {
      GotComma = FALSE;  /* Assume this character is not an inv comma   */
      
      /* Set state machine if it is a ' or "                            */
      switch(p[i])
      {
      case '\'':                 /* Got a single inverted comma         */
         switch(Commas)
         {
         case 0:                 /* Not yet in commas, set state 1      */
            Commas   = 1;
            GotComma = TRUE;
            break;
         case 1:                 /* In single commas, set state 0       */
            Commas   = 0;
            GotComma = TRUE;
            break;
         case 2:                 /* In double commas, ignore single     */
            break;
         }
         break;
      case '"':                  /* Got a double inverted comma         */
         switch(Commas)
         {
         case 0:                 /* Not yet in commas, set state 2      */
            Commas = 2;
            GotComma = TRUE;
            break;
         case 1:                 /* In single commas, ignore double     */
            break;
         case 2:                 /* In double commas, set state 0       */
            Commas = 0;
            GotComma = TRUE;
            break;
         }
         break;
      default:                   /* Got some other character            */
         break;
      }
      
      /* If we're not in commas (state 0) and we get a , or white space
         then our word has ended.
      */
      if(!Commas && (p[i]==',' || p[i]==' ' || p[i]=='\t'))
         break;

      /* If this wasn't a state-changing ' or " then store it           */
      if(CommaMode || !GotComma)
         word[j++] = p[i];
   }

   /* Terminate output string                                           */
   word[j] = '\0';

   /* Move p on to the next word                                        */
   p += i;                  /* Move p onto the next character           */
   KILLLEADSPACES(p,p);     /* Strip any leading spaces                 */
   if(*p == ',') p++;       /* Kill a comma if found                    */

   if(*p == '\0') p = NULL;
   
   return(p);
}


/************************************************************************/
/*>BOOL OpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out)
   ---------------------------------------------------------------------
   Input:   char     *infile     Input filename
            char     *outfile    Output filename
   Output:  FILE     **in        Input file pointer
            FILE     **out       Output file pointer
   Returns: BOOL                 Success?

   Open the files if specified. In and out are not modified if files
   are not specified.

   29.06.94 Original    By: ACRM
   24.08.94 Name changed from OpenFiles() and placed in gen lib.
   18.12.95 Now treats a filename of - as stdin/stdout
   05.02.96 Allows NULL pointers
*/
BOOL OpenStdFiles(char *infile, char *outfile, FILE **in, FILE **out)
{
   if(infile!=NULL && infile[0] && strcmp(infile,"-"))
   {
      if((*in = fopen(infile,"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",infile);
         return(FALSE);
      }
   }
      
   if(outfile!=NULL && outfile[0] && strcmp(outfile,"-"))
   {
      if((*out = fopen(outfile,"w"))==NULL)
      {
         fprintf(stderr,"Unable to open output file: %s\n",outfile);
         return(FALSE);
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv)
   ---------------------------------------------------------------------
   Input:     char    *filename     Filename to be opened
              char    *envvar       Unix/MS-DOS environment variable
                                    Other OS assign name (with :)
              char    *mode         Mode in which to open file (r, w, etc)
   Output:    BOOL    *noenv        Set to TRUE under Unix/MS-DOS if 
                                    the reason for failure was that the
                                    environment variable was not set.
   Returns:   FILE    *             File pointer or NULL on failure

   Attempts to open a filename as specified. Returns a file
   pointer. If this fails:

   Under UNIX/MS-DOS:
   gets the contents of the envvar environment variable and prepends
   that to the filename and tries again. If envvar was not set, noenv
   is set to TRUE and the routine returns a NULL pointer.

   Under other OSs:
   prepends the envvar string onto the filename and tries to open the
   file again.

   Returns the pointer returned by the open() command after all this.

   22.09.94 Original    By: ACRM
   11.09.94 Puts a : in for the assign type.
   24.11.94 Added __unix define. Checks for trailing / in environment
            variable
   08.03.95 Corrected basename to filename in non-unix version
   09.03.95 Checks that filename is not a NULL or blank string
*/
FILE *OpenFile(char *filename, char *envvar, char *mode, BOOL *noenv)
{
   char *datadir,
        buffer[160];
   FILE *fp;
   
   if(filename == NULL || filename[0] == '\0')
      return(NULL);

   if(noenv != NULL) *noenv = FALSE;

   /* Try to open the filename as specified                             */
   if((fp=fopen(filename,mode)) == NULL)
   {
      /* Failed, so build alternative directory/filename                */
#if (unix || __unix__ || msdos || __msdos__ || __unix)
      if((datadir = getenv(envvar)) != NULL)
      {
         if(datadir[strlen(datadir)-1] == '/')
            sprintf(buffer,"%s%s",datadir,filename);
         else
            sprintf(buffer,"%s/%s",datadir,filename);
         fp = fopen(buffer,mode);
      }
      else
      {
         if(noenv != NULL) *noenv = TRUE;
         return(NULL);
      }
#else
      sprintf(buffer,"%s:%s",envvar,filename);
      fp = fopen(buffer,mode);
#endif
   }
   
   return(fp);
}


/************************************************************************/
/*>int countchar(char *string, char ch)
   ------------------------------------
   Input:   char     *string      String to search for characters
            char     ch           Character for which to search
   Returns: int                   Number of occurrences of ch in string

   Counts occurrences of charcater ch in string, string.

   17.07.95 Original    By: ACRM
*/
int countchar(char *string, char ch)
{
   char *chp;
   int  count = 0;
   
   if(string==NULL)
      return(0);
   
   for(chp=string, count=0; *chp; chp++)
   {
      if(*chp == ch)
         count++;
   }
   return(count);
}


/************************************************************************/
/*>STRINGLIST *StoreString(STRINGLIST *StringList, char *string)
   -------------------------------------------------------------
   Input:     STRINGLIST  *StringList   The current linked list or NULL 
                                        if nothing yet allocated
              char        *string       The string to store
   Returns:   STRINGLIST  *             Start of linked list. Used on
                                        first call (when input StringList
                                        is NULL) to return the pointer to
                                        the start of the linked list.
                                        NULL if unable to allocate.

   Stores strings (of any length) in a linked list of type STRINGLIST.
   Return a pointer to the start of the linked list which is used on
   the first call to access the newly allocated memory.

   If allocation fails, memory allocated so far is freed and the routine
   returns NULL.

   06.11.95 Original    By: ACRM
*/
STRINGLIST *StoreString(STRINGLIST *StringList, char *string)
{
   STRINGLIST *p, 
              *start;
   
   if((StringList!=NULL) && ((string == NULL) || (string[0] == '\0')))
      return(StringList);
   
   /* Set p to the String List and move to the end of the linked list   */
   start = StringList;

   /* If nothing in the list, initialise it                             */
   if(start == NULL)
   {
      INIT(start,STRINGLIST);
      p=start;
   }
   else  /* Move to end of current list and add another item            */
   {
      p=start;
      LAST(p);

      /* Only allocate another slot if a string is inserted in this one */
      if(p->string != NULL)
         ALLOCNEXT(p,STRINGLIST);
   }
   
   /* Check allocation                                                  */
   if(p==NULL)
   {
      /* If failed, free the list so far and return NULL                */
      FREELIST(start, STRINGLIST);
      return(NULL);
   }
   p->string = NULL;
   
   /* Everything OK, allocate memory for the string                     */
   if((string != NULL) && (string[0] != '\0'))
   {
      if((p->string = (char *)malloc((1+strlen(string))*sizeof(char)))
         ==NULL)
      {
         /* No memory, free linked list and return                      */
         FREELIST(start, STRINGLIST);
         return(NULL);
      }
      
      /* Still OK, copy in the string and return                        */
      strcpy(p->string,string);
   }
   
   return(start);
}


/************************************************************************/
/*>BOOL InStringList(STRINGLIST *StringList, char *string)
   -------------------------------------------------------
   Input:     STRINGLIST   *StringList    Linked list of strings
              char         *string        String for which to search
   Returns:   BOOL                        Is string found?

   Searches a string list for an *exact match* with the specified string.

   06.11.95 Original    By: ACRM
*/
BOOL InStringList(STRINGLIST *StringList, char *string)
{
   STRINGLIST *p;
   
   if(string != NULL)
   {
      for(p=StringList; p!=NULL; NEXT(p))
      {
         if((p->string != NULL)  && !strcmp(p->string, string))
            return(TRUE);
      }
   }
   
   return(FALSE);
}


/************************************************************************/
/*>void FreeStringList(STRINGLIST *StringList)
   -------------------------------------------
   Input:     STRINGLIST   *StringList    Linked list of strings

   Frees memory allocated for a string list.

   06.11.95 Original    By: ACRM
*/
void FreeStringList(STRINGLIST *StringList)
{
   STRINGLIST *p;
   
   for(p=StringList; p!=NULL; NEXT(p))
   {
      if(p->string != NULL)
         free(p->string);
   }
   
   FREELIST(StringList, STRINGLIST);
}


/************************************************************************/
/*>char *QueryStrStr(char *string, char *substring)
   ------------------------------------------------
   This is like strstr() but allows a ? character in the substring
   which matches any character.

   15.12.95 Original    By: ACRM
*/
char *QueryStrStr(char *string, char *substring)
{
   int  i, j, 
        lenstr, 
        lensubstr;
   char *retval = NULL;

   /* Get lengths of the 2 strings                                      */
   lenstr    = strlen(string);
   lensubstr = strlen(substring);
   
   /* Walk through the main string                                      */
   for(i=0; i<=lenstr-lensubstr; i++)
   {
      /* Set return value to point to current position in main string   */
      retval = string + i;
      
      /* Walk through substring                                         */
      for(j=0; j<lensubstr; j++)
      {
         /* If there is a mismatch, set return to NULL and break out of
            the substring counting
         */
         if((substring[j] != '?') && 
            (substring[j] != '%') && 
            (string[i+j] != substring[j]))
         {
            retval = NULL;
            break;
         }
      }

      /* If retval hasn't been set to NULL, then we've got a match, so
         return pointer
      */
      if(retval)
         return(retval);
   }

   /* This should always be NULL...                                     */
   return(retval);
}

