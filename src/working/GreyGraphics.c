#include <stdio.h>

#include <intuition/intuition.h>
#include <proto/all.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

#include "qtree.h"
#include "graphics.p"


#define XSIZE 256
#define YSIZE 256


struct IntuitionBase *IntuitionBase = NULL;
struct GfxBase       *GfxBase       = NULL;
struct Screen        *scrn          = NULL;
struct Window        *wind          = NULL;

/************************************************************************/
BOOL InitGraphics(void)
{
   int i;
   
   /* Open intuition library                                            */
   if((IntuitionBase = (struct IntuitionBase *)
                       OpenLibrary("intuition.library",36)) == NULL)
   {
      return(FALSE);
   }
      
   /* Open graphics library                                             */
   if((GfxBase = (struct GfxBase *)
                 OpenLibrary("graphics.library",36)) == NULL)
   {
      return(FALSE);
   }
      
                   
   /* Open screen                                                       */
   if((scrn = OpenScreenTags(NULL, 
//                             SA_LikeWorkbench, TRUE,
                             SA_ShowTitle,   FALSE,
                             SA_Type,        CUSTOMSCREEN,
                             SA_Overscan,    OSCAN_STANDARD,
                             SA_Left,        0,
                             SA_Top,         0,
                             SA_Width,       XSIZE,
                             SA_Height,      YSIZE,
                             SA_Depth,       4,
                             SA_AutoScroll,  TRUE,
                             TAG_END)) == NULL)
   {
      return(FALSE);
   }
   
   /* Open window                                                       */
   if((wind = OpenWindowTags(NULL,
                             WA_Flags,          WFLG_BACKDROP |
                                                WFLG_BORDERLESS,
                             WA_CustomScreen,   scrn,
                             TAG_END)) == NULL)
   {
      return(FALSE);
   }
   
   /* Initialise colour map                                             */
   for(i=0;i<16;i++)
   {
      SetRGB4(&scrn->ViewPort,i,i,i,i);
   }
   SetRGB4(&scrn->ViewPort,1,15,15,15);

   SetAPen(wind->RPort, 1);

   return(TRUE);
}

/************************************************************************/
void EndGraphics(void)
{
   char  buffer[80];

   /* Wait for return to be pressed unless gQuickExit is set            */
   if(gQuickExit)
      printf("\n");
   else
   {
      printf("\nPress <return> to exit...");
      gets(buffer);
   }

   if(wind != NULL)           CloseWindow(wind);
   if(scrn != NULL)           CloseScreen(scrn);
   if(GfxBase != NULL)        CloseLibrary((struct Library *)
                                           GfxBase);
   if(IntuitionBase != NULL)  CloseLibrary((struct Library *)
                                           IntuitionBase);
}

/************************************************************************/
void SetPixel(int x0, int y0, REAL r, REAL g, REAL b)
{
   REAL  sum;
   int   pen;
   
   sum = r+g+b / 3;
   
   pen = (int)(16*sum);
   if(pen > 15) pen = 15;
   
   SetAPen(wind->RPort, pen);
   WritePixel(wind->RPort, x0, YSIZE - y0 - 1);
}

