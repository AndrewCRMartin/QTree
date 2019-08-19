/*************************************************************************

   Program:    QTree
   File:       writepng.h
   
   Version:    V3.0
   Date:       19.09.19
   Function:   Write a PNG image
   
   Copyright:  (c) SciTech Software 2019
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
   Based loosely on sample code from https://www.lemoda.net/c/write-png/

**************************************************************************

   Usage:
   ======

**************************************************************************

   Notes:
   ======

**************************************************************************

   Revision History:
   =================
   V3.0  19.08.19 Original

*************************************************************************/

#include "writepng.h"

/************************************************************************/
/*>blPNGPIXEL *blPNGPixelAt(blPNGIMAGE *bitmap, int x, int y)
   -----------------------------------------------------------
*//**
   Given a linearized bitmap, find the offset to the pixel at x,y

-  19.08.19 Original   By: ACRM
*/
blPNGPIXEL *blPNGPixelAt(blPNGIMAGE *bitmap, int x, int y)
{
    return bitmap->pixels + bitmap->width * y + x;
}
    
/************************************************************************/
/*>BOOL blSavePNGToFile(blPNGIMAGE *bitmap, const char *path)
   -----------------------------------------------------------
*//**
   Writes an image containing a linearized bitmap to a PNG file.
   If path is a blank string, then writes to stdout.
   Returns TRUE on success, FALSE on failure.

-  19.08.19 Original   By: ACRM
*/
BOOL blSavePNGToFile(blPNGIMAGE *image, const char *path)
{
    FILE        *fp           = stdout;
    png_structp pngPtr        = NULL;
    png_infop   infoPtr       = NULL;
    png_byte    **rowPointers = NULL;
    BOOL        retval        = FALSE;
    size_t      x, y;

    int pixelSize = 3;
    int depth     = 8;

    if(path[0])
       fp = fopen(path, "wb");

    if(fp!=NULL)
    {
       if((pngPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                            NULL, NULL, NULL))!=NULL)
       {
          if((infoPtr = png_create_info_struct(pngPtr))!=NULL)
          {
             /* Set up error handling                                    */
             if(setjmp(png_jmpbuf(pngPtr)))
                goto png_failure;

             /* Set image attributes                                     */
             png_set_IHDR(pngPtr, infoPtr,
                          image->width, image->height,
                          depth,
                          PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                          PNG_COMPRESSION_TYPE_DEFAULT,
                          PNG_FILTER_TYPE_DEFAULT);
             
             /* Initialize rows of PNG                                   */
             rowPointers = png_malloc(pngPtr,
                                      image->height *
                                      sizeof(png_byte *));

             for(y=0; y < image->height; y++)
             {
                png_byte *row = 
                   png_malloc(pngPtr,
                              sizeof(uint8_t) * image->width * pixelSize);
                rowPointers[y] = row;

                for(x=0; x < image->width; x++)
                {
                   blPNGPIXEL *pixel = blPNGPixelAt(image, x, y);
                   *row++ = pixel->red;
                   *row++ = pixel->green;
                   *row++ = pixel->blue;
                }
             }
             
             /* Write the image data to "fp"                             */
             png_init_io(pngPtr, fp);
             png_set_rows(pngPtr, infoPtr, rowPointers);
             png_write_png(pngPtr, infoPtr,
                           PNG_TRANSFORM_IDENTITY, NULL);
             
             retval = TRUE;
             for(y=0; y < image->height; y++)
                png_free(pngPtr, rowPointers[y]);

             png_free(pngPtr, rowPointers);
          }
          else
          {
             retval = FALSE;
          }
       }
       else
       {
          retval = FALSE;
       }
       
    png_failure:
       png_destroy_write_struct(&pngPtr, &infoPtr);

       if(fp != stdout)
          fclose(fp);
    }
    else
    {
       retval = FALSE;
    }
    
    return(retval);
}

#ifdef DEMO
/* Given a value and the maximum value, returns an integer
   between 0 and 255 proportional to value/maxval 
*/
static int pix(int value, int maxval)
{
    if(value < 0)
        return 0;

    return((int)(256.0 *((double)value / (double) maxval)));
}

int main()
{
    blPNGIMAGE image;
    int        x, y,
               status=0;

    /* Create an image                                                   */
    image.width  = 100;
    image.height = 100;
    image.pixels = calloc(image.width * image.height,
                          sizeof(blPNGPIXEL));

    if(!image.pixels)
       return(-1);

    for(y=0; y < image.height; y++)
    {
       for(x=0; x < image.width; x++)
       {
          blPNGPIXEL *pixel = blPNGPixelAt(&image, x, y);
          pixel->red        = pix(x, image.width);
          pixel->green      = pix(y, image.height);
       }
    }

    /* Write the image to a file 'image.png'                             */
    if(!blSavePNGToFile(&image, "image.png")) 
    {
       fprintf(stderr, "Error writing file.\n");
       status = -1;
    }

    free(image.pixels);
    return status;
}

#endif
