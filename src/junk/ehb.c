
/********************************************************************/
void EHB2(int x, int y, int rgb[3])
{
   int pen,
       diff;
       
   pen = NearestPenEHB(rgb,&diff);
   
   SetAPen(rp,pen);
   WritePixel(rp,x,y);
}
/********************************************************************/
int NearestPenEHB(int *c, int *dist)
{
   int   i,
         mindist,
         nearest,
         d1,d2;
 
   mindist=32000;
 
   for(i=0; i<nallocr; ++i)
   {
      d1=ColourDist2EHB(c,creg[i]);
      if(d1 < mindist)
      {
         mindist=d1;
         nearest=i;
      }
      d2=ColourDist2bEHB(c,creg[i]);
      if(d2 < mindist)
      {
         mindist=d2;
         nearest=i+32;
      }
   }
   if(mindist > THRESHHOLD && nallocr < 32)
   {
      if(c[0] < 8 && c[1] < 8 && c[2] < 8)
      {
         for(i=0; i<3; ++i) creg[nallocr][i]=2*c[i];
         SetRGB4(viewport,nallocr, c[0],c[1],c[2]);
         nearest=32+nallocr++;
      }
      else
      {
         for(i=0; i<3; ++i) creg[nallocr][i]=c[i];
         SetRGB4(viewport,nallocr, c[0],c[1],c[2]);
         nearest=nallocr++;
      }
   }
   *dist = mindist;

   return(nearest);
}
 
/********************************************************************/
int ColourDist2EHB(int *a, int *b)   /* the 'distance' between two colors */
{
   int   k,
         r,
         d;
 
   for(r=k=0; k<3; ++k)
   {
      d = a[k]-b[k];
      if(d < 0) d= -d;

      r += d;
   }

   return(r);
}
/********************************************************************/
int ColourDist2bEHB(int *a, int *b)   /* the 'distance' between two colors */
{
   int   k,
         r,
         d;
 
   for(r=k=0; k<3; ++k)
   {
      d = a[k]-(int)(b[k]/2);
      if(d < 0) d= -d;

      r += d;
   }

   return(r);
}
/********************************************************************/
void EHB(int i, int j,REAL brite[3])
{
   int   k,
         pix[3];

   for(k=0; k<3; ++k)
   {
      pix[k] = 16.0 * brite[k] + 0.5;
      if(pix[k] < 0)    pix[k] = 0;
      if(pix[k] > 15)   pix[k] = 15;
   }
   EHB2(i,j,pix);
}
