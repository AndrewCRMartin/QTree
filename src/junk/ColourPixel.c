

/************************************************************************/
/*>SPHERE **SortSpheresOnZ(SPHERE **Spheres, int NSphere)
   --------------------------------------------------------
   Performs a heapsort returning an array of pointers which represent the 
   input array AllSpheres sorted on x.
   
   19.07.93 Original    By: ACRM
   20.07.93 Made j and l register ints
   23.07.93 Changed to sort on zmax rather than z
*/
SPHERE **SortSpheresOnZ(SPHERE **Spheres, int NSphere)
{
   SPHERE         **sp = NULL,
                  *temp;
   int            i, ir;
   register int   j, l;
   REAL           q;
   
   /* Allocate memory for index                                         */
   if((sp = (SPHERE **)malloc(NSphere * sizeof(SPHERE *))) == NULL)
      return(NULL);

   /* Set all pointers to input ordering                                */
   for(j=0; j<NSphere; j++) 
      sp[j] = Spheres[j];

   if(NSphere > 1)
   {
      l  = NSphere/2 + 1;
      ir = NSphere;
   
      for(;;)
      {
         if(l>1)
         {
            temp = sp[--l - 1];
            q    = temp->z + temp->rad;
         }
         else
         {
            temp = sp[ir-1];
            q    = temp->z + temp->rad;
            
            sp[ir-1]=sp[0];
            if(--ir == 1)
            {
               sp[0]=temp;
               return(sp);
            }
         }
         
         i = l;
         j = l+l;
   
         while(j<=ir)
         {
            if(j<ir)
            {
               if((sp[j-1])->z + (sp[j-1])->rad < 
                  (sp[j])->z   + (sp[j])->rad) j++;
            }
            
            if(q<(sp[j-1])->z + (sp[j-1])->rad)
            {
               sp[i-1] = sp[j-1];
               i       = j;
               j      += j;
            }
            else
            {
               j       = ir+1;
            }
         }
         sp[i-1]       = temp;
      }
   }

   /* We get here if we had only one sphere                             */
   return(sp);
}


/************************************************************************/
/*>void ColourPixel(REAL x, REAL y, SPHERE **spheres, int NSphere)
   ---------------------------------------------------------------
   Sort the spheres along Z which line up with this pixel. Step back from
   the front sphere checking max z on spheres farther back until we find
   one farther away. We've then found the front-most sphere, so we can
   call the shading routine.
   19.07.93 Original    By: ACRM
   20.07.93 Made q and z register; Fixed Z search just to look at nearest
            5 centre points.
   21.07.93 Made x and y real
*/
void ColourPixel(REAL x, REAL y, SPHERE **spheres, int NSphere)
{
   SPHERE         **SrtSph = NULL;
   REAL           XOff,
                  YOff,
                  MaxZ;
   register REAL  q, z;
   int            i,
                  k,
                  FrontSphere = (-1);
   
   
static int MaxSphere = 0;

if(NSphere > MaxSphere)
{
   printf("NSphere = %d\n",NSphere);
   MaxSphere = NSphere;
}   
   
   
   
   /* Sort the sphere list for this pixel on z                          */
   if((SrtSph = SortSpheresOnZ(spheres, NSphere)) != NULL)
   {
      /* Search back through 5 spheres for the nearest z position at this
         pixel. Only look at 5 since an atom will only be bonded to 4
         others.
      */
      for(i=NSphere-1, k=0; i>=0 /* && k<5 */; i--, k++)
      {
         XOff = x - SrtSph[i]->x;
         YOff = y - SrtSph[i]->y;
         
         q = (SrtSph[i]->rad * SrtSph[i]->rad) - 
             (XOff * XOff) -
             (YOff * YOff);
         
         if(q >= 0.0)
         {
            /* Find z for this sphere on this pixel                     */
            z = sqrt(q) + SrtSph[i]->z;
            
            if(FrontSphere == (-1))
            {
               MaxZ = z;
               FrontSphere = i;
            }
            else
            {
               if(z > MaxZ)
               {
                  MaxZ = z;
                  FrontSphere = i;
               }
               
               else
               {
                  break;
               }
            }
         }
             
      }
      
      /* Shade the pixel                                                */
      if(FrontSphere != (-1))
         ShadePixel(x, y, (REAL)MaxZ, SrtSph[FrontSphere]);
      
      free(SrtSph);
   }
   
}

