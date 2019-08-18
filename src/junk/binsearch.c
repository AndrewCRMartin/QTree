// N.B. This code has a small bug somewhere; some sort of access violation
// occurs in certain conditions...


/************************************************************************/
/*>int FarLeftSearch(SPHERE **spheres, int NSphere, REAL x)
   ---------------------------------------------------------
   Finds the offset into the spheres array of the most leftward atom 
   which is (at least partically) to the right of x. i.e. Find all spheres
   with coords >= x

   19.07.93 Original    By: ACRM
*/
int FarLeftSearch(SPHERE **spheres, int NSphere, REAL x)
{
   int   bottom = 0,
         top    = NSphere,
         icut   = 0;
         
   
   /* Check that at least one value will satisfy                        */
   if((spheres[NSphere-1]->x + spheres[NSphere-1]->rad) < x)
   {
      return(-1);
   }
   
   for(;;)
   {
      if(top==bottom) return(top);
      
      icut = bottom + (int)((top - bottom) / 2);
      
      if((spheres[icut]->x + spheres[icut]->rad) >= x)
      {
         /* In range, so reset top                                      */
         top = icut;
      }
      else
      {
         /* Out of range. If we are at the end of the array or the
            position above *is* in range, this is the position we are 
            looking for.
         */
         if(icut == NSphere-1)
            return(icut);
            
         if((spheres[icut+1]->x + spheres[icut+1]->rad) >= x)
            return(icut+1);
            
         /* Otherwise, reset bottom                                     */
         bottom = icut;
      }
   }
}
   
/************************************************************************/
/*>int FarRightSearch(SPHERE **spheres, int NSphere, REAL x)
   ---------------------------------------------------------
   Finds the offset into the spheres array of the most rightward atom 
   which is (at least partically) to the left of x. i.e. Find all spheres
   with coords <= x
   
   19.07.93 Original    By: ACRM
*/
int FarRightSearch(SPHERE **spheres, int NSphere, REAL x)
{
   int   bottom = 0,
         top    = NSphere,
         icut   = 0;
         
   
   /* Check that at least one value will satisfy                        */
   if((spheres[0]->x - spheres[0]->rad) > x)
   {
      return(-1);
   }
   
   for(;;)
   {
      if(top==bottom || top-bottom == 1)
         return(bottom);
      
      icut = bottom + (top - bottom) / 2;
      
      if((spheres[icut]->x - spheres[icut]->rad) <= x)
      {
         /* In range, so reset bottom                                   */
         bottom = icut;
      }
      else
      {
         /* Out of range. If we are at the end of the array or the
            position above *is* in range, this is the position we are 
            looking for.
         */
         if(icut == 0)
            return(icut);
            
         if((spheres[icut+1]->x - spheres[icut+1]->rad) <= x)
            return(icut-1);
            
         /* Otherwise, reset top                                        */
         top = icut;
      }
   }
}
   


