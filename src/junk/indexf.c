void indexf(int n, REAL *arrin, int *indx)
{
   int j,l,ir;
   
   
   for(j=0; j<n; j++)
      indx[j] = j+1
      
   l  = n/2 + 1;
   ir = n;
   
   for(;;)
   {
      if(l>1)
      {
         l--;
         indxt=indx[l-1];
         q=arrin[indxt-1
      }
      else
      {
         indxt = indx[ir-1];
         q=arrin[indxt-1];
         indx[ir-1] = indx[0];
         ir--;
         if(ir==1)
         {
            indx[0] = indxt;
            return;
         }
      }
      
      i=l;
      j=l+l;
      
      while(j<=ir)
      {
         if(j<ir)
         {
            if(arrin[indx[j-1]-1] < arrin[indx[j]-1])
               j++;
         }
         
         if(q<arrin[indx[j-1]-1])
         {
            indx[i-1] = indx[j-1];
            i=j;
            j += j;
         }
         else
         {
            j = ir+1;
         }
      }
      indx[i-1]=indxt;
   }
}
      