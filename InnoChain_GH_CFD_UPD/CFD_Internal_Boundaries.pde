// Internal Boundary Routines - Author: Angelos Chronis

void setSolids()
{
  for (int i=1; i<nX-1; i++)
  {      
    for (int j=1; j<nY-1; j++)
    {
      for (int k=1; k<nZ-1; k++)
      {
        if (sldInd[i][j][k]==true)
        {
          int index = I(i, j, k);
          solid[index] =true;
        }
      }
    }
  }
  setBound();
}

void setBoundryIntrn(int b, float[] x)
{
  for (int i = 1; i <= nX; i++)
  {
    for (int j = 1; j <= nY; j++)
    {
      for (int k = 1; k <= nZ; k++)
      { 
        int ind  = I(i, j, k);
        if (solid[ind])
        {
          switch( boundInt[ind] )
          {
          case 000:
            x[ind] = 0.0; 
          case 001:
            x[ind] = b==3 ? - x[I( i, j, k-1)]  : x[I( i, j, k-1)];
            break;
          case 002:
            x[ind] = b==3 ? - x[I( i, j, k+1)]  : x[I( i, j, k+1)];
            break;
          case 010:
            x[ind] = b==2 ? - x[I( i, j-1, k)]  : x[I( i, j-1, k)];
            break;
          case 020:
            x[ind] = b==2 ? - x[I( i, j+1, k)]  : x[I( i, j+1, k)];
            break;
          case 100:
            x[ind] = b==1 ? - x[I( i-1, j, k)]  : x[I( i-1, j, k)];
            break;
          case 200:
            x[ind] = b==1 ? - x[I( i+1, j, k)]  : x[I( i+1, j, k)];
            break;
          case 011:
            if      (b==2) x[ind] = -x[I( i, j-1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.5* (x[I( i, j-1, k   )] + x[I( i, j, k-1 ) ]);    
            break;
          case 012:
            if      (b==2) x[ind] = -x[I( i, j-1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.5* (x[I( i, j-1, k   )] + x[I( i, j, k+1 )] );     
            break;
          case 021:
            if      (b==2) x[ind] = -x[I( i, j+1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.5* (x[I( i, j+1, k   )] + x[I( i, j, k-1 ) ]); 
            break;
          case 022:
            if      (b==2) x[ind] = -x[I( i, j+1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.5* (x[I( i, j+1, k   )] + x[I( i, j, k+1 ) ]);
            break;
          case 101:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.5* (x[I( i-1, j, k   )] + x[I( i, j, k-1 )] );   
            break;
          case 102:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.5* (x[I( i-1, j, k   )] + x[I( i, j, k+1 )] );         
            break;
          case 110:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j-1, k   )];
            else           x[ind] =  0.5* (x[I( i-1, j, k   )] + x[I( i, j-1, k   )] );
            break;
          case 120:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j+1, k   )];
            else           x[ind] =  0.5* (x[I( i-1, j, k   )] + x[I( i, j+1, k   )] );
            break;
          case 201:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.5* (x[I( i+1, j, k   )] + x[I( i, j, k-1 ) ]);         
            break;
          case 202:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.5* (x[I( i+1, j, k   )] + x[I( i, j, k+1 ) ]); 
            break;
          case 210:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j-1, k   )];
            else           x[ind] =  0.5* (x[I( i+1, j, k   )] + x[I( i, j-1, k   ) ]); 
            break;
          case 220:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j+1, k   )];
            else           x[ind] =  0.5* (x[I( i+1, j, k   )] + x[I( i, j+1, k   ) ]);               
            break;
          case 111:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j-1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.33 * (x[I( i-1, j, k   )] + x[I( i, j-1, k  )] + x[I( i, j, k-1 )] );
            break;
          case 112:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j-1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.33 * (x[I( i-1, j, k   )] + x[I( i, j-1, k  )] + x[I( i, j, k+1 )] );
            break;
          case 121:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j+1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.33 * (x[I( i-1, j, k   )] + x[I( i, j+1, k  )] + x[I( i, j, k-1 )] );
            break;
          case 122:
            if      (b==1) x[ind] = -x[I( i-1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j+1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.33 * (x[I( i-1, j, k   )] + x[I( i, j+1, k  )] + x[I( i, j, k+1 )] );
            break;
          case 211:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j-1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.33 * (x[I( i+1, j, k   )] + x[I( i, j-1, k  )] + x[I( i, j, k-1 )] );
            break;
          case 212:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j-1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.33 * (x[I( i+1, j, k   )] + x[I( i, j-1, k  )] + x[I( i, j, k+1 )] );
            break;
          case 221:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j+1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k-1 )];
            else           x[ind] =  0.33 * (x[I( i+1, j, k   )] + x[I( i, j+1, k  )] + x[I( i, j, k-1 )] );
            break;
          case 222:
            if      (b==1) x[ind] = -x[I( i+1, j, k   )];
            else if (b==2) x[ind] = -x[I( i, j+1, k   )];
            else if (b==3) x[ind] = -x[I( i, j, k+1 )];
            else           x[ind] =  0.33 * (x[I( i+1, j, k   )] + x[I( i, j+1, k  )] + x[I( i, j, k+1 )] );
            break;
          }
        }
      }
    }
  }
}

void setBound()
{
  for (int i=1; i<=nX; i++)
  {      
    for (int j=1; j<=nY; j++)
    {
      for (int k=1; k <= nZ; k++)
      {
        int sI=0;
        int sJ=0;
        int sK=0;
        int ind = I(i, j, k);
        if ( solid[ind] )
        {
          if      ( !solid[I(i-1, j, k)] )  sI = 100;
          else if ( !solid[I(i+1, j, k)] )  sI = 200;
          else                              sI = 0;
        
          if      ( !solid[I(i, j-1, k)] )  sJ = 10;
          else if ( !solid[I(i, j+1, k)] )  sJ = 20;
          else                              sJ = 0;
       
          if      ( !solid[I(i, j, k-1)] )  sK = 1;
          else if ( !solid[I(i, j, k+1)] )  sK = 2;
          else                              sK = 0;
        }
        
        boundInt[ind] = sI + sJ + sK;
      }
    }
  }
}