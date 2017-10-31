// External Boundaries - Author: Angelos Chronis

void setBoundry(int b, float[] x)
{
  for (int j = 0; j <= nY+1; j++) {
    for (int i = 0; i <= nX+1; i++) {
      //        x[I(i, j, 0    )] = 0.0;
      //        x[I(i, j, nZ+1 )] = b == 3 ? 0.0 : x[I(i, j, nZ  )];
    }
  }
  for (int k = 0; k <= nZ+1; k++) {
    for (int i = 0; i <= nX+1; i++) {
      if      (windDir==3) x[I(i, 0, k)] = x[I(i, 1, k)];  //SOUTH
      else if (windDir==2) x[I(i, nY+1, k)] = x[I(i, nY, k)];  //NORTH
      else
      {
        x[I(i, 0, k)] = b == 2 ?  0.0 : x[I(i, 1, k)];     //EAST-WEST
        x[I(i, nY+1, k)] = b == 2 ?  0.0 : x[I(i, nY, k)];
      }
    }
  }
  for (int k = 0; k <= nZ+1; k++) {
    for (int j = 0; j <= nY+1; j++) {
      if      (windDir==0) x[I(0, j, k)] = x[I(1, j, k)];  //EAST
      else if (windDir==1) x[I(nX+1, j, k)] = x[I(nX, j, k)];  //WEST
      else
      {
        x[I(0, j, k)] = b == 1 ?  0.0 : x[I( 1, j, k)];    //NORTH-SOUTH
        x[I(nX+1, j, k)] = b == 1 ?  0.0 : x[I( nX, j, k)];
      }
    }
  }
  x[I(0, 0, 0)]       = 0.33f * (x[I(1, 0, 0)]
    + x[I(0, 1, 0)]
    + x[I(0, 0, 1)]);
  x[I(0, nY+1, 0)]     = 0.33f * (x[I(1, nY+1, 0)]
    + x[I(0, nY, 0)]
    + x[I(0, nY+1, 1)]);
  x[I(0, 0, nZ+1)]     = 0.33f * (x[I(1, 0, nZ+1)]
    + x[I(0, 1, nZ+1)]
    + x[I(0, 0, nZ+1)]);
  x[I(0, nY+1, nZ+1)]   = 0.33f * (x[I(1, nY+1, nZ+1)]
    + x[I(0, nY, nZ+1)]
    + x[I(0, nY+1, nZ)]);
  x[I(nX+1, 0, 0)]     = 0.33f * (x[I(nX, 0, 0)]
    + x[I(nX+1, 1, 0)]
    + x[I(nX+1, 0, 1)]);
  x[I(nX+1, nY+1, 0)]   = 0.33f * (x[I(nX, nY+1, 0)]
    + x[I(nX+1, nY, 0)]
    + x[I(nX+1, nY+1, 1)]);
  x[I(nX+1, 0, nZ+1)]   = 0.33f * (x[I(nX, 0, nZ+1)]
    + x[I(nX+1, 1, nZ+1)]
    + x[I(nX+1, 0, nZ)]);
  x[I(nX+1, nY+1, nZ+1)] = 0.33f * (x[I(nX, nY+1, nZ+1)]
    + x[I(nX+1, nY, nZ+1)]
    + x[I(nX+1, nY+1, nZ)]);
  setBoundryIntrn(b, x);
}

void inletDensity()
{
  if (windDir<2)
  {
    for (int j= 1; j<=nY; j++ )
    {  
      for (int k= 1; k<=nZ; k++ )
      {
        if      (windDir==0) dOld[I(nX, j, k)] = 1.20f;    //____EAST
        else if (windDir==1) dOld[I(1, j, k)]  = 1.20f;    //____WEST
      }
    }
  } else
  {
    for (int i= 1; i<=nX; i++ )
    {  
      for (int k= 1; k<=nZ; k++ )
      {
        if      (windDir==2) dOld[I(i, 1, k)]  = 1.20f;   //____NORTH
        else if (windDir==3) dOld[I(i, nY, k)] = 1.20f;   //____SOUTH
      }
    }
  }
}

void inletVelocity()
{
  if (windDir<2)
  {
    for (int j= 1; j<=nY; j++ )
    {  
      for (int k= 1; k<=nZ; k++ )
      {
        if      (windDir==0) uOld[I(nX, j, k)]  = -inletV[k];    //____EAST
        else if (windDir==1) uOld[I(1, j, k)]  =  inletV[k];     //____WEST
      }
    }
  } else
  {
    for (int i= 1; i<=nX; i++ )
    {  
      for (int k= 1; k<=nZ; k++ )
      {
        if      (windDir==2) vOld[I(i, 1, k)]   =  inletV[k];     //____NORTH 
        else if (windDir==3) vOld[I(i, nY, k)]  = -inletV[k];    //____SOUTH
      }
    }
  }
}
void createVelocityProfile(int h)
{
  for (int i=0; i<inletV.length; i++) 
  {
    inletV[i] = 13.0;
  }
}

//void createVelocityProfile(int h)
//{
//  /*****************************************************************************/
//  /* UDF for specifying the wind profile using the expressions provided by     */
//  /* Richards and Hoxey   (1993)                                                 */
//  /* Assumes:                                                                  */
//  /* Z direction is upwards from the ground                                    */
//  /*****************************************************************************/
//  float ZGROUND  = 0.0;    /* Z location of the ground relative to the domain */
//  float RHT      = 1.0;    /* Aerodynamic surface roughness z                 */
//  float UREF     = 10.0;    /* Velocity at Reference height                    */
//  float ZREF     = 30.0;   /* Reference height                                */
//  float CAPPA    = 0.41;   /* von Karmans constant                            */
//  float C_MU     = 0.09;   /* Turbulence Model Parameter                      */
//  float U_STAR   = (CAPPA*UREF/log(ZREF/RHT)); /* Expression for friction velocity */

//  float   zloc; 

//  for (int i=0; i<inletV.length; i++) 
//  {
//    zloc=  i * cell - cell2 - ZGROUND;
//    inletV[i] = (U_STAR/0.41)*log(zloc/RHT);
//  }
//}