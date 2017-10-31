// Internal Boundary Routines - Author: Angelos Chronis

float avdens;  
float avvelX;
float avvelY;
float avvelZ;

boolean converged = false;

void writeData()
{
  for (int i = 0; i < nX; i++ )
  {
    for (int j = 0; j < nY; j++ )
    {
      for (int k = 0; k < nZ; k++ )
      {
        out.println( i +","+ j +","+ k +","+ dm[I(i, j, k)] +","+ um[I(i, j, k)]  +","+ vm[I(i, j, k)]  +","+ wm[I(i, j, k)]);
      }
    }
  }
  out.flush();
}


void calculateMeans()
{
  for (int i=0; i<size; i++)
  {
    um[i] += u[i]; 
    um[i] *=0.5;
    vm[i] += v[i]; 
    vm[i] *=0.5;
    wm[i] += w[i]; 
    wm[i] *=0.5;
    dm[i] += d[i]; 
    dm[i] *=0.5;
  }
}

void densityConvergence(float threshold)
{
  if (!converged)
  {
    float curravd = 0.0;
    for (int i = 1; i <= nX; i++)
    {
      for (int j = 1; j <= nY; j++)
      {
        for (int k = 1; k <= nZ; k++)
        {
          curravd += d[I(i, j, k)];
        }
      }
    }
    curravd /= nX*nY*nZ;

    float diff = abs(avdens - curravd);
    avdens = curravd;

    if ( diff < threshold)
    {
      converged = true;
      println ("Converged!");
    }
    else converged = false;
  }
}
 