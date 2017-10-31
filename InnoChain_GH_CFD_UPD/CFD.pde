/* Multithreaded solver based on:
 * FluidSolver.java by Alexander McKenzie (12 March, 2004)
 * and:
 * Jos Stam style fluid solver with vorticity confinement
 * and buoyancy force.
 * Author: Angelos Chronis
 */

class FluidSolver
{
  int nX, nY, nZ, size;

  float visc = 0.0;
  float diff = 0.0;
  float dt   = 0.1;

  float[] tmp;

  float[] d, dOld;
  float[] u, uOld;
  float[] v, vOld;
  float[] w, wOld;
  float[] curl;

  float[] dm;
  float[] um;
  float[] vm;
  float[] wm;

  boolean[] solid;
  int[] boundInt;
  float dg, dg_2;
  float cell, cell2;

  float[] inletV;
  PImage sectionX, sectionZ, sectionY;
  Mesh mesh;

  FluidSolver(PApplet p, int x, int y, int z, float cel)
  {
    nX= x;
    nY= y;
    nZ= z;

    size = (nX + 2) *  (nY + 2) * (nZ + 2);
    cell = cel;
    cell2 = cel/2;
    sectionX = createImage(nY+2, nZ+2, ARGB);
    sectionZ = createImage(nX+2, nY+2, ARGB);
    sectionY = createImage(nX+2, nZ+2, ARGB);
    mesh=new Mesh(p);
  }

  /**
   * Reset the datastructures.
   * We use 1d arrays for speed.
   **/

  void reset()
  {
    d    = new float[size];
    dOld = new float[size];
    u    = new float[size];
    uOld = new float[size];
    v    = new float[size];
    vOld = new float[size];
    w    = new float[size];
    wOld = new float[size];
    curl = new float[size];
    solid = new boolean[size];
    boundInt = new int[size];    

    dm = new float[size];
    um = new float[size];
    vm = new float[size];
    wm = new float[size];

    inletV = new float[nZ+2];

    for (int i = 0; i < size; i++)
    {
      u[i] = uOld[i] = v[i] = vOld[i] = w[i] = wOld[i] = 0.0f;
      d[i] = dOld[i] = curl[i] = 0.0f;

      um[i] = vm[i] = wm[i] = 0.0f;
      solid[i] = false;
      boundInt[i] = 0;
    }
    
    sectionX = createImage(nY+2, nZ+2, ARGB);
    sectionZ = createImage(nX+2, nY+2, ARGB);
    sectionY = createImage(nX+2, nZ+2, ARGB);
  }

  /**
   * Calculate the buoyancy force as part of the velocity solver.
   * Fbuoy = -a*d*Y + b*(T-Tamb)*Y where Y = (0,1). The constants
   * a and b are positive with appropriate (physically meaningful)
   * units. T is the temperature at the current cell, Tamb is the
   * average temperature of the fluid grid. The density d provides
   * a mass that counteracts the buoyancy force.
   *
   * In this simplified implementation, we say that the tempterature
   * is synonymous with density (since smoke is *hot*) and because
   * there are no other heat sources we can just use the density
   * field instead of a new, seperate temperature field.
   *
   * @param Fbuoy Array to store buoyancy force for each cell.
   **/
  //  void buoyancy(float[] Fbuoy)
  //  {
  //    float Tamb = 0;
  //    float a = 0.000625f;
  //    float b = 0.025f;
  //
  //    // sum all temperatures
  //    for (int i = 1; i <= nX; i++)
  //    {
  //      for (int j = 1; j <= nY; j++)
  //      {
  //        for (int k = 1; k <= nZ; k++)
  //        {
  //
  //          Tamb += d[I(i, j, k)];
  //        }
  //      }
  //    }
  //
  //    // get average temperature
  //    Tamb /= (nX * nY);
  //
  //    // for each cell compute buoyancy force
  //    for (int i = 1; i <= nX; i++)
  //    {
  //      for (int j = 1; j <= nY; j++)
  //      {
  //        for (int k = 1; k <= nZ; k++)
  //        {
  //          Fbuoy[I(i, j, k)] = a * d[I(i, j, k)] + -b * (d[I(i, j, k)] - Tamb);
  //        }
  //      }
  //    }
  //  }

  /**
   * Calculate the curl at position (i, j, k) in the fluid grid.
   * Physically this represents the vortex strength at the
   * cell. Computed as follows: w = (del x U) where U is the
   * velocity vector at (i, j, k).
   *
   * @param i The x index of the cell.
   * @param j The y index of the cell.
   **/

  float curl(int i, int j, int k)
  {
    float du_dy = (u[I(i, j + 1, k    )] - u[I(i, j - 1, k )]) * 0.50f;
    float dv_dx = (v[I(i + 1, j, k    )] - v[I(i - 1, j, k     )]) * 0.50f;
    //    float dw_dz = (w[I(i,     j, k + 1)] - w[I(i,     j, k - 1 )]) * 0.50f;
    return du_dy - dv_dx;// - dw_dz;
  }

  /**
   * Calculate the vorticity confinement force for each cell
   * in the fluid grid. At a point (i,j), Fvc = N x w where
   * w is the curl at (i,j) and N = del |w| / |del |w||.
   * N is the vector pointing to the vortex center, hence we
   * add force perpendicular to N.
   *
   * @param Fvc_x The array to store the x component of the
   *        vorticity confinement force for each cell.
   * @param Fvc_y The array to store the y component of the
   *        vorticity confinement force for each cell.
   **/

  void vorticityConfinement(float[] Fvc_x, float[] Fvc_y, float[] Fvc_z)
  {
    float dw_dx, dw_dy, dw_dz;
    float leng;
    float v;
    // Calculate magnitude of curl(u,v) for each cell. (|w|)
    for (int i = 1; i <= nX; i++)
    {
      for (int j = 1; j <= nY; j++)
      {
        for (int k = 1; k <= nZ; k++)
        {

          if ( !solid[I(i, j, k)] )
          {
            curl[I(i, j, k)] = Math.abs(curl(i, j, k));
          }
        }
      }
    }
    for (int i = 2; i < nX; i++)
    {
      for (int j = 2; j < nY; j++)
      {
        for (int k = 2; k <= nZ; k++)
        {
          if ( !solid[I(i, j, k)] )
          {
            // Find derivative of the magnitude (n = del |w|)
            dw_dx = (curl[I(i + 1, j, k     )] - curl[I(i - 1, j, k     )] ) * 0.5f;
            dw_dy = (curl[I(i, j + 1, k     )] - curl[I(i, j - 1, k     )] ) * 0.5f;
            dw_dz = (curl[I(i, j, k + 1 )] - curl[I(i, j, k - 1 )] ) * 0.5f;
            // Calculate vector length. (|n|)
            // Add small factor to prevent divide by zeros.
            leng = (float) Math.sqrt(dw_dx * dw_dx + dw_dy * dw_dy + dw_dz * dw_dz)
              + 0.000001f;
            // N = ( n/|n| )
            dw_dx /= leng;
            dw_dy /= leng;
            dw_dz /= leng;
            v = curl(i, j, k);
            // N x w
            Fvc_x[I(i, j, k)] = dw_dy * -v;
            Fvc_y[I(i, j, k)] = dw_dx *  v;
            Fvc_z[I(i, j, k)] = dw_dz *  v;
          }
        }
      }
    }
  }

  /**
   * The basic velocity solving routine as described by Stam.
   **/

  void velocitySolver()
  {

    inletVelocity();

    // add velocity that was input by mouse
    addSource(u, uOld);
    addSource(v, vOld);
    addSource(w, wOld);

    // add in vorticity confinement force
    vorticityConfinement(uOld, vOld, wOld);
    addSource(u, uOld);
    addSource(v, vOld);
    addSource(w, wOld);

    // add in buoyancy force
    //    buoyancy(vOld);
    //    addSource(v, vOld);

    // swapping arrays for economical mem use
    // and calculating diffusion in velocity.
    swapU();
    diffuse(1, u, uOld, visc);

    swapV();
    diffuse(2, v, vOld, visc);

    swapW();
    diffuse(3, w, wOld, visc);

    // we create an incompressible field
    // for more effective advection.
    project(u, v, w, uOld, vOld);

    swapU(); 
    swapV();
    swapW();

    // self advect velocities
    advect(1, u, uOld, uOld, vOld, wOld);
    advect(2, v, vOld, uOld, vOld, wOld);
    advect(3, w, wOld, uOld, vOld, wOld);

    // make an incompressible field
    project(u, v, w, uOld, vOld);

    // clear all input velocities for next frame
    for (int i = 0; i < size; i++) { 
      uOld[i] = 0.0f; 
      vOld[i] = 0.0f;
      wOld[i] = 0.0f;
    }
  }


  /**
   * The basic density solving routine.
   **/

  void densitySolver()
  {

    inletDensity();

    // add density inputted by mouse
    addSource(d, dOld);
    swapD();

    diffuse(0, d, dOld, diff);
    swapD();

    advect(0, d, dOld, u, v, w);

    // clear input density array for next frame
    for (int i = 0; i < size; i++) dOld[i] = 0.0f;
  }


  void addSource(float[] x, float[] x0)
  {
    for (int i = 0; i < size; i++)
    {
      x[i] += dt * x0[i];
    }
  }


  /**
   * Calculate the input array after advection. We start with an
   * input array from the previous timestep and an and output array.
   * For all grid cells we need to calculate for the next timestep,
   * we trace the cell's center position backwards through the
   * velocity field. Then we interpolate from the grid of the previous
   * timestep and assign this value to the current grid cell.
   *
   * @param b Flag specifying how to handle boundries.
   * @param d Array to store the advected field.
   * @param d0 The array to advect.
   * @param du The x component of the velocity field.
   * @param dv The y component of the velocity field.
   **/

  void advect(int b, float[] d, float[] d0, float[] du, float[] dv, float[] dw )
  {
    int i0, j0, i1, j1, k0, k1;
    float x, y, z, s0, t0, s1, t1, u0, u1;

    float dtx = dt * nX ;
    float dty = dt * nY ;
    float dtz = dt * nZ ;

    for (int i = 1; i <= nX; i++)
    {
      for (int j = 1; j <= nY; j++)
      {
        for (int k = 1; k <= nZ; k++)
        {        
          if ( !solid[I(i, j, k)] )
          {
            // go backwards through velocity field
            x = i - dtx    * du[I(i, j, k)];
            y = j - dty    * dv[I(i, j, k)];
            z = k - dtz    * dw[I(i, j, k)];

            // interpolate results
            if (x > nX + 0.5) x = nX + 0.5f;
            if (x < 0.5)     x = 0.5f;

            i0 = (int) x;
            i1 = i0 + 1;

            if (y > nY + 0.5) y = nY + 0.5f;
            if (y < 0.5)     y = 0.5f;

            j0 = (int) y;
            j1 = j0 + 1;

            if (z > nZ + 0.5) z = nZ + 0.5f;
            if (z < 0.5)     z = 0.5f;

            k0 = (int) z;
            k1 = k0 + 1;



            s1 = x - i0;
            s0 = 1.0 - s1;
            t1 = y - j0;
            t0 = 1.0 - t1;
            u1 = z - k0;
            u0 = 1.0 - u1;



            //                       d[I(i, j, k)] = s0 * (t0 * d0[I(i0, j0, k0)] + t1 * d0[I(i0, j1, k0)])
            //                          + s1 * (t0 * d0[I(i1, j0, k0)] + t1 * d0[I(i1, j1, k0)] );


            d[I(i, j, k)] = 

              s0 * ( t0 * (u0 * d0[I(i0, j0, k0)]  +  u1 * d0[I(i0, j0, k1)]) + (t1 * (u0 * d0[I(i0, j1, k0)]  +  u1 * d0[I(i0, j1, k1)]) ) )
              + s1 * ( t0 * (u0 * d0[I(i1, j0, k0)]  +  u1 * d0[I(i1, j0, k1)]) + (t1 * (u0 * d0[I(i1, j1, k0)]  +  u1 * d0[I(i1, j1, k1)]) ) );
          }
        }
      }
    }
    setBoundry(b, d);
  }

  /**
   * Recalculate the input array with diffusion effects.
   * Here we consider a stable method of diffusion by
   * finding the densities, which when diffused backward
   * in time yield the same densities we started with.
   * This is achieved through use of a linear solver to
   * solve the sparse matrix built from this linear system.
   *
   * @param b Flag to specify how boundries should be handled.
   * @param c The array to store the results of the diffusion
   * computation.
   * @param c0 The input array on which we should compute
   * diffusion.
   * @param diff The factor of diffusion.
   **/

  void diffuse(int b, float[] c, float[] c0, float diff)
  {
    float a = dt * diff * nX * nY;
    linearSolver(b, c, c0, a, 1 + 6 * a);
  }

  /**
   * Use project() to make the velocity a mass conserving,
   * incompressible field. Achieved through a Hodge
   * decomposition. First we calculate the divergence field
   * of our velocity using the mean finite differnce approach,
   * and apply the linear solver to compute the Poisson
   * equation and obtain a "height" field. Now we subtract
   * the gradient of this field to obtain our mass conserving
   * velocity field.
   *
   * @param x The array in which the x component of our final
   * velocity field is stored.
   * @param y The array in which the y component of our final
   * velocity field is stored.
   * @param p A temporary array we can use in the computation.
   * @param div Another temporary array we use to hold the
   * velocity divergence field.
   *
   **/

  void project(float[] x, float[] y, float []z, float[] p, float[] div)
  {
    for (int i = 1; i <= nX; i++)
    {
      for (int j = 1; j <= nY; j++)
      {
        for (int k = 1; k <= nZ; k++)
        {       
          if ( !solid[I(i, j, k)] )
          {
            div[I(i, j, k)] = - 0.5f *(   
              x[I(i+1, j, k  )] 
              - x[I(i-1, j, k  )]
              + y[I(i, j+1, k  )] 
              - y[I(i, j-1, k  )]
              + z[I(i, j, k+1)]
              - z[I(i, j, k-1)] 
              )/ nX;

            p[I(i, j, k)] = 0.0f;
          }
        }
      }
    }
    setBoundry(0, div);
    setBoundry(0, p);

    linearSolver(0, p, div, 1, 6);

    for (int i = 1; i <= nX; i++)
    {
      for (int j = 1; j <= nY; j++)
      {
        for (int k = 1; k <= nZ; k++)
        {    
          if ( !solid[I(i, j, k)] )
          {        
            x[I(i, j, k)] -= 0.5f * nX * (p[I(i+1, j, k    )] - p[I(i-1, j, k   )]);
            y[I(i, j, k)] -= 0.5f * nY * (p[I(i, j+1, k    )] - p[I(i, j-1, k   )]);
            z[I(i, j, k)] -= 0.5f * nZ * (p[I(i, j, k+1  )] - p[I(i, j, k-1 )]);
          }
        }
      }
    }
    setBoundry(1, x);
    setBoundry(2, y);
    setBoundry(3, z);
  }

  /**
   * Iterative linear system solver using the Gauss-sidel
   * relaxation technique. Room for much improvement here...
   *
   **/
  void linearSolver(int b, float[] x, float[] x0, float a, float c) {

    if (multiThreaded) {

      int div = nX/threadNum;

      Thread[] threads  = new Thread[threadNum];

      for (int i=0; i<threads.length-1; i++)
      {
        threads[i] = new RowCalc(b, x, x0, a, c, i*div+1, (i+1)*div);
      }
      threads[threads.length-1] = new  RowCalc(b, x, x0, a, c, (threads.length-1)*div+1, nX);
      for (int i=0; i<threads.length; i++)
      {
        threads[i].start();
      }
      for (int i=0; i<threads.length; i++)
      {
        try {
          threads[i].join();
        } 
        catch (InterruptedException e) {
          e.printStackTrace();
        }
      }
    } else {

      solve(b, x, x0, a, c, 1, nX);
    }
  }

  void solve(int b, float[] x, float[] x0, float a, float c, int startRow, 
    int endRow) {
    for (int n = 0; n < 20; n++) {
      for (int i = startRow; i <= endRow; i++) {
        for (int j = 1; j <= nY; j++) {
          for (int k = 1; k <= nZ; k++) {
            if (!solid[I(i, j, k)]) {
              x[I(i, j, k)] = (  x0[I(i, j, k   )] + 
                a * (  x[I(i + 1, j, k   )] + x[I(i - 1, j, k    )] + 
                x[I(i, j + 1, k   )] + x[I(i, j - 1, k    )] + 
                x[I(i, j, k+ 1)] + x[I(i, j, k - 1)]  ) )/ c;
            }
          }
        }
      }
      setBoundry(b, x);
    }
  }

  private class RowCalc extends Thread {
    int b, startRow, endRow, a, c;
    float[] x;
    float[] x0;

    private RowCalc(int b, float[] x, float[] x0, float a, float c, int startRow, int endRow) {
      this.startRow = startRow;
      this.endRow = endRow;
      this.a = (int) a;
      this.c = (int) c;
      this.b = b;
      this.x = x;
      this.x0 = x0;
    }

    public void run() {
      solve(b, x, x0, a, c, startRow, endRow);
    }
  }

  // util array swapping methods
  void swapU() { 
    tmp = u; 
    u = uOld; 
    uOld = tmp;
  }
  void swapV() { 
    tmp = v; 
    v = vOld; 
    vOld = tmp;
  }
  void swapW() { 
    tmp = w; 
    w = wOld; 
    wOld = tmp;
  }
  void swapD() { 
    tmp = d; 
    d = dOld; 
    dOld = tmp;
  }

  // util method for indexing 1d arrays
  int I(int i, int j, int k) { 
    return i + (nX + 2) * j + ((nX+2)*(nY+2))*k;
  }