// CFD Visualization - Author: Angelos Chronis

void drawSectionZ(float scalarv, float scalard, int zind, boolean stbl)
{  
  //stroke(180, 150);
  //noFill();
  //rect(0, 0, nX*cell, nY*cell);
  pushStyle();
  colorMode(HSB, 360);
  sectionZ.loadPixels();
  for (int i = 0; i < sectionZ.pixels.length; i++)
  {
    sectionZ.pixels[i] = color(0);
  }
  for (int i = 1; i <= nX; i++)
  {
    //float dx = (int)( (i - 0.5f) * cell );
    for (int j = 1; j <= nY; j++)
    {
      //float dy = (int)( (j - 0.5f) * cell );
      if (!solid[I(i, j, zind)])
      {
        float dc=0.0;
        if (!stbl) dc = map (   d[I(i, j, zind)], 0, scalard, 0, 360 )  ;
        if (stbl)  dc = map (  dm[I(i, j, zind)], 0, scalard, 0, 360 )  ;
        float vc=0.0;
        if (!stbl) vc = map ( ( abs( u[I(i, j, zind)]) + abs( v[I(i, j, zind)])) /2.0, 0, scalarv, 250, 0)  ;
        if (stbl)  vc = map ( ( abs(um[I(i, j, zind)]) + abs(vm[I(i, j, zind)])) /2.0, 0, scalarv, 250, 0)  ;

        if (dc < 0) dc = 0;
        vc = constrain(vc, 0, 360);
        //dc = constrain(vc, 0, 360);
        //fill(vc, 360, 360, dc);
        //fill(vc, 360, 360);
        //noStroke();
        //rect(dx-cell2, dy-cell2, cell, cell);
        sectionZ.pixels[i+j*sectionZ.width] = color(vc, 360, 360);
      } else
      {
        //fill(360, 0, 200);
        //noStroke();
        //rect(dx-cell2, dy-cell2, cell, cell);
        sectionZ.pixels[i+j*sectionZ.width] = color(250);
      }
      //if (solid[I(i, j, zind)])
      //{     
      //  int col = boundInt[I(i, j, zind)]%100*13 + boundInt[I(i, j, zind)]%10*13 + int(boundInt[I(i, j, zind)]/100)*13;
      //  stroke(col,360, 360);
      //  point(dx, dy);
      //}
    }
  }
  pushMatrix();
  translate(0,0,zind*CellSize);
  sectionZ.updatePixels();
  scale(CellSize);
  image(sectionZ, 0, 0);
  popMatrix();
  popStyle();
}

void drawSectionY(float scalarv, float scalard, int jind, boolean stbl)
{  
  //stroke(180, 150);
  //noFill();
  //rect(0, 0, nX*cell, nZ*cell);
  pushStyle();
  colorMode(HSB, 360); 
  sectionY.loadPixels(); 
  for (int i = 1; i <= nX; i++)
  {
    // x position of current cell
    //float dx = (int)( (i - 0.5f) * cell );
    for (int k = 1; k <= nZ; k++)
    {
      // z position of current cell
      //float dz = (int)( (nZ - k + 0.5f) * cell );
      if (!solid[I(i, jind, k)])
      {

        float dc = 0.0;
        if (!stbl) dc = map (  d[I(i, jind, k)], 0, scalard, 0, 360 )  ;
        if (stbl)  dc = map (  dm[I(i, jind, k)], 0, scalard, 0, 360 )  ;

        float vc  = 0.0;
        if (!stbl) vc  = map ( ( abs(u[I(i, jind, k)]) + abs(v[I(i, jind, k)]) ) /2.0, 0, scalarv, 250, 0)  ;
        if (stbl)  vc  = map ( ( abs(um[I(i, jind, k)]) + abs(vm[I(i, jind, k)]) ) /2.0, 0, scalarv, 250, 0)  ;

        if (dc < 0) dc = 0;
        vc = constrain(vc, 0, 360);
        //dc = constrain(vc, 0, 360);
        //fill(vc, 360, 360, dc);
        //noStroke();
        //rect(dx-cell2, dz-cell2, cell, cell);
        sectionY.pixels[i+k*sectionZ.width] = color(vc, 360, 360);
      } else if (solid[I(i, jind, k)])
      {
        //fill(360, 0, 200);
        //noStroke();
        //rect(dx-cell2, dy-cell2, cell, cell);
        sectionY.pixels[i+k*sectionZ.width] = color(250);
      }
      //if (solid[I(i, jind, k)])
      //{          
      //  stroke(360);
      //  point(dx, dz);
      //}
    }
  }
  pushMatrix();
  translate(0,jind*CellSize,0);
  sectionY.updatePixels();
  scale(CellSize);
  rotateX(HALF_PI);
  image(sectionY, 0, 0);
  popMatrix();
  popStyle();
}

void drawSectionX(float scalarv, float scalard, int iind, boolean stbl)
{  
  //stroke(180, 150);
  //noFill();
  //rect(0, 0, nX*cell, nZ*cell);
  pushStyle();
  colorMode(HSB, 360); 
  sectionX.loadPixels(); 
  for (int j = 1; j <= nY; j++)
  {
    // x position of current cell
    //float dx = (int)( (i - 0.5f) * cell );
    for (int k = 1; k <= nZ; k++)
    {
      // z position of current cell
      //float dz = (int)( (nZ - k + 0.5f) * cell );
      if (!solid[I(iind, j, k)])
      {

        float dc = 0.0;
        if (!stbl) dc = map (  d[I(iind, j, k)], 0, scalard, 0, 360 )  ;
        if (stbl)  dc = map (  dm[I(iind, j, k)], 0, scalard, 0, 360 )  ;

        float vc  = 0.0;
        if (!stbl) vc  = map ( ( abs(u[I(iind, j, k)]) + abs(v[I(iind, j, k)]) ) /2.0, 0, scalarv, 250, 0)  ;
        if (stbl)  vc  = map ( ( abs(um[I(iind, j, k)]) + abs(vm[I(iind, j, k)]) ) /2.0, 0, scalarv, 250, 0)  ;

        if (dc < 0) dc = 0;
        vc = constrain(vc, 0, 360);
        //dc = constrain(vc, 0, 360);
        //fill(vc, 360, 360, dc);
        //noStroke();
        //rect(dx-cell2, dz-cell2, cell, cell);
        sectionX.pixels[j+k*sectionX.width] = color(vc, 360, 360);
      } else if (solid[I(iind, j, k)])
      {
        //fill(360, 0, 200);
        //noStroke();
        //rect(dx-cell2, dy-cell2, cell, cell);
        sectionX.pixels[j+k*sectionX.width] = color(250);
      }
      //if (solid[I(i, jind, k)])
      //{          
      //  stroke(360);
      //  point(dx, dz);
      //}
    }
  }
  pushMatrix();
  translate(iind*CellSize,0,0);
  sectionX.updatePixels();
  scale(CellSize);
  rotateX(HALF_PI);
  rotateY(HALF_PI);
  image(sectionX, 0, 0);
  popMatrix();
  popStyle();
}

void drawIso(float scalarv)
{ 
  float[][][] isoVal = new float[nX+1][nY+1][nZ+1];
  for (int i = 0; i <= nX; i++)
  {
    for (int j = 0; j <= nX; j++)
    {
      for (int k = 0; k <= nZ; k++)
      {
        isoVal[i][j][k] = scalarv-um[I(i,j,k)]; //+ v[I(i,j,k)] + w[I(i,j,k)]
      }
    } 
  }
  mesh.makeFromVoxels(isoVal, cell, scalarv*1.15);
  pushStyle();
  fill(255);
  noStroke();
  mesh.draw();
  popStyle();
}

void draw3D(float scalarv, float scalard, float t)
{  
  pushStyle();
  colorMode(HSB, 360);  
  for (int i = 1; i <= nX; i++)
  {
    // x position of current cell
    float dx = (int)( (i - 0.5f) * cell );
    for (int j = 1; j <= nY; j++)
    {
      // y position of current cell
      float dy = (int)( (j - 0.5f) * cell );
      for (int k = 1; k <= nZ; k++)
      {
        float dz = (int)( (k - 0.5f) * cell );          
        if (!solid[I(i, j, k)] && d[I(i, j, k)] > t)
        {
          float dc = map (  d[I(i, j, k)], 0, scalard, 0, 360 )  ;
          float vc = map ( ( abs(u[I(i, j, k)]) + abs(v[I(i, j, k)])  + abs(w[I(i, j, k)]) ) /3.0, 0, scalarv, 250, 0)  ;
          if (dc < 0) dc = 0;
          fill(vc, 360, 360, dc);
          noStroke();
          pushMatrix();
          translate(dx, dy, dz);          
          box(cell);
          popMatrix();
        }
        //else if (solid[I(i, j, k)])
        //{
        //  fill(250);
        //  noStroke();
        //  //          stroke(180);
        //  //          strokeWeight(3);
        //  pushMatrix();
        //  translate(dx, dy, dz);          
        //  //point(dx, dy, dz);
        //  box(cell);
        //  popMatrix();
        //}
      }
    }
  }
  popStyle();
}

void drawSolids()
{  
  pushStyle();
  colorMode(HSB, 360);  
  for (int i = 1; i <= nX; i++)
  {
    // x position of current cell
    float dx = (int)( (i - 0.5f) * cell );
    for (int j = 1; j <= nY; j++)
    {
      // y position of current cell
      float dy = (int)( (j - 0.5f) * cell );
      for (int k = 1; k <= nZ; k++)
      {
        float dz = (int)( (k - 0.5f) * cell );          
        if (solid[I(i, j, k)])
        {
          fill(250);
          noStroke();
          //          stroke(180);
          //          strokeWeight(3);
          pushMatrix();
          translate(dx, dy, dz);          
          //point(dx, dy, dz);
          box(cell);
          popMatrix();
        }
      }
    }
  }
  popStyle();
}

void drawBox()
{
  pushMatrix();
  noFill();
  stroke(180);
  strokeWeight(0.2);
  translate(nX*cell2, nY*cell2, nZ*cell2);
  box(nX*cell, nY*cell, nZ*cell);
  popMatrix();
}

void drawLegend(float scalar)
{
  pushStyle();
  colorMode(HSB, 360); 
  noStroke();
  for (int i=0; i<16; i++)
  {
    float velf = map(i, 0, 16, 0, 250);
    fill(velf, 720, 720);
    rect(0, i*height/20, 20, height/20 );
    textFont(font, 12);
    float velt = map(i, 0, 16, scalar, 0);
    fill(360);
    text(nf(velt, 1, 2), 30, 30+i*height/20);
  }
  popStyle();
}
}