StreamLine[][] streams;

class StreamLine
{
  PVector[] keyPts;
  PVector startPt;
  PVector velocity;

  StreamLine(float x, float y, float z)
  {
    startPt = new PVector(x, y, z);
    velocity = new PVector(0, 0, 0);
    keyPts = new PVector[0];

    keyPts = (PVector[]) append(keyPts, startPt);
  }

  void step()
  {
    PVector previousPt;
    previousPt = keyPts[keyPts.length-1];

    PVector nextPt = new PVector(previousPt.x, previousPt.y, previousPt.z);

    int i = int(previousPt.x/cell);
    int j = int(previousPt.y/cell);
    int k = int(previousPt.z/cell);

    PVector fluidVelocity = new PVector();

    if (i>0 && i<nX && j>0 && j<nY && k>0 && k<nZ)
    {
      fluidVelocity.set(u[I(i, j, k)], v[I(i, j, k)], w[I(i, j, k)]);

      fluidVelocity.normalize();
      fluidVelocity.mult(cell);
      velocity.add(fluidVelocity);
      velocity.normalize();
      velocity.mult(cell);
      
      //velocity.div(2.0);
      //velocity.mult(cell);
      
      //nextPt.add(fluidVelocity);
      nextPt.add(velocity);
      
      keyPts = (PVector[]) append(keyPts, nextPt);
    }
  }

  void draw()
  {
    pushStyle();
    stroke(255);
    strokeWeight(1);
    beginShape();
    curveVertex(startPt.x, startPt.y, startPt.z);
    for (int i=0; i<keyPts.length; i++)
    {
      curveVertex(keyPts[i].x, keyPts[i].y, keyPts[i].z);
    }
    curveVertex(keyPts[keyPts.length-1].x, keyPts[keyPts.length-1].y, keyPts[keyPts.length-1].z);
    endShape();
    popStyle();
  }
}

void addStreamLinesXZ(float sx, float sz, float l, float h, int numX, int numZ)
{
  streams = new StreamLine[numX][numZ];
  for (int i=0; i<numX; i++)
  {
    for (int j=0; j<numZ; j++)
    {
      float x = map(i, 0, numX, sx, l);
      float z = map(j, 0, numZ, sz, h);

      streams[i][j] = new StreamLine(x, 20.0, z);
      //println(x+", "+z);
    }
  }
}

void drawStreamLines()
{
  for (int i=0; i<streams.length; i++)
  {
    for (int j=0; j<streams[i].length; j++)
    {
      streams[i][j].draw();
    }
  }
}

void stepStreamLines()
{
  //streams[5][9].step();
  for (int i=0; i<streams.length; i++)
  {
    for (int j=0; j<streams[i].length; j++)
    {
      streams[i][j].step();
    }
  }
}