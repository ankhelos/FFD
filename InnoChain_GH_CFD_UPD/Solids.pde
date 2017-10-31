// CFD Input - Author: Angelos Chronis

int FluidDomainLength;
int FluidDomainWidth;
int FluidDomainHeight;
boolean [][][] sldInd;

PrintWriter out;

void setDomain()
{
  //FluidDomainLength   = 50;
  //FluidDomainWidth    = 100;
  //FluidDomainHeight   = 50;

  //sldInd = new boolean[FluidDomainLength][FluidDomainWidth][FluidDomainHeight];
  //setSolidBoxInd(10,30,0,10,20,40);
  //setSolidBoxInd(30,30,0,10,20,40);

  //for (int n=0; n<20; n++)
  //{
  //  int x = int(random(2, FluidDomainLength-10));
  //  int y = int(random(20, FluidDomainWidth-20));
  //  int z = int(random(0, FluidDomainHeight-30));

  //  int xx = int(random(5, 10));
  //  int yy = int(random(5, 20));
  //  int zz = int(random(5, 30));

  //  setSolidBoxInd(x, y, z, xx, yy, zz);
  //}
  //setSolidBoxInd(20,40,0,10,20,20);
  
   setSolidIndFromFile();
}

void setSolidBoxInd(int x, int y, int z, int l, int w, int h)
{
  for (int i=x; i<x+l; i++)
  {
    for (int j=y; j<y+w; j++)
    {
      for (int k=z; k<z+h; k++)
      {
        sldInd[i][j][k] = true;
      }
    }
  }
}

void setSolidIndFromFile()
{
  String[] lines = loadStrings("CFD_Domain_Stream.txt");
  FluidDomainLength = int(lines[1]);
  FluidDomainWidth  = int(lines[0]);
  FluidDomainHeight = int(lines[2]);
  
  sldInd = new boolean[FluidDomainLength][FluidDomainWidth][FluidDomainHeight];
  
  for (int i=3; i<lines.length; i++)
  {
    String[] pieces = split(lines[i], ';');
    int x = int(pieces[1]);
    int y = int(pieces[0]);
    int z = int(pieces[2]);
    sldInd[x][y][z] = true;
  }
}