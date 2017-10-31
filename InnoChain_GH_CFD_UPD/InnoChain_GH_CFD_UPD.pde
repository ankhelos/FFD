/*
 * CFD SIMULATION MULTITHREADED FRAMEWORK
 * Author: Angelos Chronis
 */

//___________________  Input Parameters : __________________________//

// Copy the input filename here:
String inputFile = "CFD_Domain_Stream.txt";

// Specify the wind direction - Use 0 for East and 1 for West  , 2 for North and 3 for South
int    windDir   = 2;


boolean running = true;
boolean multiThreaded = true;  

int threadNum = 8;

//_________________________________________________________________8240 Converged!//


import processing.video.*;
import processing.opengl.*;
import processing.dxf.*;
import voxellib.Mesh;

FluidSolver fs;
boolean exported;
boolean drawSectionX    = false;
boolean drawSectionY    = false;
boolean drawSectionZ    = false;
boolean drawSolids      = true;
boolean draw3D          = false;
boolean drawStreams     = false;
boolean drawIso         = false;
boolean streamsRunning  = false;
boolean export          = false;

int sectionXInd = 0;
int sectionYInd = 0;
int sectionZInd = 0;
boolean drawSteady = true;  

int   CellSize = 3;
float velscale = 1.3;
int   Scale    = 4;
PFont font;

void setup()
{
  size(1600, 1000, P3D);

  setDomain();
  fs = new FluidSolver(this, FluidDomainLength, FluidDomainWidth, FluidDomainHeight, CellSize);
  fs.reset();
  fs.createVelocityProfile(FluidDomainHeight);
  fs.setSolids();
  
  sectionXInd = FluidDomainLength/2;
  sectionYInd = FluidDomainWidth/2;
  sectionZInd = FluidDomainHeight/3;

  font = createFont("Arial", 14);
  println("FluidDomainLength:" +FluidDomainLength+ ",  FluidDomainWidth:"+FluidDomainWidth+ ", FluidDomainHeight:" +FluidDomainHeight);
  setup_navigation();
  //int streamStep = 3;
  //fs.addStreamLinesXZ(streamStep, streamStep, FluidDomainWidth*CellSize, FluidDomainHeight*CellSize, int(FluidDomainWidth/streamStep), int(FluidDomainHeight/streamStep) );
}


void draw()
{
  //if (export) beginRaw(DXF, "CFD_Output_"+frameCount);
  //println(frameCount);
  //println(millis());
  background(0);
  //___________________________FLUID__________________________________//
  if (running)
  {
    fs.densitySolver();
    fs.velocitySolver();
  }
  //if (drawSectionX)
  //{
  //  pushMatrix();
  //  translate(100, 30);
  //  scale(Scale);
  //  fs.drawSectionX(velscale, 0.01, sectionYInd, drawSteady);
  //  scale(1.0/Scale );
  //  translate(-80, 0);
  //  fs.drawLegend(velscale);
  //  popMatrix();
  //}
  //if (drawSectionY)
  //{
  //  pushMatrix();
  //  translate(100, 30);
  //  scale(Scale);
  //  fs.drawSectionY(velscale, 0.01, sectionYInd, drawSteady);
  //  scale(1.0/Scale );
  //  translate(-80, 0);
  //  fs.drawLegend(velscale);
  //  popMatrix();
  //}
  
  //if (drawSectionZ)
  //{
  //  pushMatrix();
  //  translate(100, 30);
  //  scale(Scale);
  //  fs.drawSectionZ(velscale, 0.01, sectionZInd, drawSteady);
  //  scale(1.0/Scale);
  //  translate(-80, 0);
  //  fs.drawLegend(velscale);
  //  popMatrix();
  //}
  //if (draw3D)
  //{
  //  lights();
  //  pushMatrix();
  //  navigate();
  //  fs.draw3D(velscale, 0.1, 0.001);
  //  fs.drawBox();
  //  popMatrix();
  //  pushMatrix();
  //  translate(20, 30);
  //  fs.drawLegend(velscale);
  //  popMatrix();
  //}
  
  //if (draw3D)
  {
    lights();
    pushMatrix();
    navigate();
    if (draw3D) fs.draw3D(velscale, 0.1, 0.001);
    if (drawSectionX) fs.drawSectionX(velscale, 0.013, sectionXInd, drawSteady);
    if (drawSectionY) fs.drawSectionY(velscale, 0.013, sectionYInd, drawSteady);
    if (drawSectionZ) fs.drawSectionZ(velscale, 0.013, sectionZInd, drawSteady);
    if (drawSolids)   fs.drawSolids();
    if (drawIso)      fs.drawIso(velscale);
    
    fs.drawBox();
    popMatrix();
    pushMatrix();
    translate(20, 30);
    fs.drawLegend(velscale);
    popMatrix();
  }
  
  //if (drawStreams)
  //{
  //  lights();
  //  pushMatrix();
  //  if (!export) navigate();
  //  fs.drawBox();
  //  fs.drawSolids();
  //  if (streamsRunning) fs.stepStreamLines();
  //  fs.drawStreamLines();
  //  popMatrix();
  //}
  if (frameCount%1==0 && !exported && running)
  {
    fs.calculateMeans();
    fs.densityConvergence(0.001);
    //if (fs.converged ) 
    //{
    //println("export");
    //fs.writeData(); 
    //exported = true;
    //exit();
    //}
  }
  //if (export) { endRaw(); export = false;}
  //displayInfo();
}

void displayInfo()
{
  pushStyle();
  fill(255);
  textFont(font, 20);
  text("Av.dens.  = " + nf(fs.avdens, 1, 5), 100, height-30);
  text("Av.veloc. = " + nf(fs.avvelX, 1, 3), 100, height-50);
}