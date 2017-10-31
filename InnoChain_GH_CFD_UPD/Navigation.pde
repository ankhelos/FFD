// Navigation Functions - Author: Angelos Chronis  based on Przemek Jaworski

float alpha1 = -0.35;
float alpha2 = 0.51;
float distance = -150;
float translateX = 0;
float translateY = 0;
float translateZ = 0;

void setup_navigation()
{
  translateX =-121.52772; 
  translateY=-171.83762; 
  translateZ=83.08213; 
  alpha1 = 3.719959; 
  alpha2 = 0.800009; 
  distance = -114.39174;
}

void navigate()
{
  translate(width/2, height/2);
  translate(0, 0, distance);
  scale(3);
  rotateX(alpha2);
  rotateZ(alpha1);  
  translate(translateX, translateY, translateZ);
  drawCS();
}

void drawCS()
{
  pushStyle();
  colorMode(RGB, 255);
  strokeWeight(1);
  stroke(255, 0, 0);
  line(0, 0, 0, 10, 0, 0);
  stroke(0, 255, 0);
  line(0, 0, 0, 0, 10, 0);    
  stroke(0, 0, 255);
  line(0, 0, 0, 0, 0, 10); 
  popStyle();
}

void mouseDragged()
{
  if (mouseButton == LEFT)
  {
    alpha1 += (pmouseX - mouseX)/100.0;
    alpha2 += (pmouseY - mouseY)/100.0;
  }
  if (mouseButton == RIGHT)
  {
    float msx = (pmouseX - mouseX)*0.4;
    float msy = (pmouseY - mouseY)*0.4;

    translateX -= msx*cos(-alpha1);
    translateY -= msx*sin(-alpha1);

    translateX -= msy* -cos(-alpha2)*sin(-alpha1);
    translateY -= msy*cos(-alpha2)*cos(-alpha2);
    translateZ -= msy*sin(-alpha2);
  }
  if (mouseButton == CENTER)
  {
    float msz = (pmouseX -mouseX)+(pmouseY - mouseY);
    distance += msz*0.4;
  }
}
void keyPressed()
{
  if (key== 'f') { 
    saveFrame();
  }
  if (key== 'p')
  {
    println("translateX ="+translateX+"; translateY="+translateY +  "; translateZ="+translateZ+
      "; alpha1 = "+alpha1+"; alpha2 = "+alpha2+"; distance = "+distance+";");
  }
  if (key== 'r') {
    fs.reset();
    setDomain();
    fs.createVelocityProfile(FluidDomainHeight);
    fs.setSolids();
    println("Resetting...");
    fs.addStreamLinesXZ(50.0, 20.0, 450.0, 450.0, 10, 10);
  }
  if (key== 'l') {
    export = !export;
  }
  if (key== 'g') {
    running = !running;
  }
  //if (key== 'z') {
  //  streamsRunning = !streamsRunning;
  //}

  if (key == CODED) {
    if (keyCode == ESC) 
    {
      //       mm.finish();
      println("export");
      fs.writeData(); 
      out.close();
    }
  }

  if (key== '1') { 
    drawSectionX = true;
    drawSectionY = false;
    drawSectionZ = false;
    draw3D       = false;
    drawStreams  = false;
  }
  if (key== '2') { 
    drawSectionX = false;
    drawSectionY = true;
    drawSectionZ = false;
    draw3D       = false;
    drawStreams  = false;
  }  
  if (key== '3') { 
    drawSectionX = false;
    drawSectionY = false;
    drawSectionZ = true;
    draw3D       = false;
    drawStreams  = false;
  }
  if (key== '4') { 
    draw3D       = !draw3D;
  }
  if (key== '5') { 
    drawSolids   = !drawSolids;
  }
  if (key== '6') { 
    drawStreams  = false;
  }
  if (key== '7') { 
    drawIso   = !drawIso;
  }
  if (key== 'w')
  {
    if (sectionZInd<FluidDomainHeight) sectionZInd++;
    println("sectionZInd :"+sectionZInd);
  }
  if (key== 'q')
  {
    if (sectionZInd>0)                 sectionZInd--;
    println("sectionZInd :"+sectionZInd);
  }
  if (key== 's')
  {
    if (sectionYInd<FluidDomainWidth) sectionYInd++;
    println("sectionYInd :"+sectionYInd);
  }
  if (key== 'a')
  {
    if (sectionYInd>0)                 sectionYInd--;
    println("sectionYInd :"+sectionYInd);
  }
  if (key== 'x')
  {
    if (sectionXInd<FluidDomainLength) sectionXInd++;
    println("sectionXInd :"+sectionXInd);
  }
  if (key== 'z')
  {
    if (sectionXInd>0)                 sectionXInd--;
    println("sectionXInd :"+sectionXInd);
  }
  if (key== '`') { 
    drawSteady = !drawSteady;
  }
}