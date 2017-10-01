import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hvc import *

#pngDir = None
pngDir = "./png/"
seismicDir = "../../../data/seis/hvc/3d/f3d/"
gxfile = "gx"
s1 = Sampling(220,1.0,0.000)
s2 = Sampling(400,1.0,0.000)
s3 = Sampling(400,1.0,0.000)
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
k3,k2,k1 = 69,168,172; azimuth=240; elevation=20 # for 3D views
fmin,fmax = -5.5,5.5
plotOnly = False


def main(args):
  #display("tpst")
  #displayHorizons()
  #figures()
  #goSlopes()
  #goFlattenWithSlopes()
  goFlattenWithSlopesAndCorrelations()

def goSlopes():
  f = readImage(gxfile)
  f = gain(f)
  writeImage(gxfile,f)
  if not plotOnly:
    p2 = copy(f)
    p3 = copy(f)
    ep = copy(f)
    lsf = LocalSlopeFinder(8.0,2.0)
    lsf.findSlopes(f,p2,p3,ep);
    writeImage("p2",p2)
    writeImage("p3",p3)
    writeImage("ep",ep)
  else:
    p2 = readImage("p2")
    p3 = readImage("p3")
    ep = readImage("ep")
  for g in [p2,p3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
    makeFrame(world)
def goFlattenWithSlopes():
  f = readImage(gxfile)
  if not plotOnly:
    p2 = readImage("p2")
    p3 = readImage("p3")
    ep = readImage("ep")
    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    ep = pow(ep,8.0)
    fl = Flattener3Dw()
    fl.setIterations(0.01,100)
    fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
    g = fm.flatten(f)
    writeImage("ggs",g)
    s = fm.getShiftsS()
    writeImage("gss",s)
    print "s min =",min(s),"max =",max(s),"avg =",sum(s)/n1/n2/n3
  else:
    g = readImage("ggs")
  display3(f)
  display3(g)
  '''
  world = World()
  addImageToWorld(world,f)
  addImageToWorld(world,g)
  #addImage2ToWorld(world,f,s)
  makeFrame(world)
  '''

def goFlattenWithSlopesAndCorrelations():
  f = readImage(gxfile)
  if not plotOnly:
    p2 = readImage("p2")
    p3 = readImage("p3")
    ep = readImage("ep")
    p2 = mul(d1/d2,p2)
    p3 = mul(d1/d3,p3)
    wp = zerofloat(n2,n3)
    for i3 in range(n3):
      for i2 in range(n2):
        wp[i3][i2] = sum(ep[i3][i2])/n1
    fs = zerofloat(n1,n2,n3)
    rgf = RecursiveExponentialFilter(1)
    rgf.apply1(f,fs)
    gc = GlobalCorrelationFinder(-55,55)
    gc.setStrainMax(0.4)
    ks = gc.getTraceIndexes(5,5,100,0.2,wp)
    ts = gc.findCorrelations(ks,fs)
    ep = pow(ep,8.0)
    fl = Flattener3Dw()
    fl.setIterations(0.01,100)
    fm = fl.getMappingsFromSlopesAndCorrelations(s1,s2,s3,0.01,p2,p3,ep,ks,ts)
    g = fm.flatten(f)
    writeImage("gg",g)
    s = fm.getShiftsS()
    writeImage("gs",s)
  else:
    g = readImage("gg")
  display3(f)
  display3(g)
  world = World()
  addImageToWorld(world,f)
  addImageToWorld(world,g)
  #addImage2ToWorld(world,f,s)
  makeFrame(world)

def figures():
  f = readImage("tpst")
  g = readImage("tpsf")
  s = readImage("tpss")
  p2 = readImage("tpp2")
  p3 = readImage("tpp3")
  ep = readImage("tpep")
  display3(f,p3,cmin=-1.1,cmax=1.1,clabel="Inline slope",png="p3")
  display3(f,p2,cmin=-1.1,cmax=1.1,clabel="Crossline slope",png="p2")
  display3(f,ep,cmin=0.0,cmax=1.0,clabel="Planarity",png="ep")
  display3(f,s,cmin=-17,cmax=17,clabel="Time shift (samples)",png="s")
  """
  display3(f,png="f")
  display3(g,png="g")
  global k1
  for k1 in k1s:
    display3(g,png="g"+str(k1))
  global k1
  for k1 in k1s:
    display2(f,png="fh"+str(k1))
  global k1
  for k1 in k1s:
    display2(g,png="gh"+str(k1))
  """

def display(filename):
  f = readImage(filename)
  world = World()
  ipg = addImageToWorld(world,f)
  ipg.setSlices(k1,k2,k3)
  makeFrame(world)

def displayHorizons():
  f = readImage("tpst")
  s = readImage("tpss")
  r = rampfloat(0.0,1.0,0.0,0.0,n1,n2,n3)
  t = add(s1.first,mul(s1.delta,sub(r,s)))
  for h1 in k1s:
    h = slice23(h1,t)
    world = World()
    ipg = addImageToWorld(world,f)
    ipg.setSlices(k1,k2,k3)
    tg = TriangleGroup(True,s3,s2,h)
    tg.setColor(Color.YELLOW)
    world.addChild(tg)
    makeFrame(world)



#############################################################################
# read/write files

def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

def readSlice3(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(50.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s

from org.python.util import PythonObjectInputStream
def readTensors(name):
  """
  Reads tensors from file with specified basename; e.g., "tpet".
  """
  fis = FileInputStream(seismicDir+name+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  fis.close()
  return tensors
def writeTensors(name,tensors):
  """
  Writes tensors to file with specified basename; e.g., "tpet".
  """
  fos = FileOutputStream(seismicDir+name+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  fos.close()

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  return frame
"""
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
"""
 
def display2(s,png=None):
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1RIGHT_X2UP,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.setHInterval(2.0)
  pp.setVInterval(2.0)
  pp.setHLabel("Crossline (km)")
  pp.setVLabel("Inline (km)")
  pv = pp.addPixels(s2,s3,slice23(k1,s))
  pv.setClips(fmin,fmax)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  pf.setSize(926,510)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")
 
def display3(s,c=None,clabel="",cmin=0,cmax=0,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(fmin,fmax)
  if c:
    cb = pp.addColorBar(clabel)
    #cb.setInterval(1.0)
    pp.setColorBarWidthMinimum(140)
    pp.setLineColor(Color.BLACK)
  else:
    pp.setLineColor(Color.YELLOW)
    #cb = pp.addColorBar("Amplitude")
    #cb.setInterval(5.0)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if c:
    pv12 = PixelsView(s1,s2,slice12(k3,c))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,c))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,c))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  #pf.setFontSizeForSlide(1.0,1.0)
  pf.setFontSize(18)
  if c:
    pf.setSize(1036,814)
  else:
    pf.setSize(859,814)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
