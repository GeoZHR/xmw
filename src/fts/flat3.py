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

from dnp import *

#pngDir = None
pngDir = "../../../png/fts/"
seismicDir = "../../../data/seis/tpd/"
s1 = Sampling(109,0.004,0.620)
s2 = Sampling(109,0.025,0.000)
s3 = Sampling(109,0.025,1.150)

n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta

k1,k2,k3 = 99,90,90; azimuth=240; elevation=20 # for 3D views
fmin,fmax = -5.5,5.5

fxfile = "fx"
gtfile = "gt"
p2file = "p2"
p3file = "p3"
epfile = "ep"

def main(args):
  goSlope()
  goFlatten()

def goSlope():
  fx = readImage(fxfile)
  p2 = copy(fx)
  p3 = copy(fx)
  ep = copy(fx)
  lsf = LocalSlopeFinder(8.0,2.0)
  lsf.findSlopes(fx,p2,p3,ep);
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)

def goFlatten():
  findShifts = True
  fx = readImage("fx")
  p2 = readImage("p2")
  p3 = readImage("p3")
  ep = readImage("ep")
  p2 = mul(d1/d2,p2)
  p3 = mul(d1/d3,p3)
  ep = pow(ep,6.0)
  fl = Flattener3()
  fl.setIterations(0.01,100)
  fm = fl.getMappingsFromSlopes(s1,s2,s3,p2,p3,ep)
  gt = fm.flatten(fx)
  writeImage(gtfile,gt)
  gt = readImage(gtfile)
  plot3(fx,g=ep,clab="Planarity",cint=0.2)
  plot3(gt,clab="Amplitude",cint=2)
  plot3s(fx,clabel="Amplitude",png="fx")
  plot3s(gt,clabel="Amplitude",png="gt")

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

#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)

def addColorBar(frame, clab=None, cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial", Font.PLAIN, 32))  # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar, BorderLayout.EAST)
  return cbar

def plot3(f,g=None,cmin=None,cmax=None,
          cmap=None,clab=None,cint=None,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(fmin,fmax)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(fmin,fmax)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.4)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(837,700)
  else:
    sf.setSize(700,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)

  ov = sf.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.2)
  ov.setAzimuth(azimuth)
  ov.setElevation(elevation)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot3s(s,c=None,clabel="",cmin=0,cmax=0,png=None):
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
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(2.0)
  pp.setInterval1(0.2)
  pp.setInterval2(1.0)
  pp.setInterval3(1.0)
  pp.mosaic.setHeightElastic(1,100)
  #pp.mosaic.setHeightElastic(1,200)
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
  #pf.setFontSizeForSlide(1.0,0.95,16.0/9.0)
  pf.setFontSizeForSlide(0.7,0.8,16.0/9.0)
  pf.setSize(960,800)
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
