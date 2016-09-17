import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

from ssr import *
from util import *

pngDir = None
pngDir = "../../../png/ssr/semblance/"

seismicDir = "../../../data/seis/ssr/semblance/"
#seismicDir = "../../../data/seis/beg/jake/subs/"
fxfile = "gx"
snfile = "sn"
sdfile = "sd"
smfile = "sm"
epfile = "ep"
p2file = "p2"
p3file = "p3"
smhfile = "smh"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 140,880,500
#n1,n2,n3 = 426,800,830
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = False

def main(args):
  #goPlanarity()
  #goSlopes()
  goShapeSemblance()
  #goSemblanceHale()
def goPlanarity():
  sig1,sig2=4,2
  fx = readImage(fxfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(sig1,sig2,sig2)
    lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
    writeImage(epfile,ep)
  else:
    ep = readImage(epfile)
  plot3(fx)
  plot3(ep,cmin=0.2,cmax=1.0,png="ep")

def goSlopes():
  sig1,sig2=4,2
  fx = readImage(fxfile)
  if not plotOnly:
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig2,5)
    lsf.findSlopes(fx,p2,p3,ep)
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    #writeImage(epfile,ep)
  else:
    ep = readImage(epfile)
  plot3(fx)
  plot3(ep,cmin=0.2,cmax=1.0,png="ep")
  plot3(fx,g=p2,cmin=-0.5,cmax=0.5,cmap=jetFill(0.6))

def goShapeSemblance():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=4,2
    sn = readImage(snfile)
    sd = readImage(sdfile)
    ep = readImage(epfile)
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    sm = Semblance()
    #sd = mul(fx,fx)
    #sn = sm.smoothVW(2,et,fx)
    #sd = sm.smoothVW(2,et,sd)
    #sn = mul(sn,sn)
    #p2 = readImage(p2file)
    #p3 = readImage(p3file)
    #sn,sd=sm.applyForSemblanceNumDen(p2,p3,fx)
    #writeImage(snfile,sn)
    #writeImage(sdfile,sd)
    sn = readImage(snfile)
    sd = readImage(sdfile)
    wp = sub(1,ep)
    et.setEigenvalues(0.2,0.0001,1.0000)
    sem = sm.shapeSemblance(et,wp,sn,sd)
    writeImage(smfile,sem)
  else:
    sem = readImage(smfile)
  plot3(fx,png="seis")
  plot3(sem,cmin=0.2,cmax=1.0,png="semShp")

def goSemblanceHale():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,4
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    lsf = LocalSemblanceFilter(2,2)
    sem = lsf.semblance(LocalSemblanceFilter.Direction3.VW,et,fx)
    writeImage(smhfile,sem)
  else:
    sem = readImage(smhfile)
  plot3(sem,cmin=0.2,cmax=1.0,png="semh")


def normalize(ss):
  sub(ss,min(ss),ss)
  div(ss,max(ss),ss)
  
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

def readImage(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeImage(basename,image):
  """ 
  Writes an image to a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# graphics

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
def hueFill(alpha):
  return ColorMap.getHue(0.0,1.0,alpha)
def hueFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.getHue(0.0,1.0),a)

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def plot3(f,g=None,k1=51,
    cmin=None,cmax=None,cmap=None,clab=None,cint=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
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
      ipg.setClips(-1.5,1.2)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1.5,1.2)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(85)
  ipg.setSlices(k1,857,450)
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  view = sf.getOrbitView()
  zscale = 0.35*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.8)
  view.setAzimuth(235.0)
  view.setElevation(38)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.05,-0.12,0.06))
  sf.viewCanvas.setBackground(sf.getBackground())
  #sf.setSize(850,700)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

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
