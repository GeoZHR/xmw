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

from ad import *
from sso import *
from util import *

pngDir = None
pngDir = "../../../png/ad/fed3/"

seismicDir = "../../../data/seis/ad/fed/3d/"
#seismicDir = "../../../data/seis/beg/jake/subs/"
fxfile = "fx"
gxlfile = "gxl"
gxnfile = "gxn"
scnfile = "scn"
scifile = "sci"
ssnfile = "ssn"
gxsfile = "gxs"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 180,880,500
#n1,n2,n3 = 426,800,830
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = True

def main(args):
  #goLocalSmoothingFilter()
  #goLinearDiffusion()
  #goStratigraphyOrientedDiffusion()
  goNonlinearDiffusion()
  #goDiffusivity()
  #goFaultChannelSmooth()
  #goSemblance()
  #goSemblanceHale()
  #goCovariance()
  #goFastCovariance()
  #goVariance()
  #goShapeSemblance()
def goNormalPlanar():
  fx = readImage(fxfile)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  el = zerofloat(n1,n2,n3)
  sig1,sig2=4,2
  lof = LocalOrientFilter(sig1,sig2)
  lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
  lof.applyForInlineLinear(fx,u1,u2,u3,el)
  plot3(ep,cmin=0.4,cmax=1.0)
  plot3(el,cmin=0.0,cmax=0.9)

def goLinearDiffusion():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,4
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    ets.setEigenvalues(0.0001,1.0,1.0)
    sig = 10
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,ets,fx)
    writeImage(gxlfile,gx)
  else:
    gx = readImage(gxlfile)
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5,png="fgl")
  plot3(gx,cmin=-1.0,cmax=1.0,png="gxl")
  plot3(fx,cmin=-1.0,cmax=1.0,png="fx")

def goLocalSmoothingFilter():
  fx = readImage(fxfile)
  sig1,sig2=4,2
  lof = LocalOrientFilter(sig1,sig2)
  ets = lof.applyForTensors(fx)
  ets.setEigenvalues(0.0001,1.0,1.0)
  sig = 12.5
  gx = zerofloat(n1,n2,n3)
  lsf = LocalSmoothingFilter()
  lsf.apply(ets,sig,fx,gx)
  plot3(sub(fx,gx),cmin=-1.0,cmax=1.0)
  plot3(gx)
  plot3(fx)

def goStratigraphyOrientedDiffusion():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,4
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    ets.setEigenvalues(0.0001,0.0001,1.0000)
    sig = 8
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,ets,fx)
    writeImage(gxsfile,gx)
  else:
    gx = readImage(gxsfile)
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5,png="fgs")
  plot3(gx,cmin=-1.0,cmax=1.0,png="gxs")
  plot3(fx,cmin=-1.0,cmax=1.0)

def goNonlinearDiffusion():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,4
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    sig = 10
    cycle,limit=3,0.5
    lbd = 0.1
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    sc, gx = fed.apply(sig,lbd,1.0,ets,fx)
    #writeImage(gxnfile,gx)
    writeImage(scnfile,sc)
  else:
    gx = readImage(gxnfile)
    sc = readImage(scnfile)
  sc = pow(sc,2.0)
  sc = sub(sc,min(sc))
  sc = div(sc,max(sc))
  #plot3(sub(fx,gx),cmin=-0.5,cmax=0.5,clab="Amplitude",png="fgn")
  plot3(sc,cmin=0.0,cmax=1.0,clab="Diffusivity",png="scn")
  #plot3(gx,cmin=-1,cmax=1,clab="Amplitude",png="gxn59")
  #plot3(fx,cmin=-1,cmax=1)

def goDiffusivity():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,4
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    lbd = 0.1
    fed = FastExplicitDiffusion()
    sci = fed.applyForWeights(lbd,ets,fx)
    writeImage(scifile,sci)
  else:
    sci = readImage(scifile)
  plot3(sci,cmin=0.0,cmax=1.0,clab="1-diffusivity",png="sci")

def goDiffusion(sig,et,fx):
  cycle,limit=3,0.5
  fed = FastExplicitDiffusion()
  fed.setCycles(cycle,limit)
  return fed.apply(sig,et,fx)

def goDiffusionX(sig,et,fx):
  lbd = 0.1
  cycle,limit=3,0.5
  fed = FastExplicitDiffusion()
  fed.setCycles(cycle,limit)
  return fed.apply(sig,lbd,et,fx)


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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          png=None):
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
    cbar.setWidthMinimum(120)
  ipg.setSlices(59,760,450)
  ipg.setSlices(91,760,450)
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  view = sf.getOrbitView()
  #zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.72)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(238.0)
  view.setElevation(45)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.05,-0.2,0.06))
  sf.viewCanvas.setBackground(sf.getBackground())

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
