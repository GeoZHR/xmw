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

pngDir = "../../../png/sso/3d/"
pngDir = None

seismicDir = "../../../data/seis/sso/3d/real/"
seismicDir = "../../../data/seis/gbc/dat/"
#seismicDir = "../../../data/seis/beg/jake/subs/"
fxfile = "fx"
fxfile = "pp"
ellfile = "ell"
elsfile = "els"
eplfile = "epl"
epsfile = "eps"
etlfile = "etl"
etsfile = "ets"
gxlfile = "gxl"
gxsfile = "gxs"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 240,880,500
n1,n2,n3 = 2000,150,145
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = False

def main(args):
  gx = readImage(fxfile)
  print min(gx)
  plot3(gx,cmin=-1,cmax=1)
  #goLof()
  #goLoe()
  #goSmoothL()
  #goSmoothS()
  #goSemblanceHale()
def goLof():
  fx = readImage(fxfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    sig1,sig2=8,2
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
    writeImage(eplfile,ep)
    writeTensors(etlfile,ets)
  else:
    el = readImage(ellfile)
    ep = readImage(eplfile)
  plot3(fx)
  ep = pow(ep,2)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  plot3(ep,cmin=0.2,cmax=1.0)

def goLoe():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    sig1,sig2=8,2
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    loe = LocalOrientEstimator(et,20)
    loe.setEigenvalues(1.0,0.01,0.8)
    loe.applyForNormalPlanar(fx,u1,u2,u3,ep)
    writeImage(epsfile,ep)
    writeTensors(etsfile,et)
  else:
    el = readImage(elsfile)
    ep = readImage(epsfile)
    #et = readTensors(etsfile)
  plot3(fx)
  ep = pow(ep,2)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  plot3(ep,cmin=0.2,cmax=1.0)


def goSmoothL():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = readImage(eplfile)
    et = readTensors(etlfile)
    '''
    eu = zerofloat(n1,n2,n3)
    ev = zerofloat(n1,n2,n3)
    ew = zerofloat(n1,n2,n3)
    et.getEigenvalues(eu,ev,ew)
    eu = clip(0.005,1.0,eu)
    ev = clip(0.005,1.0,ev)
    ew = clip(0.005,1.0,ew)
    et.setEigenvalues(eu,ev,ew)
    et.invertStructure(1.0,4.0,4.0)
    et.getEigenvalues(eu,ev,ew)
    eu = mul(eu,0.01)
    et.setEigenvalues(eu,ev,ew)
    sig = 50
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,et,fx)
    '''
    eu = fillfloat(0.0001,n1,n2,n3)
    ev = pow(ep,30)
    ew = fillfloat(1.0000,n1,n2,n3)
    et.setEigenvalues(eu,ev,ew)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,50,fx,gx)
    writeImage(gxlfile,gx)
  else:
    gx = readImage(gxlfile)
  plot3(fx)
  plot3(gx)
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5)

def goSmoothS():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = readImage(epsfile)
    et = readTensors(etsfile)
    '''
    eu = zerofloat(n1,n2,n3)
    ev = zerofloat(n1,n2,n3)
    ew = zerofloat(n1,n2,n3)
    et.getEigenvalues(eu,ev,ew)
    eu = clip(0.005,1.0,eu)
    ev = clip(0.005,1.0,ev)
    ew = clip(0.005,1.0,ew)
    et.setEigenvalues(eu,ev,ew)
    et.invertStructure(1.0,4.0,4.0)
    et.getEigenvalues(eu,ev,ew)
    eu = mul(eu,0.01)
    et.setEigenvalues(eu,ev,ew)
    '''
    eu = fillfloat(0.0001,n1,n2,n3)
    ev = pow(ep,30)
    ew = fillfloat(1.0000,n1,n2,n3)
    et.setEigenvalues(eu,ev,ew)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,50,fx,gx)
    '''
    sig = 10
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,et,fx)
    '''
    writeImage(gxsfile,gx)
  else:
    gx = readImage(gxsfile)
  plot3(fx)
  plot3(gx)
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5)

def goLinearDiffusion():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=4,2
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    ets.setEigenvalues(0.0001,1.0,1.0)
    sig = 5
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,ets,fx)
    writeImage(gxlfile,gx)
  else:
    gx = readImage(gxlfile)
  plot3(sub(fx,gx),cmin=-1.0,cmax=1.0)
  plot3(gx)
  plot3(fx)
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
    sig1,sig2=8,2
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    ets.setEigenvalues(0.0001,0.0001,1.0000)
    sig = 10
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,ets,fx)
    writeImage(gxsfile,gx)
  else:
    gx = readImage(gxsfile)
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5)
  plot3(gx,cmin=-1.0,cmax=1.0)
  plot3(fx,cmin=-1.0,cmax=1.0)

def goStratigraphyOrientedDiffusionX():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=4,2
    ep = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    sof = StratigraphicOrientFilter(sig1,sig2)
    et = sof.applyForTensors(ets,fx,ep)
    et.setEigenvalues(0.0001,0.0001,1.0000)
    sig = 10
    cycle,limit=3,0.5
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,et,fx)
    writeImage(gxsfile,gx)
    writeImage("epx",ep)
  else:
    fx = readImage(fxfile)
    gx = readImage(gxsfile)
    ep = readImage("epx")
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5)
  plot3(gx,cmin=-1.0,cmax=1.0)
  plot3(fx,cmin=-1.0,cmax=1.0)
  plot3(ep,cmin=0.4,cmax=1.0)

def goNonlinearDiffusion():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=4,2
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    ets.setEigenvalues(0.0001,1.0,1.0)
    sig = 5
    cycle,limit=3,0.5
    lbd = 0.1
    fed = FastExplicitDiffusion()
    fed.setCycles(cycle,limit)
    gx = fed.apply(sig,lbd,ets,fx)
    writeImage(gxnfile,gx)
  else:
    gx = readImage(gxnfile)
  plot3(sub(fx,gx),cmin=-1.0,cmax=1.0)
  plot3(gx)
  plot3(fx)

def goShapeSemblance():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,6
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
    sm = Semblance()
    '''
    sd = mul(fx,fx)
    sn = sm.smoothVW(2,et,fx)
    sd = sm.smoothVW(2,et,sd)
    sn = mul(sn,sn)
    writeImage("sn",sn)
    writeImage("sd",sd)
    '''
    sn = readImage("sn")
    sd = readImage("sd")
    wp = sub(1,ep)
    #wp = fillfloat(1,n1,n2,n3)
    et.setEigenvalues(0.2,0.0001,1.0000)
    ss = sm.shapeSemblance(et,wp,sn,sd)
    writeImage("ss",ss)
  else:
    ss = readImage("ss")
  plot3(fx)
  plot3(ss,cmin=0.2,cmax=1.0)


def goSemblance():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=2,6
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    et.setEigenvalues(0.0001,0.5,1.000)
    fxs = goDiffusion(2,et,fx)
    fxs = mul(fxs,fxs)
    fss = mul(fx,fx)
    fss = goDiffusion(2,et,fss)
    et.setEigenvalues(0.2,0.001,1.0)
    fxs = goDiffusion(4,et,fxs)
    fss = goDiffusion(4,et,fss)
    sem = div(fxs,fss)
    writeImage("sem",sem)
  else:
    sem = readImage("sem")
  plot3(fx)
  plot3(sem,cmin=0.3,cmax=1.0)

def goSemblanceHale():
  fx = readImage(fxfile)
  if not plotOnly:
    et = readTensors(etlfile)
    lsf = LocalSemblanceFilter(2,2)
    sem = lsf.semblance(LocalSemblanceFilter.Direction3.VW,et,fx)
    writeImage("semHale",sem)
  else:
    sem = readImage("semHale")
  plot3(fx)
  plot3(sem,cmin=0.2,cmax=1.0)

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

def plot3(f,g=None,et=None,ep=None,k1=120,
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
  if et:
    node = TensorEllipsoids(s1,s2,s3,et,ep)
    states = StateSet.forTwoSidedShinySurface(Color.YELLOW);
    node.setStates(states)
    sf.world.addChild(node)
  if cbar:
    cbar.setWidthMinimum(120)
  #ipg.setSlices(153,760,450)
  ipg.setSlices(k1,760,450)
  #ipg.setSlices(85,5,102)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  if cbar:
    sf.setSize(837,700)
  else:
    sf.setSize(700,700)

  view = sf.getOrbitView()
  #zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  zscale = 0.6*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.72)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(225.0)
  view.setElevation(45)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.05,-0.1,0.06))
  sf.viewCanvas.setBackground(sf.getBackground())
  sf.setSize(850,700)

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
