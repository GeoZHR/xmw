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
from util import *
from sso import *

pngDir = None
pngDir = "../../../png/sso/3d/fake/"

seismicDir = "../../../data/seis/sso/3d/fake/"
#seismicDir = "../../../data/seis/beg/jake/subs/"
fxfile = "fx"
fxcfile = "fxc"
ellfile = "ell"
elsfile = "els"
eplfile = "epl"
epsfile = "eps"
etlfile = "etl"
etsfile = "ets"
gxlfile = "gxl"
gxsfile = "gxs"
p2kfile = "p2k"
p3kfile = "p3k"
p2lfile = "p2l"
p3lfile = "p3l"
p2sfile = "p2s"
p3sfile = "p3s"
gxsfile = "gxs"
gxlfile = "gxl"

hzfile = "hz"
halfile = "hal"
hasfile = "has"
hacfile = "hac"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 121,152,153
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = False

def main(args):
  #goFakeData()
  #goLof()
  #goLoe()
  #goStratigraphy()
  goChannel()
  #goSmoothS()

def goFakeData():
  sequence = 'OA' # 1 episode of folding, followed by one episode of faulting
  #sequence = 'OOOOOAAAAA' # 5 episodes of folding, then 5 of faulting
  nplanar = 1 # number of planar faults
  conjugate = False # if True, two large planar faults will intersect
  conical = False # if True, may want to set nplanar to 0 (or not!)
  impedance = False # if True, data = impedance model
  wavelet = True # if False, no wavelet will be used
  lateralViriation = True
  noise = 0.4 # (rms noise)/(rms signal) ratio
  if not plotOnly:
    fx,p2,p3,hs = FakeData.seismicAndSlopes3d2015A(sequence,
    nplanar,conjugate,conical,impedance,wavelet,lateralViriation,noise)
    hz = hs[0]
    ha = hs[1]
    writeImage(fxfile,fx)
    writeImage(hzfile,hz)
    writeImage(hacfile,ha)
    writeImage(p2kfile,p2)
    writeImage(p3kfile,p3)
  else:
    fx = readImage(fxfile)
    p2 = readImage(p2kfile)
    p3 = readImage(p3kfile)
    hz = readImage2D(hzfile)
    ha = readImage2D(hacfile)

  print "ha min =",min(ha)," max =",max(ha)
  print "fx min =",min(fx)," max =",max(fx)
  print "p2 min =",min(p2)," max =",max(p2)
  print "p3 min =",min(p3)," max =",max(p3)
  gmin,gmax,gmap = -2.0,2.0,ColorMap.GRAY
  plot3(fx,cmin=gmin,cmax=gmax,png="fx")
  plot3(fx,hz=hz,cmin=gmin,cmax=gmax,cmap=gmap,clab="Amplitude",png="fxhz")
  plot3(fx,ha=ha,cmin=gmin,cmax=gmax,png="fxha")
  plot3(fx,g=p2,cmin=-1.2,cmax=1.2,cmap=jetFill(1.0),png="p2k")
  plot3(fx,g=p3,cmin=-1.2,cmax=1.2,cmap=jetFill(1.0),png="p3k")

def goLof():
  fx = readImage(fxfile)
  hz = readImage2D(hzfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    w1 = zerofloat(n1,n2,n3)
    w2 = zerofloat(n1,n2,n3)
    w3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    sig1,sig2=4,2
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    lof.apply(fx,None,None,u1,u2,u3,None,None,None,w1,w2,w3,None,None,None,ep,None)
    hp = Helper()
    ha = hp.channelAzimuth(w2,w3,hz)
    p2 = mul(div(u2,u1),-1)
    p3 = mul(div(u3,u1),-1)
    writeImage(p2lfile,p2)
    writeImage(p3lfile,p3)
    writeImage(halfile,ha)
    writeImage(eplfile,ep)
    writeTensors(etlfile,et)
  else:
    ep = readImage(eplfile)
    et = readTensors(etlfile)
    p2 = readImage(p2lfile)
    p3 = readImage(p3lfile)
  p2k = readImage(p2kfile)
  p3k = readImage(p3kfile)
  ha = readImage2D(halfile)
  hc = readImage2D(hacfile)
  ep = pow(ep,6)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  dp3 = abs(sub(p3k,p3))
  dp2 = abs(sub(p2k,p2))
  dh = abs(sub(hc,ha))
  plot3(fx,dh=dh,cmin=-2,cmax=2,png="dhl")
  plot3(ep,hz=hz,cmin=0.2,cmax=1.0,png="epl")
  plot3(fx,ha=ha,cmin=-2,cmax=2,png="hal")
  plot3(fx,g=p2,cmin=-1.2,cmax=1.2,cmap=jetFill(1.0),png="p2l")
  plot3(fx,g=p3,cmin=-1.2,cmax=1.2,cmap=jetFill(1.0),png="p3l")
  plot3(fx,g=dp2,cmin=0.0,cmax=0.25,cmap=jetFill(1.0),cint=0.1,png="dp2l")
  plot3(fx,g=dp3,cmin=0.0,cmax=0.25,cmap=jetFill(1.0),cint=0.1,png="dp3l")


def goLoe():
  fx = readImage(fxfile)
  hz = readImage2D(hzfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    w1 = zerofloat(n1,n2,n3)
    w2 = zerofloat(n1,n2,n3)
    w3 = zerofloat(n1,n2,n3)
    et = readTensors(etlfile)
    loe = LocalOrientEstimator(et,5)
    loe.setEigenvalues(0.001,0.2,0.2)
    loe.setGradientSmoothing(3)
    loe.applyForSlopePlanar(10,fx,p2,p3,ep)
    loe.applyForInline(fx,w1,w2,w3)
    hp = Helper()
    ha = hp.channelAzimuth(w2,w3,hz)
    writeImage(p2sfile,p2)
    writeImage(p3sfile,p3)
    writeImage(hasfile,ha)
    writeImage(epsfile,ep)
  else:
    p2 = readImage(p2sfile)
    p3 = readImage(p3sfile)
    ep = readImage(epsfile)
    ha = readImage2D(hasfile)
  hz = readImage2D(hzfile)
  hc = readImage2D(hacfile)
  p2k = readImage(p2kfile)
  p3k = readImage(p3kfile)
  ep = pow(ep,6)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  dp3 = abs(sub(p3k,p3))
  dp2 = abs(sub(p2k,p2))
  dh = abs(sub(hc,ha))
  plot3(fx,cmin=-2,cmax=2)
  plot3(ep,hz=hz,cmin=0.2,cmax=1.0)
  plot3(fx,dh=dh,cmin=-2,cmax=2)
  plot3(fx,ha=ha,cmin=-2,cmax=2)
  plot3(fx,g=p2,cmin=-1.2,cmax=1.2,cmap=jetFill(1.0),clab="p2")
  plot3(fx,g=p3,cmin=-1.2,cmax=1.2,cmap=jetFill(1.0),clab="p3")
  plot3(fx,g=dp2,cmin=0.0,cmax=0.25,cmap=jetFill(1.0),cint=0.1,clab="Absolute errors")
  plot3(fx,g=dp3,cmin=0.0,cmax=0.25,cmap=jetFill(1.0),cint=0.1,clab="Absolute errors")


def goStratigraphy():
  fx = readImage(fxfile)
  hz = readImage2D(hzfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    w1 = zerofloat(n1,n2,n3)
    w2 = zerofloat(n1,n2,n3)
    w3 = zerofloat(n1,n2,n3)
    et = readTensors(etsfile)
    loe = LocalOrientEstimator(et,12)
    loe.setEigenvalues(1.00,0.05,1.0)
    loe.setGradientSmoothing(3)
    loe.applyForInline(fx,w1,w2,w3)
    hp = Helper()
    ha = hp.channelAzimuth(w2,w3,hz)
  else:
    hz = readImage2D(hzfile)
  hc = readImage2D(hacfile)
  ep = pow(ep,6)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  dh = abs(sub(hc,ha))
  plot3(fx,cmin=-2,cmax=2)
  plot3(ep,hz=hz,cmin=0.2,cmax=1.0)
  plot3(fx,dh=dh,cmin=-2,cmax=2)
  plot3(fx,ha=ha,cmin=-2,cmax=2)

def goChannel():
  fx = readImage(fxfile)
  hz = readImage2D(hzfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    w2 = zerofloat(n1,n2,n3)
    w3 = zerofloat(n1,n2,n3)
    et = readTensors(etsfile)
    loe = LocalOrientEstimator(et,5)
    loe.setEigenvalues(0.1,1.0,1.0)
    loe.setGradientSmoothing(3)
    loe.applyForStratigraphy(fx,w2,w3,ep)
    hp = Helper()
    ha = hp.channelAzimuth(w2,w3,hz)
  else:
    hz = readImage2D(hzfile)
  hc = readImage2D(hacfile)
  print min(w3)
  print max(w3)
  ep = pow(ep,2)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  dh = abs(sub(hc,ha))
  plot3(fx,cmin=-2,cmax=2)
  plot3(ep,hz=hz,cmin=0.1,cmax=1.0)
  plot3(fx,dh=dh,cmin=-2,cmax=2)
  plot3(fx,ha=ha,cmin=-2,cmax=2)

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
  hz = readImage2D(hzfile)
  plot3(fx,hz=hz,cmin=-2,cmax=2)
  plot3(gx,hz=hz,cmin=-2,cmax=2)

def goSmoothS():
  fx = readImage(fxfile)
  if not plotOnly:
    ets = readTensors(etlfile)
    '''
    loe = LocalOrientEstimator(ets,12)
    loe.setEigenvalues(1.00,0.05,1.0)
    loe.setGradientSmoothing(3)
    ets = loe.applyForTensors(fx)
    '''
    ets.setEigenvalues(0.001,0.1,1.0)
    lsf = LocalSmoothingFilter()
    gx = zerofloat(n1,n2,n3)
    lsf.apply(ets,20,fx,gx)
    writeImage(gxsfile,gx)
  else:
    gx = readImage(gxsfile)
  hz = readImage2D(hzfile)
  plot3(fx,hz=hz,cmin=-2,cmax=2)
  plot3(gx,hz=hz,cmin=-2,cmax=2)

def toDegrees(fx):
  pi = Math.PI
  scale = 180/pi
  return mul(scale,fx)

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

def readImage2D(basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

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

def plot3(f,g=None,et=None,ep=None,hz=None,ha=None,dh=None,k1=120,
    cmin=None,cmax=None,cmap=None,mk=-1,clab=None,cint=None,png=None):
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
      ipg.setClips(-2.0,2.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,2.0)
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
    hz = readImage2D(hzfile)
    node = TensorEllipsoids(s1,s2,s3,et,ep,hz)
    states = StateSet.forTwoSidedShinySurface(Color.YELLOW);
    node.setStates(states)
    sf.world.addChild(node)
  if hz:
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude(mk,[cmin,cmax],hz,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if ha:
    amin = -45 #-0.25*Math.PI
    amax =  45 #0.25*Math.PI
    sd = SurfaceDisplay()
    hz = readImage2D(hzfile)
    ha = toDegrees(ha)
    ts = sd.horizonWithChannelAzimuth([cmin,cmax],[amin,amax],hz,f,ha)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if dh:
    pi = Math.PI
    amin = 0
    amax = 20
    sd = SurfaceDisplay()
    hz = readImage2D(hzfile)
    dh = toDegrees(dh)
    ts = sd.horizonWithChannelAzimuth([cmin,cmax],[amin,amax],hz,f,dh)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)

  if cbar:
    cbar.setWidthMinimum(120)
  #ipg.setSlices(153,760,450)
  ipg.setSlices(101,138,39)
  ipg.setSlices(101,135,35)
  #ipg.setSlices(85,5,102)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  if cbar:
    sf.setSize(887,700)
  else:
    sf.setSize(750,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.48*sqrt(n1*n1+n2*n2+n3*n3)
  zscale = 0.80*max(n2*d2,n3*d3)/(n1*d1)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.4*n2,0.4*n3,radius))
  ov.setAzimuthAndElevation(140.0,40.0)
  ov.setTranslate(Vector3(-0.06,0.12,-0.27))
  ov.setScale(1.25)


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
