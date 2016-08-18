import sys

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *
from java.util import *


from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *

from hv import *
from uff import *
from util import *
from ad import *
from sso import *

pngDir = None
pngDir = "../../../png/sso/3d/nwc/"

#seismicDir = "../../../data/seis/sso/3d/pnz/"
seismicDir = "../../../data/seis/sso/3d/nwc/"
fxfile = "gn"
ellfile = "ell"
elsfile = "els"
eplfile = "epl"
p2lfile = "p2l"
p3lfile = "p3l"
p2sfile = "p2s"
p3sfile = "p3s"
epsfile = "eps"
etlfile = "etl"
etsfile = "ets"
etcfile = "etc"
gslfile = "gsl"
gssfile = "gss"
w2cfile = "w2c"
w3cfile = "w3c"
azcfile = "azc"
epcfile = "epc"
gclfile = "gcl"
gcsfile = "gcs"
hvlfile = "hvl"
hvsfile = "hvs"
hzfile = "hz"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 320,1001,701
n1,n2,n3 = 651,601,401
n1,n2,n3 = 120,350,401
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = True

def main(args):
  #goLof()
  #goLoe()
  #goChannel()
  #goSlopes()
  #goSmoothSL()
  #goSlices()
  #goFirstLook()
  #goHorizonS()
  #goHorizons()
  #goFolding()
  goHorizon()
def goHorizon():
  gx = readImage(fxfile)
  if not plotOnly:
    p2 = readImage(p2sfile)
    p3 = readImage(p3sfile)
    ep = readImage(eplfile)
    wp = pow(ep,8)
    lmt = n1-1
    k11 = [ 74, 85, 76, 73, 79, 77, 79, 77, 79, 73, 72, 71, 72, 71, 73, 73,87]
    k12 = [ 97,170,261,273,104,243,261,309,328, 72, 10, 98,287,295,292,302,62]
    k13 = [196, 70,232,279,309,338,153,202,192,280,245,251,316,338,272,259,20]
    se = SurfaceExtractorC()
    se.setWeights(0.0)
    se.setSmoothings(4.0,4.0)
    se.setCG(0.01,100)
    surf = se.surfaceInitialization(n2,n3,lmt,k11,k12,k13)
    se.surfaceUpdateFromSlopes(wp,p2,p3,k11,k12,k13,surf)
    writeImage(hzfile,surf)
  else:
    surf = readImage2D(n2,n3,hzfile)
  etc = readTensors(etcfile)
  etl = readTensors(etlfile)
  plot3(gx,hz=surf,scale=1.3,cmin=-1.5,cmax=1.5,png="seisF")
  plot3(gx,hz=surf,cmin=-1.5,cmax=1.5,png="seis")
  plot3(gx,hz=surf,et=etc,cmin=-1.5,cmax=1.5,png="etl")
  plot3(gx,hz=surf,et=etl,cmin=-1.5,cmax=1.5,png="ets")

def goFolding():
  fx = readImage("gxfull")
  fx = gain(fx)
  hp = Helper()
  gx = hp.applyFolding(fx)
  gxs = copy(200,401,n3,0,150,0,gx)
  fxs = copy(200,401,n3,0,150,0,fx)
  writeImage("gxs",gxs)
  writeImage("fxs",fxs)
  plot3(fx)
  plot3(gx)

def goHorizons():
  ns = 100
  gx = readImage(fxfile)
  hs = readHorizons(ns,hvsfile)
  et = readTensors(etcfile)
  #et = readTensors(etlfile)
  ks = [10,36,96]
  for ih in range(10,40,1):
    plot3(gx,hz=hs[ih],et=et,cmin=-1.5,cmax=1.5,clab=str(ih))
    #plot3(gx,hz=hs[ih],et=etl,cmin=-1.5,cmax=1.5,clab=str(ih))
def goChannel():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    w2 = zerofloat(n1,n2,n3)
    w3 = zerofloat(n1,n2,n3)
    et = readTensors(etsfile)
    loe = LocalOrientEstimator(et,5)
    loe.setEigenvalues(0.2,1.0,1.0)
    loe.setGradientSmoothing(3)
    loe.applyForStratigraphy(fx,w2,w3,ep)
    loe.updateTensors(et,w2,w3)
    hp = Helper()
    az = hp.channelAzimuth(w2,w3)
    writeTensors(etcfile,et)
    writeImage(w2cfile,w2)
    writeImage(w3cfile,w3)
    writeImage(epcfile,ep)
    writeImage(azcfile,az)
  else:
    ep = readImage(epcfile)
  plot3(fx)
  plot3(ep,cmin=0.1,cmax=0.9)

def goHorizonS():
  gx = readImage(fxfile)
  ns = 100
  if not plotOnly:
    p2s = readImage(p2sfile)
    p3s = readImage(p3sfile)
    ep = fillfloat(1,n1,n2,n3)
    c1 = rampfloat(10,1,119)
    c2 = fillfloat(200,ns)
    c3 = fillfloat(200,ns)
    hv = HorizonVolume()
    hv.setCG(0.01,50)
    hs = hv.applyForHorizonVolume(c1,c2,c3,ep,p2s,p3s)
    writeImage(hvsfile,hs)
  else:
    hs = readHorizons(ns,hvsfile)
  k1,k2,k3=160,390,316
  plot3p(s1,s2,s3,gx,hv=hs,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,
          clab="Amplitude",png="hvs")
  k1,k2,k3=160,386,389
  hd = HorizonDisplay()
  cs = hd.horizonCurves(k2,k3,hs)
  plot3(gx,hs=cs,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="cvs")

def goFirstLook():
  fx = readImage(fxfile)
  fx = gain(fx)
  writeImage(fxfile,fx)
  plot3(fx)

def goLof():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=4,2
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(sig1,sig2)
    lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
    writeImage(eplfile,ep)
    ets = lof.applyForTensors(fx)
    writeTensors(etlfile,ets)
  else:
    ets = readTensors(etlfile)
  k1,k2,k3=76,348,376
  plot3(fx,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="gx")

def goLoe():
  fx = readImage(fxfile)
  if not plotOnly:
    et = readTensors(etlfile)
    loe = LocalOrientEstimator(et,10)
    loe.setGradientSmoothing(0)
    loe.setEigenvalues(1.0,0.2,0.2)
    ets = loe.applyForTensorsX(fx)
    writeTensors(etsfile,ets)
  else:
    et = readTensors(etsfile)
  plot3(fx)

def goSlopes():
 etl = readTensors(etlfile)
 ets = readTensors(etsfile)
 loe = LocalOrientEstimator(etl,5)
 p2l,p3l = loe.slopesFromTensors(10,etl)
 p2s,p3s = loe.slopesFromTensors(10,ets)
 writeImage(p2lfile,p2l)
 writeImage(p3lfile,p3l)
 writeImage(p2sfile,p2s)
 writeImage(p3sfile,p3s)
 gx = readImage(fxfile)
 clab2 = "Inline slope (samples/trace)"
 clab3 = "Crossline slope (samples/trace)"
 k1,k2,k3=76,348,376
 cmap = jetFill(0.5)
 plot3(gx,p2l,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p2l")
 plot3(gx,p3l,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p3l")
 plot3(gx,p2s,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p2s")
 plot3(gx,p3s,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p3s")
 '''
 clab2 = "Inline slope (samples/trace)"
 clab3 = "Crossline slope (samples/trace)"
 cmap = jetFill(0.5)
 k1,k2,k3=40,390,316
 plot3p(s1,s2,s3,gx,g=p2l,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab2,png="p2lp")
 plot3p(s1,s2,s3,gx,g=p3l,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab3,png="p3lp")
 plot3p(s1,s2,s3,gx,g=p2s,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab2,png="p2sp")
 plot3p(s1,s2,s3,gx,g=p3s,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab3,png="p3sp")
 '''

def goSmoothSL():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=4,2
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    et.setEigenvalues(0.001,0.001,1)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,50,fx,gx)
    writeImage(gslfile,gx)
  else:
    gx = readImage(gslfile)
  plot3(fx,cmin=-0.5,cmax=0.5,clab="Amplitude",png="gx")
  plot3(gx,cmin=-0.5,cmax=0.5,clab="Amplitude",png="gsl")
  plot3(sub(fx,gx),cmin=-0.2,cmax=0.2,clab="Amplitude",png="gxl")
  '''
  k1,k2,k3=66,390,408
  plot3p(fx,k1,k2,k3,s1,s2,s3,cmin=-0.3,cmax=0.3,clab="Amplitude")
  plot3p(gx,k1,k2,k3,s1,s2,s3,cmin=-0.3,cmax=0.3,clab="Amplitude",png="gsl")
  plot3p(sub(fx,gx),k1,k2,k3,s1,s2,s3,cmin=-0.1,cmax=0.1,clab="Amplitude",png="gxl")
  '''


def goSmoothSS():
  fx = readImage(fxfile)
  if not plotOnly:
    et = readTensors(etsfile)
    et.setEigenvalues(0.001,1,1)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,20,fx,gx)
    writeImage(gssfile,gx)
  else:
    gx = readImage(gssfile)
  plot3(fx,cmin=-0.5,cmax=0.5,clab="Amplitude",png="gx")
  plot3(gx,cmin=-0.5,cmax=0.5,clab="Amplitude",png="gss")
  plot3(sub(fx,gx),cmin=-0.2,cmax=0.2,clab="Amplitude",png="gxs")

  '''
  k1,k2,k3=66,390,408
  plot3p(fx,k1,k2,k3,s1,s2,s3,cmin=-0.3,cmax=0.3,clab="Amplitude")
  plot3p(gx,k1,k2,k3,s1,s2,s3,cmin=-0.3,cmax=0.3,clab="Amplitude",png="gss")
  plot3p(sub(fx,gx),k1,k2,k3,s1,s2,s3,cmin=-0.1,cmax=0.1,clab="Amplitude",png="gxs")
  '''


def goSmoothC():
  fx = readImage(fxfile)
  if not plotOnly:
    et = readTensors(etsfile)
    w2 = readImage(w2cfile)
    w3 = readImage(w3cfile)
    loe = LocalOrientEstimator(et,10)
    loe.updateTensors(et,w2,w3)
    et.setEigenvalues(0.001,0.001,1)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,50,fx,gx)
    writeImage(gcsfile,gx)
  else:
    gx = readImage(gcsfile)
  plot3(gx,cmin=-1.5,cmax=1.5)
  #plot3(sub(fx,gx),cmin=-0.5,cmax=0.5)

def goSlices():
  ns = 290
  gx = readImage(fxfile)
  etc = readTensors(etcfile)
  etl = readTensors(etlfile)
  hs = readHorizons(ns,hvsfile)
  for ih in range(145,146,5):
    #plot3(gs,hz=hs[ih],cmin=-1.5,cmax=1.5)
    plot3(gx,hz=hs[ih],et=etc,cmin=-1.5,cmax=1.5)
    plot3(gx,hz=hs[ih],et=etl,cmin=-1.5,cmax=1.5)
def normalize(ss):
  sub(ss,min(ss),ss)
  div(ss,max(ss),ss)

def addNoise(nrms, fx):
  r = Random(1);
  gx = mul(2.0,sub(randfloat(r,n1,n2,n3),0.5));
  rgf = RecursiveGaussianFilter(2.0);
  rgf.apply2XX(gx,gx)
  gx = mul(gx,nrms*rms(fx)/rms(gx));
  return add(fx,gx,fx);

def rms(fx):
  fs = sum(mul(fx,fx))
  return sqrt(fs/n1/n2/n3)

  
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

def readImage2D(n2,n3,basename):
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

def readHorizons(ns,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n2,n3,ns)
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
  cbar.setFont(Font("Arial",Font.PLAIN,24)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def plot3(f,g=None,hz=None,hs=None,et=None,k1=290,k2=17,k3=72,
    scale=2.5,cmin=None,cmax=None,cmap=None,clab=None,cint=None,png=None):
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
      ipg.setClips(-1.5,1.5)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1.5,1.5)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if hs:
    for hi in hs:
      lg = LineGroup(hi[0],hi[1])
      ss = StateSet()
      lg.setStates(ss)
      ls = LineState()
      ls.setWidth(4)
      ls.setSmooth(False)
      ss.add(ls)
      sf.world.addChild(lg)
  if hz:
    sd = SurfaceDisplay()
    ts = sd.horizonWithAmplitude(-1,[cmin,cmax],hz,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if et:
    tv = TensorView()
    hs = tv.applyForSegments(5,et,hz)
    cp = ColorMap(0,1,ColorMap.JET)
    vi = fillfloat(0.9,6)
    cb = cp.getRgbFloats(vi)
    for hi in hs:
      lg = LineGroup(hi,cb)
      ss = StateSet()
      lg.setStates(ss)
      ls = LineState()
      ls.setWidth(8)
      ls.setSmooth(False)
      ss.add(ls)
      sf.world.addChild(lg)
  if cbar:
    cbar.setWidthMinimum(100)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(970,700)
  else:
    sf.setSize(870,700)
  view = sf.getOrbitView()
  #zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(scale)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(45.0)
  view.setElevation(50)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(-0.01,-0.01,-0.01))
  #sf.viewCanvas.setBackground(sf.getBackground())
  sf.viewCanvas.setBackground(Color.WHITE)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3p(s1,s2,s3,f,g=None,hv=None,k1=None,k2=None,k3=None,cmap=ColorMap.GRAY,
        cmin=-1,cmax=1,clab=None,cint=0.1,png=None):
  width,height,cbwm = 800,550,200
  n1,n2,n3 = s1.count,s2.count,s3.count
  orient = PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT;
  axespl = PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM
  panel = PlotPanelPixels3(orient,axespl,s1,s2,s3,f)
  #panel.mosaic.setWidthElastic(0,100)
  #panel.mosaic.setWidthElastic(1,75)
  panel.mosaic.setHeightElastic(0,150)
  #panel.mosaic.setHeightElastic(1,100)
  panel.setSlice23(k1)
  panel.setSlice13(k2)
  panel.setSlice12(k3)
  #panel.setSlice103(70)
  panel.setClips(cmin,cmax)
  if clab:
    cbar = panel.addColorBar(clab)
    cbar.setInterval(cint)
  panel.setColorBarWidthMinimum(50)
  panel.setLabel1("Samples")
  panel.setLabel2("Inline (traces)")
  panel.setLabel3("Crossline (traces)")
  panel.setInterval2(100)
  panel.setInterval3(100)
  panel.setColorModel(ColorMap.GRAY)
  panel.setLineColor(Color.WHITE)
  panel.setHLimits(0,s2.first,s2.last)
  panel.setVLimits(1,s1.first,s1.last)
  if g:
    pv12 = PixelsView(s1,s2,slice12(k3,g))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,g))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,g))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(cmap)
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    panel.pixelsView12.tile.addTiledView(pv12)
    panel.pixelsView13.tile.addTiledView(pv13)
    panel.pixelsView23.tile.addTiledView(pv23)
  if hv:
    nh = len(hv)
    hd = HorizonDisplay()
    cv12 = hd.slice12(k3,s2,hv)
    cv13 = hd.slice13(k2,s3,hv)
    cv23 = hd.slice23X(k1,s2,s3,hv)
    mp = ColorMap(0,nh,ColorMap.JET)
    print nh
    for ih in range(4,nh,1):
      pv12 = PointsView(cv12[ih][1],cv12[ih][0])
      pv13 = PointsView(cv13[ih][1],cv13[ih][0])
      pv12.setLineWidth(1.0)
      pv13.setLineWidth(1.0)
      pv12.setLineColor(mp.getColor(ih))
      pv13.setLineColor(mp.getColor(ih))
      panel.pixelsView12.tile.addTiledView(pv12)
      panel.pixelsView13.tile.addTiledView(pv13)
      nc = len(cv23[ih][0])
      for ic in range(nc):
        pv23 = PointsView(cv23[ih][0][ic],cv23[ih][1][ic])
        pv23.setLineWidth(1.0)
        pv23.setLineColor(mp.getColor(ih))
        panel.pixelsView23.tile.addTiledView(pv23)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSize(12)#ForSlide(1.0,0.8)
  frame.setSize(width,height)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,3.3,pngDir+"/"+png+".png")

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
