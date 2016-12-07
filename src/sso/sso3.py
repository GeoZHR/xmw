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
pngDir = "../../../png/sso/3d/poseidon/"

#seismicDir = "../../../data/seis/sso/3d/real/"
seismicDir = "../../../data/seis/sso/3d/poseidon/"
#seismicDir = "../../../data/seis/beg/jake/subs/"
fxfile = "fxs"
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
gclfile = "gcl"
gcsfile = "gcs"
hvlfile = "hvl"
hvsfile = "hvs"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 80,500,600
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = True

def main(args):
  #goLof()
  #goLoe()
  goSlopes()
  #goHorizonL()
  #goHorizonS()
  #goSurfaces()
  #goChannel()
  #goSmoothSL()
  #goSmoothSS()
  #goSmoothC()
  #goFirstLook()
  #gl = readImage(gslfile)
  #plot3(sub(gs,gl),cmin=-0.5,cmax=0.5)
def goHorizonL():
  gx = readImage(fxfile)
  ns = 40
  if not plotOnly:
    etl = readTensors(etlfile)
    loe = LocalOrientEstimator(etl,5)
    p2l,p3l = loe.slopesFromTensors(10,etl)
    ep = fillfloat(1,n1,n2,n3)
    c1 = rampfloat(25,1,65)
    c2 = fillfloat(390,ns)
    c3 = rampfloat(340,4,500)
    hv = HorizonVolume()
    hv.setCG(0.01,50)
    hs = hv.applyForHorizonVolume(c1,c2,c3,ep,p2l,p3l)
    writeImage(hvlfile,hs)
  else:
    hs = readHorizons(ns,hvlfile)
  k1,k2,k3=41,390,416
  plot3p(s1,s2,s3,gx,hv=hs,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,
          clab="Amplitude",png="hvl")
  k1,k2,k3=76,348,376
  k1,k2,k3=76,386,399
  hd = HorizonDisplay()
  cs = hd.horizonCurves(k2,k3,hs)
  plot3(gx,hs=cs,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="cvl")

def goHorizonS():
  gx = readImage(fxfile)
  ns = 40
  if not plotOnly:
    ets = readTensors(etsfile)
    loe = LocalOrientEstimator(ets,5)
    p2s,p3s = loe.slopesFromTensors(10,ets)
    ep = fillfloat(1,n1,n2,n3)
    c1 = rampfloat(25,1,65)
    c2 = fillfloat(390,ns)
    c3 = rampfloat(340,4,500)
    hv = HorizonVolume()
    hv.setCG(0.01,50)
    hs = hv.applyForHorizonVolume(c1,c2,c3,ep,p2s,p3s)
    writeImage(hvsfile,hs)
  else:
    hs = readHorizons(ns,hvsfile)
  k1,k2,k3=41,390,416
  plot3(gx,surf=hs[10],k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="sf10")
  plot3p(s1,s2,s3,gx,hv=hs,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,
          clab="Amplitude",png="hvs")
  k1,k2,k3=76,348,376
  k1,k2,k3=76,386,399
  hd = HorizonDisplay()
  cs = hd.horizonCurves(k2,k3,hs)
  plot3(gx,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="gx3d")
  plot3(gx,hs=cs,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="cvs")
  
def goSurfaces():
  ns = 40
  gx = readImage(fxfile)
  hs = readHorizons(ns,hvsfile)
  k1,k2,k3=79,499,599
  for k in range(ns):
    plot3s(gx,surf=hs[k],k1=k1,k2=k2,k3=k3,cmin=-0.8,cmax=0.8,png="sf"+str(k))
def goFirstLook():
  fx = readImage(fxfile)
  fx = gain(fx)
  writeImage(fxfile,fx)
  plot3(fx)
def goLof():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=8,2
    lof = LocalOrientFilter(sig1,sig2)
    ets = lof.applyForTensors(fx)
    writeTensors(etlfile,ets)
  else:
    ets = readTensors(etlfile)
  k1,k2,k3=41,390,416
  plot3p(s1,s2,s3,fx,k1=k1,k2=k2,k3=k3,
        cmap=ColorMap.GRAY,cmin=-0.5,cmax=0.5,clab="Amplitude",png="gxp")
  k1,k2,k3=76,348,376
  plot3(fx,k1=k1,k2=k2,k3=k3,cmin=-0.5,cmax=0.5,png="gx")

def goLoe():
  fx = readImage(fxfile)
  if not plotOnly:
    et = readTensors(etlfile)
    loe = LocalOrientEstimator(et,10)
    loe.setGradientSmoothing(0)
    loe.setEigenvalues(1.0,0.1,0.1)
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
 gx = readImage(fxfile)
 clab2 = "Inline slope (samples/trace)"
 clab3 = "Crossline slope (samples/trace)"
 k1,k2,k3=76,348,376
 k1,k2,k3=76,386,399
 cmap = jetFill(0.5)
 plot3(gx,p2l,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p2l")
 plot3(gx,p3l,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p3l")
 plot3(gx,p2s,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p2s")
 plot3(gx,p3s,k1=k1,k2=k2,k3=k3,cmap=cmap,cmin=-0.6,cmax=0.6,png="p3s")
 clab2 = "Inline slope (samples/trace)"
 clab3 = "Crossline slope (samples/trace)"
 cmap = jetFill(0.5)
 k1,k2,k3=40,390,416
 plot3p(s1,s2,s3,gx,g=p2l,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab2,png="p2lp")
 plot3p(s1,s2,s3,gx,g=p3l,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab3,png="p3lp")
 plot3p(s1,s2,s3,gx,g=p2s,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab2,png="p2sp")
 plot3p(s1,s2,s3,gx,g=p3s,k1=k1,k2=k2,k3=k3,
        cmap=cmap,cmin=-0.6,cmax=0.6,clab=clab3,png="p3sp")
def goChannel():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    w2 = zerofloat(n1,n2,n3)
    w3 = zerofloat(n1,n2,n3)
    et = readTensors(etsfile)
    loe = LocalOrientEstimator(et,5)
    loe.setEigenvalues(0.1,1.0,1.0)
    loe.setGradientSmoothing(3)
    loe.applyForStratigraphy(fx,w2,w3,ep)
    loe.updateTensors(et,w2,w3)
    writeTensors(etcfile,et)
  else:
    ep = readImage(epfile)
  plot3(fx,cmin=-2,cmax=2)
  plot3(ep,hz=hz,cmin=0.1,cmax=1.0)


def goSmoothSL():
  fx = readImage(fxfile)
  if not plotOnly:
    et = readTensors(etlfile)
    et.setEigenvalues(0.001,1,1)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,20,fx,gx)
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
    et = readTensors(etcfile)
    et.setEigenvalues(0.001,0.001,1)
    gx = zerofloat(n1,n2,n3)
    lsf = LocalSmoothingFilter()
    lsf.apply(et,50,fx,gx)
    writeImage(gxcfile,gx)
  else:
    gx = readImage(gxcfile)
  plot3(fx)
  plot3(gx)
  plot3(sub(fx,gx),cmin=-0.5,cmax=0.5)
  plot3p(fx,80,390,408,s1,s2,s3,cmin=-1,cmax=1,clab="Amplitude")
  plot3p(gx,80,390,408,s1,s2,s3,cmin=-1,cmax=1,clab="Amplitude")
  plot3p(sub(fx,gx),80,390,408,s1,s2,s3,cmin=-0.5,cmax=0.5,clab="Amplitude")


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

def plot3s(f,surf=None,k1=None,k2=None,k3=None,cmin=None,cmax=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  ipg = sf.addImagePanels(s1,s2,s3,f)
  if cmin!=None and cmax!=None:
    ipg.setClips(cmin,cmax)
  else:
    ipg.setClips(-1,1)
  if surf:
    sd = SurfaceDisplay()
    xyz,rgb = sd.horizonWithAmplitude([cmin,cmax],surf,f)
    tgs = TriangleGroup(True,xyz,rgb)
    sf.world.addChild(tgs)
  ipg.setSlices(k1,k2,k3)
  sf.setSize(870,700)
  view = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.2)
  view.setAzimuth(210.0)
  view.setElevation(40)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(-0.18,0.30,-0.25))
  sf.viewCanvas.setBackground(Color.WHITE)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3(f,g=None,hs=None,surf=None,k1=None,k2=None,k3=None,
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
      ipg.setClips(-1,1)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1,1)
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
      ls.setWidth(8)
      ls.setSmooth(False)
      ss.add(ls)
      sf.world.addChild(lg)
  if surf:
    sd = SurfaceDisplay()
    xyz,rgb = sd.horizonWithAmplitude([cmin,cmax],surf,f)
    tgs = TriangleGroup(True,xyz,rgb)
    sf.world.addChild(tgs)
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
  view.setScale(2.8)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(145.0)
  view.setElevation(40)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(-0.18,-0.06,-0.25))
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
