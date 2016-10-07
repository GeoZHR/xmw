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
pngDir = "../../../png/sso/3d/nwc/"

seismicDir = "../../../data/seis/sso/3d/nwc/"
fxfile = "gn"
ellfile = "ell"
elsfile = "els"
eplfile = "epl"
epsfile = "eps"
etlfile = "etl"
etsfile = "ets"
hzfile = "hz"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 120,350,401
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = False
k1 = 51
k1 = 56

def main(args):
  #goLof()
  goSta()
def goLof():
  fx = readImage(fxfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    el = zerofloat(n1,n2,n3)
    sig1,sig2=4,2
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
    writeImage(eplfile,ep)
    writeTensors(etlfile,et)
  else:
    ep = readImage(eplfile)
  ep = pow(ep,2.5)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  surf = readImage2D(n2,n3,hzfile)
  plot3(ep,hz=surf,cmin=0.2,cmax=1.0,png="epl")
def goSta():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    el = zerofloat(n1,n2,n3)
    et = readTensors(etlfile)
    sta = StructureTensorAttribute(et,20)
    sta.setEigenvalues(1.0,0.01,0.6)
    #sta.updateTensors(4,fx)
    sta.applyForPlanarLinear(fx,ep,el)
    writeImage(epsfile,ep)
    writeImage(elsfile,el)
    writeTensors(etsfile,et)
  else:
    ep = readImage(epsfile)
    el = readImage(elsfile)
    #et = readTensors(etsfile)
  ep = pow(ep,2.5)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  surf = readImage2D(n2,n3,hzfile)
  plot3(ep,hz=surf,cmin=0.2,cmax=1.0,png="eps")

def goSemblance():
  fx = readImage(fxfile)
  if not plotOnly:
    sig1,sig2=8,2
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig2,5)
    lsf.findSlopes(fx,p2,p3,ep)
    cov = Covariance()
    em,es=cov.covarianceEigen(8,p2,p3,fx)
    print "done..."
    sem = div(em,es)
    writeImage(semfile,sem)
  else:
    sem = readImage(semfile)
  sem = pow(sem,2.0)
  sem = sub(sem,min(sem))
  sem = div(sem,max(sem))
  plot3(sem,k1=k1,cmin=0.2,cmax=1.0,clab="Semblance",cint=0.1,png="sem"+str(k1))

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
    ts = sd.horizonWithAmplitude([cmin,cmax],hz,f)
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
