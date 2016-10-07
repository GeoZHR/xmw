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

from ipfx import *
from ad import *
from util import *
from sso import *


pngDir = None
pngDir = "../../../png/sso/3d/fake/"

seismicDir = "../../../data/seis/sso/3d/fake/"
fxfile = "fx"
hzfile = "hz"
epsfile = "eps"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 121,152,153
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = True

def main(args):
  #goLof()
  goSta()
  #goSemblance()
def goSta():
  fx = readImage(fxfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    el = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(4,2)
    et = lof.applyForTensors(fx)
    sta = StructureTensorAttribute(et,30)
    sta.setGradientSmoothing(2)
    sta.setEigenvalues(1.0,0.1,0.6)
    sta.applyForPlanarLinear(fx,ep,el)
    writeImage(epsfile,ep)
  else:
    ep = readImage(epsfile)
  ep = pow(ep,6)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  hz = readImage2D(n2,n3,hzfile)
  plot3(ep,hz=hz,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="eps")

def goVW():
  sk = readSkins(fskbase)
  gx = readImage(gxfile)
  hz = readImage2D(n2,n3,hzfile)
  lof = LocalOrientFilter(4,3)
  et = lof.applyForTensors(gx)
  plot3(gx,cmin=-2,cmax=1.5,clab="Amplitude",hz=hz,mk=-1,png="gnhz")
  plot3(gx,cmin=-2,cmax=1.5,hz=hz,mk=-1,et=et,png="gnhzw")
  plot3(gx,cmin=-2,cmax=1.5,skins=sk,hz=hz,mk=-1,et=et,w=False,png="gnhzv")


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

def readImage2D(n1,n2,name):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2)
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

def skinName(basename,index):
  return basename+("%03i"%(index))
def skinIndex(basename,fileName):
  assert fileName.startswith(basename)
  i = len(basename)
  return int(fileName[i:i+3])

def listAllSkinFiles(basename):
  """ Lists all skins with specified basename, sorted by index. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  return fileNames

def removeAllSkinFiles(basename):
  """ Removes all skins with specified basename. """
  fileNames = listAllSkinFiles(basename)
  for fileName in fileNames:
    File(seismicDir+fileName).delete()

def readSkin(basename,index):
  """ Reads one skin with specified basename and index. """
  return FaultSkin.readFromFile(seismicDir+skinName(basename,index)+".dat")

def readSkins(basename):
  """ Reads all skins with specified basename. """
  fileNames = []
  for fileName in File(seismicDir).list():
    if fileName.startswith(basename):
      fileNames.append(fileName)
  fileNames.sort()
  skins = []
  for iskin,fileName in enumerate(fileNames):
    index = skinIndex(basename,fileName)
    skin = readSkin(basename,index)
    skins.append(skin)
  return skins

def writeSkin(basename,index,skin):
  """ Writes one skin with specified basename and index. """
  FaultSkin.writeToFile(seismicDir+skinName(basename,index)+".dat",skin)

def writeSkins(basename,skins):
  """ Writes all skins with specified basename. """
  for index,skin in enumerate(skins):
    writeSkin(basename,index,skin)


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
    print min(ha)
    print max(ha)
    ts = sd.horizonWithChannelAzimuth([cmin,cmax],[amin,amax],hz,f,ha)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if dh:
    pi = Math.PI
    amin = 0
    amax = 15
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
