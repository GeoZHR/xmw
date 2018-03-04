import sys

from java.awt import *
from java.io import *
from java.nio import *
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
pngDir = "../../../png/nam/"

seismicDir = "../../../data/seis/nam/"
fxfile = "gxSub"
ellfile = "ell"
elsfile = "els"
eplfile = "epl"
epsfile = "eps"
etlfile = "etl"
etsfile = "ets"
gxlfile = "gxl"
gxsfile = "gxs"
semfile = "sem"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 100,280,100
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = True
k1 = 51
k1 = 58

def main(args):

  #goLof()
  #goSta()
  goSemblance()
  goDlResult()

def goDlResult():
  k1 = 58
  fx = readImageL(n1,n2,n3,fxfile)
  ct = readImageL(n1,n3,n2,"dlCh")
  ch = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
        ch[i3][i2] = ct[i2][i3]
  plot3(fx,k1=k1,clab="Amplitude",cint=0.2,png="seis")
  plot3(fx,g=ch,k1=k1,cmap=jetRamp(1.0),cmin=0.2,cmax=1,
          clab="DL channel scores",cint=0.2,png="dlCh")
def goLof():
  #fx = readImage(fxfile)
  fx = readImageL(n1,n2,n3,fxfile)
  if not plotOnly:
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    el = zerofloat(n1,n2,n3)
    sig1,sig2=8,2
    lof = LocalOrientFilter(sig1,sig2)
    et = lof.applyForTensors(fx)
    lof.applyForNormalPlanar(fx,u1,u2,u3,ep)
    lof.applyForInlineLinear(fx,u1,u2,u3,el)
    writeImage(eplfile,ep)
    writeImage(ellfile,el)
    writeTensors(etlfile,et)
  else:
    ep = readImage(eplfile)
    el = readImage(ellfile)
  ep = pow(ep,2.5)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  plot3(fx,k1=k1,clab="Amplitude",cint=0.4,png="fx"+str(k1))
  plot3(ep,k1=k1,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="epl"+str(k1))
  plot3(el,k1=k1,cmin=0.0,cmax=0.4,clab="Linearity",cint=0.1,png="ell"+str(k1))
def goSta():
  #fx = readImage(fxfile)
  fx = readImageL(n1,n2,n3,fxfile)
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
  writeImage(epsfile,ep)
  plot3(ep,k1=k1,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="eps"+str(k1))
  #plot3(el,k1=k1,cmin=0.0,cmax=0.4,clab="Linearity",cint=0.1,png="els"+str(k1))

def goSemblance():
  #fx = readImage(fxfile)
  fx = readImageL(n1,n2,n3,fxfile)
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
  sem = pow(sem,2.5)
  sem = sub(sem,min(sem))
  sem = div(sem,max(sem))
  #writeImage(semfile,sem)
  plot3(sem,k1=k1,cmin=0.2,cmax=1.0,clab="Coherence",cint=0.2,png="sem"+str(k1))

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

def readImageL(n1,n2,n3,basename):
  """ 
  Reads an image from a file with specified basename
  """
  fileName = seismicDir+basename+".dat"
  image = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,ByteOrder.LITTLE_ENDIAN)
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
    cbar.setWidthMinimum(85)
  ipg.setSlices(k1,36,10)
  if cbar:
    sf.setSize(987,650)
  else:
    sf.setSize(850,700)
  view = sf.getOrbitView()
  zscale = 0.35*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.68)
  view.setAzimuth(55.0)
  view.setElevation(48)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.05,-0.01,0.00))
  sf.viewCanvas.setBackground(sf.getBackground())
  #sf.setSize(850,700)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot2(s1,s2,f,cmin=None,cmax=None,cint=None,clab=None,png=None): 
  f1 = s1.getFirst()
  f2 = s2.getFirst()
  d1 = s1.getDelta()
  d2 = s2.getDelta()
  n1 = s1.getCount()
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation,PlotPanel.AxesPlacement.NONE)
  #panel.setVInterval(0.1)
  #panel.setHInterval(1.0)
  panel.setHLabel("Crossline (traces)")
  panel.setVLabel("Samples")
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  #pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #cb = panel.addColorBar();
  if cint:
    cb.setInterval(cint)
  if clab:
    cb.setLabel(clab)
  #panel.setColorBarWidthMinimum(50)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  #frame.setTitle("normal vectors")
  frame.setVisible(True);
  #frame.setSize(1020,700) #for f3d
  frame.setSize(500,325) #for poseidon
  #frame.setFontSize(13)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")

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
