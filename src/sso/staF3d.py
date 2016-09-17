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

pngDir = "../../../png/sso/3d/crf/"
pngDir = None

seismicDir = "../../../data/seis/sso/3d/sta/f3d/"
gxfile = "gx"
epfile = "ep"
p2file = "p2"
p3file = "p3"
ellfile = "ell"
elsfile = "els"
eplfile = "epl"
epsfile = "eps"
etlfile = "etl"
etsfile = "ets"
gxlfile = "gxl"
gxsfile = "gxs"
semfile = "smSub"
f1,f2,f3 = 0,0,0
d1,d2,d3 = 1,1,1
n1,n2,n3 = 210,920,825
n1,n2,n3 = 65,380,591

s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
plotOnly = False
k1 = 38
k1 = 51
k1 = 110

def main(args):
  #goSeis()
  goSta()
  #goSlope()
  #goSemblance()
  #goSub()
  #go2dPlots()
  #goSubPlots()
def goSubPlots():
  fx = readImage("fxSub")
  ep = readImage("epSub")
  sm = readImage("smSub")
  fs  = copy(n1,450,n3,0,0,0,fx)
  es  = copy(n1,450,n3,0,0,0,ep)
  ss  = copy(n1,450,n3,0,0,0,sm)
  plot3s(fs,k1=k1,clab="Amplitude",cint=0.2,png="fxSub"+str(k1))
  plot3s(es,k1=k1,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="epSub"+str(k1))
  plot3s(ss,k1=k1,cmin=0.2,cmax=1.0,clab="Coherence",cint=0.1,png="smSub"+str(k1))
def go2dPlots():
  fx = readImage("fxSub")
  ep = readImage("epSub")
  sm = readImage("smSub")
  fs = slice(k1,fx)
  es = slice(k1,ep)
  ss = slice(k1,sm)
  plot2(s2,s3,fs,cmin=-2.0,cmax=2.0,png="seis2d"+str(k1))
  plot2(s2,s3,es,cmin=0.2,cmax=1.0)
  plot2(s2,s3,ss,cmin=0.2,cmax=1.0)
  plot3(ep,k1=k1,cmin=0.2,cmax=1.0)
  plot3(sm,k1=k1,cmin=0.2,cmax=1.0)
def goSub():
  fx = readImage(fxfile)
  ep = readImage(epfile)
  sem = readImage(semfile)
  fx  = copy(210,n2,n3,30,0,0,fx)
  ep  = copy(210,n2,n3,30,0,0,ep)
  sem = copy(210,n2,n3,30,0,0,sem)
  writeImage(fxfile,fx)
  writeImage(epfile,ep)
  writeImage(semfile,sem)
def goSeis():
  fx = readImage(fxfile)
  ep = readImage(epfile)
  plot3(fx,k1=k1,clab="Amplitude",png="fx"+str(k1))
  plot3(ep,k1=k1,cmin=0.1,cmax=1.0,clab="Planarity",png="eps"+str(k1))
def goSta():
  gx = readImage(gxfile)
  if not plotOnly:
    ep = zerofloat(n1,n2,n3)
    el = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(6,4)
    et = lof.applyForTensors(gx)
    sta = StructureTensorAttribute(et,20)
    sta.setEigenvalues(1.0,0.05,0.4)
    sta.applyForPlanarLinear(gx,ep,el)
    writeImage(epsfile,ep)
    writeImage(elsfile,el)
  else:
    ep = readImage(epsfile)
    #et = readTensors(etsfile)
  plot3f(gx)
  plot3f(ep,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="eps"+str(k1))

def goSlope():
  print "goSlope..."
  if not plotOnly:
    sig1,sig2=8,2
    fx = readImage(gxfile)
    p2 = zerofloat(n1,n2,n3)
    p3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lsf = LocalSlopeFinder(sig1,sig2,sig2,5)
    lsf.findSlopes(fx,p2,p3,ep)
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(eplfile,ep)
  else:
    ep = readImage(eplfile)
  plot3f(ep,cmin=0.2,cmax=1.0)
def goSemblance():
  print "goSemblance..."
  fx = readImage(fxfile)
  if not plotOnly:
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    cov = Covariance()
    em,es=cov.covarianceEigen(12,p2,p3,fx)
    sem = div(em,es)
    #writeImage(semfile,sem)
    writeImage("sem12",sem)
  else:
    sem = readImage(semfile)
  plot3(sem,k1=k1,cmin=0.1,cmax=1.0,clab="Semblance",cint=0.1,png="sem"+str(k1))

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

def slice(k1,fx):
  fs = zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      fs[i3][i2] = fx[i3][i2][k1]
  return fs
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

def plot3f(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,htgs=None,fbs=None,png=None):
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
      #ipg.setClips(-2.0,2.0)
      ipg.setClips(-2.0,1.5) # use for subset plots
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,1.5)
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
  if xyz:
    pg = PointGroup(0.2,xyz)
    ss = StateSet()
    cs = ColorState()
    cs.setColor(Color.YELLOW)
    ss.add(cs)
    pg.setStates(ss)
    #ss = StateSet()
    #ps = PointState()
    #ps.setSize(5.0)
    #ss.add(ps)
    #pg.setStates(ss)
    sf.world.addChild(pg)
  if cells:
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    cmap = ColorMap(0.0,1.0,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.7,cmap,cells,False)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if htgs:
    for htg in htgs:
      sf.world.addChild(htg)
  if fbs:
    mc = MarchingCubes(s1,s2,s3,fbs)
    ct = mc.getContour(0.0)
    tg = TriangleGroup(ct.i,ct.x,ct.u)
    states = StateSet()
    cs = ColorState()
    cs.setColor(Color.MAGENTA)
    states.add(cs)
    lms = LightModelState()
    lms.setTwoSide(True)
    states.add(lms)
    ms = MaterialState()
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setSpecular(Color.WHITE)
    ms.setShininess(100.0)
    states.add(ms)
    tg.setStates(states);
    sf.world.addChild(tg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setLocalViewer(True)
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    if links:
      size = 0.65 
      ls = LineState()
      ls.setWidth(1.5)
      ls.setSmooth(True)
      ss.add(ls)
    ct = 0
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(-smax,smax,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForThrow(size,cmap,False)
      else: # show fault likelihood
        cmap = ColorMap(0.15,0.8,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,False)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
      if curve or trace:
        cell = skin.getCellNearestCentroid()
        if curve:
          xyz = cell.getFaultCurveXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
        if trace:
          xyz = cell.getFaultTraceXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
      if links:
        if ct==0:
          r,g,b=0,0,0
        if ct==1:
          r,g,b=0,0,1
        if ct==2:
          r,g,b=0,1,1
        if ct==3:
          #r,g,b=0.627451,0.12549,0.941176
          r,g,b=1,1,1
        r,g,b=0,0,1
        xyz = skin.getCellLinksXyz()
        #rgb = skin.getCellLinksRgb(r,g,b,xyz)
        #lg = LineGroup(xyz,rgb)
        lg = LineGroup(xyz)
        sg.addChild(lg)
        #ct = ct+1
    sf.world.addChild(sg)
  ipg.setSlices(106,80,207)
  ipg.setSlices(56,25,436)
  #ipg.setSlices(115,25,167)
  if cbar:
    sf.setSize(987,720)
  else:
    sf.setSize(850,720)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.3*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.4)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.20,0.12))
  ov.setAzimuthAndElevation(-56.0,35.0)
  #ov.setAzimuthAndElevation(-56.0,40.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3s(f,g=None,et=None,ep=None,k1=120,
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
  ipg.setSlices(k1,425,5)
  if cbar:
    sf.setSize(1351,600)
  else:
    sf.setSize(950,650)
  view = sf.getOrbitView()
  zscale = 0.25*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(2.7)
  view.setAzimuth(160.0)
  view.setElevation(35)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.02,-0.10,0.1))
  sf.viewCanvas.setBackground(sf.getBackground())
  #sf.setSize(850,700)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

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
  ipg.setSlices(k1,755,736)
  if cbar:
    sf.setSize(1051,750)
  else:
    sf.setSize(950,750)
  view = sf.getOrbitView()
  zscale = 0.25*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(2.0)
  view.setAzimuth(248.0)
  view.setElevation(35)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  view.setTranslate(Vector3(0.12,-0.02,0.1))
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
  #orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  orientation = PlotPanel.Orientation.X1RIGHT_X2UP;
  panel = PlotPanel(1,1,orientation,PlotPanel.AxesPlacement.NONE)
  #panel.setVInterval(0.1)
  #panel.setHInterval(1.0)
  #panel.setHLabel("Crossline (traces)")
  #panel.setVLabel("Samples")
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
  frame.setSize(700,325) #for poseidon
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
