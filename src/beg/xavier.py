from utils import *
#setupForSubset("bpSub1")
setupForSubset("bahamas")
#setupForSubset("australia")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
# Names and descriptions of image files used below.
pngDir = "../../../png/sso/3d/sta/"
pngDir = None

gxfile = "gx"
p2file = "p2"
p3file = "p3"
u1file = "u1"
u2file = "u2"
u3file = "u3"
epfile = "ep"
epsfile = "eps"
hvsfile = "hvs"
plotOnly = True
k1 = 51
k1 = 59

def main(args):
  #goSta()
  #goNormals()
  #goSlopes()
  goHorizonS()

def goSta():
  if not plotOnly:
    gx = readImage(gxfile)
    ep = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(3,4,4)
    et = lof.applyForTensors(gx)
    sta = StructureTensorAttribute(et,15)
    sta.setEigenvalues(0.6,0.01,1.0)
    sta.applyForPlanar(gx,ep)
    writeImage(epsfile,ep)
  else:
    ep = readImage(epsfile)
  '''
  ep = pow(ep,2)
  ep = sub(ep,min(ep))
  ep = div(ep,max(ep))
  '''
  ep = mask(ep)
  writeImage(epsfile,ep)
  print min(ep)
  print max(ep)
  plot3(ep,k1=k1,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1,png="eps"+str(k1))
  
def goNormals():
  gx = readImage(gxfile)
  if not plotOnly:
    lof = LocalOrientFilter(8,2)
    etl = lof.applyForTensors(gx)
    loe = LocalOrientEstimator(etl,10)
    loe.setGradientSmoothing(0)
    loe.setEigenvalues(1.0,0.2,0.2)
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    loe.applyX(gx,None,None,u1,u2,u3,None,None,None,None,None,None,None,None,None,None,None)
    writeImage(u1file,u1)
    writeImage(u2file,u2)
    writeImage(u3file,u3)
  plot3(gx)

def goSlopes():
  gx = readImage(gxfile)
  u1 = readImage(u1file)
  u2 = readImage(u2file)
  u3 = readImage(u3file)
  hp = Helper()
  p2,p3=hp.slopesFromNormals(10,u1,u2,u3)
  zm = ZeroMask(0.1,4,1,1,gx)
  zm.setValue(0.0,p2)
  zm.setValue(0.0,p3)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  plot3(p2,cmin=0.2,cmax=1)

def goHorizonS():
  ns = 2
  ep = readImage(epsfile)
  if not plotOnly:
    ep = fillfloat(1,n1,n2,n3)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    c1 = rampfloat(45,5,ns)
    c2 = fillfloat(1000,ns)
    c3 = fillfloat(2000,ns)
    hv = HorizonVolume()
    hv.setCG(0.01,50)
    hs = hv.applyForHorizonVolume(c1,c2,c3,ep,p2,p3)
    writeImage(hvsfile,hs)
  else:
    hs = readHorizons(ns,hvsfile)
  plot3(ep,surf=hs[0],cmin=0.2,cmax=1.0,png="sf0")
  plot3(ep,surf=hs[1],cmin=0.2,cmax=1.0,png="sf1")

def mask(ep,mv):
  gx = readImage(gxfile)
  zm = ZeroMask(0.1,4,1,1,gx)
  zm.setValue(mv,ep)
  return ep

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

def plot3(f,g=None,et=None,ep=None,surf=None,k1=120,
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
  if surf:
    sd = SurfaceDisplay()
    xyz,rgb = sd.horizonWithAmplitude([cmin,cmax],surf,f)
    tgs = TriangleGroup(True,xyz,rgb)
    sf.world.addChild(tgs)
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
