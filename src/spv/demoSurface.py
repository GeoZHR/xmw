#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("surf")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
#############################################################################
gxfile = "gx" # input semblance image
epfile = "ep"  # planarity
effile = "ef"  # 1-planarity
fefile = "fe"  # 1-planarity
flfile = "fl"  # fault likelihood
fpfile = "fp"  # fault strike;
ftfile = "ft"  # fault dip;
fvfile = "fv"  # fault dip;
vpfile = "vp"  # fault dip;
vtfile = "vt"  # fault dip;
fetfile = "fet" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fskfile = "skin"

pngDir = None
pngDir = getPngDir()
plotOnly = False
# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 65,80
minPhi,maxPhi = 0,360
sigmaPhi,sigmaTheta=3,6

def main(args):
  #goUvxBox()
  #goPlanar()
  goSurface()
  #goTranspose()
def goUvxBox():
  gx = readImage3D("f3draw")
  c1,c2,c3=230,632,162
  u = [0,-0.5*sqrt(2),-0.5*sqrt(2)]
  v = [1,0,0]
  w = [0, 0.5*sqrt(2),-0.5*sqrt(2)]
  os = OptimalSurfacer(50,150,50)
  gb = os.getUvwBox(c1,c2,c3,u,v,w,gx)
  gx = os.transpose21(gb)
  gx = gain(gx)
  gx = copy(155,101,101,145,0,0,gx)
  writeImage(gxfile,gx)
  plot3(gx)
def goPlanar():
  gx = readImage3D(gxfile)
  gx = gain(gx)
  if not plotOnly:
    lof = LocalOrientFilter(2,1,1)
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
    writeImage(epfile,ep)
  else:
    ep = readImage3D(epfile)
  plot3(gx)
  ep = pow(ep,8)
  writeImage("gx1",gx[1])
  plot3(gx,sub(1,ep),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-planarity",png="ep")
def goSurface():
  gx = readImage3D(gxfile)
  ep = readImage3D(epfile)
  ep = pow(ep,6)
  gx = gain(gx)
  os = OptimalSurfacer(50,77,50)
  os.setStrainMax(0.35,0.35)
  os.setAttributeSmoothing(1)
  os.setSurfaceSmoothing(2,1)
  gxt = os.transpose21(gx)
  ept = os.transpose21(ep)
  c1,c2,c3=55,74,33
  os.setControlPoint(c1,c2,c3,ept)
  epc = copy(ept)
  plot3(gxt,sub(1,ept),k1=100,k2=c2,k3=c3,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-planarity",png="ept")
  sft = os.findSurface(c1,c2,c3,ept)
  xyz = os.getSurfaceTriangles(s3,s1,sft)
  xzy = os.getSurfaceTrianglesT(s3,s1,sft)
  plot3(gxt,sub(1,epc),k1=100,k2=c2,k3=c3,surf=xyz,cmin=0.25,cmax=1.0,
          cmap=jetRamp(1.0), clab="1-planarity",png="epts")
  m3 = len(epc)
  m2 = len(epc[0])
  m1 = len(epc[0][0])
  epa1 = zerofloat(m1,m2,m3)
  epa2 = zerofloat(m1,m2,m3)
  os.smoothAttribute1(epc, epa1)
  os.setControlPoint(c1,c2,c3,epa1)
  os.smoothAttribute2(epa1,epa2)
  os.setControlPoint(c1,c2,c3,epa2)
  plot3(gxt,sub(1,epa1),k1=100,k2=c2,k3=c3,cmin=0.25,cmax=1.0,
          cmap=jetRamp(1.0), clab="1-planarity",png="epa1")
  plot3(gxt,sub(1,epa2),k1=100,k2=c2,k3=c3,cmin=0.25,cmax=1.0,
          cmap=jetRamp(1.0), clab="1-planarity",png="epa2")
  plot3(gxt,sub(1,epa2), surf=xyz,k1=100,k2=c2,k3=c3,cmin=0.25,cmax=1,
          cmap=jetRamp(1.0),clab="1-planarity",png="epa2s")

  k1,k2,k3=258,48,72
  plot3x(gx,sub(1,ep), surf=xzy,k1=k1,k2=k2,k3=k3,cmin=0.25,cmax=1,
          cmap=jetRamp(1.0),clab="1-planarity",png="eps")
  plot3x(gx,surf=xzy,k1=k1,k2=k2,k3=k3,png="gxs")
  plot3x(gx,sub(1,ep), k1=k1,k2=k2,k3=k3,cmin=0.25,cmax=1,
          cmap=jetRamp(1.0),clab="1-planarity",png="ep")
  plot3x(gx,k1=k1,k2=k2,k3=k3,png="gx")
  epf = os.planarityOnFault(sub(1,epc),sft)
  xzy,rgb = os.getSurfaceTrianglesT(0.25,1.0,s3,s1,sft,epf)
  tg = TriangleGroup(True,xzy,rgb)
  plot3x(gx,tg=tg,k1=k1,k2=k2,k3=k3,png="epps")
  rgf = RecursiveGaussianFilterP(6)
  epfs = zerofloat(m2,m3)
  rgf.apply00(epf,epfs)
  xzy,rgbs = os.getSurfaceTrianglesT(0.25,1.0,s3,s1,sft,epfs)
  tgs = TriangleGroup(True,xzy,rgbs)
  plot3x(gx,tg=tgs,k1=k1,k2=k2,k3=k3,png="eppss")
  epft = zerofloat(m1,m2,m3)
  os.surfaceToImage(sft,epfs,epft)
  epf = os.transpose21(epft)
  plot3x(gx,epf, k1=k1,k2=k2,k3=k3,cmin=0.25,cmax=1,
          cmap=jetRamp(1.0),clab="1-planarity",png="epf")

def gain(x):
  m3 = len(x)
  m2 = len(x[0])
  m1 = len(x[0][0])
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = zerofloat(m1,m2,m3)
  div(x,sqrt(g),y)
  return y

def smooth(sig,u):
  v = copy(u)
  rgf = RecursiveGaussianFilterP(sig)
  rgf.apply0(u,v)
  return v

def smooth2(sig1,sig2,u):
  v = copy(u)
  rgf1 = RecursiveGaussianFilterP(sig1)
  rgf2 = RecursiveGaussianFilterP(sig2)
  rgf1.apply0X(u,v)
  rgf2.applyX0(v,v)
  return v


def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

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

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

def plot3(f,g=None,surf=None,cmin=-2,cmax=2,cmap=None,clab=None,cint=None,
          k1=50,k2=150,k3=99,png=None):
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
    ipg.setClips1(-2,2)
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
  if surf:
    tg = TriangleGroup(True,surf)
    tg.setColor(Color.CYAN)
    sf.world.addChild(tg)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(987,720)
  else:
    sf.setSize(850,720)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 1
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.55)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.06,-0.05))
  ov.setAzimuthAndElevation(45,36)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")

def plot3x(f,g=None,surf=None,tg=None,cmin=-2,cmax=2,cmap=None,clab=None,cint=None,
          k1=50,k2=150,k3=99,png=None):
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
    ipg.setClips1(-2,2)
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
  if surf:
    tg = TriangleGroup(True,surf)
    tg.setColor(Color.CYAN)
    sf.world.addChild(tg)
  if tg:
    sf.world.addChild(tg)
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(687,800)
  else:
    sf.setSize(550,800)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 1
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.06,-0.05))
  ov.setAzimuthAndElevation(-45,36)
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
