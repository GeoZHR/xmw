#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("f3d")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
#############################################################################
gxfile = "xs" # input semblance image
epfile = "ep"  # planarity
effile = "ef"  # 1-planarity
fefile = "fe"  # 1-planarity
fpfile = "fp"  # fault strike;
ftfile = "ft"  # fault dip;
fvfile = "fv"  # fault dip;
ftfile  = "ft" # fault dip (theta)
fetfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned

pngDir = getPngDir()
pngDir = None
plotOnly = True
# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 65,80
minPhi,maxPhi = 0,360
sigmaPhi,sigmaTheta=3,6

def main(args):
  #goPlanar()
  #goFaultOrientScan()
  #goSurfaceVoting()

def goPlanar():
  gx = readImage3D(gxfile)
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
  plot3(gx,sub(1,ep),cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-planarity",png="fl")

def goFaultOrientScan():
  gx = readImage3D(gxfile)
  ep = readImage3D(epfile)
  fos = FaultOrientScanner3(sigmaPhi,sigmaTheta)
  if not plotOnly:
    fe,fp,ft = fos.scan(minPhi,maxPhi,minTheta,maxTheta,ep)
    fet,fpt,ftt=fos.thin([fe,fp,ft])
    writeImage(fefile,fe)
    writeImage(fetfile,fet)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    fe = readImage3D(fefile)
    fp = readImage3D(fpfile)
  print min(fe) 
  print max(fe) 
  ep = sub(1,pow(ep,8))
  plot3(gx,ep,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Planarity",png="ep")
  plot3(gx,fe,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Ray tracing",png="fl")
  plot3(gx,fp,cmin=1,cmax=360,cmap=jetFill(1.0),
      clab="Fault strike (degrees)",png="ph")

def goSurfaceVoting():
  gx = readImage3D(gxfile)
  gx = gain(gx)
  if not plotOnly:
    fet = readImage3D(fetfile)
    fpt = readImage3D(fptfile)
    ftt = readImage3D(fttfile)
    osv = OptimalSurfaceVoterP(10,20,30)
    osv.setStrainMax(0.25,0.25)
    osv.setSurfaceSmoothing(2,2)
    fv = osv.applyVoting(4,0.3,fet,fpt,ftt)
    writeImage(fvfile,fv)
  else:
    fv = readImage3D(fvfile)
  ep = readImage3D(epfile)
  ep = sub(1,pow(ep,8))
  plot3(gx)
  plot3(gx,ep,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Planarity",png="ep")
  plot3(gx,fv,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Surface voting",png="fv")
  '''
  plot3(gx,fvp,cmin=1,cmax=180,cmap=jetFill(1.0),
      clab="Fault strike (degrees)",png="ph")
  '''

def goRayTracing2():
  gx = readImage3D(gxfile)
  fp = readImage3D(fpfile)
  if not plotOnly:
    fe = FaultEnhance(4,0.2)
    se = fe.applyTracing(80,20,minTheta,maxTheta,fp)
    writeImage(sefile,se)
  else:
    se = readImage3D(sefile)
  plot3(gx,fp,cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")
  plot3(gx,se,cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="Fault likelihood",png="fl")

def goRayTracing21():
  gx = readImage3D(gxfile)
  fp = readImage3D(sefile)
  ep = readImage3D(fpfile)
  if not plotOnly:
    fp = mul(fp,sub(2,fp))
    fp = sub(fp,min(fp))
    fp = div(fp,max(fp))
    fe = FaultEnhance(10,0.5)
    se = fe.applyTracing1(80,15,minTheta,maxTheta,fp)
    writeImage(sefile+"21",se)
  else:
    se = readImage3D(sefile+"21")
  plot3(gx,ep,cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="Plaanrity",png="fl")
  plot3(gx,fp,cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="Ray tracing",png="fl")
  plot3(gx,se,cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="Ray tracing",png="fl")
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
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

def plot3(f,g=None,cmin=-2,cmax=2,cmap=None,clab=None,cint=None,
          sx=None,sy=None,surf=None,fd=None,cells=None,skins=None,png=None):
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
  if surf:
    tg = TriangleGroup(True,surf)
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
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    sf.world.addChild(sg)
  ipg.setSlices(106,80,207)
  ipg.setSlices(93,25,32)
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
  ov.setScale(1.55)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,0.06,-0.05))
  ov.setAzimuthAndElevation(45.0,36.0)
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
