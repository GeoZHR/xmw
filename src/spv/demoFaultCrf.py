#############################################################################
"""
Demo of dynamic warping for automatic picking
Author: Xinming Wu, University of Texas at Austin
Version: 2016.06.01
"""


from utils import * 
setupForSubset("crf")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
f1,f2,f3 = s1.getFirst(),s2.getFirst(),s3.getFirst()
d1,d2,d3 = s1.getDelta(),s2.getDelta(),s3.getDelta()
#############################################################################
gxfile = "gx" # input semblance image
epfile = "ep"  # planarity
effile = "ef"  # 1-planarity
fefile = "fe"  # 1-planarity
fpfile = "fp"  # fault strike;
ftfile = "ft"  # fault dip;
fvfile = "fv"  # fault dip;
ftfile  = "ft" # fault dip (theta)
fetfile = "fet" # fault likelihood thinned
fptfile = "fpt" # fault likelihood thinned
fttfile = "ftt" # fault dip thinned


pngDir = getPngDir()
pngDir = None
plotOnly = False
# These parameters control the scan over fault strikes and dips.
# See the class FaultScanner for more information.
minTheta,maxTheta = 65,80
minPhi,maxPhi = 0,360
sigmaPhi,sigmaTheta=4,8


def main(args):
  #goPlanarX()
  goFaultOrientScan()
  #goSurfaceVoting()
def goPlanarX():
  gx = readImage3D(gxfile)
  if not plotOnly:
    lof = LocalOrientFilter(16,4)
    et3 = lof.applyForTensors(gx)
    et3.setEigenvalues(1.0,0.01,0.5)
    fer = FaultEnhancer(sigmaPhi,sigmaTheta)
    ep = fer.applyForPlanar(10,et3,gx)
    writeImage(epfile,ep)
    print min(ep)
    print max(ep)
  else:
    ep = readImage3D(epfile)
  #plot3(gx,cmin=-3,cmax=3)
  plot3(ep,cmin=0.2,cmax=1.0,clab="Planarity",cint=0.1)
  #plot3(gx,sub(1,ep),cmin=0.1,cmax=0.8,cmap=jetRamp(1.0),
  #    clab="1-planarity",png="fl")

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
  plot3(gx,sub(1,ep),cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="1-planarity",png="ep")

def goFaultOrientScan():
  gx = readImage3D(gxfile)
  ep = readImage3D(epfile)
  fos = FaultOrientScanner3(sigmaPhi,sigmaTheta)
  if not plotOnly:
    fe,fp,ft = fos.scan(minPhi,maxPhi,minTheta,maxTheta,ep)
    fet,fpt,ftt=fos.thin([fe,fp,ft])
    writeImage(fefile,fe)
    writeImage(fpfile,fp)
    writeImage(fetfile,fet)
    writeImage(fptfile,fpt)
    writeImage(fttfile,ftt)
  else:
    fp = readImage3D(fpfile)
    fe = readImage3D(fefile)
  print min(fe) 
  print max(fe) 
  plot3(gx,ep,cmin=0.1,cmax=0.7,cmap=jetRamp(1.0),
      clab="1-planarity",png="ep")
  plot3(gx,fe,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Enhanced",png="fe")
  plot3(gx,fp,cmin=1,cmax=360,cmap=jetFill(1.0),
      clab="Fault strike (degrees)",png="fp")

def goSurfaceVoting():
  gx = readImage3D(gxfile)
  if not plotOnly:
    fet = readImage3D(fetfile)
    fpt = readImage3D(fptfile)
    ftt = readImage3D(fttfile)
    osv = OptimalSurfaceVoterP(10,30,20)
    osv.setStrainMax(0.2,0.2)
    osv.setSurfaceSmoothing(2,2)
    fv = osv.applyVoting(4,0.3,fet,fpt,ftt)
  else:
    fv = readImage3D(fvfile)
  ep = readImage3D(epfile)
  ep = sub(1,pow(ep,8))
  plot3(gx,ep,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="1-planarity",png="ep")
  plot3(gx,fv,cmin=0.25,cmax=1.0,cmap=jetRamp(1.0),
      clab="Surface voting",png="sv")

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2)
  div(x,sqrt(g),y)
  return y
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
  cbar.setFont(Font("Arial",Font.PLAIN,24)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

def plot3(f,g=None,k1=120,cmin=-2,cmax=2,cmap=None,clab=None,cint=None,
          sx=None,sy=None,horizon=None,fd=None,cells=None,skins=None,png=None):
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
    cbar.setWidthMinimum(85)
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
  if horizon and not fd:
    hp = Helper()
    print "fvalues"
    print min(f)
    print max(f)
    ts = hp.horizonWithAmplitude(s1,s2,s3,s1,s2,s3,[cmin,cmax],horizon,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if horizon and fd:
    hp = Helper()
    ts = hp.horizonWithFaultDensity(s1,s2,s3,[0.0,0.15],horizon,fd)
    tg = TriangleGroup(True,ts[0],ts[1])
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
