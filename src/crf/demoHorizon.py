"""
Demonstrate 3D seismic image processing for faults
Author: Xinming Wu, Colorado School of Mines
Version: 2016.01.22
"""

from utils import *
setupForSubset("nathanSub8")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count
# Names and descriptions of image files used below.
gxfile  = "gx" # input image (maybe after bilateral filtering)
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
hl1file = "hl1"

pngDir = "../../../png/beg/nathan/sub8/"
pngDir = None
plotOnly = True

# Processing begins here. When experimenting with one part of this demo, we
# can comment out earlier parts that have already written results to files.
def main(args):
  goPickedSurfaces()
  #goSlopes()
  #goHorizonL1()

def goPickedSurfaces():
  fn = "su1"
  gx = readImage(gxfile)
  if not plotOnly:
    nl1 = 1646796
    nu1 = 1330218
    nm1 =  888165
    #sp = readImage2D(nl1,3,"hzs/L1")
    sp = readImage2D(nu1,3,"hzs/U1");  dmax=100
    #sp = readImage2D(nm1,3,"hzs/M1"); dmax=200
    print min(sp[0])
    print max(sp[0])
    print min(sp[2])
    print max(sp[2])
    hp = Helper()
    ndfs = zerofloat(3,2)
    sf = hp.surfaceResample(s2,s3,1.5,dmax,sp,ndfs)
    sy = Sampling(round(ndfs[0][0]),ndfs[0][1],ndfs[0][2])
    sx = Sampling(round(ndfs[1][0]),ndfs[1][1],ndfs[1][2])
    writeImage(fn,sf)
    writeImage(fn+"ndfs",ndfs)
  else:
    ndfs = readImage2D(3,2,fn+"ndfs")
    ny = round(ndfs[0][0])
    nx = round(ndfs[1][0])
    sf = readImage2D(ny,nx,fn)
    sy = Sampling(round(ndfs[0][0]),ndfs[0][1],ndfs[0][2])
    sx = Sampling(round(ndfs[1][0]),ndfs[1][1],ndfs[1][2])
    hp = Helper()
    x1,x2,x3=hp.controlPointsFromSurface(s1,s2,s3,sy,sx,sf)
  plot3(gx,cmin=-3,cmax=3,sx=sx,sy=sy,horizon=sf)

def goSlopes():
  print "goSlopes ..."
  gx = readImage(gxfile)
  sigma1,sigma2,sigma3,pmax = 8.0,3.0,3.0,5.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  writeImage(epfile,ep)
  print "p2  min =",min(p2)," max =",max(p2)
  print "p3  min =",min(p3)," max =",max(p3)
  print "ep min =",min(ep)," max =",max(ep)
  plot3(gx,p2, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Inline slope (sample/sample)",png="p2")
  plot3(gx,p3, cmin=-1,cmax=1,cmap=bwrNotch(1.0),
        clab="Crossline slope (sample/sample)",png="p3")
  plot3(gx,sub(1,ep),cmin=0,cmax=1,cmap=jetRamp(1.0),
        clab="Planarity")

def goHorizonL1():
  d1 = s1.getDelta();
  f1 = s1.getFirst();
  d2 = s2.getDelta();
  f2 = s2.getFirst();
  d3 = s3.getDelta();
  f3 = s3.getFirst();
  gx = readImage(gxfile)
  if not plotOnly:
    fn = "sl1"
    ndfs = readImage2D(3,2,fn+"ndfs")
    ny = round(ndfs[0][0])
    nx = round(ndfs[1][0])
    sf = readImage2D(ny,nx,fn)
    sy = Sampling(round(ndfs[0][0]),ndfs[0][1],ndfs[0][2])
    sx = Sampling(round(ndfs[1][0]),ndfs[1][1],ndfs[1][2])
    hp = Helper()
    x1,x2,x3=hp.controlPointsFromSurface(s1,s2,s3,sy,sx,sf)
    k1 = [1795,1795,1795,2495,1965,2193,2002,2298,2210,1924,2120,2001,2105,2417,2204]
    k2 = [6969,6969,6964,4122,4029,3839,3782,3821,3414,3631,3428,3344,3398,3707,3431]
    k3 = [2109,2195,2870,2765,2765,2771,2787,2201,2221,2718,2594,2594,2550,2128,2128]
    np = len(x1)
    nk = len(k1)
    c1 = zerofloat(np+nk)
    c2 = zerofloat(np+nk)
    c3 = zerofloat(np+nk)
    for ip in range(np):
      c1[ip] = x1[ip]
      c2[ip] = x2[ip]
      c3[ip] = x3[ip]
    for ik in range(nk):
      c1[ik+np] = (k1[ik]-f1)/d1
      c2[ik+np] = (k2[ik]-f2)/d2
      c3[ik+np] = (k3[ik]-f3)/d3
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    ep = readImage(epfile)

    ep = pow(ep,6)
    se = SurfaceExtractorC()
    se.setWeights(0.0)
    se.setCG(0.01,100)
    sf = se.surfaceInitializationFast(n2,n3,n1-1,c1,c2,c3)
    se.surfaceUpdateFromSlopes(ep,p2,p3,c1,c2,c3,sf)
    sf = add(sf,f1)
    sf = mul(sf,d1)
    writeImage(hl1file,sf)
  else:
    sf = readImage2D(n2,n3,hl1file)
    hp = Helper()
    hx = zerofloat(n1,n2,n3)
    hp.horizonToImage(s1,sf,hx)
  #plot3(gx,cmin=-3,cmax=3,sx=s3,sy=s2,horizon=sf)
  plot3(gx,hx,cmin=0,cmax=1,cmap=jetFillExceptMin(1.0))

def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(100.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y
def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))




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

def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
        sx=None,sy=None,horizon=None,fd=None,cells=None,skins=None,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  #s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
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
  if horizon and not fd:
    hp = Helper()
    print "fvalues"
    print min(f)
    print max(f)
    ts = hp.horizonWithAmplitude(s1,s2,s3,s1,sy,sx,[cmin,cmax],horizon,f)
    tg = TriangleGroup(True,ts[0],ts[1])
    sf.world.addChild(tg)
  if horizon and fd:
    hp = Helper()
    ts = hp.horizonWithFaultDensity(n1-2,[0.0,0.15],horizon,fd)
    tg = TriangleGroup(True,ts[0],ts[1])
    states = StateSet()
    cs = ColorState()
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
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    for skin in skins:
      cmap = ColorMap(0.2,0.8,ColorMap.JET)
      xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,False)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
    sf.world.addChild(sg)
  ipg.setSlices(150,5,56)
  #ipg.setSlices(85,5,43)
  #ipg.setSlices(85,5,102)
  #ipg.setSlices(n1,0,n3) # use only for subset plots
  ipg.setSlices(n1,376,308)
  if cbar:
    sf.setSize(1037,700)
  else:
    sf.setSize(900,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  #zscale = 1.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.5)
  #ov.setScale(2.5)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(0.0,-0.00,-0.05))
  ov.setAzimuthAndElevation(45.0,35.0)
  #ov.setAzimuthAndElevation(-55.0,35.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(720,1,pngDir+png+"cbar.png")


#############################################################################
run(main)
