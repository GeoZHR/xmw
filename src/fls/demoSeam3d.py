"""
Demo of salt likelihoods and salt boundaries
Author: Xinming Wu, Colorado School of Mines
Version: 2015.12.19
"""

from utils import *

setupForSubset("seamSub")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # input image 
gefile  = "ge" # envelope image
slfile  = "sl" # salt likelihood image
gsfile  = "gs" # combined envelope and salt likelihood
sffile  = "sf" # salt indicator function
wx1file = "wx1"
wx2file = "wx2"
ws1file = "ws1"
ws2file = "ws2"
fp1file = "fp1"
fp2file = "fp2"
dxsfile = "dxs"
fsfile = "fs"
psfile = "ps"
psdfile = "psd"

pngDir = "../../../png/fls/seam/sub3d/"
pngDir = False

plotOnly = False

def main(args):
  #goEnvelope()
  #goSaltLike()
  #goEnvAndSaltLike()
  goPicker3()
def goPicker3():
  gx = readImage(gxfile)
  plot3(gx,png="seis")
  gs = readImage(gsfile)
  gm = max(gs)*0.4
  for i2 in range(n2):
    for i1 in range(n1):
      gs[0][i2][i1] = gm
  for i3 in range(n3):
    for i2 in range(n2):
      gs[i3][i2][n1-1] = gm
  zs,ys,xs=getPiks()
  fs = zerofloat(n1,n2,n3)
  sp3 = SaltPicker3()
  pks = sp3.pick3(25,xs,ys,zs,gs,fs)
  lgs = getLineGroups(1,zs,ys,xs)
  rgf1 = RecursiveGaussianFilterP(3)
  rgf2 = RecursiveGaussianFilterP(6)
  rgf3 = RecursiveGaussianFilterP(3)
  rgf1.apply0XX(fs,fs)
  rgf2.applyX0X(fs,fs)
  rgf3.applyXX0(fs,fs)
  plot3(gx,g=fs,cmin=-0.5,cmax=0.5,png="saltBody")
  plot3(gx,fbs=fs,png="saltBound")
  '''
  gx2 = zerofloat(n1,n3)
  fs2 = zerofloat(n1,n3)
  for i2 in range(0,n2,50):
    for i3 in range(n3):
      gx2[i3] = gx[i3][i2]
      fs2[i3] = fs[i3][i2]
    plot2(gx2)
    plot2(fs2)
    writeImage("salt"+str(i2),fs2)

  for i2 in range(399,n2,1):
    for i3 in range(n3):
      gx2[i3] = gx[i3][i2]
      fs2[i3] = fs[i3][i2]
    plot2(gx2)
    plot2(fs2)
    writeImage("salt"+str(i2),fs2)
    #plot2(gx2,xp=pks[i2],clab=str(i2))
  '''
  '''
  plot3(gs,cmin=0.2,cmax=0.6,lgs=lgs,png="slicesInitial")
  xps,yps,zps=[],[],[]
  xrs,yrs,zrs=[],[],[]
  for i2 in range(0,n2,50):
    yrs.append(i2)
    xrs.append(pks[i2][1])
    zrs.append(pks[i2][0])
  yrs.append(n2-1)
  xrs.append(pks[n2-1][1])
  zrs.append(pks[n2-1][0])
  lrs = getLineGroups(1,zrs,yrs,xrs)
  plot3(gx,lgs=lrs,png="slicesFinal")
  '''


def goSaltLike():
  gx = readImage(gxfile)
  lof = LocalOrientFilter(4,3,3)
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  e1 = zerofloat(n1,n2,n3)
  e2 = zerofloat(n1,n2,n3)
  e3 = zerofloat(n1,n2,n3)
  lof.applyForNormalPlanar(gx,u1,u2,u3,ep)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        if (isNaN(ep[i3][i2][i1])):
          ep[i3][i2][i1] = 1
  rgf = RecursiveGaussianFilterP(4)
  rgf.apply100(ep,e1)
  rgf.apply010(ep,e2)
  rgf.apply001(ep,e3)
  e1 = mul(e1,e1)
  e2 = mul(e2,e2)
  e3 = mul(e3,e3)
  es = add(e1,e2)
  es = add(es,e3)
  es = normalize(es)
  writeImage(slfile,es)
  plot3(es,cmin=0.2,cmax=0.6)

def goEnvAndSaltLike():
  ge = readImage(gefile)
  sl = readImage(slfile)
  ge = mul(ge,0.6)
  sl = mul(sl,0.4)
  gs = add(ge,sl)
  gs = normalize(gs)
  writeImage(gsfile,gs)
  plot3(ge,cmin=0.2,cmax=0.6)
  plot3(sl,cmin=0.2,cmax=0.6)
  plot3(gs,cmin=0.2,cmax=0.6)

def goEnvelope():
  gx = readImage(gxfile)
  rgf = RecursiveGaussianFilter(1)
  rgf.apply000(gx,gx)
  ge = zerofloat(n1,n2,n3)
  FastLevelSet3.applyForInsAmp(gx,ge)
  ge = normalize(ge)
  writeImage(gefile,ge)
  plot3(gx,cmin=min(gx)/2,cmax=max(gx)/2)
  plot3(ge,cmin=0.2,cmax=0.6)


def isNaN(num):
    return num != num

def normalize(x):
  x = sub(x,min(x))
  return div(x,max(x))


def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

def getPiks():
  y1 = [  0, 42, 22,52,281,149,122, 77,  0,  0]
  z1 = [302,197,155,72,231,359,512,353,359,302]

  y2 = [  0,82,222,380,196,130, 35,  0,  0] 
  z2 = [316,30,109,108,330,510,325,355,316]

  y3 = [ 15, 69,49,390,213,154, 98, 41, 15]
  z3 = [319,159,28,129,330,549,349,314,319]

  y4 = [ 28, 68,47,390,304,203,149, 86,60, 28]
  z4 = [306,155,22,113,311,326,549,296,288,306]

  y5 = [ 26, 76,55,48,94,272,379,238,174,139, 74, 26]
  z5 = [298,150,60,42,20,101, 44,361,374,416,262,298]

  y6 = [ 38, 81,65,72,268,319,370,351,169, 88, 61, 38]
  z6 = [271,151,71,33, 75, 21, 91,239,288,280,263,271]

  y7 = [ 47, 90,86,284,306,364,230, 47]
  z7 = [294,110,38, 82, 54,285,253,294]

  y8 = [  0, 91,117,105,133,280,312,358,227, 89,  0,  0]
  z8 = [376,140,107, 49, 35, 81, 49,308,241,312,400,376]

  y9 = [  0,114,111,293,322,338,305,386,224,149, 91,  0,  0]
  z9 = [366,129, 34, 62, 10, 24, 88,315,255,355,340,453,366]
  zs = [z1,z2,z3,z4,z5,z6,z7,z8,z9]
  ys = [0,50,100,150,200,250,300,350,399]
  xs = [y1,y2,y3,y4,y5,y6,y7,y8,y9]
  return zs,ys,xs

def getLineGroups(dx,zs,ys,xs):
  lgs = []
  for ic in range(0,len(xs),dx):
    xyz = []
    rgb = []
    for ip in range(len(xs[ic])):
      xyz.append(xs[ic][ip])
      xyz.append(ys[ic])
      xyz.append(zs[ic][ip])
      rgb.append(1)
      rgb.append(0)
      rgb.append(0)
    lg = LineGroup(xyz,rgb)
    lgs.append(lg)
  return lgs

#############################################################################
# graphics
gray = ColorMap.GRAY
jet = ColorMap.JET

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def rwbFill(alpha):
  return ColorMap.setAlpha(ColorMap.RED_WHITE_BLUE,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)

def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrRamp(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,rampfloat(0.0,alpha/256,256))
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

def plot2(f,xp=None,ps=None,xs=None,xu=None,nr=50,phi=None,v1=None,v2=None,
        cmin=None,cmax=None,w1=None,w2=None,clab=None,png=None): 
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  #panel.setVInterval(0.2)
  n2 = len(f)
  n1 = len(f[0])
  s2 = Sampling(n2)
  s1 = Sampling(n1)
  panel.setHLabel("Inline (traces)")
  panel.setVLabel("Time (samples)")

  panel.setHLimits(0,0,n2-1)
  panel.setVLimits(0,0,n1-1)
  pxv = panel.addPixels(0,0,s1,s2,f);
  #pxv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pxv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pxv.setColorModel(ColorMap.GRAY)
  if cmin and cmax:
    pxv.setClips(cmin,cmax)
  else:
    pxv.setClips(min(f),max(f))
  #panel.setTitle("normal vectors")
  if phi:
    cv = panel.addContours(phi)
    cv.setContours([0])
    cv.setLineColor(Color.RED)
    cv.setLineWidth(1.0)
  if xp:
    ptv = panel.addPoints(0,0,xp[0],xp[1])
    ptv.setLineColor(Color.RED)
    ptv.setLineWidth(3.0)
  if ps:
    ptv1 = panel.addPoints(0,0,ps[0][0],ps[0][1])
    ptv2 = panel.addPoints(0,0,ps[1][0],ps[1][1])
    ptv3 = panel.addPoints(0,0,ps[2][0],ps[2][1])
    ptv1.setLineStyle(PointsView.Line.NONE)
    ptv2.setLineStyle(PointsView.Line.NONE)
    ptv3.setLineStyle(PointsView.Line.NONE)
    ptv1.setMarkStyle(PointsView.Mark.PLUS)
    ptv2.setMarkStyle(PointsView.Mark.PLUS)
    ptv3.setMarkStyle(PointsView.Mark.PLUS)
    ptv1.setMarkSize(5)
    ptv2.setMarkSize(5)
    ptv3.setMarkSize(5)
    ptv1.setMarkColor(Color.WHITE)
    ptv2.setMarkColor(Color.RED)
    ptv3.setMarkColor(Color.BLUE)
  if xu:
    np = len(xu[0])
    for ip in range(np):
      x1c = xu[0][ip]
      x2c = xu[1][ip]
      u1c = xu[2][ip]
      u2c = xu[3][ip]
      x1m = x1c-u1c*nr
      x2m = x2c-u2c*nr
      x1p = x1c+u1c*nr
      x2p = x2c+u2c*nr
      x1s = [x1m,x1c,x1p]
      x2s = [x2m,x2c,x2p]
      ptv2 = panel.addPoints(x1s,x2s)
      ptv2.setLineColor(Color.YELLOW)
      ptv2.setLineWidth(1.0)
    if xp:
      ptv = panel.addPoints(0,0,xp[0],xp[1])
      ptv.setLineColor(Color.RED)
      ptv.setLineWidth(3.0)
  if xs:
    for ip in range(len(xs)):
      ptv = panel.addPoints(xs[ip][0],xs[ip][1])
      ptv.setLineStyle(PointsView.Line.NONE)
      ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      ptv.setMarkColor(Color.RED)
      ptv.setMarkSize(2.0)
  if(clab):
    cb = panel.addColorBar();
    cb.setLabel(clab)
  panel.setColorBarWidthMinimum(130)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True);
  if w1 and w2:
    frame.setSize(w2,w1)
  else:
    frame.setSize(round(n2*0.8),round(n1*0.8))
  #frame.setSize(1190,760)
  frame.setFontSize(14)
  if pngDir and png:
    frame.paintToPng(720,3.333,pngDir+png+".png")


def plot3(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          lgs=None,cells=None,tg=None, fbs=None,png=None):
  n1 = len(f[0][0])
  n2 = len(f[0])
  n3 = len(f)
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  cbar = None
  k1 = [280]
  k3 = [62,320]
  k2 = [9,399]
  if g==None:
    ipg = sf.addImagePanelsX(s1,s2,s3,f,k1,k2,k3)
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
    axes = [Axis.X,Axis.Y]
    ipg = ImagePanelGroup2X(s1,s2,s3,f,g,k1,k2,k3)
    ipg.setClips1(-1.5,1.5)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = rwbFill(0.4)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if lgs:
    for lg in lgs:
      ss = StateSet()
      lg.setStates(ss)
      ls = LineState()
      ls.setWidth(8)
      ls.setSmooth(False)
      ss.add(ls)
      sf.world.addChild(lg)
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
  if tg:
    states = StateSet()
    cs = ColorState()
    cs.setColor(Color.MAGENTA)
    #cs.setColor(Color.ORANGE)
    #cs.setColor(Color.CYAN)
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
  if fbs:
    mc = MarchingCubes(s1,s2,s3,fbs)
    ct = mc.getContour(0.0)
    xyz = ct.x
    np = len(xyz)/3
    zs = zerofloat(np)
    for ip in range(np):
      zs[ip] = -xyz[ip*3+2]
    cp = ColorMap(-n1-1,0,ColorMap.JET);
    rgb = cp.getRgbFloats(zs)
    tg = TriangleGroup(ct.i,ct.x,ct.u,rgb)
    states = StateSet()
    cs = ColorState()
    #cs.setColor(Color.MAGENTA)
    #cs.setColor(Color.ORANGE)
    #cs.setColor(Color.CYAN)
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
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.8*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.02,-0.10,-0.05))
  ov.setAzimuthAndElevation(140.0,45.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
