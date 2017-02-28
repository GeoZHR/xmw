"""
Demo of salt likelihoods and salt boundaries
Author: Xinming Wu, Colorado School of Mines
Version: 2015.12.19
"""

from utils import *

setupForSubset("bag")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of image files used below.
gxfile  = "gx" # bag input image 
pafile  = "pa" # bag input image 
p2file  = "p2" # eigenvalue-derived planarity
p3file  = "p3" # eigenvalue-derived planarity
u1file  = "u1" # eigenvalue-derived planarity
u2file  = "u2" # eigenvalue-derived planarity
u3file  = "u3" # eigenvalue-derived planarity
gefile  = "ge" # eigenvalue-derived planarity
dpfile  = "dp" # eigenvalue-derived planarity
epfile  = "ep" # eigenvalue-derived planarity
slfile  = "sl" # eigenvalue-derived planarity
slefile  = "sle" # eigenvalue-derived planarity
sffile  = "sf" # salt indicator function
mkfile  = "mk" # mask file
phfile  = "ph"
phdxfile = "phdx"
sfcfile  = "sfc" # salt indicator function with constraints
flfile  = "fl" # fault likelihood
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
wxfile = "wx"
fsfile = "fs"

pngDir = False
pngDir = "../../../png/fls/bag/3d/"

plotOnly = True

def main(args):
  #goPikSlices()
  goSaltSurface()
def goSaltSurface():
  gx = readImage(gxfile)
  fs = zerofloat(n1,n2,n3)
  if not plotOnly:
    wx = readImage(wxfile)
    isr = ImplicitSurfaceReconstructor()
    isr.signAsignment(wx,fs)
    rgf1 = RecursiveGaussianFilter(3)
    rgf2 = RecursiveGaussianFilter(6)
    rgf3 = RecursiveGaussianFilter(12)
    rgf1.apply0XX(fs,fs)
    rgf2.applyX0X(fs,fs)
    rgf3.applyXX0(fs,fs)
    writeImage(fsfile,fs)
  else:
    fs = readImage(fsfile)
  plot3(gx,g=fs,cmin=-0.5,cmax=0.5,png="saltBody")
  plot3(gx,fbs=fs,png="saltBound")
def goPikSlices():
  fx = readImage(gxfile)
  plot3(fx,png="seis")
  pa = readImage(pafile)
  pm = max(pa)*0.5
  for i3 in range(n3):
    for i1 in range(n1):
      pa[i3][0][i1] = pm
      pa[i3][n2-1][i1] = pm
  zs,ys,xs=getPiks()
  lgs = getLineGroups(2,zs,ys,xs)
  plot3(pa,cmin=0.5,cmax=max(pa)*0.8,lgs=lgs,png="slicesInitial")
  zus = []
  yus = []
  xus = []
  xrs = []
  yrs = []
  zrs = []
  for ix in range(len(xs)):
    k3 = xs[ix]
    c2 = ys[ix]
    c1 = zs[ix]
    gx = fx[k3]
    sp = SaltPicker2()
    xu = sp.initialBoundary(1,c1,c2)
    bs = sp.refine(50,1,10,1,xu,pa[k3])
    xp = copy(xu)
    xm = copy(xu)
    xt = sp.regridBoundary(0.5,xu[0],xu[1])
    xus.append(k3)
    zus.append(xt[0])
    yus.append(xt[1])
    xrs.append(k3)
    zrs.append(xt[0])
    yrs.append(xt[1])
    e3 = min(k3+25,n3)
    for i3 in range(k3+1,e3,1):
      xp = sp.pickNext(5,1,3,1.0,xp[0],xp[1],pa[i3])
      xt = sp.regridBoundary(0.5,xp[0],xp[1])
      xus.append(i3)
      yus.append(xt[1])
      zus.append(xt[0])
    e3 = max(k3-25,-1)
    for i3 in range(k3-1,e3,-1):
      xm = sp.pickNext(5,1,3,1.0,xm[0],xm[1],pa[i3])
      xt = sp.regridBoundary(0.5,xm[0],xm[1])
      xus.append(i3)
      zus.append(xt[0])
      yus.append(xt[1])
  lgr = getLineGroups(2,zrs,yrs,xrs)
  lgu = getLineGroups(1,zus,yus,xus)
  plot3(fx,lgs=lgr,png="slicesFinal")
  plot3(fx,lgs=lgu,png="piks")
  '''
  wx = zerofloat(n1,n2,n3)
  isr = ImplicitSurfaceReconstructor()
  isr.pointsToImage(xus,yus,zus,pa,wx)
  writeImage(wxfile,wx)
  '''

def getLineGroups(dx,zs,ys,xs):
  lgs = []
  for ic in range(0,len(ys),dx):
    xyz = []
    rgb = []
    for ip in range(len(ys[ic])):
      xyz.append(xs[ic])
      xyz.append(ys[ic][ip])
      xyz.append(zs[ic][ip])
      rgb.append(1)
      rgb.append(0)
      rgb.append(0)
    lg = LineGroup(xyz,rgb)
    lgs.append(lg)
  return lgs

def goTF():
  fx = readImage(gxfile)
  sp = SaltPicker2()
  ks = rampint(60,1,30)
  ft = copy(fx)
  sp.applyTF(90,10,100,ks,fx[550],ft[550])
  plot3(fx)
  plot3(ft)
def goEnvAndSaltLike():
  pa = readImage(pafile)
  sl = readImage(slfile)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  '''
  p2 = zerofloat(n1,n2,n3)
  p3 = zerofloat(n1,n2,n3)
  ep = zerofloat(n1,n2,n3)
  lsf = LocalSlopeFinder(8,2,10)
  lsf.findSlopes(sl,p2,p3,ep)
  print min(p2)
  print min(p3)
  writeImage(p2file,p2)
  writeImage(p3file,p3)
  '''
  sp = SaltPicker2()
  pmin = 2.5
  sp.combineEnvAndSaltLike(pmin,p2,p3,pa,sl)
  writeImage(slefile,sl)
  plot3(pa)
  plot3(sl)

def goPik():
  fx = readImage(gxfile)
  pa = readImage(slefile)
  k3 = 442
  p3 = pa[k3]
  # pick the first slice
  gx = fx[k3]
  c1 = [680,144, 96,526,530,580,690,680]
  c2 = [410,326,120, 60,202,259,390,410]
  sp = SaltPicker2()
  #pa = sp.applyForInsAmp(gx)
  #pa = mul(pa,sl[k3])
  pm = max(p3)*0.5
  for i1 in range(n1):
    p3[0   ][i1] = pm
    p3[n2-1][i1] = pm
  xu = sp.initialBoundary(1,c1,c2)
  #xu = sp.regridBoundary(1,[c1,c2])
  plot(gx,cmin=-2,cmax=2)
  plot(gx,cmin=-2,cmax=2,pp=[c1,c2])
  bs = sp.refine(95,1,20,2,xu,p3)
  plot(gx,cmin=-2,cmax=2,xp=[xu[0],xu[1]])
  for i3 in range(443,200,-1):
    gn = fx[i3]
    p3 = pa[i3]
    for i1 in range(n1):
      p3[0   ][i1] = pm
      p3[n2-1][i1] = pm
    xu = sp.pickNext(5,1,5,1.0,xu[0],xu[1],p3)
    if(i3%10==0):
      plot(p3,cmin=0.001,cmax=0.5,xp=[xu[0],xu[1]],title="Slice"+str(i3))

def goEnvelope():
  gx = readImage(gxfile)
  rgf = RecursiveGaussianFilter(1)
  rgf.apply000(gx,gx)
  ge = zerofloat(n1,n2,n3)
  FastLevelSet3.applyForInsAmp(gx,ge)
  #writeImage(gefile,ge)
  plot3(gx,cmin=min(gx)/2,cmax=max(gx)/2)
  plot3(ge,cmin=min(ge),cmax=max(ge)/2)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

def getPiks():
  y1 = [  0,105,264,337,360,385,429,489,549,549,144,110, 72, 35,  0,  0]
  z1 = [410,145,318,236,226,334,320,315,338,733,445,341,368,631,641,410]
  y2 = [  0, 43,135,241,334,436,498,549,549,288,112, 47,  0,  0]
  z2 = [408,392,126,305,172,352,300,310,757,545,375,442,513,408]
  y3 = [  0,125,235,325,417,549,549,115,  0,  0]
  z3 = [440,121,303,166,404,318,785,386,507,440]
  y4 = [  0, 37,123,194,241,340,380,479,549,549,442,319,215,102,  0,  0]
  z4 = [189,258,136,172,295,169,467,577,383,696,787,577,486,367,453,189]
  y5 = [  0, 70,110,275,296,319,362,415,484,549,549,414,273,102,  0,  0]
  z5 = [157,138,158,193,275,306,235,502,527,347,751,750,559,415,468,157]
  y6 = [  0,259,293,351,490,549,549,415,272,132,  0,  0]
  z6 = [ 90,213,270,292,388,316,760,718,604,469,498, 90]
  y7 = [  0,304,421,484,549,549,227,218,139,  0,  0]
  z7 = [ 89,246,350,544,526,775,650,650,490,530, 89]
  y8 = [  0,365,393,395,443,549,549,463,272,169,134,  0,  0]
  z8 = [ 60,192,420,529,688,708,784,825,688,634,520,561, 60]
  y9 = [  0,294,314,353,379,439,549,549,402,384,225,153, 38, 42,  0,  0]
  z9 = [ 48,122,104,246,623,747,745,800,773,743,683,511,584,141,140, 48]
  y10 = [ 96,142,314,430,549,549,420,320,182, 96]
  z10 = [486, 80,128,710,735,802,773,619,546,486]
  y11= [141,172,332,439,482,549,549,434,324,187,141]
  z11= [408, 89,170,634,733,710,782,740,558,566,408]
  y12= [170,151,196,384,423,483,419,333,170]
  z12= [588,406, 71,153,609,735,703,559,588]
  y13= [196,163,194,365,432,420,493,549,549,420,196]
  z13= [614,561, 63,159,148,487,647,643,718,586,614]
  y14= [154,224,385,495,474,500,360,300,154]
  z14= [647, 87,145, 66,467,617,572,566,647]
  y15= [187,296,549,549,481,321,187]
  z15= [632,103, 84,678,685,556,632]
  y16= [243,290,253,253,324,380,433,496,549,549,440,356,243]
  z16= [630,291,139, 71, 88, 66, 99, 56, 68,697,699,614,630]

  y17= [232,266,253,207,282,315,381,489,549,549,232]
  z17= [600,250,174, 58, 30, 94, 41, 36, 60,679,600]

  y18= [139,193,166,197,279,329,389,432,501,549,549,139]
  z18= [616,204, 66, 21, 24, 64, 26, 87, 67, 37,616,616]
  zs = [z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18]
  ys = [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15,y16,y17,y18]
  xs = [  0, 50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,849]
  return zs,ys,xs

#############################################################################
# graphics

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
  k1 = [408]
  k2 = [240,470]
  k3 = [0,100,200,300,400,500,600,700,800]
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
    tg = TriangleGroup(ct.i,ct.x,ct.u)
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
  if cbar:
    sf.setSize(987,700)
  else:
    sf.setSize(850,700)
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  ov.setAxesScale(1.0,1.0,zscale)
  ov.setScale(1.6)
  ov.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  ov.setTranslate(Vector3(-0.02,0.13,-0.11))
  ov.setAzimuthAndElevation(-148.0,45.0)
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

#############################################################################
run(main)
